# ── Trigger types ─────────────────────────────────────────────────────────────

abstract type AbstractTrigger end

"""
Fire when a `simulate!` call exits (stop condition met or extinction) — once per call,
so an accumulator carried across chained blocks gets one snapshot per block.
"""
struct AtEnd <: AbstractTrigger end

"""
Fire once with the population state exactly at simulation time `t`, i.e. after every
event at or before `t` and before any later one. The snapshot is labelled `t`.
"""
struct AtTime <: AbstractTrigger
    t::Float64
    function AtTime(t::Real)
        isnan(t) && throw(ArgumentError("AtTime: t must not be NaN"))
        return new(t)
    end
end

"""Fire once, at the first moment the population size is at least `N`."""
struct AtPopSize <: AbstractTrigger
    N::Int
    function AtPopSize(N::Integer)
        N >= 1 || throw(ArgumentError("AtPopSize: N must be >= 1, got $N"))
        return new(N)
    end
end

# ── Statistic types ───────────────────────────────────────────────────────────

abstract type AbstractStatistic end

"""Driver mutation site-frequency spectrum: `sfs[k]` = number of driver mutations in exactly k cells."""
struct SFS <: AbstractStatistic end

"""Full fitness distribution: fitness value of every alive cell."""
struct FitnessDistribution <: AbstractStatistic end

"""Total driver mutations per alive cell, summed along each cell's lineage."""
struct DriversPerCell <: AbstractStatistic end

# ── User-facing specification ─────────────────────────────────────────────────

"""
    MeasurementSpec(; trajectory_dt, snapshot_triggers, snapshot_stats)

Declares what to record during `simulate!`.

# Keyword arguments
- `trajectory_dt::Float64` — time between trajectory points (must be `> 0`); `Inf`
  disables trajectory recording.
- `snapshot_triggers` — vector of `AbstractTrigger`s specifying when to take snapshots.
- `snapshot_stats` — vector of `AbstractStatistic`s specifying what to compute at each snapshot.

# Example
```julia
spec = MeasurementSpec(
    trajectory_dt     = 0.5,
    snapshot_triggers = [AtEnd()],
    snapshot_stats    = [SFS(), FitnessDistribution(), DriversPerCell()],
)
acc  = MeasurementAccumulator(spec)
simulate!(pop, block, rng; accumulator = acc)
m    = finalize_measurements(acc)
```
"""
struct MeasurementSpec
    trajectory_dt::Float64
    snapshot_triggers::Vector{AbstractTrigger}
    snapshot_stats::Vector{AbstractStatistic}

    function MeasurementSpec(trajectory_dt::Real, snapshot_triggers, snapshot_stats)
        # `> 0` also rejects NaN. A zero step would make trajectory recording loop forever.
        trajectory_dt > 0 || throw(ArgumentError(
            "MeasurementSpec: trajectory_dt must be > 0 (Inf disables it), got $trajectory_dt"))
        return new(trajectory_dt, snapshot_triggers, snapshot_stats)
    end
end

function MeasurementSpec(;
    trajectory_dt     = Inf,
    snapshot_triggers = [AtEnd()],
    snapshot_stats    = [FitnessDistribution()],
)
    return MeasurementSpec(
        Float64(trajectory_dt),
        Vector{AbstractTrigger}(snapshot_triggers),
        Vector{AbstractStatistic}(snapshot_stats),
    )
end

# ── Output types ──────────────────────────────────────────────────────────────

"""
    TrajectoryPoint

One sample in a continuously recorded population trajectory, taken every
`MeasurementSpec.trajectory_dt` simulation-time units. Each point is the exact state at
its grid time `t`. `k` is a cell's total driver count; the variances are `NaN` while only
one cell is alive.
"""
struct TrajectoryPoint
    t::Float64
    N_total::Int
    mean_fitness::Float64
    var_fitness::Float64
    mean_k::Float64
    var_k::Float64
end

"""
    SnapshotData

Full statistical snapshot taken when a trigger fires.
Fields are `nothing` when the corresponding statistic was not requested.
"""
struct SnapshotData
    t::Float64
    trigger::AbstractTrigger
    sfs::Union{Vector{Int64}, Nothing}
    fitness_distribution::Union{Vector{Float64}, Nothing}
    drivers_per_cell::Union{Vector{Int64}, Nothing}
end

"""
    Measurements

All data collected during one `simulate!` call:
- `trajectory` — time-series of population state (empty if `trajectory_dt = Inf`)
- `snapshots`  — full statistical snapshots at each trigger event
"""
struct Measurements
    trajectory::Vector{TrajectoryPoint}
    snapshots::Vector{SnapshotData}
end

# ── Internal accumulator ──────────────────────────────────────────────────────

"""
    MeasurementAccumulator(spec::MeasurementSpec)

Mutable collector that `simulate!` writes into. Build one from a
[`MeasurementSpec`](@ref), pass it as the `accumulator` keyword, and call
[`finalize_measurements`](@ref) afterwards to get an immutable [`Measurements`](@ref).

```julia
acc = MeasurementAccumulator(spec)
simulate!(pop, block, rng; accumulator = acc)
m   = finalize_measurements(acc)
```

One accumulator can be carried across chained `simulate!` calls: trajectory recording
resumes at the current `population.t` rather than back-filling from zero, `AtTime` and
`AtPopSize` triggers that already fired do not fire again, and `AtEnd` fires at the end
of every call. An extinction restart discards only what was recorded during the current
call, so earlier blocks' records survive and the current block's describe its successful
attempt.

Its fields are internal bookkeeping — read results from `finalize_measurements`, which
copies, leaving the accumulator usable afterwards.
"""
mutable struct MeasurementAccumulator
    spec::MeasurementSpec
    next_trajectory_t::Float64
    trajectory_points::Vector{TrajectoryPoint}
    snapshots::Vector{SnapshotData}
    fired_triggers::Set{Int}
    # Time at which the current `simulate!` call started. An `AtTime` earlier than this
    # cannot be recorded exactly any more, so it never fires.
    call_start_t::Float64
end

MeasurementAccumulator(spec::MeasurementSpec) =
    MeasurementAccumulator(spec, 0.0, TrajectoryPoint[], SnapshotData[], Set{Int}(), 0.0)

# Called by `simulate!` on entry.
function _start_call!(acc::MeasurementAccumulator, pop::Population)
    acc.call_start_t = pop.t
    # Start recording at the current time: a fresh accumulator handed to a chained call
    # would otherwise back-fill points from t = 0 up to pop.t.
    acc.next_trajectory_t = max(acc.next_trajectory_t, pop.t)
    return acc
end

# ── Internal helpers ──────────────────────────────────────────────────────────
#
# Timing convention: the population state is right-continuous — the state "at time t"
# includes every event at or before t. `simulate!` calls `_record_until!` with the
# time of the next event *before* applying it, so every grid time and every `AtTime`
# strictly earlier than that event is recorded with the state that held there exactly.

function _compute_snapshot(
    stats::Vector{AbstractStatistic},
    trigger::AbstractTrigger,
    pop::Population,
    t::Float64,
)
    sfs_val     = nothing
    fitness_val = nothing
    drivers_val = nothing
    for stat in stats
        if stat isa SFS
            sfs_val = sitefrequencyspectrum(pop)
        elseif stat isa FitnessDistribution
            fitness_val = fitness_per_cell(pop)
        elseif stat isa DriversPerCell
            drivers_val = mutations_per_cell(pop)
        end
    end
    return SnapshotData(t, trigger, sfs_val, fitness_val, drivers_val)
end

function _push_snapshot!(acc::MeasurementAccumulator, i::Int, trigger, pop, t)
    push!(acc.fired_triggers, i)
    push!(acc.snapshots, _compute_snapshot(acc.spec.snapshot_stats, trigger, pop, t))
end

# Record every trajectory point and fire every `AtTime` trigger strictly before `t_next`
# (or at or before it when `inclusive`), using the current state.
function _record_until!(acc::MeasurementAccumulator, pop::Population, t_next::Float64;
                        inclusive::Bool = false)
    due(t) = inclusive ? t <= t_next : t < t_next
    dt = acc.spec.trajectory_dt
    if !isinf(dt) && due(acc.next_trajectory_t) && popsize(pop) > 0
        fitnesses = fitness_per_cell(pop)
        ks        = Float64.(mutations_per_cell(pop))
        mf, vf    = mean(fitnesses), var(fitnesses)
        mk, vk    = mean(ks), var(ks)
        while due(acc.next_trajectory_t)
            push!(acc.trajectory_points,
                  TrajectoryPoint(acc.next_trajectory_t, popsize(pop), mf, vf, mk, vk))
            acc.next_trajectory_t += dt
        end
    end
    for (i, trigger) in enumerate(acc.spec.snapshot_triggers)
        trigger isa AtTime && !(i in acc.fired_triggers) && due(trigger.t) &&
            trigger.t >= acc.call_start_t &&
            _push_snapshot!(acc, i, trigger, pop, trigger.t)
    end
end

function _check_popsize_triggers!(acc::MeasurementAccumulator, pop::Population)
    N = popsize(pop)
    for (i, trigger) in enumerate(acc.spec.snapshot_triggers)
        trigger isa AtPopSize && !(i in acc.fired_triggers) && N >= trigger.N &&
            _push_snapshot!(acc, i, trigger, pop, pop.t)
    end
end

# At exit: the final state holds at `pop.t`, so grid points and `AtTime`s at exactly
# `pop.t` are still due; then every `AtEnd` fires (once per call, never marked fired).
function _finish_call!(acc::MeasurementAccumulator, pop::Population)
    _record_until!(acc, pop, pop.t; inclusive = true)
    for trigger in acc.spec.snapshot_triggers
        trigger isa AtEnd && push!(acc.snapshots,
            _compute_snapshot(acc.spec.snapshot_stats, trigger, pop, pop.t))
    end
end

# What the accumulator held when a `simulate!` call began, so that an extinction
# restart can discard exactly what that call recorded.
_checkpoint(acc::MeasurementAccumulator) = (
    n_points  = length(acc.trajectory_points),
    n_snaps   = length(acc.snapshots),
    fired     = copy(acc.fired_triggers),
    next_t    = acc.next_trajectory_t,
)

function _rollback!(acc::MeasurementAccumulator, cp)
    resize!(acc.trajectory_points, cp.n_points)
    resize!(acc.snapshots, cp.n_snaps)
    acc.fired_triggers    = copy(cp.fired)
    acc.next_trajectory_t = cp.next_t
    return acc
end

# ── Public API ────────────────────────────────────────────────────────────────

"""
    finalize_measurements(acc::MeasurementAccumulator) -> Measurements

Package collected trajectory points and snapshots into a `Measurements` object.
Call once after `simulate!` has returned.
"""
finalize_measurements(acc::MeasurementAccumulator) =
    Measurements(copy(acc.trajectory_points), copy(acc.snapshots))
