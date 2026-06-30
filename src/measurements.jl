# ── Trigger types ─────────────────────────────────────────────────────────────

abstract type AbstractTrigger end

"""Fire once when `simulate!` exits (stop condition met or N = 0)."""
struct AtEnd <: AbstractTrigger end

"""Fire once when simulation time first reaches or exceeds `t`."""
struct AtTime <: AbstractTrigger
    t::Float64
end

"""Fire once when population size first reaches or exceeds `N`."""
struct AtPopSize <: AbstractTrigger
    N::Int
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
- `trajectory_dt::Float64` — time between trajectory points; `Inf` disables trajectory recording.
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

One sample in a continuously recorded population trajectory.
Collected every `MeasurementSpec.trajectory_dt` simulation-time units.
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

mutable struct MeasurementAccumulator
    spec::MeasurementSpec
    next_trajectory_t::Float64
    trajectory_points::Vector{TrajectoryPoint}
    snapshots::Vector{SnapshotData}
    fired_triggers::Set{Int}
end

MeasurementAccumulator(spec::MeasurementSpec) =
    MeasurementAccumulator(spec, 0.0, TrajectoryPoint[], SnapshotData[], Set{Int}())

# ── Internal helpers ──────────────────────────────────────────────────────────

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

function record_trajectory_if_due!(acc::MeasurementAccumulator, pop::Population)
    isinf(acc.spec.trajectory_dt) && return
    t = pop.t
    N = popsize(pop)
    N == 0 && return
    while t >= acc.next_trajectory_t
        fitnesses = fitness_per_cell(pop)
        ks        = Float64.(mutations_per_cell(pop))
        push!(acc.trajectory_points, TrajectoryPoint(
            acc.next_trajectory_t, N,
            mean(fitnesses), var(fitnesses),
            mean(ks), var(ks),
        ))
        acc.next_trajectory_t += acc.spec.trajectory_dt
    end
end

function check_timed_triggers!(acc::MeasurementAccumulator, pop::Population)
    t = pop.t
    N = popsize(pop)
    for (i, trigger) in enumerate(acc.spec.snapshot_triggers)
        i in acc.fired_triggers && continue
        fire = (trigger isa AtTime    && t >= trigger.t) ||
               (trigger isa AtPopSize && N >= trigger.N)
        if fire
            push!(acc.fired_triggers, i)
            push!(acc.snapshots,
                _compute_snapshot(acc.spec.snapshot_stats, trigger, pop, Float64(t)))
        end
    end
end

function _fire_end_triggers!(acc::MeasurementAccumulator, pop::Population)
    for (i, trigger) in enumerate(acc.spec.snapshot_triggers)
        if trigger isa AtEnd && !(i in acc.fired_triggers)
            push!(acc.fired_triggers, i)
            push!(acc.snapshots,
                _compute_snapshot(acc.spec.snapshot_stats, trigger, pop, pop.t))
        end
    end
end

function _reset_accumulator!(acc::MeasurementAccumulator)
    empty!(acc.trajectory_points)
    empty!(acc.snapshots)
    empty!(acc.fired_triggers)
    acc.next_trajectory_t = 0.0
end

# ── Public API ────────────────────────────────────────────────────────────────

"""
    finalize_measurements(acc::MeasurementAccumulator) -> Measurements

Package collected trajectory points and snapshots into a `Measurements` object.
Call once after `simulate!` has returned.
"""
finalize_measurements(acc::MeasurementAccumulator) =
    Measurements(copy(acc.trajectory_points), copy(acc.snapshots))
