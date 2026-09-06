# Output

A finished run leaves you three things: the [`Population`](@ref) (what is alive now), the
lineage tree reachable from it (everything that led there), and — if you asked for one —
a [`Measurements`](@ref) record of what happened along the way.

```julia
simulate!(pop, block, rng)                        # no recording
simulate!(pop, block, rng; accumulator = acc)     # with recording
```

`simulate!` mutates `pop` in place and returns it, so the two forms differ only in
whether anything was written down as it ran.

## The population

```julia
popsize(pop)              # number of living cells
allcells(pop)             # Vector{BinaryNode{NonMarkovCell}} of the living cells
pop.t                     # time of the most recently processed event
age(pop)                  # === pop.t
pop                       # Population: 10000 cells (t = 14.212)
```

`pop.cells` is a `Dict` keyed by cell id. Two things follow from that:

- **Iteration order is arbitrary.** `allcells(pop)`, `fitness_per_cell(pop)` and
  `mutations_per_cell(pop)` all iterate that same `Dict`, so they are co-indexed with
  each other — entry `i` is the same cell in all three — but the order itself is not
  meaningful and is not stable across Julia versions or insertion histories. If you need
  a stable order, work from the tree instead; see
  [Which results are co-indexed](statistics.md#Which-results-are-co-indexed).
- **You can inspect it, but editing it breaks the event queue.** Adding or removing
  entries by hand leaves the heap holding a different number of events than there are
  cells, and the next `simulate!` throws naming both counts.
  [`reset_schedule!`](@ref) is the documented recovery.

## The tree

The tree is not stored separately — it is reachable from any living cell by following
`parent`:

```julia
root = getsingleroot(allcells(pop))     # nothing if the population is a forest
```

From there, `AbstractTrees` works directly (`Leaves`, `PreOrderDFS`, `nodevalue`,
`children`, `print_tree`), and so does every root-based statistic on the
[Tree statistics](statistics.md) page.

```julia
node.data          # the NonMarkovCell
node.parent        # nothing at the root
node.left, node.right   # nothing at a leaf; one may be nothing at a unary node

endtime(node)      # when this cell divided; nothing if it is still alive
celllifetime(node, pop.t)   # birth to division, or birth to now
celllifetimes(root)         # every completed cell lifetime in the tree
getalivecells(root)         # the living leaves
popsize(root)               # how many there are
age(root)                   # birthtime of the most recently born leaf
```

!!! warning "`age(root)` is not `pop.t`"
    `age(::BinaryNode)` returns the birthtime of the **last-born leaf**, which is the time
    of the last *division*. `pop.t` is the time of the last event of any kind, so a run
    that ended on a death has `age(root) < pop.t`. `celllifetime(node)` with no second
    argument uses `age(getroot(node))`, so pass `pop.t` explicitly when you mean "how old
    is this cell right now":

    ```julia
    celllifetime(node, pop.t)
    ```

## Recording while it runs

[`MeasurementSpec`](@ref) declares what to collect;
`MeasurementAccumulator` collects it; [`finalize_measurements`](@ref) packages it.

```julia
spec = MeasurementSpec(
    trajectory_dt     = 0.5,
    snapshot_triggers = [AtTime(5.0), AtPopSize(1_000), AtEnd()],
    snapshot_stats    = [SFS(), FitnessDistribution(), DriversPerCell()],
)
acc = MeasurementAccumulator(spec)

simulate!(pop, block, rng; accumulator = acc)
m = finalize_measurements(acc)

m.trajectory     # Vector{TrajectoryPoint}
m.snapshots      # Vector{SnapshotData}
```

Defaults are `trajectory_dt = Inf` (no trajectory), `snapshot_triggers = [AtEnd()]` and
`snapshot_stats = [FitnessDistribution()]`.

Call `finalize_measurements` once, after `simulate!` has returned. It copies, so the
accumulator can keep being used afterwards.

### The trajectory

A [`TrajectoryPoint`](@ref) is recorded every `trajectory_dt` units of simulation time:

| Field | Meaning |
|:---|:---|
| `t` | the grid time this point is labelled with |
| `N_total` | population size |
| `mean_fitness`, `var_fitness` | over living cells |
| `mean_k`, `var_k` | total driver burden per living cell |

!!! note "`t` is the grid time; the state is the first one at or after it"
    Recording is driven by events, not by a clock. When an event lands past one or more
    grid points, each is emitted with its own grid `t` but with the population state as
    of that event. Points are therefore exactly `trajectory_dt` apart in `t`, while the
    state attached to them lags by up to one inter-event interval — negligible once the
    population is large, visible at the very start of a run from a single cell.

    Nothing is recorded once the population is empty, so an extinct run's trajectory
    simply stops.

!!! warning "Trajectory recording is the expensive part of a run"
    Every point recomputes `mutations_per_cell`, which walks each living cell's whole
    ancestral path — ``O(N \times \text{depth})`` per point. At ``N = 10^5`` and depth
    ``\approx 20`` that is a few million pointer hops for one point. A `trajectory_dt`
    fine enough to resolve the early growth phase will dominate the runtime of the late
    one. If you need both, record the early phase with a fine `dt` in a first block and
    the late phase with a coarse one in a chained second block.

### Snapshot triggers

| Trigger | Fires |
|:---|:---|
| [`AtEnd`](@ref)`()` | once, when `simulate!` returns |
| [`AtTime`](@ref)`(t)` | the first event at or after simulation time `t` |
| [`AtPopSize`](@ref)`(N)` | the first event at which population size reaches `N` |

Each trigger fires **at most once** per run. `AtTime` and `AtPopSize` are tested after
every event; `AtEnd` fires when the block exits, whether it exited on the stop condition
or on extinction.

Two details worth knowing:

- Triggers are checked *after* an event, never before the first one. `AtPopSize(N₀)` set
  to the starting size therefore fires on the first event, when the size is already
  ``N_0 \pm 1``.
- If the block restarts on extinction, the accumulator is reset — trajectory, snapshots
  and fired-trigger record all cleared — so what you get back describes the **successful
  attempt only**, with times measured from the restart point.

`AtTime` is also the right tool when you want the state at a specific time without ending
the run there, since a time-based `stopfunction` overshoots by one event (see
[Growth modes](blocks.md#Growth-modes:-the-stop-condition)).

### Snapshot statistics

| Statistic | Field on [`SnapshotData`](@ref) | Contents |
|:---|:---|:---|
| [`SFS`](@ref)`()` | `.sfs` | `Vector{Int64}`; `sfs[k]` = driver mutation events carried by exactly `k` living cells |
| [`FitnessDistribution`](@ref)`()` | `.fitness_distribution` | `Vector{Float64}`, one entry per living cell |
| [`DriversPerCell`](@ref)`()` | `.drivers_per_cell` | `Vector{Int64}`, total burden per living cell |

Fields you did not request are `nothing`, so check before use:

```julia
for snap in m.snapshots
    println(snap.trigger, " at t = ", snap.t)
    isnothing(snap.sfs) || println("  SFS: ", snap.sfs[1:min(end, 10)])
end
```

`snap.trigger` is the trigger object itself, so a run with several triggers stays
attributable:

```julia
endsnap = only(s for s in m.snapshots if s.trigger isa AtEnd)
```

## Getting more than the three built-in statistics

The snapshot statistics are a fixed set. Anything else is computed after the fact from
the tree, which loses nothing, because the tree is a complete record of the survivors:

```julia
root = getsingleroot(allcells(pop))

branch_spectrum(root)                  # topology, free of the mutation rate
leaf_depths(root)                      # divisions per surviving cell
leaf_fitness(root)                     # fitness per surviving leaf, in tree order
filtered_mutations_per_cell(root, 0.1) # burden, high-frequency variants dropped
coalescence_times(pop)                 # pairwise times to MRCA
celllifetimes(root)                    # realised cell-cycle durations
```

The one thing you cannot recover afterwards is a quantity at an **intermediate** time —
the tree only holds the present. If you need a mid-run observable that is not one of the
three, take it with an `on_division` hook, which sees the population at every division
(see [Hooks](blocks.md#on_division)):

```julia
sizes = Tuple{Float64, Int}[]
on_division = function (pop, parent, d1, d2)
    push!(sizes, (pop.t, popsize(pop)))
    return nothing
end
```

A hook that does not draw from `rng` leaves the random stream untouched, so adding one to
an existing run reproduces it exactly.

## Persisting a run

`Population`, `BinaryNode` and `NonMarkovCell` are plain Julia types with no external
resources, so `Serialization` round-trips them:

```julia
using Serialization
serialize("run.jls", pop)
pop = deserialize("run.jls")
```

Two caveats. A serialised `Population` carries its whole tree — every internal node back
to the founder — so the file scales with the total number of divisions, not with
`popsize`. And `_pending` holds live event objects, which deserialize fine but describe a
schedule drawn under the block you were running; call [`reset_schedule!`](@ref) before
continuing a deserialised population under a different block.

To store the observables rather than the run, serialise the statistics — or a
[`LeafSample`](@ref), which is a real tree at a fraction of the size and replays exactly
from its recorded seed. See [Sampling](sampling.md).
