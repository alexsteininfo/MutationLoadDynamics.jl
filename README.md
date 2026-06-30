# MutationLoadDynamics.jl

Julia package for simulating somatic evolution under **non-Markovian** (arbitrary waiting-time distribution) birth-death dynamics with per-cell fitness tracking.

Unlike Gillespie-based simulators where waiting times are exponentially distributed, this package allows division and death times to follow any distribution (Gamma, Weibull, log-normal, …). Each driver mutation draws a random fitness increment, so individual cells accumulate fitness independently — there are no shared fitness subclones.

## Key features

- **Non-Markovian dynamics** — division and death waiting times drawn from arbitrary user-supplied distributions (functions of per-cell fitness)
- **Per-cell fitness** — each driver mutation draws `δ ~ driver_dist`; fitness updated via a user-supplied rule (`f → fitness_update(f, δ)`)
- **Global min-heap event queue** — O(log N) per event, no acceptance-rejection step
- **Full binary tree preserved** — entire lineage tree stored as `BinaryNode{NonMarkovCell}`; supports SFS, pairwise distances, coalescence times
- **Competing-risks birth/death** — both waiting times sampled at birth; the earlier fires

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/alexanderstein/MutationLoadDynamics.jl")
```

## Quick start

```julia
using MutationLoadDynamics, Distributions, Random

rng = MersenneTwister(42)

pop = initialize_population(fitness_init = 1.0)

block = NonMarkovBlock(
    birth_dist     = f -> Gamma(2.0, 1.0 / f),     # mean division time = 2/f (fitter → faster)
    death_dist     = f -> Gamma(2.0, 5.0),          # mean death time = 10 (fitness-independent)
    stopfunction   = pop -> popsize(pop) >= 5_000,
    driver_dist    = Exponential(0.05),             # fitness increment per driver mutation
    fitness_update = (f, δ) -> f + δ,              # additive fitness model
    ν              = 0.5,                           # mean drivers per daughter per division
    restart_on_extinction = true,
)

simulate!(pop, block, rng)

println("Population size : ", popsize(pop))
println("Mean fitness    : ", mean(fitness_per_cell(pop)))
println("Mean drivers/cell: ", mean_k(pop))
```

## How it works

### Algorithm

At each birth event the two daughter cells immediately **pre-schedule** their next event via competing risks:

```julia
t_div  = birthtime + rand(rng, block.birth_dist(fitness))
t_die  = birthtime + rand(rng, block.death_dist(fitness))
event_time = min(t_div, t_die)          # whichever comes first
event_type = t_div ≤ t_die ? :birth : :death
```

All living cells' events are stored in a **global min-heap**. The simulation loop pops the cell with the smallest event time, processes its event (division or death), and pushes new events for daughter cells. This guarantees events are processed in globally correct order — no Gillespie maximum-rate computation or acceptance-rejection needed.

### Mutation model

At each division, each daughter cell independently:
1. Draws `j ~ Poisson(ν)` — number of new driver mutations
2. For each mutation: draws `δ ~ driver_dist`, then updates `f ← fitness_update(f, δ)`
3. Stores `j` as `NonMarkovCell.mutations` and the final `f` as `NonMarkovCell.fitness`

The `fitness_update` function is applied **once per mutation event**, not to the sum:

```julia
# Additive model
fitness_update = (f, δ) -> f + δ

# Multiplicative model
fitness_update = (f, δ) -> f * (1 + δ)
```

## API

### `NonMarkovBlock`

```julia
NonMarkovBlock(;
    birth_dist     = f -> Gamma(2.0, 1.0 / f),   # (fitness::Float64) -> Distribution
    death_dist     = f -> Gamma(2.0, 5.0),        # (fitness::Float64) -> Distribution
    stopfunction   = pop -> popsize(pop) >= 1000, # (pop::Population) -> Bool
    driver_dist    = Exponential(0.05),           # Distribution for Δ per driver mutation
    fitness_update = (f, δ) -> f + δ,            # (f::Float64, δ::Float64) -> Float64
    ν              = 0.5,                         # Poisson mean drivers per daughter per division
    restart_on_extinction = false,                # restart from initial state on extinction
)
```

| Parameter | Type | Description |
|---|---|---|
| `birth_dist` | `f -> Distribution` | Waiting-time distribution for division; receives the cell's current fitness |
| `death_dist` | `f -> Distribution` | Waiting-time distribution for death |
| `stopfunction` | `pop -> Bool` | Simulation stop criterion (time, population size, or arbitrary) |
| `driver_dist` | `Distribution` | Distribution from which each driver mutation's fitness increment `δ` is drawn |
| `fitness_update` | `(f, δ) -> f` | How a single driver mutation changes fitness; applied once per mutation event |
| `ν` | `Float64` | Mean number of driver mutations per daughter cell per division (Poisson) |
| `restart_on_extinction` | `Bool` | If `true`, restart from the initial population on extinction (default `false`) |

### `simulate!`

```julia
simulate!(pop, block, rng)
simulate!(pop, block, rng; accumulator = acc)   # with trajectory/snapshot recording
```

Runs the simulation in-place and returns `pop`. The RNG defaults to `Random.GLOBAL_RNG`.

### `initialize_population`

```julia
# Single founding cell
pop = initialize_population(fitness_init = 1.0)

# N identical cells (useful for starting from a pre-existing pool)
pop = initialize_population(100; fitness_init = 1.5)
```

## Cell representation

Cells are `BinaryNode{NonMarkovCell}` nodes forming a binary tree. Dead cells are pruned
from the tree; internal nodes (cells that divided) are preserved for ancestry tracking.

```julia
struct NonMarkovCell
    id::Int64         # unique identifier
    birthtime::Float64
    mutations::Int64  # new driver mutations at this cell's birth (Poisson draw)
    fitness::Float64  # cumulative fitness (all ancestors + self)
end
```

The `Population` keeps **only alive cells** in a `Dict{Int64, BinaryNode{NonMarkovCell}}`
for O(1) insertion and removal:

```julia
popsize(pop)          # number of alive cells
allcells(pop)         # Vector of all alive BinaryNodes
pop.t                 # current simulation time
```

## Measurements and trajectory recording

```julia
spec = MeasurementSpec(
    trajectory_dt     = 0.5,           # record a TrajectoryPoint every 0.5 time units
    snapshot_triggers = [AtEnd()],     # also take a full snapshot at simulation end
    snapshot_stats    = [SFS(), FitnessDistribution(), DriversPerCell()],
)
acc = MeasurementAccumulator(spec)
simulate!(pop, block, rng; accumulator = acc)
m = finalize_measurements(acc)

m.trajectory   # Vector{TrajectoryPoint}
m.snapshots    # Vector{SnapshotData}
```

### Trajectory

Each `TrajectoryPoint` (recorded every `trajectory_dt` time units) contains:

| Field | Description |
|---|---|
| `t` | Simulation time |
| `N_total` | Population size |
| `mean_fitness` | Mean fitness across alive cells |
| `var_fitness` | Variance of fitness |
| `mean_k` | Mean total driver count per cell |
| `var_k` | Variance of driver count |

### Snapshot triggers

| Trigger | When it fires |
|---|---|
| `AtEnd()` | Once when the simulation exits |
| `AtTime(t)` | First time simulation time reaches `t` |
| `AtPopSize(N)` | First time population size reaches `N` |

### Snapshot statistics

| Statistic | `SnapshotData` field | Description |
|---|---|---|
| `SFS()` | `.sfs` | Driver mutation site-frequency spectrum: `sfs[k]` = number of driver mutations present in exactly `k` alive cells |
| `FitnessDistribution()` | `.fitness_distribution` | Fitness values of all alive cells |
| `DriversPerCell()` | `.drivers_per_cell` | Total driver count per alive cell (summed along lineage) |

## Analysis functions

### Fitness

```julia
fitness_per_cell(pop)    # Vector{Float64} — fitness of each alive cell
fitness_distribution(pop) # alias for fitness_per_cell
mean_k(pop)              # mean total driver count per alive cell
var_k(pop)               # variance of driver count
```

### Mutation burden

```julia
mutations_per_cell(pop)  # Vector{Int64} — total driver mutations per cell (summed to root)
average_mutations(pop)   # Float64 — mean across population
clonal_mutations(pop)    # Int64 — driver mutations at the MRCA (shared by all cells)
```

### Site-frequency spectrum

```julia
sitefrequencyspectrum(pop)   # Vector{Int64} — driver mutation SFS
                              # sfs[k] = number of driver mutation events shared by exactly k cells
```

### Pairwise distances

```julia
pairwisedistance(cell1, cell2)     # Int64 — driver mutation differences between two cells
pairwisedistances(pop[, idx])      # Vector{Int64} — all pairwise distances
pairwise_differences(pop[, idx])   # Dict{Int64,Int64} — histogram of distances
```

### Coalescence

```julia
coalescence_times(pop[, idx])         # Vector{Float64} — time to MRCA for every cell pair
coalescence_times(root::BinaryNode)   # same, from a specific subtree root
```

### Tree utilities

```julia
findMRCA(pop)              # BinaryNode — most recent common ancestor of all alive cells
findMRCA(node1, node2)     # MRCA of two specific cells
getsingleroot(cells)       # unique root node, or nothing if multiple trees
endtime(node)              # Float64 — when the cell divided (nothing if still alive)
celllifetime(node)         # Float64 — cell age from birth to division (or current time)
celllifetimes(root)        # Vector{Float64} — all division lifetimes in a subtree
age(pop)                   # Float64 — simulation time (pop.t)
getalivecells(root)        # alive leaf nodes under a given root
```

## Example

A worked example is in [`examples/GrowthWithDrivers.jl`](examples/GrowthWithDrivers.jl):

```julia
julia --project=. examples/GrowthWithDrivers.jl
```

The example simulates growth to 5 000 cells with additive fitness effects, records a
full trajectory, and prints summary statistics.

## Background

This package targets the **high driver-mutation-rate** regime where cells accumulate many
drivers per generation (ν ≥ 0.1). In this regime:

- Subclone identities are not meaningful (every cell has a unique fitness history)
- The Gillespie algorithm is inappropriate because division times are not exponential
- Per-cell fitness tracking is required

The **min-heap event queue** (Next Reaction Method for non-Markovian processes) is the natural algorithm: each cell pre-samples its own next-event time from its fitness-dependent distribution, and the simulation always processes the globally earliest event. This is exact (no approximation) and scales as O(N log N) per simulation run.

For the low driver-mutation-rate regime (ν ≪ 1) with clone-level dynamics and neutral
passenger mutations, see [`BirthDeathMutation.jl`](https://github.com/alexanderstein/BirthDeathMutation.jl).
