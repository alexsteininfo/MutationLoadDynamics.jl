"""
    NonMarkovBlock{F1,F2,F3,F4,F5,F6,D}

Defines a non-Markovian birth-death simulation. Division and death waiting times are
drawn from arbitrary distributions (functions of the cell's fitness). The block runs
until `stopfunction(pop)` returns `true`.

# Fields
- `birth_dist::F1` — `(fitness::Float64) -> Distribution` — waiting-time distribution
  for cell division; mean = expected time from birth to next division.
- `death_dist::F2` — `(fitness::Float64) -> Distribution` — waiting-time distribution
  for cell death.
- `stopfunction::F3` — `(pop::Population) -> Bool` — simulation stop criterion.
- `driver_dist::D` — distribution from which each driver mutation's fitness increment
  `δ` is drawn (e.g. `Exponential(0.05)`).
- `fitness_update::F4` — `(parent_fitness::Float64, δ::Float64) -> Float64` — how each
  individual driver mutation changes the cell's fitness; applied once per mutation event.
- `ν::Float64` — mean number of driver mutations per daughter cell per division
  (Poisson distributed).
- `restart_on_extinction::Bool` — if `true`, restart from the initial state whenever
  the population goes extinct (default: `false`).
- `on_division::F5` — optional `(pop, parent, d1, d2) -> Nothing` callback fired at every
  division, or `nothing` (default) to disable. See "Division hook" below.
- `on_restart::F6` — optional `(pop) -> Nothing` callback fired after the population is
  restored following an extinction, or `nothing` (default). Only ever called when
  `restart_on_extinction = true`. See "Restart hook" below.

# Division hook

`on_division` lets you apply a deterministic, state-triggered fitness change to a single
cell at a single moment — something the Poisson-`ν` driver channel cannot express. It is
called once per birth event, *after* the two daughters exist and *before* either is
scheduled, so a fitness change applies to the boosted daughter's own first division
rather than only to its descendants.

Because `NonMarkovCell` is immutable while `BinaryNode` is mutable, the supported way to
change a daughter is to replace its `data`:

```julia
c = d1.data
d1.data = NonMarkovCell(c.id, c.birthtime, c.mutations, 1.0 + s)
```

The new fitness propagates to all future descendants automatically, since `_make_daughter`
reads `parent_node.data.fitness`. Changing `id` is **not** supported — ids key
`population.cells`.

Note that `celldivision!` removes the parent and inserts both daughters before the hook
runs, so inside the callback `popsize(pop)` is the pre-division size **+ 1**. To inject at
`N_critic` cells, test `popsize(pop) == N_critic + 1`.

Keep all hook state in the closure (no package-level mutable state), and do not consume
the simulation's `rng` if seed reproducibility matters.

# Restart hook

On restart the package resets `cells`, `t`, `_next_id` and the pending event queue, but it
cannot reset state owned by your closure. A hook holding an `injected::Ref{Bool}` flag
would therefore silently fail to inject on the second attempt. Use `on_restart` to reset
that state:

```julia
injected  = Ref(false)
on_restart = pop -> (injected[] = false)
```

Alternatively set `restart_on_extinction = false` and drive your own retry loop with a
fresh closure per attempt, which is preferable when you also need to retry on other
criteria (e.g. loss of the driver clone to drift).

# Example
```julia
block = NonMarkovBlock(
    birth_dist  = f -> Gamma(2.0, 1.0 / f),
    death_dist  = f -> Gamma(2.0, 5.0),
    stopfunction = pop -> popsize(pop) >= 10_000,
    driver_dist = Exponential(0.05),
    fitness_update = (f, δ) -> f + δ,
    ν = 0.5,
)
```

Injecting a single driver into one daughter of the division that first takes the
population past `N_critic` cells:

```julia
injected = Ref(false)
block = NonMarkovBlock(
    birth_dist  = f -> Gamma(5.0, 1.0 / (5.0 * f)),
    death_dist  = f -> Gamma(5.0, 1.0 / (5.0 * 0.5)),
    stopfunction = pop -> popsize(pop) >= 10_000,
    driver_dist = Dirac(0.0),
    fitness_update = (f, δ) -> f,
    ν = 2.0,
    on_division = function (pop, parent, d1, d2)
        injected[] && return nothing
        popsize(pop) == N_critic + 1 || return nothing
        injected[] = true
        c = d1.data
        d1.data = NonMarkovCell(c.id, c.birthtime, c.mutations, 1.0 + s)
        return nothing
    end,
)
```
"""
@kwdef struct NonMarkovBlock{F1, F2, F3, F4, F5, F6, D}
    birth_dist::F1
    death_dist::F2
    stopfunction::F3
    driver_dist::D
    fitness_update::F4
    ν::Float64
    restart_on_extinction::Bool = false
    on_division::F5 = nothing
    on_restart::F6  = nothing
end
