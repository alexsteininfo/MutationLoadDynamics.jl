"""
    NonMarkovBlock{F1,F2,F3,F4,D}

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
"""
@kwdef struct NonMarkovBlock{F1, F2, F3, F4, D}
    birth_dist::F1
    death_dist::F2
    stopfunction::F3
    driver_dist::D
    fitness_update::F4
    ν::Float64
    restart_on_extinction::Bool = false
end
