# Mutations and selection

Three fields of [`NonMarkovBlock`](@ref) describe how cells acquire drivers and what a
driver does. They are separate on purpose: `ν` is *how many*, `driver_dist` is *how big*,
and `fitness_update` is *how it combines with what the cell already has*. Selection
models differ almost entirely in the third.

## The driver channel

At every division, **each daughter independently**:

1. draws a number of new drivers, ``j \sim \mathrm{Poisson}(ν)``;
2. for each of them draws an effect size ``δ \sim`` `driver_dist` and applies
   ``f \leftarrow`` `fitness_update(f, δ)`, starting from the **parent's** fitness;
3. stores ``j`` in its own `mutations` field and the final ``f`` in `fitness`.

```math
j \sim \mathrm{Poisson}(ν), \qquad
f_\text{daughter} = \bigl(\texttt{fitness\_update}(\cdot,\, δ_j) \circ \cdots \circ
\texttt{fitness\_update}(\cdot,\, δ_1)\bigr)(f_\text{parent})
```

Two properties follow, and both matter:

**`fitness_update` is applied once per mutation, not once to their sum.** Two drivers of
effect ``δ`` are `update(update(f, δ), δ)`, which for anything non-linear — a cap, a
maximum, a saturating factor — is not `update(f, 2δ)`. This is what lets the third
argument express epistasis at all.

**The two daughters are drawn independently.** A division does not produce one mutant and
one wild-type daughter; it produces two draws from the same law. Both can mutate, neither
can, or one can. The mean number of new drivers per division is therefore ``2ν``, not
``ν``.

`ν` is the mean drivers **per daughter per division**. The package targets the regime
``ν \gtrsim 0.1``, where every cell has a distinct fitness history and a subclone label
would be meaningless. For ``ν \ll 1`` with clone-level dynamics, use
[`BirthDeathMutation`](https://github.com/alexsteininfo/BirthDeathMutation) instead.

!!! note "Drivers only, and only at division"
    `mutations` counts driver events — mutations that pass through `fitness_update`.
    Neutral passengers are not simulated, because they need no state: their expected
    spectrum is ``m`` times the tree's [`branch_spectrum`](@ref) and their per-cell
    burden is `Poisson(m * depth)` drawn against [`leaf_depths`](@ref). See
    [Neutral evolution](#Neutral-evolution) below.

    Mutations also occur only at division. A cell's fitness never changes between its
    birth and its own division, which is precisely what makes the competing-risks
    scheduling exact — see [Competing risks](concepts.md#Competing-risks).

## Selection modes

Every mode below is a two-line change: pick `driver_dist`, pick `fitness_update`. `s` is
the selection coefficient throughout, and `f₀ = 1.0` the founding fitness.

| Mode | `driver_dist` | `fitness_update` |
|:---|:---|:---|
| [Neutral](#Neutral-evolution) | anything (`ν = 0`) | anything |
| [Additive, fixed](#Additive,-fixed-effect) | `Dirac(s)` | `(f, δ) -> f + δ` |
| [Additive, random](#Additive,-random-effect) | `Exponential(s)` | `(f, δ) -> f + δ` |
| [Multiplicative](#Multiplicative) | `Dirac(s)` or `Exponential(s)` | `(f, δ) -> f * (1 + δ)` |
| [Winner-takes-all](#Winner-takes-all-(max-random)) | `Exponential(s)` | `(f, δ) -> max(f, 1 + δ)` |
| [Capped / diminishing returns](#Capped-selection-(diminishing-returns)) | `Dirac(1.0)` or `Exponential(1.0)` | `(f, δ) -> min(f * (1 + s*δ*(1 - f/M)), M)` |
| [Deleterious](#Deleterious-mutations-and-mutational-load) | `Normal(-μ, σ)` | `(f, δ) -> max(f + δ, f_min)` |

### Neutral evolution

```julia
block = NonMarkovBlock(
    birth_dist     = f -> Gamma(k, 1 / (k * b * f)),
    death_dist     = f -> Gamma(k, 1 / (k * d)),
    stopfunction   = pop -> popsize(pop) >= N,
    driver_dist    = Dirac(0.0),      # unused
    fitness_update = (f, δ) -> f,      # unused
    ν              = 0.0,
)
```

With `ν = 0` no effect sizes are ever drawn, every cell keeps `f₀`, and all cells share
one waiting-time distribution. `driver_dist` and `fitness_update` are still required
fields — the block has no defaults — but neither is consulted.

This is the null model, and it is also where the package earns its keep: the tree is
still non-trivial, because non-exponential timing changes its shape even without
selection. The observables to compare against are [`leaf_depths`](@ref) and
[`branch_spectrum`](@ref).

!!! tip "Neutral passengers come afterwards, not during"
    For a neutral per-division mutation rate ``m``, do not simulate the mutations —
    derive them:

    ```julia
    root   = getsingleroot(allcells(pop))
    burden = [rand(rng, Poisson(m * d)) for d in leaf_depths(root)]  # per-cell burden
    sfs    = m .* branch_spectrum(root)                              # expected SFS, k ≥ 2
    ```

    This is exact in distribution, costs nothing, and lets one simulated tree serve every
    value of ``m``. Note the ``k \ge 2`` caveat on the spectrum identity —
    see [`branch_spectrum`](@ref).

!!! warning "`ν = 0` still consumes a random draw per daughter"
    `rand(rng, Poisson(0.0))` returns `0` but does advance the stream, so a neutral run
    does not share a random sequence with a run that has no mutation step at all. This
    only matters when comparing seeds across models.

### Additive, fixed effect

```julia
driver_dist    = Dirac(s)
fitness_update = (f, δ) -> f + δ
```

Every driver adds the same increment, so after ``n`` drivers ``f = 1 + ns``. There is no
ceiling; fitness grows without bound as drivers accumulate.

The natural reading is a set of independent gain-of-function events in one pathway, each
shortening the mean division time by the same absolute amount. Because ``n`` is a sum of
Poisson draws along a lineage, the endpoint fitness distribution is close to normal, and
at ``s = 0`` the model degenerates exactly to neutral.

Use it as the reference selection model: it is the simplest thing that is not neutral,
and the one whose analytics are most tractable.

### Additive, random effect

```julia
driver_dist    = Exponential(s)      # mean s, heavy right tail
fitness_update = (f, δ) -> f + δ
```

Same accumulation rule, but drivers now vary in impact. Mean fitness tracks the
fixed-effect model at the same ``s``; the variance is larger, and rare large-effect
drivers produce lineages that sweep much faster than the mean would suggest.

The exponential is the maximum-entropy choice for a positive quantity with a known mean,
and it matches measured distributions of beneficial effect sizes in microbial evolution
experiments. Comparing this mode against the fixed-effect one at matched ``s`` isolates
the consequences of effect-size heterogeneity alone.

!!! warning "`Exponential(0.0)` throws"
    Sweeping ``s`` down to zero needs a special case, because `Exponential` requires a
    strictly positive scale:

    ```julia
    driver_dist = s > 0 ? Exponential(s) : Dirac(0.0)
    ```

    Note that this changes the random stream as well as the model — `Dirac` consumes no
    draws and `Exponential` consumes one.

### Multiplicative

```julia
driver_dist    = Dirac(s)            # or Exponential(s) for random effects
fitness_update = (f, δ) -> f * (1 + δ)
```

After ``n`` drivers ``f = (1+s)^n``, so fitness compounds and log-fitness is additive.
This is the right rule when drivers act on independent multiplicative components of the
division rate rather than on a shared additive budget.

Growth is faster than additive at the same ``s``, and much faster in the tail; combined
with a fitness-dependent birth rate it can produce division times short enough that the
population becomes numerically dominated by one lineage well before the size target. If
that is not what you want, cap it — see the next mode.

### Winner-takes-all (max-random)

```julia
driver_dist    = s > 0 ? Exponential(s) : Dirac(0.0)
fitness_update = (f, δ) -> max(f, 1 + δ)
```

Each driver proposes a new fitness ``1 + δ`` and the cell keeps the better of that and
what it already had. Fitness is non-decreasing and, once high, permanent.

This is the sharpest available statement of epistasis in the package: a second driver in
an already-activated pathway does nothing unless it exceeds the current level. Only the
running maximum of the lineage's draws matters, so after ``n`` drivers the fitness is the
order statistic ``\max\{δ_1, \ldots, δ_n\}``, whose expectation for exponential draws
grows like ``s \ln n`` — logarithmically, not linearly. The result is strong diminishing
returns and a characteristic plateau, with a sweep signature in the SFS more pronounced
than the additive model shows at the same mean effect.

Note that this rule reads `f` but discards it whenever the proposal wins, so unlike every
other mode it is not a function of the accumulated history beyond its maximum.

### Capped selection (diminishing returns)

```julia
M              = 10.0                 # ceiling on fitness
driver_dist    = Dirac(1.0)           # fixed magnitude …
driver_dist    = Exponential(1.0)     # … or random, same mean
fitness_update = (f, δ) -> min(f * (1 + s * δ * (1 - f / M)), M)
```

A multiplicative rule with a logistic brake: the factor ``(1 - f/M)`` shrinks the gain as
fitness approaches the ceiling ``M``, and the outer `min` guarantees it. This is the
discrete analogue of logistic growth in fitness space, with a stable attractor at
``f = M``.

The ceiling stands for a maximum achievable division rate — nutrient supply, checkpoint
constraints, physical crowding — rather than for anything about the mutations
themselves. Early drivers have close to their full multiplicative effect; near ``M`` the
gain per driver goes to zero and further drivers are nearly neutral. Mutational load
therefore matters non-linearly, which is the qualitative difference from the additive
modes.

Writing the magnitude as a separate ``δ`` with mean 1 is what makes the fixed and random
variants a one-line switch at matched mean effect. The random variant is overdispersed
at the endpoint, and for large ``s`` a single right-tail draw can carry a cell to ``M`` in
one step, producing a bimodal final fitness distribution — cells near ``M`` alongside
cells still near 1. At ``s = 0`` both variants are neutral regardless of ``δ``.

### Deleterious mutations and mutational load

Nothing requires ``δ > 0``. The package is named for the regime where most mutations are
mildly harmful and load accumulates against selection:

```julia
f_min          = 0.05                       # hard floor — see below
driver_dist    = Normal(-0.02, 0.05)        # mostly deleterious, occasionally beneficial
fitness_update = (f, δ) -> max(f + δ, f_min)
```

!!! danger "A floor is mandatory, not stylistic"
    Fitness divides into a scale parameter. The moment a cell's fitness reaches ``0`` or
    goes negative, its next scheduling raises

    ```
    DomainError: Gamma: the condition θ > zero(θ) is not satisfied
    ```

    partway through the run rather than at setup. Any `fitness_update` whose effect
    distribution has support below zero **must** clamp — `max(f + δ, f_min)` — with
    `f_min` strictly positive. A multiplicative rule `f * (1 + δ)` needs the same care
    whenever ``δ`` can reach ``-1``.

A floored additive rule makes ``f = f_\text{min}`` absorbing in effect: a cell there
divides slowly and is likely to be out-competed before it recovers. If you want load to
be lethal rather than merely crippling, express that through death instead — give
`death_dist` a mean that shortens as ``f`` falls — which is both more biological and
free of the floor problem.

## Choosing ν, s and k together

The three interact, and it is easy to pick a corner of parameter space where nothing
interesting happens.

- **``ν`` sets how many drivers a cell has by the time the run ends**, which is roughly
  ``2ν`` times its divisional depth. Growing to ``N = 10^4`` from one cell means a depth
  of order ``\log_2 N \approx 13``, so ``ν = 0.2`` gives about 5 drivers per cell. Below
  ``ν \approx 0.01`` most cells carry none and the run is neutral in all but name.
- **``s`` and ``ν`` trade off.** What drives the fitness distribution is the product ``νs``
  (per-division expected gain), so a sweep over ``s`` at fixed ``ν`` and a sweep over
  ``ν`` at fixed ``s`` explore much the same axis — until a non-linear
  `fitness_update` breaks the symmetry, which is exactly what the capped and
  winner-takes-all modes are for.
- **``k`` changes the growth rate at fixed mean division time** (see
  [What non-exponential timing changes](concepts.md#What-non-exponential-timing-changes)),
  so a selection coefficient calibrated at ``k = 1`` does not mean the same thing at
  ``k = 5``. Fix ``k`` first, then calibrate ``s``.

A useful diagnostic is the variance of the fitness distribution against its mean: under
additive selection both grow linearly in ``νs``; under a capped rule the variance turns
over and falls as the population approaches ``M``; under winner-takes-all the mean
saturates logarithmically while the variance collapses.

```julia
fitness_per_cell(pop)     # the full distribution
mean_k(pop), var_k(pop)   # driver counts, the load itself
```

## Injecting a single driver at a chosen moment

The Poisson channel is stochastic and population-wide. To place **one** driver in **one**
cell at **one** moment — the standard setup for asking how a single advantageous clone
fares against a growing background — turn the channel off and use the `on_division` hook
instead. See [Hooks](blocks.md#on_division).
