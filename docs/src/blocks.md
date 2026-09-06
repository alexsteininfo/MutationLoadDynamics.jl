# The simulation block

[`NonMarkovBlock`](@ref) is the complete description of what cells do. Every field is a
function or a distribution you supply, so each is independently replaceable and cheap to
sweep — which matters, because fitting these is the point.

```julia
block = NonMarkovBlock(
    birth_dist     = f -> Gamma(5.0, 1 / (5 * f)),
    death_dist     = f -> Gamma(5.0, 1 / (5 * 0.3)),
    stopfunction   = pop -> popsize(pop) >= 10_000,
    driver_dist    = Exponential(0.05),
    fitness_update = (f, δ) -> f + δ,
    ν              = 0.2,
    restart_on_extinction = false,   # optional
    on_division    = nothing,        # optional
    on_restart     = nothing,        # optional
)
```

| Field | Signature | Role |
|:---|:---|:---|
| `birth_dist` | `f -> Distribution` | waiting time from birth to division, given this cell's fitness |
| `death_dist` | `f -> Distribution` | waiting time from birth to death |
| `stopfunction` | `pop -> Bool` | when the block is finished |
| `driver_dist` | `Distribution` | effect size ``δ`` of one driver mutation |
| `fitness_update` | `(f, δ) -> Float64` | how one driver changes fitness |
| `ν` | `Float64` | mean drivers per daughter per division (Poisson) |
| `restart_on_extinction` | `Bool` | retry from the initial state if the population dies |
| `on_division` | `(pop, parent, d1, d2) -> Nothing` | hook at every division |
| `on_restart` | `(pop) -> Nothing` | hook after an extinction restart |

The first two are covered below; `driver_dist`, `fitness_update` and `ν` have their own
page, [Mutations and selection](selection.md).

There are no defaults for the six required fields. That is deliberate: a default
waiting-time distribution would be a scientific claim smuggled in as a convenience.

## Waiting-time modes

`birth_dist` and `death_dist` are functions from a cell's fitness to *any*
`Distributions.jl` distribution over non-negative reals. Everything downstream —
the heap, the tree, the statistics — is indifferent to which you pick.

The useful way to parameterise them is to fix the **mean** at ``1/(b f)`` and let the
shape control only the variability. Then ``b`` is a rate you can compare across models
and ``k`` is a pure noise parameter.

| Mode | `birth_dist` | Mean | CV |
|:---|:---|:---|:---|
| Deterministic | `f -> Dirac(1 / (b * f))` | ``1/(bf)`` | ``0`` |
| Exponential (Markov) | `f -> Exponential(1 / (b * f))` | ``1/(bf)`` | ``1`` |
| Gamma | `f -> Gamma(k, 1 / (k * b * f))` | ``1/(bf)`` | ``1/\sqrt{k}`` |
| Weibull | `f -> Weibull(α, 1 / (b * f * Γ₁))`, `Γ₁ = mean(Weibull(α, 1.0))` | ``1/(bf)`` | ``\sqrt{\Gamma(1+2/α)/\Gamma(1+1/α)^2 - 1}`` |
| Log-normal | `f -> LogNormal(-log(b * f) - σ^2/2, σ)` | ``1/(bf)`` | ``\sqrt{e^{σ^2}-1}`` |

The Weibull row needs ``\Gamma(1 + 1/α)`` to normalise the mean. Rather than pull in
`SpecialFunctions`, get it from `Distributions` itself: `mean(Weibull(α, 1.0))` **is**
``\Gamma(1 + 1/α)``, since `mean(Weibull(α, θ)) = θ·Γ(1 + 1/α)`.

Every mean and CV in this table is asserted against `Distributions.mean` and
`Distributions.std` in `test/events.jl`, so these recipes cannot drift from what the
package actually recommends.

The Gamma row is the one to reach for first, and note the parameterisation:
`Gamma(k, 1/(k*b*f))`, **not** `Gamma(k, 1/(b*f))`. `Distributions.Gamma(α, θ)` has mean
``αθ``, so the second form has mean ``k/(bf)`` — changing ``k`` would silently rescale
time as well as changing the noise, and no comparison across ``k`` would mean anything.

### Deterministic — a fixed clock

```julia
birth_dist = f -> Dirac(1 / f)      # every cell divides exactly 1/f after birth
death_dist = f -> Dirac(Inf)        # and never dies
```

`Dirac` is the ``k \to \infty`` limit. From a synchronous start the population doubles in
lockstep: after ``m`` rounds there are ``2^m`` cells, all at divisional depth ``m``, and
`leaf_depths` returns a single repeated value. Exact ties in the heap are fine — tied
events simply fire consecutively, and `t_div <= t_die` breaks a birth/death tie in favour
of division.

This is the cleanest way to sanity-check a pipeline end to end, because every observable
has a closed form.

!!! tip "`Dirac(Inf)` is the idiom for \"no death\""
    `rand` on it returns `Inf`, which always loses the competing-risks comparison, and it
    consumes **no** random numbers — so a pure-birth run has the same stream as the same
    model written any other way. Prefer it to a large finite scale like
    `Exponential(1e8)`, which is only approximately deathless and does burn a draw.
    `Exponential(Inf)` and `Gamma(k, Inf)` also construct and return `Inf`, but they
    consume a draw, so `Dirac(Inf)` is the better spelling.

### Exponential — the Markov special case

```julia
b, d = 1.0, 0.3
birth_dist = f -> Exponential(1 / (b * f))
death_dist = f -> Exponential(1 / d)
```

This reproduces the ordinary Markovian birth–death process **exactly**, not
approximately. Pre-drawing each cell's event time and taking the global minimum is
equivalent to Gillespie precisely because the exponential is memoryless: an untouched
cell's residual waiting time has the same law as a fresh one, so never revising a
pending event costs nothing.

Use it when you want a baseline to compare non-exponential runs against, or when you
want the classical results to hold: net Malthusian rate ``r = b - d``, extinction
probability ``d/b`` from a single cell.

!!! warning "`Exponential(θ)` requires ``θ > 0``"
    `Exponential(0.0)` throws a `DomainError`. Writing `Exponential(1/d)` with `d = 0.0`
    does *not* throw — it gives `Exponential(Inf)`, whose `rand` returns `Inf` and so
    behaves as "no death" — but it still consumes a random draw and reads as an
    accident. For no death, write `Dirac(Inf)` and mean it.

### Gamma — the recommended default

```julia
k, b, d = 5.0, 1.0, 0.3
birth_dist = f -> Gamma(k, 1 / (k * b * f))
death_dist = f -> Gamma(k, 1 / (k * d))
```

``k`` is the only dimensionless timing parameter and CV ``= 1/\sqrt{k}``. See
[Why not the exponential?](concepts.md#Why-not-the-exponential?) for the empirical CV
ranges that motivate ``k = 5`` as a tumour-cell default, and for what changing ``k`` does
to the growth rate at fixed mean.

Death is often better left exponential even when division is not: apoptosis has no
comparable refractory structure, and a population-level exponential time-to-commitment is
the more defensible model (Fennell et al. 2005). Mixing is fine —
`birth_dist = f -> Gamma(5.0, 1/(5f))` with `death_dist = f -> Exponential(1/d)` is a
perfectly coherent block.

### Coupling fitness to death instead of division

Nothing forces fitness into `birth_dist`. A cell that survives longer is exactly as fit
as one that divides faster, and which one your model means is a biological claim.

```julia
# faster division (the usual choice)
birth_dist = f -> Gamma(k, 1 / (k * b * f));  death_dist = f -> Gamma(k, 1 / (k * d))

# longer survival, same division speed
birth_dist = f -> Gamma(k, 1 / (k * b));      death_dist = f -> Gamma(k, f / (k * d))

# both
birth_dist = f -> Gamma(k, 1 / (k * b * f));  death_dist = f -> Gamma(k, f / (k * d))
```

These are genuinely different processes, not reparameterisations of one another: the same
value of ``f`` buys a different growth rate under each, and even tuned to a common ``r``
they leave different trees behind. A birth advantage makes the fit lineage divide more
often, deepening it and compressing its coalescence times; a survival advantage makes it
lose fewer branches, widening it at the same depth. [`leaf_depths`](@ref) and
[`coalescence_times`](@ref) together distinguish the two.

!!! danger "Fitness must stay strictly positive"
    Fitness divides into a scale parameter, so a rule that lets ``f`` reach ``0`` or go
    negative fails at the next scheduling with a `DomainError` from the distribution
    constructor — *`Gamma: the condition θ > zero(θ) is not satisfied`* — partway
    through the run, not at setup. Deleterious effects need an explicit floor; see
    [Deleterious mutations and mutational load](selection.md#Deleterious-mutations-and-mutational-load).

## Growth modes: the stop condition

`stopfunction` receives the whole [`Population`](@ref) and returns a `Bool`. It is
evaluated **once per event**, and also once before the first event — so a block whose
criterion is already met does nothing at all, which is what makes chained blocks
composable.

```julia
stopfunction = pop -> popsize(pop) >= 10_000                    # fixed size
stopfunction = pop -> pop.t >= 20.0                             # fixed time
stopfunction = pop -> popsize(pop) >= 10_000 || pop.t >= 20.0   # whichever comes first
stopfunction = pop -> false                                     # run until extinction
```

**A size target is hit exactly.** Population size moves by ``\pm 1`` per event, so a
target above the current size is reached on the nose: ask for 10 000 and you get 10 000.

**A time target overshoots, by design.** The condition is tested after the event that
crossed it has already fired, so the run exits at ``t \gtrsim T`` — the state is the
first one at or after ``T``, not the state exactly at ``T``. The overshoot is one
inter-event interval and shrinks like ``1/N``. If you need the state exactly at ``T``,
record it with an [`AtTime`](@ref) snapshot trigger instead of stopping there.

**`pop -> false` runs until the heap empties**, which happens only when the last cell
dies. Under a supercritical process that is a long wait, and under pure birth
(`death_dist = f -> Dirac(Inf)`) it never happens at all — the loop will not terminate.
Pair it with a subcritical process, or add a ceiling:
`pop -> popsize(pop) > 10^6`.

!!! warning "Keep the stop condition ``O(1)``"
    It runs once per event, so anything that scans the population makes the whole
    simulation quadratic. `pop -> mean(fitness_per_cell(pop)) > 2.0` costs ``O(N)`` per
    event and therefore ``O(N^2)`` per run — at ``N = 10^5`` that is ``10^{10}`` field
    reads. If you need a criterion like that, track it incrementally in an
    `on_division` closure and have `stopfunction` read the closure variable.

### Growth regimes

With `p = P(T_div < T_die)` the probability that a cell divides rather than dies, the
lineage is a branching process with mean offspring number ``2p``:

| Regime | Condition | Extinction probability from one cell |
|:---|:---|:---|
| Supercritical | ``p > 1/2`` | ``(1-p)/p`` |
| Critical | ``p = 1/2`` | ``1`` |
| Subcritical | ``p < 1/2`` | ``1`` |

For exponentials this is the textbook result: ``p = b/(b+d)`` and the extinction
probability is ``d/b``. For any other pair of distributions ``p`` is a one-dimensional
integral you can evaluate numerically, or just estimate by running the block a few
hundred times.

This is the number to consult when a single-founder run keeps dying: at ``b = 1``,
``d = 0.5`` roughly half of all attempts go extinct, which is not a bug but the model.
`restart_on_extinction = true` handles it — with the caveat in the next section.

### Density-dependent and homeostatic growth

`birth_dist` receives only fitness, but a closure can capture the population and read its
size, which is enough to write logistic or homeostatic rules:

```julia
K   = 10_000
pop = initialize_population(fitness_init = 1.0)

block = NonMarkovBlock(
    birth_dist   = f -> Gamma(k, 1 / (k * b * f * max(1e-6, 1 - popsize(pop) / K))),
    death_dist   = f -> Gamma(k, 1 / (k * d)),
    stopfunction = p -> p.t >= 100.0,
    driver_dist = Dirac(0.0), fitness_update = (f, δ) -> f, ν = 0.0,
)
```

(The `max(1e-6, …)` floor is not optional: at ``N = K`` the bare factor is ``0`` and the
scale becomes `Inf`.)

!!! warning "Density dependence is frozen at each cell's birth"
    A cell's waiting time is drawn once, at birth, and the simulator has no mechanism to
    invalidate a pending event when the population changes. So the rule above evaluates
    the density **as of each cell's birth**, and a cell born at ``N = K/2`` keeps that
    slow-growth clock even after the population reaches ``K``.

    Over a regime where ``N`` changes little within one cell cycle this is a reasonable
    approximation. Approaching a homeostatic ceiling — exactly where you would want it —
    it is not, and the population will overshoot ``K``. This is true even for exponential
    waiting times: a real Gillespie implementation would resample on every rate change,
    and this one does not.

    The honest alternative is to change the *regime* rather than the rate, by chaining:
    grow to ``K`` under one block, then continue under a second block with ``b = d``.
    Turnover at constant size is then exact, because neither rate depends on ``N``.

```julia
grow = NonMarkovBlock(birth_dist = f -> Gamma(k, 1/(k*b*f)), death_dist = f -> Dirac(Inf),
                      stopfunction = p -> popsize(p) >= K, ...)
hold = NonMarkovBlock(birth_dist = f -> Gamma(k, 1/(k*b*f)), death_dist = f -> Gamma(k, 1/(k*b)),
                      stopfunction = p -> p.t >= t_end, ...)

simulate!(pop, grow, rng)
simulate!(pop, hold, rng)     # same tree, critical turnover, N fluctuates around K
```

## Extinction and restarting

`restart_on_extinction = true` makes `simulate!` snapshot the population as it is when
the call begins, and restore that snapshot whenever the population hits zero, retrying
until a lineage survives to meet the stop condition.

What a restart resets: `cells`, `t`, `_next_id`, the pending event queue, and the
[`MeasurementAccumulator`](@ref) if one was
passed — so the trajectory and snapshots describe the successful attempt only. What it
cannot reset is state owned by your own closures, which is what `on_restart` is for.

Three things to know:

- **The retry is unbounded.** A subcritical or critical block with
  `restart_on_extinction = true` never terminates, because every attempt is doomed.
  Check ``p > 1/2`` first.
- **Ids are reused across attempts.** `_next_id` is restored too, so cells from a failed
  attempt and cells from the successful one can share ids. The failed tree is dropped, so
  nothing observable collides — but do not cache node ids across a `simulate!` call that
  might restart.
- **It does not combine with chaining.** The snapshot is taken at the *start of this
  call*, so restarting a chained second block restores cells that were born under the
  first block and rewinds them to birthtimes that predate the boundary; their
  rescheduling then has exactly the age-conditioning defect that
  [Chaining blocks](#Chaining-blocks) exists to avoid. The snapshot is also a `deepcopy`
  of a live tree, so on a large chained population it copies the whole history. Use
  `restart_on_extinction` only on the first block, where the population is small and its
  birthtimes are the initial ones.

For anything more selective than "retry on extinction" — retry because the driver clone
was lost to drift, say — leave the flag off and write your own loop with a fresh closure
per attempt:

```julia
function attempt(seed)
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, block, MersenneTwister(seed))
    return popsize(pop) >= N_target ? pop : nothing
end

pop = nothing
for seed in 1:1000
    pop = attempt(seed)
    isnothing(pop) || break
end
```

## Hooks

### `on_division`

`on_division(pop, parent, d1, d2)` fires once per birth event. It exists for
deterministic, state-triggered interventions that the stochastic Poisson-``ν`` channel
cannot express: injecting one driver into one cell at one moment, tallying a statistic
incrementally, watching the trajectory of a specific lineage.

It runs **after** both daughters exist and **before** either is scheduled. That ordering
is load-bearing: `schedule_cell!` reads `node.data.fitness` at push time, so a hook that
ran after scheduling would leave a boosted daughter's own first division drawn at its
pre-boost fitness, and the change would only take effect one generation late.

`NonMarkovCell` is immutable and `BinaryNode` is mutable, so the way to change a daughter
is to replace its `data` wholesale:

```julia
c = d1.data
d1.data = NonMarkovCell(c.id, c.birthtime, c.mutations, 1.0 + s)
```

The new fitness propagates to every descendant automatically, because `_make_daughter`
reads the parent node's fitness at division time. **Do not change `id`** — ids key
`population.cells`.

!!! note "`popsize` inside the hook is the post-division size"
    `celldivision!` has already removed the parent and inserted both daughters, so a
    population that had `M` cells reads as `M + 1` inside the hook. To act on the
    division that takes the population past `N_critic`, test
    `popsize(pop) == N_critic + 1`.

A complete example — inject a single driver of effect `s` into one daughter of the first
division after the population reaches `N_critic`, in an otherwise neutral run:

```julia
injected = Ref(false)

block = NonMarkovBlock(
    birth_dist     = f -> Gamma(5.0, 1 / (5 * f)),
    death_dist     = f -> Gamma(5.0, 1 / (5 * 0.5)),
    stopfunction   = pop -> popsize(pop) >= 10_000,
    driver_dist    = Dirac(0.0),          # the Poisson channel is off …
    fitness_update = (f, δ) -> f,          # … in both respects
    ν              = 2.0,
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

Keep hook state in the closure — the package holds no global mutable state, so one
population, block and closure per simulation stays thread-safe. A hook that does not draw
from `rng` leaves the random stream untouched, which is what lets you add a watcher to an
existing run without changing its result.

### `on_restart`

`on_restart(pop)` fires after the population has been restored following an extinction.
It exists because a restart cannot reset state your closure owns:

```julia
injected   = Ref(false)
on_restart = pop -> (injected[] = false)
```

Without it, a hook holding a one-shot flag would spend its single injection on an attempt
that later went extinct, and every subsequent attempt would silently run neutral. It is
only ever called when `restart_on_extinction = true`.

## Chaining blocks

Calling [`simulate!`](@ref) again on the same population continues the **same lineage
tree** under new rules. This is how you express a change of regime: growth then
homeostasis, mutagen on then off, treatment then relapse.

```julia
pop = initialize_population(fitness_init = 1.0)
simulate!(pop, growth_block, rng)      # ν = 1.0, grow to 1 000
simulate!(pop, neutral_block, rng)     # ν = 0.0, grow to 5 000 — same tree
```

The queue of already-drawn, not-yet-fired events is carried on the population and
**reused**, not redrawn.

That is the whole subtlety, and it matters for any non-exponential waiting time. Cells
alive at the boundary have already survived part of a cell cycle without an event.
Redrawing would anchor each of them at its `birthtime` and ignore that accumulated age;
the correct law is the conditional ``T \mid T > \text{age}``, and drawing
unconditionally both loses the conditioning and places some events *before* the current
`pop.t`, running the clock backwards. Carrying the queue sidesteps the problem exactly:
each cell simply keeps the event it had already committed to. A chained run is then
identical draw-for-draw to the equivalent uninterrupted one, which the test suite asserts
bit for bit.

The consequence to plan around:

!!! note "A carried event was drawn under the previous block's distributions"
    A cell that was already scheduled at the boundary has committed to that division or
    death time. A second block that changes `birth_dist` or `death_dist` therefore only
    affects cells scheduled *after* the boundary — the change phases in over roughly one
    cell cycle rather than taking effect instantly.

    This is usually what you want. `driver_dist`, `fitness_update` and `ν` are read at
    division time, not at scheduling time, so those take effect immediately.

[`reset_schedule!`](@ref) discards the queue so the next `simulate!` redraws everything:

```julia
reset_schedule!(pop)
simulate!(pop, second_block, rng)
```

But this reintroduces the age-conditioning error above, so it is only sound when every
living cell has just been born — a freshly initialised population, essentially. Its real
purpose is recovery: if you add or remove entries in `pop.cells` by hand, the queue no
longer matches the population, and `simulate!` throws naming both counts rather than
running with a corrupt heap. `reset_schedule!` is the documented way back.
