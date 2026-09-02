# PLAN — division hook + correct block chaining

**Status:** specification / handoff. Written 2026-09-01 from the consumer side
(`gITH-nonMarkovian`). Not yet planned in detail, not implemented, not tested.

**How to use this file:** open a Claude Code session in this repo
(`/Users/alexanderstein/Documents/GitHub/MutationLoadDynamics.jl`), tell it to read
`PLAN-block-hooks.md`, then plan → implement → test. The "Open design decisions"
section lists what still needs a call; everything else is settled requirements.

---

## 1. Why

The consumer repo `gITH-nonMarkovian` (same author) is building a new family of
selection scenarios under `analysis/sim_runs/selection_1/`. **Scenario 1** is:

> Take the neutral simulation settings exactly as they are, and add **exactly one**
> fitness-enhancing mutation per simulation. It is placed in one daughter of the
> division that takes the population from `N_critic` to `N_critic + 1` cells (the
> first such birth event). Selection strength `s ∈ {0.1, 0.2, …, 2.0}` (20 values,
> additive: `f ← 1 + s`). `N_critic` is log-spaced from 1 to `N_target/2`
> (10 values, so injection times are roughly equally spaced, since `t ≈ ln(N)/r`
> during exponential growth). Full 20 × 10 factorial grid × 5 replicates = 1000
> simulations per parameter set.

Parameter sets mirror the neutral runs: deterministic (`Dirac` timing, `N ∈ {1024,
16384}`, `d = 0`), gamma `k = 5` and markov (`N ∈ {1000, 10000}`, `d ∈ {0, 0.5, 0.9}`)
— 14 sets in total. Neutral mutations keep being generated the existing way, with
`ν = 2.0` and `driver_dist = Dirac(0.0)`; the injected driver is a *separate* event
that must not disturb that stream.

**This package currently offers no way to do that**, and the one API that looks like
it should work (chaining two `simulate!` calls) is silently broken for non-exponential
waiting times. Two work items follow. **Item A is required for scenario 1. Item B is a
latent-bug fix that is independently worth doing** and unblocks later scenarios — do
not let B block A.

---

## 2. Item A — `on_division` hook on `NonMarkovBlock` (required)

### A.1 The gap

`NonMarkovBlock` (`src/blocks.jl:35-43`) exposes only a Poisson-`ν` driver channel:
every daughter draws `j ~ Poisson(ν)` increments from `driver_dist` and folds each
through `fitness_update` (`src/cellupdates.jl:31-45`). There is no way to say "apply
this fitness change to this one cell, at this one moment". Scenario 1 needs a
deterministic, single, state-triggered fitness event.

### A.2 Proposed API

Add one optional field:

```julia
@kwdef struct NonMarkovBlock{F1, F2, F3, F4, F5, D}
    birth_dist::F1
    death_dist::F2
    stopfunction::F3
    driver_dist::D
    fitness_update::F4
    ν::Float64
    restart_on_extinction::Bool = false
    on_division::F5 = nothing        # (pop, parent, d1, d2) -> Nothing
end
```

and call it from `simulate!` (`src/simulations.jl:35-38`), **between `celldivision!`
and the two `schedule_cell!` calls**:

```julia
if event.event_type == :birth
    d1, d2 = celldivision!(population, event.node, event.time, block, rng)
    isnothing(block.on_division) || block.on_division(population, event.node, d1, d2)
    schedule_cell!(heap, d1, block, rng)
    schedule_cell!(heap, d2, block, rng)
else
    ...
```

### A.3 Why that exact position

This is the whole point of the hook and the thing most likely to be got wrong.

`schedule_cell!` (`src/events.jl:26-41`) reads `f = node.data.fitness` and draws the
cell's next event time from `block.birth_dist(f)`. If the hook fires **after**
scheduling, the boosted daughter's *own first division* is still drawn at `f = 1` and
only its descendants divide faster — a systematic one-cell-cycle delay in every sweep.
Firing before scheduling makes the injection exact.

Note also that `celldivision!` (`src/cellupdates.jl:24-26`) deletes the parent and
inserts both daughters before returning, so inside the callback
`popsize(pop) == (pre-division size) + 1`. A consumer injecting at `N_critic` therefore
tests `popsize(pop) == N_critic + 1`.

### A.4 What the callback is allowed to do

`BinaryNode` is mutable but `NonMarkovCell` is immutable, so the documented way to
change a daughter is to replace its `data`:

```julia
c = d1.data
d1.data = NonMarkovCell(c.id, c.birthtime, c.mutations, 1.0 + s)
```

The fitness then propagates to descendants automatically — `_make_daughter` reads
`parent_node.data.fitness` (`src/cellupdates.jl:15,33`). Whether the consumer also
increments `c.mutations` (making the driver visible in the SFS alongside the neutral
mutations) is a consumer-side choice; the hook must permit either. Document that
replacing `.data` is supported and that changing `id` is not (ids key
`population.cells`).

### A.5 Hard requirements

1. **`on_division === nothing` must be bit-for-bit behaviour-preserving.** The consumer
   holds ~1.8 GB of serialized neutral results that must stay reproducible from their
   seeds. No new `rand` call may be introduced on the default path, and the branch must
   not perturb the rng stream. A golden-value regression test is required (see A.6.1).
2. **No global or package-level mutable state.** The consumer runs simulations under
   `Threads.@threads`, constructing one `Population`, one block, and one closure per
   simulation. All hook state lives in the consumer's closure.
3. **Backward-compatible construction.** Existing code builds `NonMarkovBlock` without
   `on_division`; the `@kwdef` default covers this. But adding the `F5` type parameter
   changes the struct's type signature — grep for explicitly parameterized
   `NonMarkovBlock{...}` annotations before committing. Verified 2026-09-01: the only
   two occurrences are the docstring header (`src/blocks.jl:2`, which spells out the
   type parameters and must be updated) and the definition itself — every internal
   signature uses the bare `block::NonMarkovBlock`, so this is clean.
4. **`restart_on_extinction` interaction must be documented.** On restart
   (`src/simulations.jl:50-55`) the package resets `population.cells`, `t`, `_next_id`
   and the accumulator, but it cannot reset a user closure's state — so a
   consumer whose closure has an `injected::Ref{Bool}` flag would silently fail to
   inject on the second attempt. Either document loudly that consumers must use
   `restart_on_extinction = false` with an explicit retry loop, or add an optional
   `on_restart(pop)` callback fired at the same point. See open decision D3.

### A.6 Tests to add (`test/simulations.jl`)

1. **Regression / no-op:** fixed seed, `on_division = nothing`, grow to a target N.
   Assert exact equality of `popsize`, the sorted `fitness_per_cell` vector, and
   `mutations_per_cell` against values recorded from the pre-change code. This is the
   test that protects the existing neutral data.
2. **Hook contract:** a callback that appends `(popsize(pop), parent.data.id,
   d1.data.id, d2.data.id)` to a vector. Assert it fires once per birth event,
   that `d1.parent === parent && d2.parent === parent`, and that the recorded popsize
   increases by exactly 1 per division in a `d = 0` run.
3. **Timing is applied immediately (the A.3 test):** `birth_dist = f -> Dirac(1.0/f)`,
   negligible death, callback sets `d1` fitness to `2.0` at the very first division.
   Assert the boosted daughter's lifetime is `0.5`, not `1.0`. A hook placed after
   `schedule_cell!` fails this and nothing else catches it.
4. **Single clade:** inject once at `popsize == N_critic + 1` in a `ν = 0` run; at the
   end, the set of cells with `fitness > 1` is exactly the alive leaves under the
   injected node (use `AbstractTrees.Leaves` from that node), and every other cell has
   `fitness == 1.0`.
5. **Monotonicity (statistical, several seeds):** driver clone size at the target N
   increases with `s`. Keep tolerances loose.
6. **Constructor:** a `NonMarkovBlock` built without `on_division` still runs.

---

## 3. Item B — chaining two `simulate!` calls is broken (independent fix)

### B.1 The bug

`simulate!`'s docstring says it returns the population "so blocks can be chained", but
chaining is unsound for any non-exponential waiting-time distribution — which is the
entire point of this package. No script in either repo has ever exercised it, so it has
gone unnoticed.

`simulate!` rebuilds the event queue from scratch on entry (`src/simulations.jl:24-27`)
and `schedule_cell!` (`src/events.jl:33-35`) anchors the redraw at the cell's
**birthtime**:

```julia
t0 = node.data.birthtime
t_div = t0 + rand(rng, block.birth_dist(f))
t_die = t0 + rand(rng, block.death_dist(f))
```

For a fresh `initialize_population` every cell has `birthtime == population.t`, so this
is correct. On a **chained** call, alive cells have `birthtime < population.t` and have
already survived `age = population.t - birthtime` without an event. The correct redraw
is from the conditional law `T | T > age`; the code draws unconditionally. Two
consequences:

1. **Time runs backwards.** A fresh draw shorter than `age` yields an event time before
   the current `population.t`, and line 33 assigns `population.t = event.time`
   unconditionally. Daughters are then born "before" the present; `celllifetime` can go
   negative, and `record_trajectory_if_due!`'s `while t >= acc.next_trajectory_t` loop
   silently stops advancing.
2. **The age-conditioning is dropped**, so the process across the boundary is not the
   intended non-Markovian process — a spurious partial re-synchronization of the whole
   population's cell cycles.

The severity scales with the number of alive cells at the boundary, so it is invisible
in a 1-cell smoke test and severe at N = 5000.

### B.2 Candidate fixes

**Fix 1 — carry the pending event queue across blocks (recommended).** Keep the heap on
`Population` (or thread it through `simulate!`) so a chained call reuses the event times
already drawn rather than redrawing. Exact, assumption-free, and cheap. The semantic
question it raises: if block 2 changes `birth_dist`/`death_dist`, carried-over events
were drawn under block 1's law. That is defensible ("the cell already committed to this
division time") but must be *documented*, and `restart_on_extinction` must clear the
carried queue.

**Fix 2 — conditional resampling.** Redraw from `T | T > age` by inverse-CDF:
`T = quantile(D, u₀ + rand(rng)*(1 - u₀))` with `u₀ = cdf(D, age)`. Exact in law and it
correctly applies the *new* block's distributions, so it composes with Fix 1 as the
opt-in path when block 2 genuinely changes the timing law. Watch out for `Dirac`
(degenerate `cdf`/`quantile`, and it is used by the deterministic neutral model) and for
`u₀ → 1` on old cells, which is numerically nasty.

**Fix 3 — anchor at `population.t` instead of `birthtime` — REJECTED.** It removes the
negative-time symptom but resets every cell's age to zero, which is exactly the
artificial synchronization that makes the model wrong. Recording it here so it does not
get rediscovered as "the simple fix".

### B.3 Tests to add

7. **Monotone time:** chain two blocks (stop at N=100, continue to N=200) with a
   `Gamma(5, ...)` birth distribution and assert `population.t` never decreases —
   easiest via an accumulator trajectory, or a callback that records event times.
8. **Tree sanity:** for every non-root node, `birthtime >= parent.data.birthtime`, and
   `celllifetime >= 0` for all completed cells.
9. **Statistical equivalence (the real test):** growing 1 → 200 in a single block versus
   1 → 100 → 200 in two identical chained blocks must give the same distribution of
   (a) total time to reach 200 and (b) mean leaf depth, over ~200 seeds each, within a
   loose tolerance. This currently fails; it should pass after the fix.

---

## 4. Open design decisions

- **D1 — callback signature.** `on_division(pop, parent, d1, d2)` as specified, or also
  pass `rng` (`on_division(pop, parent, d1, d2, rng)`)? Passing `rng` lets a consumer
  randomize which daughter is hit or draw a random effect size, at the cost of a wider
  contract; a callback that ignores it leaves the stream untouched either way.
  Scenario 1 does not need it.
- **D2 — symmetric `on_death` hook?** Not needed by scenario 1. Deliberately out of
  scope unless it falls out for free.
- **D3 — `restart_on_extinction` and closure state.** Document-only ("use
  `restart_on_extinction = false` plus an explicit retry loop"), or add an optional
  `on_restart(pop)` callback? The consumer needs an explicit retry loop regardless,
  because scenario 1 also retries when the driver clone is lost to drift (see §5), so
  document-only is probably sufficient.
- **D4 — scope of Item B.** Fix 1 alone, or Fix 1 + Fix 2? Fix 1 alone is enough to make
  chaining sound for identical-law blocks, which is all anything needs today.
- **D5 — version bump.** `Project.toml` is at `0.1.0`. A struct field addition warrants
  `0.2.0`. The consumer `Pkg.develop`s this package by absolute path
  (`/Users/alexanderstein/Documents/GitHub/MutationLoadDynamics.jl`), so changes are
  picked up with no manifest edit, but bump it anyway for honesty.

---

## 5. Consumer-side context (for reference — do not implement here)

Once the hook exists, `gITH-nonMarkovian` implements scenario 1 entirely in its own
scripts, roughly:

```julia
injected = Ref(false)
on_division = function (pop, parent, d1, d2)
    injected[] && return nothing
    popsize(pop) == N_critic + 1 || return nothing
    injected[] = true
    c = d1.data
    d1.data = NonMarkovCell(c.id, c.birthtime, c.mutations, 1.0 + s)
    return nothing
end
```

with `ν = 2.0`, `driver_dist = Dirac(0.0)`, `fitness_update = (f, δ) -> f` unchanged
from the neutral runs, `restart_on_extinction = false`, and an outer retry loop that
re-runs with a fresh seed and a fresh closure until (a) the population reaches
`N_target` and (b) the driver clone is still alive at `N_target` — recording the attempt
count, since `1/E[attempts]` is the driver establishment probability and its
Gamma-vs-exponential contrast is itself a result. (Single-driver loss probability is
≈ `d/(b(1+s))`: 45% at `d = 0.5, s = 0.1`, 82% at `d = 0.9, s = 0.1`, so this matters.)

Two approaches were considered and rejected on the consumer side; recording them so the
hook is not re-litigated:

- **Stateful `fitness_update` closure** capturing `pop` and firing when
  `popsize(pop) == N_critic`. Needs no package change and would work — but the boost
  rides on the neutral-mutation Poisson stream, so when both daughters draw `j = 0`
  (probability `e^{-4} ≈ 1.8%`) the injection slips to a later division. Also entangles
  two conceptually separate mutation channels.
- **Two chained blocks** (grow to `N_critic`, boost a cell, grow to `N_target`). This is
  Item B's bug, and at `N_critic = 5000` it is fatal.

---

## 6. Suggested order of work

1. Item A: add the field, wire the call site, write tests A.6.1–A.6.6, update the
   `NonMarkovBlock` docstring (`src/blocks.jl:1-34`) and its example. Green tests here
   unblock the consumer immediately.
2. Item B: pick a fix from B.2, implement, write tests B.3.7–B.3.9, and correct or
   qualify the "so blocks can be chained" claim in the `simulate!` docstring
   (`src/simulations.jl:1-9`).
3. Bump `version` in `Project.toml`; run the full suite
   (`julia --project=. -e 'using Pkg; Pkg.test()'`).
