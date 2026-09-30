# Changelog

## 0.4.0 (unreleased)

This release is a deliberate clean break. Each item says what downstream code
(for example `gITH-nonMarkovian`) has to change.

### Breaking: data layout

- **`NonMarkovCell` has a fifth field.** The layout is now
  `NonMarkovCell(id, birthtime, mutations, total_mutations, fitness)`, where
  `total_mutations` is the cell's whole driver burden (parent's total plus its own
  `mutations`). The 4-argument constructor no longer exists.
  - Hooks that did `d1.data = NonMarkovCell(c.id, c.birthtime, c.mutations, f)` must use
    **`set_fitness!(d1, f)`**.
  - Hand-built trees must supply a consistent `total_mutations`.
  - Trees and populations serialised with 0.3 cannot be deserialised by 0.4. Re-run, or
    convert them with 0.3 loaded.

### Breaking: results that change for the same seed

- **Leaf sampling uses a new, version-stable recipe.** `sample_leaves` now draws with
  `StableRNG(seed)` and a partial Fisher–Yates shuffle, instead of
  `MersenneTwister(seed)` and `randperm`. The same seed gives a different draw than 0.3,
  but the same draw on every Julia version from now on.
  - `sample_trees` derives per-draw seeds with an explicit splitmix64 mix instead of
    `Base.hash`, so a base seed also gives different per-draw seeds than 0.3.
  - Samples stored by 0.3 can only be replayed with 0.3 on the Julia line that made them.
- **`mutations_per_cell(root; includeclonal)` changed meaning** (same name, same default).
  - `false` (default) now counts only mutations acquired strictly *below* `root`. Before,
    it included `root`'s own mutations.
  - `true` now gives each leaf's full burden back to the top of its tree.
  - For the founder of a simulated tree nothing changes, because it carries 0 mutations.
    Results change only for subtree roots, and for hand-built roots with mutations.
  - `mutations_per_cell(s.root; includeclonal = true)` is the explicit way to get a
    sampled cell's full-tree burden.
- **`AtTime(t)` is exact.** It records the state exactly at `t`, labelled `t`. Before, it
  recorded the state after the first event at or after `t`, labelled with that event's time.
- **Trajectory points are exact.** Each point is the state at its grid time. Before, it
  was the state after the first event past it. Stored trajectories differ slightly,
  mostly early in a run.
- **`AtEnd` fires at the end of every `simulate!` call.** An accumulator carried across
  a chain now gets one end snapshot per block, instead of one in total.
- **`AtPopSize` is also checked before the first event.** A threshold already met at the
  start fires at once, with the starting state.
- **An extinction restart discards only the current call's records.** Before, it wiped
  the whole accumulator, including earlier chained blocks.
- **Fresh schedules visit cells in id order and are conditioned on age.** This applies to
  a new population, `reset_schedule!` and restarts.
  - Each cell's next event is drawn conditioned on nothing having happened to it before
    `pop.t`.
  - Runs from `initialize_population(N)` with `N > 1`, runs after `reset_schedule!` on
    an aged population, and restarts of chained blocks now consume the stream
    differently. They are now statistically correct.
  - Single-founder runs, with or without chaining and hooks, are **bit-identical to 0.3**
    (verified on Julia 1.10).
- **The default rng of `simulate!`** is `Random.default_rng()` instead of the legacy
  `Random.GLOBAL_RNG`.

### Breaking: validation and requirements

- `MeasurementSpec` rejects `trajectory_dt <= 0` or `NaN`; before, `0` hung `simulate!`.
- `AtPopSize(N)` requires `N >= 1`, and `AtTime(NaN)` is rejected.
- `initialize_population(N)` requires `N >= 1`. `fitness_init` and `time` accept any
  `Real`.
- Julia 1.10 or newer is required (was 1.9).
- New dependency: StableRNGs.jl.

### Fixed

- `filtered_mutations_per_cell` on a subtree no longer throws `KeyError`.
- `findMRCA`, `pairwisedistance` and `coalescence_times` no longer crash on forests
  whose root ids are not the smallest. The coalescence climb is iterative, so it cannot
  overflow the stack.
- `simulate!` raises an error instead of silently running time backwards if an event
  precedes the current time. It also raises, instead of looping forever, if a cell is
  too old to be conditioned under the block's waiting times.
- `restart_on_extinction` now works on chained blocks, because restored cells are
  conditioned on their age. The documented limitation is gone.

### Added

- `set_fitness!(node, f)`.
- `getsingleroot(population)`, much faster than `getsingleroot(allcells(population))`.
- Test suites: `test/measurements.jl`; `test/validation.jl`, which checks extinction
  probabilities, the Gamma Malthusian rate and Poisson driver counts against theory; and
  exact-value tests for distances, MRCA, coalescence and forests.

### Performance

Measured on a 2×10⁵-cell tree on Julia 1.10; `simulate!` alone takes about 0.4 s.

| Call | 0.3 | 0.4 |
|---|---|---|
| `mutations_per_cell(root)` | 1 736 ms | 44 ms |
| `mutations_per_cell(pop)` | 164 ms | 7 ms |
| `branch_spectrum(root)` | 1 118 ms | 164 ms |
| `sitefrequencyspectrum(root)` | 1 042 ms | 152 ms |
| `sitefrequencyspectrum(pop)` | 196 ms | 154 ms |
| `sample_leaves(root, 1000)` | 404 ms | 45 ms |
| `simulate!` with `trajectory_dt = 0.1` | 2.11 s | 0.86 s |

### Unchanged (verified)

- `simulate!` random-stream consumption for single-founder runs, including chained
  runs, hooks and restarts from one founder.
- `sitefrequencyspectrum`, `branch_spectrum`, `leaf_depths` (including its frozen
  order), `fitness_per_cell`, `pairwisedistance` and `coalescence_times` values.
