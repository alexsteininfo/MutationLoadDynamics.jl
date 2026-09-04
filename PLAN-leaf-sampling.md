# PLAN — leaf sampling of lineage trees

**Status:** implemented and tested, on branch `feat/leaf-sampling`, released as
v0.3.0. Written 2026-09-04 from the consumer side (`gITH-nonMarkovian`).

**How to use this file:** open a Claude Code session in this repo
(`/Users/alexanderstein/Documents/GitHub/MutationLoadDynamics.jl`), tell it to read
`PLAN-leaf-sampling.md`, then plan → implement → test. Section 8 lists what still
needs a call; everything else is settled requirements.

**This file is self-contained.** The receiving session has no access to the
conversation it came from.

---

## 1. Why

A working implementation of uniform leaf sampling already exists in the consumer repo
`gITH-nonMarkovian` at `analysis/helpers/subsampling.jl`, has generated 4.2 GB of
output across 525 shards, and is covered by a 598-assertion suite
(`analysis/subsampling/verify_subsampling.jl`). It is being promoted into this
package because three separate consumers now need it:

1. **`gITH-nonMarkovian`** — the existing SFS / scMB study (paper 1). Already uses it.
2. **`CopyNumberEvolution.jl`** — a new package that turns a lineage tree into
   per-cell copy-number profiles. Its first step is literally "derive a sampled tree
   from the full tree carrying the same information".
3. **A second study repo** for the copy-number paper (paper 2), not yet created.

Sampling is a pure tree → tree operation with no dependency on any study's parameter
grid, filenames or scenario types, so it belongs here rather than being copied three
times. The wider repo architecture this is part of: this package generates trees;
`CopyNumberEvolution.jl` is the observation model; `EvoTracer.jl` holds
estimators and must never depend on a simulator; study repos own parameter grids,
file naming and figures.

**The governing rule for what lands here:** *packages own algorithms and
observation-model types; study repos own parameter records, file naming, grids and
figures.* Section 6 spells out what that rule keeps out of this package, and why
following it means **zero migration cost for the 4.2 GB already on disk**.

---

## 2. What sampling means here

Draw `n` of a tree's `N` leaves uniformly without replacement, and return the
**induced** lineage tree: the sampled leaves plus *every ancestor of a sampled leaf*,
with the resulting unary nodes **retained, never collapsed**.

### 2.1 Prune, do not collapse — the load-bearing decision

This is the one design choice the whole feature rests on, and it is easy to get
backwards. Because every division ancestral to a sampled cell remains a node, a
sampled cell's root-to-leaf path is unchanged. Therefore:

- `mutations_per_cell(root)` returns each sampled cell's **full-tree** mutational
  burden, not a burden restricted to surviving branches.
- Divisional depth of a sampled leaf is its **full-tree** divisional depth.
- Each retained edge is still exactly one division, so `node.data.mutations` remains
  a per-edge mutation count and `birthtime` differences remain real-time edge lengths.

Collapsing unary nodes would turn depth into "a count of bifurcations that survived
sampling" — a property of the sample rather than of the cell — and would silently
break every burden-based statistic. It would also break `CopyNumberEvolution.jl`,
which needs one CNA-drawing opportunity per real division along the path.

It happens to be free: this is exactly the shape `prune_tree!`
(`src/simulation_trees.jl:7`) already leaves behind when a lineage dies out, so a
sampled tree is a `BinaryNode{NonMarkovCell}` indistinguishable in kind from a full
tree, and every existing tree statistic applies to it unchanged.

### 2.2 The founder stays the root

The founder is retained as the root even when it ends up with a single child. This is
deliberate: it makes `sitefrequencyspectrum` accumulate the founder's mutations into
`sfs[n]` for a sampled tree exactly as it accumulates them into `sfs[N]` for a full
tree.

Consequence to be aware of and to document: **the root of a sampled tree is the
original founder, not the MRCA of the sample.** Consumers that want the sample's MRCA
must call `findMRCA` on the sampled leaves.

---

## 3. Proposed API

Three layers. Only layer 1 is strictly required; layers 2 and 3 are what make the
"full only / sample only / both" choice ergonomic, which is the consumer-side request
driving this work.

### 3.1 Layer 1 — the single draw (required, and the reproducibility anchor)

```julia
"""
    sample_leaves(root, n; seed, replicate=1) -> LeafSample
    sample_leaves(population, n; seed, replicate=1) -> LeafSample
"""
struct LeafSample
    root::BinaryNode{NonMarkovCell}   # induced tree: sampled leaves + all ancestors
    n::Int                            # cells drawn
    N_full::Int                       # leaf count of the source tree
    seed::UInt64                      # rng seed of this draw
    replicate::Int                    # which independent draw at this (tree, n)
    sampled_ids::Vector{Int64}        # NonMarkovCell.id of each drawn cell, draw order
end
```

The reference implementation to port is `subsample_tree` + `_copy_marked` in
`gITH-nonMarkovian/analysis/helpers/subsampling.jl`. It is correct and reviewed; port
it rather than rewriting it, then extend.

`seed` is a **required keyword taking a caller-supplied `UInt64`**, not an optional
one. See HR-1: this is what preserves bit-exact reproduction of the 525 existing
shards, whose seeds were derived from a study-side filename stem the package cannot
know.

### 3.2 Layer 2 — what to produce (the consumer-side request)

The consumer wants one declaration covering: *full tree only*, *one sample at a
defined size*, or *the full tree plus one or several sample sizes*.

```julia
@kwdef struct SamplingSpec
    sizes::Vector{Int}   = Int[]    # sample sizes to draw; empty = draw nothing
    replicates::Int      = 1        # independent draws per size
    retain_full::Bool    = true     # keep the full tree in the returned bundle
end

struct SampledTrees
    full::Union{BinaryNode{NonMarkovCell}, Nothing}
    samples::Vector{LeafSample}
end

sample_trees(root_or_population, spec::SamplingSpec; seed::UInt64) -> SampledTrees
```

The three requested modes:

| intent | spec |
|---|---|
| full data only | `SamplingSpec()` |
| one sample of size `n` only | `SamplingSpec(sizes=[n], retain_full=false)` |
| full data plus several sizes | `SamplingSpec(sizes=[1000, 100])` |

Add convenience constructors `SamplingSpec(n::Int)` and
`SamplingSpec(sizes::Vector{Int})` if they read well; keep `@kwdef` as the primary
form.

**Be honest in the docstring about what `retain_full = false` does and does not do.**
It controls what the returned `SampledTrees` holds. It cannot free the caller's own
reference, and it cannot free the tree reachable from `Population.cells` — every
alive leaf holds a `parent` chain back to the founder, so as long as the `Population`
is alive the full tree is alive. The memory win requires the caller to drop both. The
right thing is to document the pattern, not to add package-side cleverness:

```julia
spec = SamplingSpec(sizes=[1000], retain_full=false)
out  = sample_trees(pop, spec; seed = myseed)
pop  = nothing            # caller must do this
GC.gc()
```

`sizes` is validated (`1 <= n <= N_full`, no duplicates unless `replicates > 1`), and
sampling is ordered largest-first so a too-large `n` fails before work is done.

### 3.3 Layer 3 — the `simulate!` convenience (optional, decide in §8)

A one-call wrapper `simulate_and_sample!(pop, block, rng, spec)` is tempting but adds
little: `simulate!` already returns the population and `sample_trees` takes one.
Recommend **not** adding it in the first pass.

Do **not** attach sampling to `NonMarkovBlock` or to `MeasurementSpec`. Sampling is
post-hoc, applied once to a finished tree; the block describes dynamics and the
measurement machinery fires *during* the run. Putting a sample-size list on either
would misrepresent when it happens and would interact badly with
`restart_on_extinction` (which resets the accumulator) and with chained blocks.

### 3.4 `replicates` — decide it now, because later is expensive

The consumer's current design is one draw per `(simulation, n)`, and its
`SubsampleResult` struct deliberately has **no replicate field**: adding one would
make its 525 serialized shards unreadable. That constraint is a property of the
*consumer's* record type, not of sampling.

This package is defining `LeafSample` fresh, with no serialized data behind it, so
`replicate` costs one `Int` now and is genuinely awkward to retrofit later. Include
it. Per-draw seeds must then depend on `replicate` (HR-2), and the consumer can adopt
replicates whenever it introduces a v2 record type.

### 3.5 Companion tree statistics (in scope, decided)

The consumer's `analysis/helpers/tree_analysis.jl` holds a family of tree statistics
that are all pure functions of a `BinaryNode` root with no study parameters in them.
By the ownership rule (§6.1) they belong here, and the author has asked for the whole
family to move together rather than only the one the tests need. Moving them also
removes a genuine duplication: `_fill_sfs!` in that file is **character-for-character
the same algorithm** as `_sfs_fill!` in `src/statistics.jl:233`, reached by a
different entry point.

| consumer helper | disposition here |
|---|---|
| `compute_mut_per_cell(root)` | already exists as `mutations_per_cell(root)`; the helper is a one-line wrapper today. Nothing to do. |
| `compute_sfs(root, N)` | **add** a root-based `sitefrequencyspectrum(root, N)`. `sitefrequencyspectrum(population)` already exists and already calls `_sfs_fill!`; this is a second entry point taking an explicit length, which is what a sampled tree needs (`N = n`, not the population size). |
| `compute_leaf_depths(root)` | **add** as `leaf_depths(root)`. |
| `compute_branch_spectrum(root, N)` | **add** as `branch_spectrum(root, N)` — internal nodes subtending exactly `k` leaves; the topological SFS. |
| `compute_filtered_mut_per_cell(root, threshold)` | **add** as `filtered_mutations_per_cell(root, threshold)`. |
| `compute_leaf_fitness(root)` | **add** as `leaf_fitness(root)`. |

`leaf_depths` in particular is what makes the depth-invariance test of §5 expressible
here at all, which is the test that catches an accidental unary collapse.

**HR-10 applies to every one of these** (§4): the move must be output-preserving to
the element, including vector order, because 4.3 GB of stage-2 arrays were produced by
the current implementations. Two order facts to carry across deliberately rather than
rediscover:

- `mutations_per_cell` and `leaf_fitness` both iterate `getalivecells(root)`, so they
  are **co-indexed**: entry `i` is the same cell in both.
- `leaf_depths` walks its own explicit stack in a **different** order and is
  therefore **not** co-indexed with either. It is valid only as a pooled
  distribution. Do not "fix" this by switching it to `Leaves` order — that would
  change the stored arrays. If a co-indexed variant is wanted later, add it as a
  separate keyword or function, and say so in the docstring.

Two API details worth settling while moving:

- `sitefrequencyspectrum(root, N)` and `branch_spectrum(root, N)` should **error**
  with a clear message when the tree has more than `N` alive leaves, rather than
  throwing a raw `BoundsError` from `sfs[count]`. The consumer always passes the right
  `N`, so this is a guard against a future caller, not a bug fix.
- Consider defaulting `N` to the tree's alive-leaf count, so
  `sitefrequencyspectrum(root)` works without the caller counting leaves. Keep the
  explicit two-argument form, since the consumer relies on it.

---

## 4. Hard requirements

**HR-1 — bit-exact reproduction of existing draws.** The consumer has 4.2 GB of
sampled trees on disk whose draws must remain reproducible: its verifier's strongest
test asserts that re-drawing reproduces `data/processed/` arrays exactly, and a forced
regeneration must not silently produce *different* samples than the ones already
analysed. The draw is therefore pinned to this exact recipe, and the port must not
"improve" any part of it:

```julia
leaves = collect(Leaves(root))      # AbstractTrees order; children returns (left, right)
N_full = length(leaves)
rng    = MersenneTwister(seed)      # exactly this RNG, seeded with exactly this UInt64
idx    = randperm(rng, N_full)[1:n] # exactly this — not randperm!, not sample(), not shuffle
```

Any of the following changes the draw and is forbidden: a different RNG type,
`StableRNG`, `Random.seed!` on a global RNG, `StatsBase.sample`, drawing `n` indices
directly instead of taking the first `n` of a full permutation, or iterating leaves in
any order other than `Leaves(root)`.

Note `getalivecells(root)` (`src/simulation_trees.jl:117`) is a comprehension over
`Leaves(root)` filtered by `isalive`, and `isalive(::NonMarkovCell) = true` — dead
cells are pruned out of the tree entirely by `celldeath!`. So `getalivecells` and
`collect(Leaves(root))` are the same vector in the same order for any tree this
package produces, and either is acceptable. Use one, state in a comment why the order
is what makes the draw reproducible, and add the regression test in §5.

There is a residual risk this cannot defend against: `randperm`'s output for a given
`MersenneTwister` seed is not a documented stability guarantee across Julia
versions. Record the Julia version used to generate the existing data in the test
file, and make the golden-fixture test (§5) the thing that detects a change.

**HR-2 — seed derivation is explicit and per-draw.** Layer 1 takes a caller-supplied
`UInt64`. Layer 2 derives per-draw seeds from its own base seed as
`hash((base_seed, n, replicate))`, and the derived value is stored in
`LeafSample.seed` so any single draw replays in isolation from what the output file
records. Two different `(n, replicate)` pairs must not share a seed.

The consumer keeps its own `sample_seed(stem, sim_index, n) = hash((stem, sim_index,
n))` and passes the result into layer 1 — this is exactly why layer 1 must not own
seed derivation.

**HR-3 — non-destructive.** `root` is never mutated. The same tree is drawn from again
at other sample sizes and other replicates, and those draws must be independent. The
existing suite has a dedicated test for this; keep it.

**HR-4 — faithful sub-shape.** Left/right slots are preserved: a node whose left
lineage was dropped keeps `left = nothing`. The tree must not be silently
re-balanced. Parent links in the rebuilt tree must be consistent
(`child.parent === node` for every retained edge, `root.parent === nothing`).

**HR-5 — share `data`, don't copy it.** `NonMarkovCell` is immutable, so the sampled
tree shares `node.data` with the source tree. This keeps sampling cheap and keeps ids,
birthtimes, mutation counts and fitness identical by construction rather than by
copy. Document that a consumer must not attempt to mutate `data` on a sampled tree
(it would be visible in the source tree too, if the source is still alive) — the
supported mutation idiom, as with `on_division`, is to replace `node.data` wholesale,
which is safe because it rebinds only the sampled tree's node.

**HR-6 — cost is proportional to the sample, not the tree.** Ancestor marking stops
at the first already-marked node, so total work is O(retained nodes), not O(N). This
matters: the consumer draws from trees of 16384 leaves inside shards of up to 400 MB.

**HR-7 — no threading inside the package.** The consumer threads over simulations
within a shard (`Threads.@threads` in `subsample_shard`), and that works precisely
because each draw builds its own `MersenneTwister` from its own seed. Keep the
package's functions single-draw and thread-free so they compose; document that
independent draws are safe to run concurrently.

**HR-8 — errors, not warnings, on bad input.** `n > N_full`, `n < 1`, an empty tree,
and a forest (no single root) all error with a message naming the offending values.
Silent degradation here would corrupt a downstream sweep in a way that only shows up
as a wrong figure.

**HR-9 — no new dependencies, and no file IO.** Sampling needs only
`AbstractTrees`, `Random` and `DataStructures`, all already dependencies. This
package must not gain `Serialization`: serialization format, atomic writes, filenames
and resumability are study-repo concerns (§6).

**HR-10 — the statistics move is output-preserving.** Every function promoted in
§3.5 must return exactly what the consumer's version returns, element for element and
in the same order, for the same tree. 4.3 GB of stage-2 arrays under
`data/processed/` were produced by the current implementations, and stage 2 is
resumable-by-`isfile`, so a subtly reordered output would mix old and new conventions
inside one dataset without any error. Port the bodies as they are; do not tidy them.
The consumer will verify this against real shards before deleting its helpers (§6.2).

---

## 5. Tests

New file `test/sampling.jl`, added to the `tests` vector in `test/runtests.jl`.

The consumer's `analysis/subsampling/verify_subsampling.jl` is the source to port
from. **Take only the tree-level assertions**; the rest are study-specific and stay
where they are (§6.2). Specifically, port:

- **Golden fixture** — a small deterministic tree built in the test file; assert the
  exact `sampled_ids` for a fixed `(n, seed)`. This is the HR-1 regression test.
  Record the Julia version in a comment.
- `n = N_full` reproduces the source tree exactly — same leaf set, same ids, same
  per-cell burdens, same depths, same SFS. The single most valuable test: it is the
  only one that checks the induced tree against a *known* answer.
- `n = 1` over many seeds — the induced tree is a path; its single leaf's burden and
  depth match the source tree's for that cell; every leaf is eventually drawn.
- `n = 2` keeps exactly the branching node where the two lineages meet, and nothing
  below it on the dropped side.
- Ancestor retention — every retained node has a sampled leaf below it, and every
  ancestor of every sampled leaf is retained (both directions).
- Unary nodes retained, left/right slots preserved (HR-4).
- Parent links consistent throughout the rebuilt tree (HR-4).
- Determinism — same `(tree, n, seed)` gives the same draw; different seeds give
  different draws; different `n` at the same seed are independent draws.
- Source tree unmutated after several draws (HR-3).
- `data` is shared, not copied — `sampled.data === source.data` (HR-5).
- Bounds and error messages (HR-8): `n = 0`, `n = N_full + 1`, single-node tree.
- **Invariance** — `mutations_per_cell` and `leaf_depths` of the sampled cells equal
  their values in the full tree (§2.1, §3.5). This is the test that catches an
  accidental collapse.
- `SamplingSpec` behaviour — `retain_full` both ways; several sizes in one call;
  `replicates > 1` gives distinct draws with distinct stored seeds; largest-first
  validation errors before any work.
- A run through `simulate!` end to end: grow a population, sample it, check
  `N_full == popsize(pop)` and that the sampled ids are a subset of
  `keys(pop.cells)`.

Two tests the consumer has that should **not** move: the `Serialization`
round-trip of its record type, and `serialize_atomic`. Both are about its file
layout.

---

## 6. What stays in the consumer repo — and why this migration is free

### 6.1 The ownership split

| stays in `gITH-nonMarkovian` | reason |
|---|---|
| `SubsampleResult{P}` (`analysis/helpers/types_subsampled.jl`) | carries `SimParams`/`Sel1Params`/`Sel2Params`, `sim_index`, and a filename convention — all study-scoped |
| `sample_sizes(N_target)` fixed table (`1000 -> [100]`, `16384 -> [1638, 164]`, …) | a parameter-grid decision, deliberately a table so a new `N_target` fails loudly |
| `sample_seed(stem, sim_index, n)` | depends on the study's filename stem |
| `serialize_atomic`, `subsample_file`, `subsample_shard` | file IO, naming, resumability, threading policy |
| co-indexing with `data/processed/` | a property of the study's array layout |

**This is what makes the migration cost zero for the 4.2 GB already on disk.** Julia's
`Serialization` resolves a struct by module path and name. `SubsampleResult` is
currently defined in `Main` via `include`, and the 525 existing shards record it that
way. Moving that struct into this package would make every one of them unreadable.
Leaving the *record* in the consumer and promoting only the *algorithm* avoids the
problem entirely — the consumer constructs its unchanged `SubsampleResult` from a
`LeafSample`'s fields.

Do not propose moving `SubsampleResult` here. If a future session is tempted, the
answer is a new type with a new filename convention, never a change to that one.

### 6.2 What the consumer session will do afterwards (not part of this work)

Listed so the implementing session knows the contract it is being held to, and does
not try to do it from here:

1. `analysis/helpers/subsampling.jl` keeps `serialize_atomic`, `sample_sizes`,
   `sample_seed`, `subsample_shard`, `subsample_file`; `subsample_tree` and
   `_copy_marked` are deleted and `subsample_shard` calls
   `MutationLoadDynamics.sample_leaves(sim.tree_root, n; seed = seed)`.
2. `analysis/helpers/tree_analysis.jl` is deleted, and its six call sites across
   `analysis/processing/`, `analysis/processing_subsampled/` and `analysis/plots/`
   switch to the package names from §3.5.
3. `verify_subsampling.jl` drops the assertions that moved here and keeps the rest.
4. **Two equivalence checks before anything is deleted.** (a) HR-1: re-draw a handful
   of simulations with the package function and compare `sampled_ids` against the
   existing shards on disk. (b) HR-10: recompute all six statistics on real shards
   with both the old helpers and the new package functions and assert
   element-for-element equality. If either differs, stop — do not regenerate data.
5. Re-`Pkg.develop` and re-instantiate so `Manifest.toml` picks up the new version.

---

## 7. Versioning

New feature, no breaking change to any existing API: bump `Project.toml` to
**0.3.0**. The consumer's `Manifest.toml` pins this package by local dev path, so it
tracks the working copy — but tag the version anyway, because a third and fourth
consumer are about to appear and "which version of the sampler drew this data" needs
an answer. See the note in §8 on a local registry.

---

## 8. Open design decisions

1. **`replicates` in `LeafSample` — in or out?** Recommendation: **in** (§3.4). Costs
   one field now, genuinely awkward later, and the consumer's inability to persist it
   today is a property of its own frozen record type, not of this package.
   **Resolved: in.** `LeafSample.replicate::Int` shipped as specified, defaulting to
   1 in `sample_leaves` and derived per-draw by `sample_trees`/`SamplingSpec`.
2. **Add `leaf_depths` to `src/statistics.jl`?** Recommendation: **yes** (§3.5) — the
   depth-invariance test cannot be written here without it.
   **Resolved: yes.** `leaf_depths(root::BinaryNode)` shipped in `src/statistics.jl`
   and is exported.
3. **`sample_trees` vs `subsample_trees` vs `sample_population` for the layer-2
   name.** No strong view; pick one and be consistent. Avoid `subsample_tree` for the
   layer-1 function even though that is the consumer's current name, since layer 1
   returns a `LeafSample` rather than the old three-tuple — a different name makes the
   consumer-side migration a compile error rather than a silent shape change.
   **Resolved: `sample_trees`.** Layer 1 shipped as `sample_leaves` (distinct from the
   consumer's old `subsample_tree` name, as required); layer 2 as `sample_trees`.
4. **`simulate_and_sample!` convenience wrapper?** Recommendation: **no** for now
   (§3.3).
   **Resolved: no.** Not implemented; `simulate!` and `sample_trees` are called
   separately, as recommended.
5. **Sampling schemes other than uniform-without-replacement?** Out of scope. The
   hypergeometric projection the consumer's theory relies on (`theory/sfs.md`) assumes
   exactly this scheme; a `SamplingScheme` abstraction with one implementation would
   be speculative. Note it as a future extension point and move on.
   **Resolved: stayed out of scope.** No `SamplingScheme` abstraction was added;
   `sample_leaves`/`sample_trees` implement uniform-without-replacement only.
6. **Local registry (repo-family concern, mentioned here for context).** This package
   is about to have three dev-path consumers, and `Manifest.toml` portability is
   already a known problem (an absolute `/Users/alexanderstein/...` path). A private
   `LocalRegistry.jl` registry would let consumers depend on versioned releases
   instead of dev paths. Not part of this work; raised because tagging 0.3.0 is the
   moment it becomes worth doing.
   **Resolved: deferred, as scoped.** Not part of this work; consumers still depend
   on this package via `Pkg.develop` local paths. Remains open for a future session.

---

## 9. Suggested order of work

1. Port `sample_leaves` + `LeafSample` (layer 1) with the HR-1 recipe intact.
2. Write `test/sampling.jl`, golden fixture first, then `n = N_full`, then the rest.
3. Promote the tree statistics of §3.5 into `src/statistics.jl` under HR-10, with
   their own tests (including the `leaf_depths` / `mutations_per_cell` invariance
   test, which needs step 1 to be done).
4. Add `SamplingSpec` / `SampledTrees` / `sample_trees` (layer 2) and its tests.
5. Export the new names from `src/MutationLoadDynamics.jl`; document sampling in
   `README.md` including the "prune, do not collapse" rationale and the
   `retain_full` memory caveat.
6. Bump to 0.3.0.
