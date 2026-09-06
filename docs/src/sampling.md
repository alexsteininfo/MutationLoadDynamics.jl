# Sampling

Real experiments sequence a few hundred or a few thousand cells out of a population of
millions. A simulated observable is only comparable to data once it has been through the
same bottleneck, so sampling is a first-class operation here rather than an afterthought.

```julia
out = sample_trees(pop, SamplingSpec([1000, 100]); seed = UInt64(0xBEEF))

out.full                              # the whole tree, still there
s = out.samples[1]                    # a LeafSample: n = 1000
sitefrequencyspectrum(s.root, s.n)    # the sample's SFS
mutations_per_cell(s.root)            # each sampled cell's FULL-TREE burden
leaf_depths(s.root)                   # each sampled cell's FULL-TREE depth
```

Sampling is **post-hoc**: it applies to a tree that has finished growing. It is
deliberately not a field of [`NonMarkovBlock`](@ref) or [`MeasurementSpec`](@ref), both of
which describe things that happen *during* `simulate!`.

## One draw

```julia
s = sample_leaves(root, n; seed = UInt64(1))
s = sample_leaves(pop,  n; seed = UInt64(1), replicate = 2)
```

[`sample_leaves`](@ref) draws `n` of the tree's leaves uniformly without replacement and
returns the induced tree as a [`LeafSample`](@ref):

| Field | Contents |
|:---|:---|
| `root` | the induced tree |
| `n` | cells drawn |
| `N_full` | leaf count of the source tree |
| `seed` | the seed this draw used |
| `replicate` | which independent draw at this `(tree, n)` this is |
| `sampled_ids` | `NonMarkovCell.id` of each drawn cell, in draw order |

It is **non-destructive** — the source tree is untouched, because the same tree is
normally drawn from again at other sizes and those draws must be independent. Cost is
proportional to the number of *retained* nodes, not to the size of the source tree.

## Prune, do not collapse

The induced tree keeps the sampled leaves **plus every ancestor of a sampled leaf**, and
retains the resulting unary nodes rather than merging them into their children.

This is the load-bearing decision of the whole module. Because every division ancestral to
a sampled cell is still a node, a sampled cell's root-to-leaf path is unchanged — so
[`mutations_per_cell`](@ref) and [`leaf_depths`](@ref) on the sampled tree return exactly
that cell's **full-tree** burden and divisional depth. Collapsing would turn depth into a
count of bifurcations that happened to survive sampling, which is a property of the
sample rather than of the cell.

It is also exactly the shape [`prune_tree!`](@ref MutationLoadDynamics.prune_tree!)
already leaves behind when a lineage dies out during a run. A sampled tree is therefore a
`BinaryNode{NonMarkovCell}` indistinguishable in kind from a full one, and every function
on the [Tree statistics](statistics.md) page applies to it unchanged.

Two consequences worth internalising:

!!! warning "The root of a sampled tree is the founder, not the MRCA of the sample"
    Ancestry is retained all the way up, so `s.root` is the original founding cell even
    when it has a single child. That is what makes `sitefrequencyspectrum(s.root, s.n)`
    put the founder's mutations in `sfs[n]` exactly as the full tree puts them in
    `sfs[N]`. If you want the sample's own MRCA, ask for it:

    ```julia
    findMRCA(getalivecells(s.root))
    ```

!!! note "The sampled tree shares `data` with the source tree"
    `NonMarkovCell` is immutable, so ids, birthtimes, mutation counts and fitness are
    identical by construction rather than by copy — sampling allocates new `BinaryNode`s
    but no new cells. To change something in a sampled tree, replace a node's `data`
    wholesale (which rebinds only the sample); never mutate a cell in place, because
    there is no such thing as a cell that belongs to only one of the trees.

## Declaring what to produce

[`SamplingSpec`](@ref) says what [`sample_trees`](@ref) should build from one finished
tree.

| Intent | Spec |
|:---|:---|
| full tree only | `SamplingSpec()` |
| one sample of size `n`, nothing else | `SamplingSpec(sizes = [n], retain_full = false)` |
| full tree plus several sizes | `SamplingSpec([1000, 100])` |
| several independent draws per size | `SamplingSpec(sizes = [1000], replicates = 20)` |

- `sizes` — the sample sizes to draw. Duplicates are rejected: use `replicates` for
  repeated draws at one size.
- `replicates` — independent draws per size, each with its own derived seed.
- `retain_full` — whether the returned [`SampledTrees`](@ref) carries the full tree.

`SamplingSpec(n)` and `SamplingSpec([n1, n2])` are shorthands for the keyword form.
Validation happens in the constructor, which is the only way to build one — there is no
unvalidated positional path in.

[`sample_trees`](@ref) checks every size against the tree **before** making any draw, so
an oversized request fails immediately rather than after minutes of work. Results are
ordered by the spec's `sizes`, and within a size by `replicate`.

```julia
out = sample_trees(pop, SamplingSpec(sizes = [1000, 100], replicates = 5);
                   seed = UInt64(20260906))

length(out.samples)                    # 10 = 2 sizes × 5 replicates
[(s.n, s.replicate) for s in out.samples]
```

!!! note "`retain_full = false` does not free memory by itself"
    It controls only what the returned bundle holds. It cannot release your own reference
    to the tree, and it cannot release the tree reachable from a `Population` — every
    living leaf holds a `parent` chain back to the founder, so while the population is
    alive the whole history is alive. To actually bound memory across a sweep, drop both:

    ```julia
    out = sample_trees(pop, SamplingSpec(sizes = [1000], retain_full = false);
                       seed = myseed)
    pop = nothing
    GC.gc()
    ```

## Reproducibility

`seed` is required, not optional. The draw is a pure function of `(tree, n, seed)`, and
each `LeafSample` records the seed it used, so any single draw replays in isolation from
that record alone — you do not need the spec, the batch, or the order the batch ran in.

Callers supply the seed because they are the ones who can derive it from their own
provenance: a filename stem, a simulation index, a sweep coordinate. Two draws that
should differ must be given different seeds. `sample_trees` handles this for you inside
one batch, deriving a distinct per-draw seed from the base seed, the size and the
replicate index.

Each call builds its own `MersenneTwister` from the seed and touches no shared state, so
independent draws are safe to run concurrently from your own threads.

!!! warning "The draw recipe is frozen"
    `MersenneTwister(seed)` and `randperm(rng, N_full)[1:n]` over `Leaves(root)` order are
    load-bearing, not incidental: serialised sampled trees produced by earlier runs must
    stay reproducible. A golden test pins them. Do not substitute a different rng,
    `StatsBase.sample`, or a direct `n`-index draw — each would silently invalidate
    stored results.

## Sampling needs a single-rooted tree

`sample_leaves` and `sample_trees` throw an `ArgumentError` naming the number of roots
when handed a forest — a population from `initialize_population(N)` with `N > 1`, where
the founders share no ancestry. There is no sensible uniform draw across independent
trees that also preserves "the root is the founder", so sample each root's tree
separately:

```julia
using AbstractTrees

roots   = AbstractTrees.getroot(allcells(pop))   # the distinct roots, one per tree
samples = [sample_leaves(r, min(100, popsize(r)); seed = UInt64(i))
           for (i, r) in enumerate(roots)]
```

## A typical sweep

```julia
using MutationLoadDynamics, Distributions, Random, Serialization

for (i, s) in enumerate(0.0:0.1:0.5)
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, block(s), MersenneTwister(i))

    out = sample_trees(pop, SamplingSpec(sizes = [1000, 100], replicates = 10),
                       seed = hash(("sweep-2026-09", i)))

    serialize("s=$(s).jls", out.samples)   # small: the samples, not the population
    pop = nothing
    GC.gc()
end
```

Serialising the samples rather than the population is the point: a `LeafSample` at
``n = 1000`` is a tree of a few thousand nodes, against a population of ``10^6`` cells
whose full history is orders of magnitude larger — and it replays exactly from its
recorded seed if you ever need to prove where it came from.
