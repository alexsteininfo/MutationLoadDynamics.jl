# Tree statistics

Everything here is computed after a run, from the lineage tree. Two families of methods
exist, and picking the right one is mostly about whether you want a *stable* ordering.

**Population methods** take a [`Population`](@ref) and iterate its `Dict` of living
cells. Convenient, but the order is arbitrary.

**Root methods** take a `BinaryNode` and traverse the tree. The order is deterministic
(`Leaves` order, which follows `(left, right)` at every node), and they work on any
subtree — including the induced tree from [`sample_leaves`](@ref), which is why every
statistic applies unchanged to a sample.

```julia
root = getsingleroot(allcells(pop))
```

## Mutational burden

```julia
mutations_per_cell(pop)            # Vector{Int64}, one per living cell
mutations_per_cell(root; includeclonal = true)   # same, in Leaves(root) order
average_mutations(pop)             # mean burden
mean_k(pop), var_k(pop)            # mean and variance of the burden
clonal_mutations(pop)              # drivers carried by every living cell
```

A cell's burden is the sum of `mutations` over its whole ancestral path. It is stored
on the cell as `total_mutations`, so these calls are ``O(N)``.

For a root method, `includeclonal` says whether mutations carried by *every* leaf under
`root` count:

```julia
sub = findMRCA(some_cells)
mutations_per_cell(sub)                        # acquired strictly below sub (default)
mutations_per_cell(sub; includeclonal = true)  # full burden, back to the top of the tree
```

With the default `false`, `root`'s own mutations and those of its ancestors are left
out: they are clonal in that subtree. For the founder of a simulated tree the two agree,
because the founder carries no mutations.

[`clonal_mutations`](@ref) is the burden at the MRCA of all living cells — the mutations
no observable variation can distinguish, since every cell has them. It returns `0` when
there is no single MRCA.

### Dropping high-frequency variants

```julia
filtered_mutations_per_cell(root, 0.1)   # drop nodes subtending > 10% of leaves
filtered_mutations_per_cell(root, 1.0)   # drop nothing
```

[`filtered_mutations_per_cell`](@ref) excludes mutations from any ancestral node whose
living-descendant count exceeds `floor(threshold * N)`. This is the tree-side analogue of
filtering out high-frequency variants before estimating a mutation rate: a node
subtending a large fraction of the population contributes identically to every cell below
it, so it carries no information about within-clone divergence. It counts from each leaf
up to and including `root`, and with `threshold = 1.0` on the top of a tree equals
`mutations_per_cell(root; includeclonal = true)`. Returned in
`Leaves(root)` order.

## Spectra

```julia
sitefrequencyspectrum(pop)          # length popsize(pop)
sitefrequencyspectrum(root)         # length = the tree's living-leaf count
sitefrequencyspectrum(root, N)      # length N
branch_spectrum(root)               # topological SFS
branch_spectrum(root, N)
```

`sfs[k]` is the number of **driver mutation events** carried by exactly `k` of the tree's
leaves; `bs[k]` is the number of **internal nodes** subtending exactly `k` leaves. One
post-order traversal each.

The explicit `N` exists for sampled trees, where the meaningful length is the sample size
rather than whatever the induced tree happens to contain. `N` larger than the leaf count
pads with zeros; `N` smaller is an `ArgumentError` rather than a `BoundsError`.

```julia
s = out.samples[1]
sitefrequencyspectrum(s.root, s.n)   # a length-n spectrum, comparable across draws
```

Mutations on the root's own edge accumulate into `sfs[n]`, where `n` is the tree's actual
leaf count — they are clonal *in this tree*. Note `n`, not `N`: a zero-padded spectrum
keeps its clonal bin at `n`.

### Separating topology from mutation rate

`branch_spectrum` is the SFS a tree would have at a neutral per-division rate of one. For
any neutral rate ``m``:

```math
\mathbb{E}[\texttt{sfs}[k]] = m \cdot \texttt{bs}[k] \quad (k \ge 2), \qquad
\mathbb{E}[\texttt{sfs}[1]] = m \cdot (\texttt{bs}[1] + n).
```

The ``k = 1`` case needs the extra ``n`` because each leaf's own edge also carries
mutations into `sfs[1]`, while `bs` counts only internal nodes — so `bs[1]` is the number
of **unary** nodes (divisions whose other daughter's lineage died out), not the number of
leaves. Getting this wrong understates the singleton class by a factor of roughly ``N/bs[1]``.

The practical value: one simulated tree serves every ``m``, and a fitted ``m`` never
absorbs topology.

## Depths, fitness and lifetimes

```julia
leaf_depths(root)          # divisions from root to each leaf
leaf_fitness(root)         # fitness of each living leaf
celllifetimes(root)        # completed cell lifetimes (dividers only, by default)
celllifetimes(root; excludeliving = false)   # include cells still alive
celllifetime(node, pop.t)  # one cell
endtime(node)              # when this cell divided; nothing if alive
```

[`leaf_depths`](@ref) is the primary observable for a neutral run: it is the count of
real divisions, so unary nodes are included and correctly so. Combined with a rate ``m``
it gives per-cell burdens directly, without simulating passengers:

```julia
burden = [rand(rng, Poisson(m * d)) for d in leaf_depths(root)]
```

`celllifetimes` is the realised cell-cycle distribution, and it is **not** an unbiased
sample from `birth_dist`. Two separate effects shift it earlier:

- **Competing risks.** A cell that dies never completes a division, so what you see is the
  division time *conditioned on dividing first*.
- **Stopping at a fixed size.** By default only cells that have already divided
  contribute, and at any moment those are biased toward the short cycles — the long ones
  are still in progress. This bias is present even with no death at all: a pure-birth,
  neutral run with a true mean cycle of 1.0 returns a mean completed lifetime near 0.87.

`excludeliving = false` adds the cells still alive, but their "lifetime" is then
right-censored at `age(root)` rather than completed, so it does not fix the bias either.
Treat `celllifetimes` as the distribution an experiment measuring interdivision times
would actually observe — which is its point — and not as a way to confirm `birth_dist`
was wired up correctly. For that, sample the distribution directly.

## Which results are co-indexed

This is the one thing to get right, and it is not symmetric.

| Function | Order | Co-indexed with |
|:---|:---|:---|
| `mutations_per_cell(pop)` | `Dict` order (arbitrary) | `fitness_per_cell(pop)`, `allcells(pop)` |
| `fitness_per_cell(pop)` | `Dict` order (arbitrary) | `mutations_per_cell(pop)`, `allcells(pop)` |
| `mutations_per_cell(root)` | `Leaves(root)` | `leaf_fitness`, `filtered_mutations_per_cell` |
| `leaf_fitness(root)` | `Leaves(root)` | `mutations_per_cell(root)`, `filtered_mutations_per_cell` |
| `filtered_mutations_per_cell(root, θ)` | `Leaves(root)` | the two above |
| `leaf_depths(root)` | **its own stack order** | **nothing** |

!!! warning "`leaf_depths` is not co-indexed with anything"
    It uses an explicit stack whose order is neither `Leaves` order nor `Dict` order. It
    is valid as a **pooled distribution** — a histogram of divisional depths — and not as
    a per-cell vector to be paired with a burden or a fitness.

    The order is kept fixed on purpose: large volumes of stored results were
    produced with it, so it will not be changed to match `Leaves`. If you need depth
    aligned with the other per-leaf quantities, compute it yourself over
    `getalivecells(root)`.

The population-order guarantee is narrower than it looks: `allcells`, `fitness_per_cell`
and `mutations_per_cell` iterate the same `Dict` and so agree with each other *within one
unmodified population*, but the order itself is not meaningful, not sorted, and not
stable across Julia versions. Anything you want to reproduce should go through the tree.

## Distances and coalescence

```julia
pairwisedistance(node1, node2)     # drivers differing between two cells
pairwisedistances(pop)             # every pair, as a flat Vector{Int64}
pairwisedistances(pop, idx)        # restricted to allcells(pop)[idx]
pairwise_differences(pop[, idx])   # the same, as a countmap histogram
coalescence_times(pop[, idx])      # time to MRCA for every pair
coalescence_times(root[, idx]; t)
```

`pairwisedistance` climbs from both cells to their MRCA using the id ordering (an
ancestor always has a smaller id) and sums the mutations on both paths.

!!! warning "These are quadratic — subsample deliberately"
    All the pairwise functions build every one of the ``\binom{N}{2}`` pairs, each at
    ``O(\text{depth})``. At ``N = 10^4`` that is ``5 \times 10^7`` pairs. The `idx`
    argument indexes into `allcells(pop)`, which is `Dict` order — so `idx = 1:100` is
    **100 arbitrary cells, not a reproducible random sample**. For a sample you can
    record and replay, draw the tree instead and pair up its leaves:

    ```julia
    s     = sample_leaves(pop, 100; seed = UInt64(1))
    cells = getalivecells(s.root)
    dists = [pairwisedistance(cells[i], cells[j])
             for i in eachindex(cells) for j in i+1:length(cells)]
    ```

    The sampled tree retains every ancestor of a sampled cell, so these are the cells'
    true full-population distances, not distances within the sample.

The two `coalescence_times` methods differ in their default reference time: the
population method measures back from `pop.t`, the root method from `age(root)` — the
last-born leaf, which is earlier if the run ended on a death. Pass `t` explicitly to
compare across runs.

## Populations with more than one root

`initialize_population(N)` seeds `N` independent founders, so the population is a forest.
Behaviour is deliberate but varies by function:

| Function | On a forest |
|:---|:---|
| `getsingleroot(allcells(pop))` | `nothing` |
| `findMRCA(pop)` | `nothing` |
| `clonal_mutations(pop)` | `0` |
| `sitefrequencyspectrum(pop)` | correct — traverses every root |
| `mutations_per_cell(pop)`, `fitness_per_cell(pop)` | correct — per-cell, root-agnostic |
| `sample_leaves`, `sample_trees` | throw, naming the number of roots |
| `coalescence_times(pop)` | returns values, but see below |

!!! note "Coalescence times across independent founders"
    Two cells descending from different founders have no common ancestor. The
    implementation returns the time back to the seeding moment for such a pair — that is,
    it treats independent founders as if they coalesced when the population was
    initialised. That is a defensible convention (it is a lower bound on the true
    divergence), but it is a convention, and a histogram of coalescence times from a
    forest will show a spike at ``t - t_0`` that is an artefact of it, not a feature of
    the process.

    Sample each founder's tree separately if you want per-lineage coalescence.

To work root-by-root on a forest:

```julia
using AbstractTrees

for r in AbstractTrees.getroot(allcells(pop))   # the distinct roots, one per tree
    n = popsize(r)
    n >= 2 && println(n, " cells, SFS ", sitefrequencyspectrum(r, n))
end
```

`AbstractTrees.getroot` has a method here that takes a *vector* of nodes and returns the
distinct roots they belong to — which is exactly what `getsingleroot` checks for having
length one.

## Tree utilities

```julia
findMRCA(pop)                # MRCA of all living cells; nothing on a forest
findMRCA(node1, node2)       # of two cells
findMRCA(nodes)              # of a vector of cells
getsingleroot(cells)         # the unique root, or nothing
getalivecells(root)          # living leaves under a node
popsize(root)                # how many
age(root)                    # birthtime of the last-born leaf — not pop.t
leftchild!(parent, data)     # tree construction, for fixtures and tests
rightchild!(parent, data)
```

`findMRCA` on a vector reduces pairwise and returns `nothing` as soon as any pair has no
common ancestor, so it degrades correctly on a forest rather than erroring.

## Cost summary

``N`` is the number of living cells, ``D`` the typical depth, ``T`` the total number of
nodes in the tree.

| Function | Cost |
|:---|:---|
| `fitness_per_cell`, `leaf_fitness` | ``O(N)`` |
| `mutations_per_cell`, `mean_k`, `var_k` | ``O(N)`` |
| `sitefrequencyspectrum`, `branch_spectrum`, `leaf_depths` | ``O(T)`` |
| `filtered_mutations_per_cell` | ``O(T)`` |
| `celllifetimes` | ``O(T)`` |
| `findMRCA(pop)`, `getsingleroot(pop)` | ``O(T)`` |
| `pairwisedistances`, `pairwise_differences`, `coalescence_times` | ``O(N^2 D)`` |

On a single tree, `getsingleroot(pop)` is much faster than
`getsingleroot(allcells(pop))`, which has to climb from every cell.

The last row is the only one that will surprise you at scale, and it is why
[Sampling](sampling.md) exists.
