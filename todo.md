# TODO — split of analysis into `gITH-nonMarkovian`

Context: `analysis/`, `data/`, `plots/`, `theory/` were moved out of this repo into
`/Users/alexanderstein/Documents/GitHub/gITH-nonMarkovian`. This package now only
contains the simulation code (`src/`, `test/`, `examples/`).

## This repo (`MutationLoadDynamics.jl`)

- [ ] Run the full test suite (`julia --project=. -e 'using Pkg; Pkg.test()'`) to confirm the CairoMakie removal did not break anything downstream.
- [ ] Run each `examples/*.jl` file end-to-end to verify no example silently depended on CairoMakie or on files under the moved folders.
- [ ] Review the `README.md` — the "Analysis functions" section still describes in-package utilities in `src/`, but double-check none of the code snippets reference paths that used to live in `analysis/`, `data/`, or `plots/`.
- [ ] Commit the cleanup: `.gitignore`, `Project.toml`, `Manifest.toml` (message suggestion: "Split analysis into gITH-nonMarkovian; drop CairoMakie dep").
- [ ] Decide whether `Manifest.toml` should still be committed for a package (currently it is — the `# Manifest.toml` line in `.gitignore` is commented out). Standard practice for a registered package is to *not* commit it; for a research/reproducibility repo, committing is fine. Leave as-is unless you change your mind.

## New repo (`gITH-nonMarkovian`)

- [ ] Write a `.gitignore` **before** the first `git add`. At minimum ignore `data/` — it is 3.1 GB and will not fit on GitHub without LFS. Consider also ignoring large files under `plots/` (2.6 MB total, probably fine but check for individual large PDFs/PNGs).
- [ ] Decide on data handling:
  - Option A: keep `data/` fully local and gitignored (simplest).
  - Option B: use Git LFS for a curated subset (`git lfs track "data/**/*.csv"` etc.).
  - Option C: store raw data elsewhere (Zenodo, institutional storage) and keep only pointers/README under `data/`.
- [ ] Create a proper `Project.toml` for the analysis environment. It should include the plotting/analysis deps that this repo actually uses — at minimum `CairoMakie` (removed from the simulation package) plus whatever the scripts import. Add `MutationLoadDynamics` as a dependency, either via `Pkg.develop(path="../MutationLoadDynamics.jl")` for local work or via the GitHub URL once published.
- [ ] Flesh out the `README.md` (currently 19 bytes). Describe: purpose of the repo, how it relates to `MutationLoadDynamics.jl`, how to reproduce figures/analyses, where the data lives.
- [ ] Walk through `analysis/**/*.jl` (and `theory/`) and fix any paths that assumed the old layout — e.g. `../data/…`, `../src/…`, `include("../src/...")`, or `using MutationLoadDynamics` from a relative path. The simulation package is no longer a sibling folder in this repo.
- [ ] Confirm the LICENSE matches what you want. It was copied from a fresh `git init` so it may not match `MutationLoadDynamics.jl`'s LICENSE — align them if they should be consistent.
- [ ] First commit + push to GitHub. Before pushing: `du -sh .git` to make sure no oversized blobs slipped in via an earlier `git add`.

## Cross-repo hygiene

- [ ] Add a short paragraph to `MutationLoadDynamics.jl/README.md` pointing at `gITH-nonMarkovian` for the analyses/figures/theory that used to live here (mirrors the existing pointer to `BirthDeathMutation.jl` at the bottom).
- [ ] Add the reciprocal pointer in `gITH-nonMarkovian/README.md` back to the simulation package.
