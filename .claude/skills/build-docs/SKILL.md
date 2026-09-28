---
name: build-docs
description: Build OceanBioME documentation locally, with or without running the Literate examples
---

# Build Documentation

## Steps

1. Ask the user if they want a **full build** (runs all examples; slow, `eady`/`kelp`/`oae_experiment`
   are heavy) or a **fast build** (docs only, optionally a subset of examples)
2. Set up the docs environment if needed:
   ```sh
   julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
   ```
3. **Full build**:
   ```sh
   julia docs/prepare_examples.jl                           # runs every example via make_examples.jl
   JULIA_DEBUG=Documenter julia --project=docs/ docs/make.jl
   ```
4. **Fast build**:
   - Build only the examples you need: `julia --project=docs/ docs/make_examples.jl box`
   - Temporarily remove the unbuilt entries from `examples` in `docs/make.jl`
     (pages point at `generated/<name>.md`, which only exists once built)
   - If needed, temporarily add `warnonly = [:cross_references, :missing_docs]` and/or
     `draft = true` to `makedocs`
   - Run `julia --project=docs/ docs/make.jl`
   - **Revert all changes** to `docs/make.jl` afterwards (`git diff docs/make.jl` should be empty)
5. Preview: open `docs/build/index.html`, or
   ```julia
   using LiveServer
   serve(dir="docs/build")
   ```
6. Report warnings and errors from the build

## Notes

- `docs/make.jl` has `doctest = false`; doctests are checked by `Pkg.test()` instead
- `docs/make.jl` also generates parameter pages from `docs/display_parameters.jl`
- Build outputs (`docs/build`, `docs/src/generated`, `*.jld2`, `*.png`, ...) are gitignored
- On CI (Buildkite) examples are built in a matrix on the GPU agent, then docs, then deploy
