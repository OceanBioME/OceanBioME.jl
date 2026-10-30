---
name: babysit-ci
description: Monitor OceanBioME's Buildkite CI through GitHub commit statuses, fix small issues, pause on bigger problems
---

# Babysit CI

OceanBioME's tests, examples and docs run on **Buildkite** (`.buildkite/pipeline.yml`), on a
self-hosted agent, not GitHub Actions. GitHub Actions only runs CompatHelper, TagBot and doc
preview cleanup. So:

- `gh run ...` commands **do not** see the test jobs
- Status is visible on the PR/commit via `gh pr checks` / the commit status API
- Buildkite **logs are not reachable through `gh`**

## Step 1: Find the status

```sh
git branch --show-current
gh pr view --json number,url,headRefOid
gh pr checks <PR_NUMBER>
# or, for a commit:
gh api repos/OceanBioME/OceanBioME.jl/commits/<SHA>/status --jq '.statuses[] | [.context, .state, .target_url] | @tsv'
```

`gh pr checks` lists two sets of statuses:
- Named contexts from `notify:` in the pipeline: `Initialise environment`, `CPU tests`,
  `GPU tests`, `Documentation` (attached to the **deploy** step), `Clean up`
- Per-step `buildkite/oceanbiome/pr/<step>` statuses, e.g. `rowboat-cpu-unit-tests`,
  `speedboat-gpu-unit-tests`, `books-building-examples`, `docusaurus-documentation`,
  `rocket-deploy-documentation`, plus the overall `buildkite/oceanbiome/pr`

The URL column links to the Buildkite job. CPU tests take ~10–15 min, GPU tests ~25 min, so poll
no more often than every ~10 min.

## Step 2: Get the failure output

Pick the first of these that works:
1. If a Buildkite CLI (`bk`) or a `BUILDKITE_API_TOKEN` is configured, fetch the job log with it
2. Otherwise give the user the `target_url` and ask them to paste the failing section of the log
3. Reproduce locally: `/run-tests` for CPU test failures, `docs/make_examples.jl <name>` for an
   example failure

Don't guess from the status name alone.

## Step 3: Triage

### Auto-fix (commit and push without asking)

| Failure | Fix |
|---|---|
| Doctest output mismatch | Update the expected output after confirming the new output is correct |
| Typo in docstring / error message | Fix it |
| Missing `import` of an extended method, missing export, missing docstring for an export (`checkdocs = :exports`) | Add it |
| Broken `@ref` cross-reference | Fix the reference |

```sh
git add <specific files>
git commit -m "Fix <description>"
git push
```

### Retrigger (likely flaky / infrastructure)

- `Initialise environment` failure from `Pkg.instantiate`/network/registry errors
- Timeout with no test failure, agent lost, disk full on the agent
- DataDeps download failure

Retrigger via the **Retry** button on the Buildkite job (ask the user if you can't reach
Buildkite) or with an empty commit:

```sh
git commit --allow-empty -m "Retrigger CI"
git push
```

### Pause and describe (needs judgment)

- A test assertion fails (numerical mismatch, conservation failure)
- GPU-only failure ("dynamic invocation", "unsupported call", adapt errors) — likely a real bug
- Package fails to precompile/load, or the resolver fails after a compat change
- An example errors or blows up (NaNs)
- Several statuses fail with a common root cause
- Failures unrelated to the PR (upstream Oceananigans change? explain both possibilities)

Report: which context(s) failed, the key lines of output, your assessment of the cause, and fix options.

## Step 4: Confirm green

```sh
gh pr checks <PR_NUMBER>
```

## Notes

- CI runs a single Julia version (`JULIA_VERSION` in `.buildkite/pipeline.yml`) on CPU and CUDA GPU
- Docs are only deployed after `Documentation` passes; previews are cleaned up by the GitHub
  Action when the PR closes
- Never force-push or rewrite history to fix CI — always add new commits
