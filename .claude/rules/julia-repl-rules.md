# Julia Execution Rules

## Prefer MCP REPL over Bash

- When an MCP REPL server is available (check for MCP Julia/REPL tools), use it to run Julia
  code instead of launching Julia via the Bash tool — it has packages loaded and precompiled.
- Otherwise run Julia via Bash with `julia --project` from the repository root. Expect first-run
  precompilation of Oceananigans to take several minutes; use a long timeout or run in the background.
