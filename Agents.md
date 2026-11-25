
# Project instructions for Codex

This repository simulates scalar-field evolution around a black-hole binary.

- Prioritize physical correctness and numerical reproducibility.
- Read the existing implementation before modifying scientific code.
- Preserve existing equations, discretization, update order and diagnostics unless explicitly instructed otherwise.
- Keep circular inspiral behavior compatible with existing results.
- Fixed eccentric orbits and circular inspiral are distinct orbital models.
- Do not implement eccentric inspiral unless explicitly requested.
- Separate structural refactoring from numerical-method changes.
- Reuse existing sparse-matrix, Crank–Nicolson and checkpoint implementations.
- Prefer minimal changes and avoid unnecessary abstractions.
- For scientific changes, run relevant regression tests and compare against the pre-change baseline.
- Report numerical differences, even if tests pass.
- Keep implementation summaries concise.
