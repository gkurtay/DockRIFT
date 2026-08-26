# Contributing to DockRIFT

Contributions that improve numerical correctness, documentation, interoperability, tests, or scientifically justified diagnostics are welcome.

1. Open an issue describing the proposed change.
2. Create a focused branch and add tests for changes to numerical behavior.
3. Run `python -m pytest` before submitting a pull request.
4. Do not change the observed-support/no-extrapolation/censoring semantics of RRL without an explicit methodological proposal and regression fixtures.
5. Keep GUI presentation separate from the reusable numerical core.

Bug reports should include the DockRIFT version, Python/platform information from `dockrift doctor`, a minimal input example when shareable, and the exact error message.
