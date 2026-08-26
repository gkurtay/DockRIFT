# DockRIFT v1.0.0

DockRIFT v1.0.0 is the first public release of the Docking Rank Instability & Fragility Toolkit.

## Scientific scope

DockRIFT quantifies receptor-conditioned ligand-ranking stability while keeping seed repeatability, native-pose competence, receptor-to-receptor transfer, and retrospective activity-label performance conceptually separate. It includes directional observed-support RRL analysis and retains right-censoring whenever the selected reversal criterion is unsupported by the observed score-gap range.

## Included interfaces

- reusable Python numerical core;
- `dockrift` command-line interface;
- custom local-first DockRIFT Studio browser workbench;
- reproducibility/provenance export;
- packaged VDR reference data and regression tests.

## Verification

- 13 automated tests pass;
- VDR Vina reference: receptor Spearman ρ = 0.9585856068, conservative RRL = 0.875;
- VDR Vinardo reference: receptor Spearman ρ = 0.9517698269, conservative RRL = 0.884;
- API smoke and `dockrift doctor` records are included in `docs/release_v1.0.0/`;
- release artifacts are accompanied by SHA256 checksums.

## License

MIT License.
