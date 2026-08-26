# DockRIFT 1.0.0 feature matrix

| Area | v1.0.0 capability |
|---|---|
| Data ingestion | CSV/TSV drag-and-drop, common alias mapping, numeric coercion, duplicate-key rejection |
| QC | Balanced-design warning, receptor×seed coverage matrix, row/ligand/receptor/seed census |
| Provenance | SHA256 source fingerprint, package/environment versions, methods capsule |
| Pair transfer | Spearman, Kendall, inversion, Δrank, Top-k survival/Jaccard, receptor>seed |
| RRL | Directional isotonic observed-support model, conservative pair propagation, censor-aware display |
| RRL sensitivity | α sweep, held-out calibration diagnostic, ligand-universe subsampling |
| Bootstrap | Ligand-identity bootstrap, finite/censored counts, finite-subset interval status |
| Ensemble view | Pair heatmaps and seed-repeatability versus receptor-transferability views |
| Seed robustness | Deterministic multi-seed convergence lens for 3/5/10+ seed datasets |
| Validity | ROC-AUC, PR-AUC, EF@k, active Top-k Jaccard, per-receptor screening landscape |
| Ligand trace | Score/rank trajectory across receptor × seed branches |
| Cross scoring | Side-by-side scoring-function pair summaries when multiple score families are loaded |
| Export | ZIP with normalized scores, pair metrics, methods, provenance, interactive report |
| Privacy | Local-only HTTP server, no telemetry or remote upload service |
| Reproducibility | Built-in VDR regression demo, automated tests, `dockrift doctor` |
