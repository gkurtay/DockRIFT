# DockRIFT 1.0.0

**Docking Rank Instability & Fragility Toolkit**

DockRIFT is an open-source Python toolkit and local-first visual workbench for quantifying how molecular-docking ligand rankings change across receptor conformations. It separates three validation questions that are often conflated:

1. **within-receptor stochastic repeatability** across seeds;
2. **native-pose competence** of receptor/preparation workflows; and
3. **receptor-to-receptor rank transferability** for a frozen ligand universe.

DockRIFT also implements the **observed-support Rank Resolution Limit (RRL)**: the smallest observed score separation at which a monotone estimate of receptor-induced pairwise reversal probability reaches a chosen criterion. If the criterion is not reached within the observed score-gap range, DockRIFT reports a right-censored lower bound rather than extrapolating.

## Main capabilities

- CSV/TSV ingestion with alias normalization, numeric coercion, duplicate-key checks, balanced-design diagnostics, and receptor×seed coverage census.
- SHA256 source fingerprinting, methods capsules, and reproducibility exports.
- Spearman and Kendall rank transfer, exact pairwise inversion, absolute rank displacement, Top-k survival/Jaccard, and receptor-versus-seed attribution.
- Directional, observed-support RRL with conservative pair propagation and explicit censoring.
- RRL α sensitivity, held-out calibration diagnostics, ligand-universe subsampling, and ligand-identity bootstrap.
- Receptor-ensemble matrices, seed-convergence views, ligand-level traces, and Vina/Vinardo cross-scoring comparison.
- Optional activity-label analysis with ROC-AUC, PR-AUC, EF@k, and active Top-k stability.
- A custom **DockRIFT Studio** browser interface served locally from the same numerical core as the CLI.
- Reproducibility ZIP export containing normalized scores, pair metrics, methods, provenance, and a self-contained interactive report.

## Installation

DockRIFT requires Python 3.10 or newer.

### From a GitHub release wheel

```bash
python -m venv .venv
# Windows PowerShell: .venv\Scripts\Activate.ps1
# macOS/Linux: source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install dockrift-1.0.0-py3-none-any.whl
```

### From source

```bash
git clone https://github.com/gkurtay/DockRIFT.git
cd DockRIFT
python -m pip install .
```

For development and tests:

```bash
python -m pip install -e ".[test]"
python -m pytest
```

## Quick start

Check the installation:

```bash
dockrift --version
dockrift doctor
```

Launch the local Studio interface:

```bash
dockrift studio
```

or:

```bash
dockrift-studio
```

The default address is `http://127.0.0.1:8765/`. Uploaded score tables remain in the local Python process; DockRIFT does not require a cloud database, telemetry endpoint, or remote analysis service.

## Input contract

DockRIFT expects a long-form CSV/TSV with one row per ligand × receptor × seed:

- `ligand_id` — aliases include `dockrift_id`, `ligand`, `compound_id`, and `id`;
- `receptor` — aliases include `pdb_id`, `pdb`, `structure`, and `conformation`;
- `seed`;
- `score` — aliases include `mode1_score`, `vina_score`, and `affinity`;
- optional `scoring_function`;
- optional `class` with values such as `active` and `inactive`.

By default, lower numerical scores are treated as better rankings.

## CLI examples

Analyze a balanced score table:

```bash
dockrift analyze scores.csv --out analysis_output --alpha 0.05
```

Extract mode-1 scores from a task manifest with accessible PDBQT paths:

```bash
dockrift extract task_manifest.csv --out normalized_scores.csv
```

Freeze and verify exact inputs:

```bash
dockrift freeze scores.csv analysis.py config.conf --out analysis.sha256
dockrift verify analysis.sha256
```

## Scientific interpretation

DockRIFT measures **rank stability and transferability**, not experimental binding correctness. Activity labels form a separate practical-validity layer.

RRL values:

- are conditional on the receptor direction, scoring model, docking protocol, and ligand universe;
- are evaluated only on observed score-gap support;
- remain censored if the selected reversal criterion is not reached;
- use the numerical units of the supplied docking score; and
- are **not binding-free-energy error bars**.

A ligand ordering separated by less than a finite RRL should be interpreted as unresolved for that receptor transition under the analyzed protocol and chemical universe. A censored RRL is a lower bound, not a finite estimate.

## Verification

The v1.0.0 release includes:

- a 13-test automated suite;
- a packaged 400-ligand VDR Vina/Vinardo reference fixture;
- an API smoke workflow reproducing receptor Spearman ρ = 0.958586 and a finite Vina RRL of 0.875;
- packaged demo data and Studio assets; and
- SHA256 manifests for release artifacts.

Software verification demonstrates implementation consistency; it is not independent evidence for the scientific conclusions of the accompanying study.

## Citation

Please cite the archived DockRIFT v1.0.0 software record and the accompanying manuscript. GitHub renders the software citation metadata from [`CITATION.cff`](CITATION.cff). The version-specific Zenodo DOI is created from the tagged GitHub release.

## License

DockRIFT is released under the [MIT License](LICENSE).
