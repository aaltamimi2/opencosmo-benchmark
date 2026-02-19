# opencosmo-LLE

Systematic comparison of **openCOSMO-RS** (`opencosmorspy`, open-source) vs **COSMOtherm** (commercial, BP_TZVP_19) activity coefficient predictions for 465 binary solvent pairs at five temperatures (25–125 °C).

## Key Results

- **MAE(GE/RT) = 0.048** across 435 pairs (excluding chloroform, a known H-bond donor outlier)
- **R² = 0.931** for parity between openCOSMO-RS and COSMOtherm predictions
- 85% of pairs agree within GE/RT = 0.1; agreement improves at elevated temperatures
- COSMOtherm finds LLE for 65 pairs; openCOSMO-RS captures the same qualitative trends
- See `results/comparison/paper_summary.png` for a 4-panel overview

## Directory Structure

```
opencosmo-LLE/
├── scripts/                         # All Python scripts
│   ├── calculate_opencosmorspy_activity_coefficients.py  # Main openCOSMO-RS calculation
│   ├── calculate_incremental.py     # Resumable version with progress tracking
│   ├── compare_cosmotherm_vs_opencosmo.py  # Comparison & 10 visualization plots
│   └── visualize_results.py         # GE/RT heatmaps, phase diagrams, lit. comparison
├── data/
│   ├── calculated/opencosmo_activity_coefficients/  # openCOSMO-RS output (48,825 points)
│   └── experimental/                # IUPAC-NIST literature LLE data (24 points, 6 systems)
├── cosmotherm-reference-data/       # COSMOtherm .inp input and .tab output
├── results/comparison/              # Comparison statistics, plots, and analysis
├── plots/                           # openCOSMO-RS standalone visualizations
└── docs/                            # Migration guide
```

## Requirements

```
opencosmorspy
numpy
pandas
scipy
matplotlib
```

COSMOtherm `.cosmo` files (from COSMObase, not included due to licensing) are required to rerun the openCOSMO-RS calculations. The pre-computed results in `data/calculated/` can be used directly for comparison and visualization.

## Reproducing

1. **openCOSMO-RS calculations** (requires `.cosmo` files):
   ```bash
   python scripts/calculate_opencosmorspy_activity_coefficients.py
   ```

2. **COSMOtherm vs openCOSMO-RS comparison** (requires COSMOtherm `.tab` output):
   ```bash
   python scripts/compare_cosmotherm_vs_opencosmo.py
   ```

3. **Standalone visualizations** (openCOSMO-RS data only):
   ```bash
   python scripts/visualize_results.py
   ```

## Solvents (31)

Alkanes (6), aromatics (3), alcohols (7), polar aprotics (5), halogenated (3), ethers (2), esters (2), acetic acid, water.

Full list and .cosmo file mapping in `scripts/calculate_opencosmorspy_activity_coefficients.py`.

## Notes

- **Chloroform outlier**: openCOSMO-RS `default_turbomole` underestimates the C-H hydrogen bond donor strength of chloroform. See `results/comparison/chloroform_outlier_analysis.md`.
- **Upgrade path**: The `openCOSMORS24a` parameterization (2024) may resolve this but requires ORCA `.orcacosmo` files. See `results/comparison/openCOSMORS24a_upgrade_path.md`.
- **Migration guide**: `docs/OPENCOSMO_MIGRATION_GUIDE.md` covers the COSMOtherm → openCOSMO-RS workflow.
