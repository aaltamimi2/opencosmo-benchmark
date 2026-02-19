# Scripts

All scripts use openCOSMO-RS (`opencosmorspy`) or standard Python libraries (pandas, matplotlib). No cCOSMO or cosmopharm dependencies.

## Calculation

- **`calculate_opencosmorspy_activity_coefficients.py`** -- Main calculation script. Computes activity coefficients (gamma1, gamma2), ln(gamma), and GE/RT for all 465 binary pairs of 31 common solvents. Sweeps 21 compositions x 5 temperatures (25-125C). Also calculates infinite-dilution activity coefficients. Outputs per-temperature CSVs, a summary CSV, and metadata JSON.

- **`calculate_incremental.py`** -- Same calculation as above but with resume capability. Saves results after each pair and tracks progress in a JSON file, so interrupted runs pick up where they left off.

## Comparison

- **`compare_cosmotherm_vs_opencosmo.py`** -- Parses COSMOtherm `.tab` output and compares against openCOSMO-RS CSV data. Interpolates both onto a common composition grid using cubic splines, computes per-pair MAE/RMSE/bias for GE/RT and ln(gamma), and generates 10 publication-quality plots (parity, heatmap, class breakdown, solvent ranking, temperature dependence, LLE summary, 4-panel paper figure). Supports excluding specific solvents (e.g. chloroform).

## Visualization

- **`visualize_results.py`** -- Reads output from the calculation scripts and generates: GE/RT heatmaps (31x31 solvent matrix at each temperature), 3-panel phase diagrams (gamma, ln(gamma), GE/RT vs composition) for 6 key systems, a literature comparison figure, and a text summary of most/least ideal pairs.
