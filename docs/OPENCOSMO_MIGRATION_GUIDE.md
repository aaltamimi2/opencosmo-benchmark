# Migrating from COSMOtherm to openCOSMO-RS: Implementation Guide

This guide explains how to migrate your existing COSMOtherm `.cosmo` files to
openCOSMO-RS (`opencosmorspy`) for activity coefficient calculations, excess
Gibbs energy, and phase diagram generation.

## What is openCOSMO-RS?

openCOSMO-RS is an open-source Python implementation of the COSMO-RS
(COnductor-like Screening MOdel for Real Solvents) thermodynamic model. It
reads Turbomole-format `.cosmo` files -- the same files produced by
COSMOtherm/COSMObase -- and calculates activity coefficients from first
principles using sigma profiles.

- **Package:** `opencosmorspy` (`pip install opencosmorspy`)
- **Input:** Turbomole `.cosmo` files (no conversion needed from COSMObase)
- **Parameterization:** `default_turbomole`
- **Output:** Activity coefficients (ln gamma), excess Gibbs energy (GE/RT)

## .cosmo File Format

These are the native output from COSMOtherm/COSMObase (BP-TZVP-COSMO level).
They contain **full segment-level data**:

- Atomic positions and elements
- Segment positions (x, y, z) on the molecular cavity surface
- Segment charges, areas, and screening potentials
- Total area, volume, and energy

**Source:** COSMObase-1901 (`BP-TZVP-COSMO/` directory). The `.cosmo` files
are not included in this repo due to licensing — they must be obtained from
a COSMObase installation.

Your existing COSMOtherm `.cosmo` files work directly with openCOSMO-RS --
**no file conversion is required**.

## Basic Usage

```python
from opencosmorspy import COSMORS
from opencosmorspy.parameterization import Parameterization
import numpy as np

# Initialize with Turbomole parameterization
par = Parameterization("default_turbomole")
cosmo = COSMORS(par)

# Load molecules directly from .cosmo files
cosmo.add_molecule(["all-cosmotherm-solvents/methanol_c0.cosmo"])
cosmo.add_molecule(["all-cosmotherm-solvents/hexane_c0.cosmo"])

# Set composition and temperature
x = np.array([0.3, 0.7])
T = 298.15  # Kelvin

# Add calculation job
cosmo.add_job(x=x, T=T, refst="pure_component")

# Run
result = cosmo.calculate()

# Extract activity coefficients
ln_gamma = result['tot']['lng'][0]  # [ln(gamma1), ln(gamma2)]
gamma1 = np.exp(ln_gamma[0])
gamma2 = np.exp(ln_gamma[1])

# Excess Gibbs energy
GE_RT = x[0] * ln_gamma[0] + x[1] * ln_gamma[1]
```

**Script reference:** `scripts/calculate_opencosmorspy_activity_coefficients.py`

## Composition Sweep

To calculate activity coefficients across a full composition range:

```python
import numpy as np
from opencosmorspy import COSMORS
from opencosmorspy.parameterization import Parameterization

COSMO_DIR = "all-cosmotherm-solvents"
COMPOSITIONS = np.linspace(0.001, 0.999, 21)  # avoid exact 0 and 1
TEMPERATURES = [298.15, 323.15, 348.15, 373.15, 398.15]  # 25-125C

cosmo1_path = f"{COSMO_DIR}/methanol_c0.cosmo"
cosmo2_path = f"{COSMO_DIR}/hexane_c0.cosmo"

results = []
for T in TEMPERATURES:
    for x1 in COMPOSITIONS:
        x = np.array([x1, 1.0 - x1])

        par = Parameterization("default_turbomole")
        cosmo = COSMORS(par)
        cosmo.add_molecule([cosmo1_path])
        cosmo.add_molecule([cosmo2_path])
        cosmo.add_job(x=x, T=T, refst="pure_component")
        result = cosmo.calculate()

        ln_gamma = result['tot']['lng'][0]
        results.append({
            'T_K': T,
            'x1': x1,
            'gamma1': np.exp(ln_gamma[0]),
            'gamma2': np.exp(ln_gamma[1]),
            'ln_gamma1': ln_gamma[0],
            'ln_gamma2': ln_gamma[1],
            'GE_RT': x1 * ln_gamma[0] + (1 - x1) * ln_gamma[1]
        })
```

**Script reference:** `scripts/calculate_incremental.py` (resumable version
with progress tracking)

## All-Pairs Binary Survey

The script `scripts/calculate_opencosmorspy_activity_coefficients.py`
calculates activity coefficients for **all 465 binary pairs** of 31 common
solvents. The solvents span:

| Category | Solvents |
|----------|----------|
| Alkanes | Pentane, Hexane, Heptane, Octane, Decane, Cyclohexane |
| Aromatics | Benzene, Toluene, Ethylbenzene |
| Alcohols | Methanol, Ethanol, Propanol, 2-Propanol, 1-Butanol, 1-Hexanol, 1-Octanol |
| Polar aprotic | Acetone, Acetonitrile, DMF, DMSO, NMP, THF |
| Halogenated | Dichloromethane, Chloroform, Carbon tetrachloride |
| Ethers | Diethyl ether, Dioxane |
| Esters | Ethyl acetate, Methyl acetate |
| Other | Acetic acid, Water |

Each pair is calculated at 5 temperatures (25, 50, 75, 100, 125 C) and 21
compositions (x1 from 0.001 to 0.999), producing **48,825 total data points**.

Results are saved as CSV files per temperature and a summary with GE/RT at
x=0.5 and infinite-dilution activity coefficients for each pair.

## Experimental Reference Data

Experimental LLE data from IUPAC-NIST is included for 6 binary solvent
systems as a reference for evaluating openCOSMO-RS predictions:

| System | Lit. UCST (K) | Data Source |
|--------|--------------|-------------|
| Methanol + n-Hexane | 306.8 | IUPAC-NIST |
| Acetonitrile + n-Hexane | 350.2 | Sinegubova 1978, IUPAC-NIST |
| Water + 1-Butanol | 393 (UCST) | IUPAC-NIST |
| Water + n-Hexane | Immiscible | IUPAC-NIST |
| Water + Ethyl acetate | No CST | IUPAC-NIST |
| Water + THF | LCST 345 K, UCST 410 K | Matous 1972 |

Stored in `data/experimental/`:
- `literature_lle_numerical.csv` -- 24 data points across 6 systems
  (columns: system, solvent1, solvent2, T_K, x1_phase1, x1_phase2, source)
- `literature_lle_data.md` -- full references and notes

## Visualization

The script `scripts/visualize_results.py` reads the activity coefficient
output and generates:

- **GE/RT heatmaps** -- 31x31 solvent matrix at each temperature, showing
  which pairs are most non-ideal (see `plots/ge_heatmap_25C_complete.png`)
- **Phase diagrams** -- 3-panel plots (gamma, ln gamma, GE/RT vs composition)
  for 6 key systems (Methanol+Hexane, Acetonitrile+Hexane, Water+1-Butanol,
  Water+Hexane, Water+Ethyl acetate, Water+THF)
- **Literature comparison** -- GE/RT curves overlaid with experimental phase
  boundaries from IUPAC-NIST (see `plots/literature_lle_comparison.png`)
- **Statistics summary** -- text file listing most/least ideal pairs

## Quick-Start Migration Checklist

1. **Install opencosmorspy:**
   ```bash
   pip install opencosmorspy
   ```

2. **Locate your `.cosmo` files** from COSMOtherm (typically in
   `COSMObase-1901/BP-TZVP-COSMO/`). These work directly -- no conversion.

3. **Initialize the model:**
   ```python
   from opencosmorspy import COSMORS
   from opencosmorspy.parameterization import Parameterization

   par = Parameterization("default_turbomole")
   cosmo = COSMORS(par)
   ```

4. **Add molecules** (one `add_molecule` call per component):
   ```python
   cosmo.add_molecule(["path/to/solvent1_c0.cosmo"])
   cosmo.add_molecule(["path/to/solvent2_c0.cosmo"])
   ```

5. **Add a job and calculate:**
   ```python
   cosmo.add_job(x=np.array([x1, x2]), T=T_kelvin, refst="pure_component")
   result = cosmo.calculate()
   ln_gamma = result['tot']['lng'][0]
   ```

6. **Extract results:**
   - `ln_gamma[0]` = ln(gamma) for molecule 1
   - `ln_gamma[1]` = ln(gamma) for molecule 2
   - GE/RT = sum(x_i * ln_gamma_i)

## Key Differences: COSMOtherm vs openCOSMO-RS

| Aspect | COSMOtherm | openCOSMO-RS |
|--------|-----------|--------------|
| License | Commercial | Open-source |
| Package | Desktop application | `pip install opencosmorspy` |
| Parameterization | BP-TZVP-COSMO (proprietary fit) | `default_turbomole` |
| Input format | `.cosmo` (native) | `.cosmo` (same files, no conversion) |
| Output | GUI + text files | Python dict (`result['tot']['lng']`) |
| Reference state | Multiple options | `refst="pure_component"` |
| Batch support | Built-in | Loop over compositions/temperatures |

## Known Gotchas and Limitations

1. **Reinitializing per calculation:** The scripts in this repo create a fresh
   `COSMORS` object for every composition point. This works correctly but is
   slow. For production use, batch multiple compositions if the API supports it.

2. **Avoid exact 0 and 1 compositions:** Use `np.linspace(0.001, 0.999, N)`
   rather than `np.linspace(0, 1, N)`. Pure-component endpoints can cause
   numerical issues in the activity coefficient calculation.

3. **File naming convention:** COSMObase files typically use the pattern
   `<name>_c0.cosmo` (conformer 0). If your molecule has multiple conformers
   (`_c0`, `_c1`, etc.), you can pass multiple paths to `add_molecule`:
   ```python
   cosmo.add_molecule(["mol_c0.cosmo", "mol_c1.cosmo"])
   ```

4. **Result structure:** The result dict has nested structure. Activity
   coefficients are at `result['tot']['lng'][0]` where `[0]` indexes the
   first (and usually only) job. The `'tot'` key gives total ln(gamma);
   individual contributions are under `'comb'`, `'res'`, etc.

5. **COSMOtherm vs openCOSMO-RS agreement:** For small molecules, agreement
   is generally good. A comparison for cyclohexanone showed exact surface area
   match (138.27 A^2). However, the parameterization differs from COSMOtherm's
   proprietary fit, so expect quantitative differences -- trends and rankings
   should be consistent.

6. **Surface area convention:** openCOSMO-RS uses the **sum of segment areas**
   (not the cavity area from the `.cosmo` file header) as the effective
   molecular surface area. This matches COSMOtherm's convention.

7. **Temperature units:** All temperatures must be in **Kelvin**. Common
   values: 298.15 K (25 C), 323.15 K (50 C), 348.15 K (75 C).
