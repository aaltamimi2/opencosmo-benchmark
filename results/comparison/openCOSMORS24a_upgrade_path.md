# openCOSMORS24a: Potential Upgrade Path

## Current Setup

This project uses `Parameterization("default_turbomole")` from opencosmorspy
with Turbomole `.cosmo` files from COSMObase-1901. This parameterization
dates to Eckert & Klamt (2002) and shows systematic deviations from
COSMOtherm BP_TZVP_19, particularly for chloroform H-bond donor
interactions (see `chloroform_outlier_analysis.md`).

## openCOSMORS24a

The `openCOSMORS24a` parameterization was published in 2024 and is the most
accurate open-source COSMO-RS parameterization available. It was refit to
match COSMOtherm 24 accuracy for solvation free energies and activity
coefficients.

**Key requirement:** It is parameterized for **ORCA 6.0** DFT calculations
(BP86/def2-TZVP and def2-TZVPD), so it requires `.orcacosmo` files — not
the Turbomole `.cosmo` files from COSMObase that we currently have.

## What Would Be Needed

1. Install ORCA 6.0 (free for academic use)
2. Re-run DFT geometry optimization + COSMO surface calculation on all 31
   solvents at the BP86/def2-TZVP level using ORCA
3. Use the resulting `.orcacosmo` files with `openCOSMORS24a()`
4. Re-run all 465 binary pair calculations at 5 temperatures

## Usage

```python
from opencosmorspy.parameterization import openCOSMORS24a
from opencosmorspy import COSMORS

par = openCOSMORS24a()
cosmo = COSMORS(par)
cosmo.add_molecule(["solvent1.orcacosmo"])
cosmo.add_molecule(["solvent2.orcacosmo"])
```

Note: `openCOSMORS24a` is a class, not a string — it must be instantiated
directly rather than passed as a name to `Parameterization()`.

## References

- [openCOSMO-RS_py GitHub](https://github.com/TUHH-TVT/openCOSMO-RS_py)
- [Predicting solvation free energies with openCOSMO-RS (arXiv 2407.03434)](https://arxiv.org/html/2407.03434v1)
