# Chloroform Outlier Analysis: COSMOtherm vs openCOSMO-RS

## Finding

Chloroform is the dominant source of disagreement between COSMOtherm
(BP_TZVP_19 parameterization) and openCOSMO-RS (`default_turbomole`
parameterization) across all 465 binary solvent pairs at 5 temperatures.

| Subset | MAE(GE/RT) | Median | Max | % pairs < 0.1 |
|--------|-----------|--------|-----|---------------|
| All 465 pairs | 0.079 | 0.039 | 1.033 | 81% |
| Excluding CHCl3 (435 pairs) | 0.059 | 0.035 | 0.462 | 85% |
| Excluding all halogenated (378) | 0.049 | 0.028 | 0.458 | 89% |
| CHCl3 pairs only (30) | 0.361 | 0.451 | 1.033 | 10% |

## Evidence

1. **Every chloroform pair has positive bias** (COSMOtherm predicts more
   negative GE/RT than openCOSMO-RS), meaning COSMOtherm captures stronger
   favorable cross-interactions that openCOSMO-RS misses.

2. **Error correlates with H-bond acceptor strength of the partner:**
   - DMSO (strong acceptor): MAE = 1.03
   - NMP, DMF: MAE = 0.81-0.90
   - Ethers (THF, dioxane, diethyl ether): MAE = 0.56-0.70
   - Alcohols: MAE = 0.42-0.48
   - Alkanes (no acceptor): MAE = 0.08-0.14
   - CCl4 (no H-bond at all): MAE = 0.001

3. **Not a generic halogenated-solvent problem.** CCl4 has no C-H bond and
   sits right at the overall average (MAE = 0.079). DCM (CH2Cl2) is
   intermediate (MAE = 0.166), consistent with its weaker H-bond donor
   ability relative to CHCl3.

## Postulated Cause

Chloroform (CHCl3) is an unusually strong **hydrogen-bond donor** due to
the three electron-withdrawing chlorine atoms that activate its single C-H
bond. In COSMO-RS, hydrogen-bonding interactions are captured through the
**misfit energy** between surface segments: a positive sigma segment (H-bond
donor, e.g., the C-H of chloroform) interacts favorably with a negative
sigma segment (H-bond acceptor, e.g., the S=O of DMSO).

The magnitude of this interaction depends on two parameterization choices:

1. **The hydrogen-bond interaction threshold and strength constant** — these
   determine when two surface segments are considered to form a hydrogen
   bond and how much stabilization energy that bond provides.

2. **The sigma profile of chloroform itself** — while both implementations
   read the same `.cosmo` file (identical quantum-chemical surface), the
   binning, averaging, and effective-area weighting of segments differs
   between COSMOtherm's proprietary parameterization and the open-source
   `default_turbomole` parameters.

COSMOtherm's BP_TZVP_19 parameterization was fitted to a large experimental
database that likely includes chloroform-containing systems. The stronger
H-bond donor treatment in COSMOtherm produces more negative GE/RT (more
favorable mixing) for CHCl3 + acceptor pairs, which matches the known
experimental behavior (e.g., chloroform + acetone forms a negative-deviation
azeotrope with GE < 0).

The `default_turbomole` parameterization in openCOSMO-RS, by contrast,
**underestimates** this C-H donor interaction, leading to systematically
less negative GE/RT predictions. The error is proportional to the H-bond
acceptor strength of the partner, producing the observed monotonic trend.

## Implication

For any application using openCOSMO-RS that involves chloroform (or other
C-H hydrogen-bond donors like dichloromethane), users should be aware of
a systematic underestimation of favorable cross-interactions. For the
remaining 30 solvents without strong C-H donor character, the two
implementations agree well (MAE = 0.059, 85% of pairs within 0.1).
