#!/usr/bin/env python3
"""
Binary Solvent Phase Diagram Calculations using OpenCOSMO-RS

Uses opencosmorspy with original .cosmo files from COSMObase to calculate
activity coefficients and excess Gibbs energy for all binary pairs of ~30
common solvents.

This is the same method as the openCOSMO blue line in comparison plots.

Output:
- Activity coefficients (gamma1, gamma2)
- Excess Gibbs energy (GE/RT)
- Full composition range (0 to 1 in 21 points)
- Temperatures: 25, 50, 75, 100, 125°C
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations
import warnings
import json
from datetime import datetime

warnings.filterwarnings('ignore')

from opencosmorspy import COSMORS
from opencosmorspy.parameterization import Parameterization

# ============================================================================
# Configuration
# ============================================================================

OUTPUT_DIR = Path(__file__).parent / "results"
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Path to original COSMObase .cosmo files
COSMO_DIR = Path(__file__).parent.parent / "all-cosmotherm-solvents"

# Temperature range: 25, 50, 75, 100, 125°C in Kelvin
TEMPERATURES = [298.15, 323.15, 348.15, 373.15, 398.15]
TEMP_LABELS = ['25C', '50C', '75C', '100C', '125C']

# Composition points (21 points from 0 to 1)
N_COMPOSITIONS = 21
COMPOSITIONS = np.linspace(0.001, 0.999, N_COMPOSITIONS)  # Avoid exact 0 and 1

# ============================================================================
# Common Solvents Definition (~32 solvents)
# ============================================================================

# Format: name, cosmo_file (with _c0.cosmo), MW (g/mol), V298 (cm3/mol)
COMMON_SOLVENTS = [
    # Alkanes
    {'name': 'Pentane', 'cosmo': 'pentane_c0.cosmo', 'MW': 72.15, 'V298': 116.1},
    {'name': 'Hexane', 'cosmo': 'hexane_c0.cosmo', 'MW': 86.18, 'V298': 131.6},
    {'name': 'Heptane', 'cosmo': 'n-heptane_c0.cosmo', 'MW': 100.20, 'V298': 147.4},
    {'name': 'Octane', 'cosmo': 'octane_c0.cosmo', 'MW': 114.23, 'V298': 163.5},
    {'name': 'Decane', 'cosmo': 'n-decane_c0.cosmo', 'MW': 142.28, 'V298': 195.9},
    {'name': 'Cyclohexane', 'cosmo': 'cyclohexane_c0.cosmo', 'MW': 84.16, 'V298': 108.7},

    # Aromatics
    {'name': 'Benzene', 'cosmo': 'benzene_c0.cosmo', 'MW': 78.11, 'V298': 89.4},
    {'name': 'Toluene', 'cosmo': 'toluene_c0.cosmo', 'MW': 92.14, 'V298': 106.8},
    {'name': 'Ethylbenzene', 'cosmo': 'ethylbenzene_c0.cosmo', 'MW': 106.17, 'V298': 123.1},

    # Alcohols
    {'name': 'Methanol', 'cosmo': 'methanol_c0.cosmo', 'MW': 32.04, 'V298': 40.7},
    {'name': 'Ethanol', 'cosmo': 'ethanol_c0.cosmo', 'MW': 46.07, 'V298': 58.5},
    {'name': 'Propanol', 'cosmo': 'propanol_c0.cosmo', 'MW': 60.10, 'V298': 75.2},
    {'name': '2-Propanol', 'cosmo': '2-propanol_c0.cosmo', 'MW': 60.10, 'V298': 76.8},
    {'name': '1-Butanol', 'cosmo': '1-butanol_c0.cosmo', 'MW': 74.12, 'V298': 91.5},
    {'name': '1-Hexanol', 'cosmo': '1-hexanol_c0.cosmo', 'MW': 102.17, 'V298': 125.2},
    {'name': '1-Octanol', 'cosmo': '1-octanol_c0.cosmo', 'MW': 130.23, 'V298': 158.4},

    # Polar aprotic
    {'name': 'Acetone', 'cosmo': 'propanone_c0.cosmo', 'MW': 58.08, 'V298': 74.0},
    {'name': 'Acetonitrile', 'cosmo': 'acetonitrile_c0.cosmo', 'MW': 41.05, 'V298': 52.6},
    {'name': 'DMF', 'cosmo': 'dimethylformamide_c0.cosmo', 'MW': 73.09, 'V298': 77.0},
    {'name': 'DMSO', 'cosmo': 'dimethylsulfoxide_c0.cosmo', 'MW': 78.13, 'V298': 71.3},
    {'name': 'NMP', 'cosmo': 'n-methyl-2-pyrrolidinone_c0.cosmo', 'MW': 99.13, 'V298': 96.5},
    {'name': 'THF', 'cosmo': 'thf_c0.cosmo', 'MW': 72.11, 'V298': 81.7},

    # Halogenated
    {'name': 'Dichloromethane', 'cosmo': 'ch2cl2_c0.cosmo', 'MW': 84.93, 'V298': 63.9},
    {'name': 'Chloroform', 'cosmo': 'chcl3_c0.cosmo', 'MW': 119.38, 'V298': 80.7},
    {'name': 'Carbon tetrachloride', 'cosmo': 'ccl4_c0.cosmo', 'MW': 153.82, 'V298': 97.1},

    # Ethers
    {'name': 'Diethyl ether', 'cosmo': 'diethylether_c0.cosmo', 'MW': 74.12, 'V298': 104.8},
    {'name': 'Dioxane', 'cosmo': 'dioxane_c0.cosmo', 'MW': 88.11, 'V298': 85.7},

    # Esters
    {'name': 'Ethyl acetate', 'cosmo': 'ethylacetate_c0.cosmo', 'MW': 88.11, 'V298': 98.5},
    {'name': 'Methyl acetate', 'cosmo': 'methylacetate_c0.cosmo', 'MW': 74.08, 'V298': 79.7},

    # Other
    {'name': 'Acetic acid', 'cosmo': 'aceticacid_c0.cosmo', 'MW': 60.05, 'V298': 57.6},
    {'name': 'Water', 'cosmo': 'h2o_c0.cosmo', 'MW': 18.02, 'V298': 18.0},
]


# ============================================================================
# Helper Functions
# ============================================================================

def validate_solvents():
    """Check which solvent files exist."""
    valid_solvents = []
    missing = []

    for solvent in COMMON_SOLVENTS:
        cosmo_path = COSMO_DIR / solvent['cosmo']
        if cosmo_path.exists():
            valid_solvents.append(solvent)
        else:
            missing.append(solvent['name'])

    if missing:
        print(f"Warning: Missing COSMO files for: {missing}")

    return valid_solvents


def calculate_pair_activity(solvent1, solvent2, T, compositions):
    """
    Calculate activity coefficients for a binary pair using openCOSMO-RS.

    Returns:
        list of dicts with x1, gamma1, gamma2, GE_RT for each composition
    """
    results = []

    cosmo1_path = str(COSMO_DIR / solvent1['cosmo'])
    cosmo2_path = str(COSMO_DIR / solvent2['cosmo'])

    for x1 in compositions:
        x2 = 1.0 - x1
        x = np.array([x1, x2])

        try:
            # Initialize fresh COSMO-RS for each calculation
            par = Parameterization("default_turbomole")
            cosmo = COSMORS(par)

            # Add molecules
            cosmo.add_molecule([cosmo1_path])
            cosmo.add_molecule([cosmo2_path])

            # Add job
            cosmo.add_job(x=x, T=T, refst="pure_component")

            # Calculate - returns results directly
            result = cosmo.calculate()

            # Get activity coefficients from result dict
            # 'tot' contains total ln(gamma), 'lng' is the array
            ln_gamma = result['tot']['lng'][0]  # First job result
            gamma1 = np.exp(ln_gamma[0])
            gamma2 = np.exp(ln_gamma[1])

            # Calculate GE/RT = sum(x_i * ln(gamma_i))
            GE_RT = x1 * ln_gamma[0] + x2 * ln_gamma[1]

            results.append({
                'x1': x1,
                'x2': x2,
                'gamma1': gamma1,
                'gamma2': gamma2,
                'ln_gamma1': ln_gamma[0],
                'ln_gamma2': ln_gamma[1],
                'GE_RT': GE_RT
            })

        except Exception as e:
            results.append({
                'x1': x1,
                'x2': x2,
                'gamma1': np.nan,
                'gamma2': np.nan,
                'ln_gamma1': np.nan,
                'ln_gamma2': np.nan,
                'GE_RT': np.nan
            })

    return results


def calculate_infinite_dilution(solvent1, solvent2, T):
    """
    Calculate infinite dilution activity coefficients.

    Returns gamma1_inf (solvent1 in solvent2) and gamma2_inf (solvent2 in solvent1)
    """
    cosmo1_path = str(COSMO_DIR / solvent1['cosmo'])
    cosmo2_path = str(COSMO_DIR / solvent2['cosmo'])

    gamma1_inf, gamma2_inf = np.nan, np.nan

    # Solvent1 at infinite dilution in solvent2
    try:
        par = Parameterization("default_turbomole")
        cosmo = COSMORS(par)
        cosmo.add_molecule([cosmo1_path])
        cosmo.add_molecule([cosmo2_path])
        cosmo.add_job(x=np.array([0.001, 0.999]), T=T, refst="pure_component")
        result = cosmo.calculate()
        gamma1_inf = np.exp(result['tot']['lng'][0][0])
    except:
        pass

    # Solvent2 at infinite dilution in solvent1
    try:
        par = Parameterization("default_turbomole")
        cosmo = COSMORS(par)
        cosmo.add_molecule([cosmo1_path])
        cosmo.add_molecule([cosmo2_path])
        cosmo.add_job(x=np.array([0.999, 0.001]), T=T, refst="pure_component")
        result = cosmo.calculate()
        gamma2_inf = np.exp(result['tot']['lng'][0][1])
    except:
        pass

    return gamma1_inf, gamma2_inf


# ============================================================================
# Main Calculation
# ============================================================================

def main():
    print("=" * 70)
    print("BINARY SOLVENT PHASE DIAGRAMS - OpenCOSMO-RS")
    print("=" * 70)
    print()
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"COSMO files directory: {COSMO_DIR}")
    print()

    # Validate solvents
    valid_solvents = validate_solvents()
    n_solvents = len(valid_solvents)
    n_pairs = n_solvents * (n_solvents - 1) // 2

    print(f"Valid solvents: {n_solvents}")
    print(f"Number of pairs: {n_pairs}")
    print(f"Temperatures: {[f'{T-273.15:.0f}°C' for T in TEMPERATURES]}")
    print(f"Composition points: {N_COMPOSITIONS}")
    print(f"Total data points: {n_pairs * len(TEMPERATURES) * N_COMPOSITIONS}")
    print()

    # Generate all pairs
    pairs = list(combinations(range(n_solvents), 2))

    # Storage for all results
    all_results = []
    summary_data = []

    print(f"Calculating {len(pairs)} solvent pairs...")
    print("-" * 70)

    for pair_idx, (i, j) in enumerate(pairs):
        s1 = valid_solvents[i]
        s2 = valid_solvents[j]
        pair_name = f"{s1['name']}_{s2['name']}"

        print(f"\r[{pair_idx+1}/{n_pairs}] {s1['name']} + {s2['name']}...", end="", flush=True)

        for temp_idx, T in enumerate(TEMPERATURES):
            temp_label = TEMP_LABELS[temp_idx]

            # Calculate activity coefficients across compositions
            pair_results = calculate_pair_activity(s1, s2, T, COMPOSITIONS)

            # Calculate infinite dilution
            gamma1_inf, gamma2_inf = calculate_infinite_dilution(s1, s2, T)

            # Store detailed results
            for result in pair_results:
                all_results.append({
                    'solvent1': s1['name'],
                    'solvent2': s2['name'],
                    'cosmo1': s1['cosmo'],
                    'cosmo2': s2['cosmo'],
                    'T_K': T,
                    'T_C': T - 273.15,
                    **result
                })

            # Find GE/RT at x=0.5
            mid_idx = len(COMPOSITIONS) // 2
            GE_RT_x05 = pair_results[mid_idx]['GE_RT'] if not np.isnan(pair_results[mid_idx]['GE_RT']) else np.nan

            summary_data.append({
                'solvent1': s1['name'],
                'solvent2': s2['name'],
                'T_K': T,
                'T_C': T - 273.15,
                'gamma1_inf': gamma1_inf,
                'gamma2_inf': gamma2_inf,
                'GE_RT_x05': GE_RT_x05
            })

    print()
    print("-" * 70)

    # Convert to DataFrames
    df_all = pd.DataFrame(all_results)
    df_summary = pd.DataFrame(summary_data)

    # Save results
    df_all.to_csv(OUTPUT_DIR / "activity_coefficients_all.csv", index=False)
    df_summary.to_csv(OUTPUT_DIR / "summary_statistics.csv", index=False)

    # Save per-temperature files
    for T, temp_label in zip(TEMPERATURES, TEMP_LABELS):
        df_temp = df_all[df_all['T_K'] == T]
        df_temp.to_csv(OUTPUT_DIR / f"activity_coefficients_{temp_label}.csv", index=False)

    # Save metadata
    metadata = {
        'timestamp': datetime.now().isoformat(),
        'method': 'openCOSMO-RS (opencosmorspy)',
        'parameterization': 'default_turbomole',
        'n_solvents': n_solvents,
        'n_pairs': n_pairs,
        'temperatures_K': TEMPERATURES,
        'compositions': list(COMPOSITIONS),
        'solvents': [s['name'] for s in valid_solvents],
        'total_data_points': len(df_all)
    }

    with open(OUTPUT_DIR / 'calculation_metadata.json', 'w') as f:
        json.dump(metadata, f, indent=2)

    # Print summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total data points: {len(df_all)}")
    print(f"Valid calculations: {df_all['gamma1'].notna().sum()}")
    print(f"Failed calculations: {df_all['gamma1'].isna().sum()}")
    print()

    # Show some statistics at 25°C
    df_25 = df_summary[df_summary['T_K'] == 298.15]
    print("Most non-ideal pairs at 25°C (highest |GE/RT| at x=0.5):")
    df_25_sorted = df_25.dropna(subset=['GE_RT_x05'])
    df_25_sorted = df_25_sorted.reindex(df_25_sorted['GE_RT_x05'].abs().sort_values(ascending=False).index)
    print(df_25_sorted[['solvent1', 'solvent2', 'GE_RT_x05', 'gamma1_inf', 'gamma2_inf']].head(10).to_string(index=False))

    print()
    print(f"Results saved to: {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
