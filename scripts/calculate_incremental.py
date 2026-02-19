#!/usr/bin/env python3
"""
Incremental Binary Solvent Phase Diagram Calculations using OpenCOSMO-RS

Saves results after each pair to avoid losing progress.
Can be resumed by checking which pairs are already calculated.
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

COSMO_DIR = Path(__file__).parent.parent / "all-cosmotherm-solvents"

TEMPERATURES = [298.15, 323.15, 348.15, 373.15, 398.15]
TEMP_LABELS = ['25C', '50C', '75C', '100C', '125C']
N_COMPOSITIONS = 21
COMPOSITIONS = np.linspace(0.001, 0.999, N_COMPOSITIONS)

# Resume file
PROGRESS_FILE = OUTPUT_DIR / "progress.json"
RESULTS_FILE = OUTPUT_DIR / "activity_coefficients_incremental.csv"

# ============================================================================
# Solvents
# ============================================================================

COMMON_SOLVENTS = [
    {'name': 'Pentane', 'cosmo': 'pentane_c0.cosmo'},
    {'name': 'Hexane', 'cosmo': 'hexane_c0.cosmo'},
    {'name': 'Heptane', 'cosmo': 'n-heptane_c0.cosmo'},
    {'name': 'Octane', 'cosmo': 'octane_c0.cosmo'},
    {'name': 'Decane', 'cosmo': 'n-decane_c0.cosmo'},
    {'name': 'Cyclohexane', 'cosmo': 'cyclohexane_c0.cosmo'},
    {'name': 'Benzene', 'cosmo': 'benzene_c0.cosmo'},
    {'name': 'Toluene', 'cosmo': 'toluene_c0.cosmo'},
    {'name': 'Ethylbenzene', 'cosmo': 'ethylbenzene_c0.cosmo'},
    {'name': 'Methanol', 'cosmo': 'methanol_c0.cosmo'},
    {'name': 'Ethanol', 'cosmo': 'ethanol_c0.cosmo'},
    {'name': 'Propanol', 'cosmo': 'propanol_c0.cosmo'},
    {'name': '2-Propanol', 'cosmo': '2-propanol_c0.cosmo'},
    {'name': '1-Butanol', 'cosmo': '1-butanol_c0.cosmo'},
    {'name': '1-Hexanol', 'cosmo': '1-hexanol_c0.cosmo'},
    {'name': '1-Octanol', 'cosmo': '1-octanol_c0.cosmo'},
    {'name': 'Acetone', 'cosmo': 'propanone_c0.cosmo'},
    {'name': 'Acetonitrile', 'cosmo': 'acetonitrile_c0.cosmo'},
    {'name': 'DMF', 'cosmo': 'dimethylformamide_c0.cosmo'},
    {'name': 'DMSO', 'cosmo': 'dimethylsulfoxide_c0.cosmo'},
    {'name': 'NMP', 'cosmo': 'n-methyl-2-pyrrolidinone_c0.cosmo'},
    {'name': 'THF', 'cosmo': 'thf_c0.cosmo'},
    {'name': 'Dichloromethane', 'cosmo': 'ch2cl2_c0.cosmo'},
    {'name': 'Chloroform', 'cosmo': 'chcl3_c0.cosmo'},
    {'name': 'Carbon tetrachloride', 'cosmo': 'ccl4_c0.cosmo'},
    {'name': 'Diethyl ether', 'cosmo': 'diethylether_c0.cosmo'},
    {'name': 'Dioxane', 'cosmo': 'dioxane_c0.cosmo'},
    {'name': 'Ethyl acetate', 'cosmo': 'ethylacetate_c0.cosmo'},
    {'name': 'Methyl acetate', 'cosmo': 'methylacetate_c0.cosmo'},
    {'name': 'Acetic acid', 'cosmo': 'aceticacid_c0.cosmo'},
    {'name': 'Water', 'cosmo': 'h2o_c0.cosmo'},
]


def validate_solvents():
    valid = []
    for s in COMMON_SOLVENTS:
        if (COSMO_DIR / s['cosmo']).exists():
            valid.append(s)
    return valid


def calculate_pair(s1, s2):
    """Calculate all temperatures and compositions for one pair."""
    results = []
    cosmo1_path = str(COSMO_DIR / s1['cosmo'])
    cosmo2_path = str(COSMO_DIR / s2['cosmo'])

    for T in TEMPERATURES:
        for x1 in COMPOSITIONS:
            x2 = 1.0 - x1
            x = np.array([x1, x2])

            try:
                par = Parameterization("default_turbomole")
                cosmo = COSMORS(par)
                cosmo.add_molecule([cosmo1_path])
                cosmo.add_molecule([cosmo2_path])
                cosmo.add_job(x=x, T=T, refst="pure_component")
                result = cosmo.calculate()

                ln_gamma = result['tot']['lng'][0]
                gamma1 = np.exp(ln_gamma[0])
                gamma2 = np.exp(ln_gamma[1])
                GE_RT = x1 * ln_gamma[0] + x2 * ln_gamma[1]

                results.append({
                    'solvent1': s1['name'],
                    'solvent2': s2['name'],
                    'T_K': T,
                    'T_C': T - 273.15,
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
                    'solvent1': s1['name'],
                    'solvent2': s2['name'],
                    'T_K': T,
                    'T_C': T - 273.15,
                    'x1': x1,
                    'x2': x2,
                    'gamma1': np.nan,
                    'gamma2': np.nan,
                    'ln_gamma1': np.nan,
                    'ln_gamma2': np.nan,
                    'GE_RT': np.nan
                })

    return results


def load_progress():
    if PROGRESS_FILE.exists():
        with open(PROGRESS_FILE) as f:
            return json.load(f)
    return {'completed_pairs': [], 'last_pair_idx': -1}


def save_progress(completed_pairs, last_pair_idx):
    with open(PROGRESS_FILE, 'w') as f:
        json.dump({'completed_pairs': completed_pairs, 'last_pair_idx': last_pair_idx}, f)


def main():
    print("=" * 70)
    print("OpenCOSMO-RS INCREMENTAL CALCULATIONS")
    print("=" * 70)

    valid_solvents = validate_solvents()
    n_solvents = len(valid_solvents)
    pairs = list(combinations(range(n_solvents), 2))
    n_pairs = len(pairs)

    print(f"Solvents: {n_solvents}, Pairs: {n_pairs}")
    print(f"Results file: {RESULTS_FILE}")
    print(f"Progress file: {PROGRESS_FILE}")

    # Load progress
    progress = load_progress()
    completed = set(tuple(p) for p in progress['completed_pairs'])
    start_idx = progress['last_pair_idx'] + 1

    print(f"\nResuming from pair {start_idx}/{n_pairs}")
    print(f"Already completed: {len(completed)} pairs")
    print("-" * 70)

    # Open results file in append mode
    write_header = not RESULTS_FILE.exists() or RESULTS_FILE.stat().st_size == 0

    for pair_idx, (i, j) in enumerate(pairs):
        if pair_idx < start_idx:
            continue

        if (i, j) in completed:
            continue

        s1 = valid_solvents[i]
        s2 = valid_solvents[j]

        print(f"[{pair_idx+1}/{n_pairs}] {s1['name']} + {s2['name']}...", end="", flush=True)

        # Calculate
        results = calculate_pair(s1, s2)

        # Save immediately
        df = pd.DataFrame(results)
        df.to_csv(RESULTS_FILE, mode='a', header=write_header, index=False)
        write_header = False

        # Update progress
        completed.add((i, j))
        save_progress(list(completed), pair_idx)

        # Count valid results
        valid = sum(1 for r in results if not np.isnan(r['gamma1']))
        print(f" {valid}/{len(results)} valid")

    print("\n" + "=" * 70)
    print("COMPLETE!")
    print(f"Results: {RESULTS_FILE}")
    print("=" * 70)


if __name__ == '__main__':
    main()
