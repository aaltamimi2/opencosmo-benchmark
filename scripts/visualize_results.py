#!/usr/bin/env python3
"""
Visualize OpenCOSMO-RS Binary Solvent Phase Diagram Results

Creates:
1. GE/RT heatmaps at each temperature
2. Activity coefficient vs composition plots for selected pairs
3. Temperature dependence plots
4. Comparison with literature LLE data
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import interp1d

# ============================================================================
# Configuration
# ============================================================================

REPO = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO / "data" / "calculated" / "opencosmo_activity_coefficients"
PLOTS_DIR = REPO / "plots"
PLOTS_DIR.mkdir(exist_ok=True, parents=True)

# Literature data file
LIT_FILE = REPO / "data" / "experimental" / "literature_lle_numerical.csv"

# Key systems for detailed plots
KEY_SYSTEMS = [
    ('Methanol', 'Hexane'),
    ('Acetonitrile', 'Hexane'),
    ('Water', '1-Butanol'),
    ('Water', 'Hexane'),
    ('Water', 'Ethyl acetate'),
    ('Water', 'THF'),
]

TEMP_LABELS = ['25C', '50C', '75C', '100C', '125C']
TEMPERATURES = [298.15, 323.15, 348.15, 373.15, 398.15]


# ============================================================================
# Visualization Functions
# ============================================================================

def create_ge_heatmap(df_summary, T_K, output_path):
    """Create GE/RT heatmap at a specific temperature."""
    df_temp = df_summary[df_summary['T_K'] == T_K].copy()

    if len(df_temp) == 0:
        print(f"  No data at {T_K} K")
        return

    # Get unique solvents
    solvents = sorted(set(df_temp['solvent1'].tolist() + df_temp['solvent2'].tolist()))
    n_solvents = len(solvents)
    solvent_idx = {s: i for i, s in enumerate(solvents)}

    # Create matrix
    matrix = np.zeros((n_solvents, n_solvents))
    matrix[:] = np.nan

    for _, row in df_temp.iterrows():
        i = solvent_idx[row['solvent1']]
        j = solvent_idx[row['solvent2']]
        val = row['GE_RT_x05']
        matrix[i, j] = val
        matrix[j, i] = val

    # Plot
    fig, ax = plt.subplots(figsize=(14, 12))

    # Use diverging colormap centered at 0
    vmax = max(abs(np.nanmin(matrix)), abs(np.nanmax(matrix)))
    im = ax.imshow(matrix, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal')

    # Labels
    ax.set_xticks(range(n_solvents))
    ax.set_yticks(range(n_solvents))
    ax.set_xticklabels(solvents, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(solvents, fontsize=8)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('GE/RT at x=0.5', fontsize=12)

    T_C = T_K - 273.15
    ax.set_title(f'OpenCOSMO-RS: Excess Gibbs Energy (GE/RT) at {T_C:.0f}°C\n'
                 f'Red = positive (repulsive), Blue = negative (attractive)', fontsize=14)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def plot_phase_diagram(df_all, solvent1, solvent2, output_path):
    """Create phase diagram for a specific solvent pair."""
    # Get data for this pair
    mask = ((df_all['solvent1'] == solvent1) & (df_all['solvent2'] == solvent2)) | \
           ((df_all['solvent1'] == solvent2) & (df_all['solvent2'] == solvent1))
    df_pair = df_all[mask].copy()

    if len(df_pair) == 0:
        print(f"  No data for {solvent1} + {solvent2}")
        return

    # Ensure consistent ordering
    swap_needed = df_pair['solvent1'].iloc[0] != solvent1
    if swap_needed:
        df_pair['x1'], df_pair['x2'] = df_pair['x2'], df_pair['x1']
        df_pair['gamma1'], df_pair['gamma2'] = df_pair['gamma2'], df_pair['gamma1']

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Color map for temperatures
    colors = plt.cm.plasma(np.linspace(0, 0.8, len(TEMPERATURES)))

    # Plot 1: Activity coefficients vs composition
    ax1 = axes[0]
    for T, color, label in zip(TEMPERATURES, colors, TEMP_LABELS):
        df_t = df_pair[df_pair['T_K'] == T]
        if len(df_t) > 0:
            ax1.plot(df_t['x1'], df_t['gamma1'], '-', color=color, label=f'{label} γ₁')
            ax1.plot(df_t['x1'], df_t['gamma2'], '--', color=color, label=f'{label} γ₂')

    ax1.set_xlabel(f'{solvent1} mole fraction', fontsize=11)
    ax1.set_ylabel('Activity coefficient γ', fontsize=11)
    ax1.set_title('Activity Coefficients')
    ax1.legend(fontsize=8, ncol=2)
    ax1.set_xlim(0, 1)
    ax1.grid(True, alpha=0.3)

    # Plot 2: ln(gamma) vs composition
    ax2 = axes[1]
    for T, color, label in zip(TEMPERATURES, colors, TEMP_LABELS):
        df_t = df_pair[df_pair['T_K'] == T]
        if len(df_t) > 0:
            ax2.plot(df_t['x1'], df_t['ln_gamma1'], '-', color=color, label=f'{label}')
            ax2.plot(df_t['x1'], df_t['ln_gamma2'], '--', color=color)

    ax2.set_xlabel(f'{solvent1} mole fraction', fontsize=11)
    ax2.set_ylabel('ln(γ)', fontsize=11)
    ax2.set_title('Logarithmic Activity Coefficients')
    ax2.legend(fontsize=9)
    ax2.set_xlim(0, 1)
    ax2.grid(True, alpha=0.3)

    # Plot 3: GE/RT vs composition
    ax3 = axes[2]
    for T, color, label in zip(TEMPERATURES, colors, TEMP_LABELS):
        df_t = df_pair[df_pair['T_K'] == T]
        if len(df_t) > 0:
            ax3.plot(df_t['x1'], df_t['GE_RT'], '-o', color=color, label=label, markersize=3)

    ax3.set_xlabel(f'{solvent1} mole fraction', fontsize=11)
    ax3.set_ylabel('GE/RT', fontsize=11)
    ax3.set_title('Excess Gibbs Energy')
    ax3.legend(fontsize=9)
    ax3.set_xlim(0, 1)
    ax3.grid(True, alpha=0.3)
    ax3.axhline(0, color='k', linewidth=0.5)

    plt.suptitle(f'OpenCOSMO-RS Phase Diagram: {solvent1} + {solvent2}', fontsize=14, y=1.02)
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def load_literature_data():
    """Load literature LLE data for comparison."""
    if LIT_FILE.exists():
        return pd.read_csv(LIT_FILE)
    return None


def create_literature_comparison(df_all, df_summary, lit_df, output_path):
    """Create comparison plot with literature data."""
    if lit_df is None or len(lit_df) == 0:
        print("  No literature data available")
        return

    # Systems to compare
    systems = [
        ('methanol_hexane', 'Methanol', 'Hexane', 'Methanol'),
        ('acetonitrile_hexane', 'Acetonitrile', 'Hexane', 'Acetonitrile'),
        ('water_butanol', 'Water', '1-Butanol', 'Water'),
        ('water_hexane', 'Water', 'Hexane', 'Water'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for ax, (sys_name, s1, s2, x1_name) in zip(axes, systems):
        # Get calculated data
        mask = ((df_all['solvent1'] == s1) & (df_all['solvent2'] == s2)) | \
               ((df_all['solvent1'] == s2) & (df_all['solvent2'] == s1))
        df_calc = df_all[mask].copy()

        # Get literature data
        lit_sys = lit_df[lit_df['system'] == sys_name]

        if len(df_calc) == 0:
            ax.text(0.5, 0.5, f'No data for\n{s1} + {s2}', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{s1} + {s2}')
            continue

        # Plot calculated GE/RT at 25°C
        df_25 = df_calc[df_calc['T_K'] == 298.15]
        if len(df_25) > 0:
            # Ensure correct x1 ordering
            if df_25['solvent1'].iloc[0] != s1:
                x1 = df_25['x2'].values
                ge = df_25['GE_RT'].values
            else:
                x1 = df_25['x1'].values
                ge = df_25['GE_RT'].values

            ax.plot(x1, ge, 'b-', linewidth=2, label='OpenCOSMO-RS (25°C)')

        # Plot literature points if available (show as annotations)
        if len(lit_sys) > 0:
            for _, row in lit_sys.iterrows():
                if row['property'] in ['mutual_solubility', 'binodal']:
                    # Mark phase boundaries
                    ax.axvline(row['x1_phase1'], color='r', linestyle='--', alpha=0.5)
                    ax.axvline(row['x1_phase2'], color='r', linestyle='--', alpha=0.5)
            ax.axvline(np.nan, color='r', linestyle='--', alpha=0.5, label='Lit. phase boundary')

        ax.set_xlabel(f'{x1_name} mole fraction')
        ax.set_ylabel('GE/RT')
        ax.set_title(f'{s1} + {s2}')
        ax.set_xlim(0, 1)
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='k', linewidth=0.5)
        ax.legend(fontsize=9)

    plt.suptitle('OpenCOSMO-RS vs Literature: Key LLE Systems', fontsize=14, y=1.02)
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def create_statistics_summary(df_summary, output_path):
    """Create summary statistics text file."""
    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("OPENCOSMO-RS BINARY SOLVENT PHASE DIAGRAM SUMMARY\n")
        f.write("=" * 80 + "\n\n")

        df_25 = df_summary[df_summary['T_K'] == 298.15].copy()

        # Most non-ideal
        f.write("Most Non-Ideal Pairs (highest |GE/RT| at x=0.5, 25°C):\n")
        f.write("-" * 80 + "\n")
        df_sorted = df_25.dropna(subset=['GE_RT_x05'])
        df_sorted = df_sorted.reindex(df_sorted['GE_RT_x05'].abs().sort_values(ascending=False).index)
        f.write(df_sorted[['solvent1', 'solvent2', 'GE_RT_x05', 'gamma1_inf', 'gamma2_inf']].head(20).to_string(index=False))
        f.write("\n\n")

        # Most ideal
        f.write("Most Ideal Pairs (lowest |GE/RT| at x=0.5, 25°C):\n")
        f.write("-" * 80 + "\n")
        df_sorted = df_25.dropna(subset=['GE_RT_x05'])
        df_sorted = df_sorted.reindex(df_sorted['GE_RT_x05'].abs().sort_values(ascending=True).index)
        f.write(df_sorted[['solvent1', 'solvent2', 'GE_RT_x05', 'gamma1_inf', 'gamma2_inf']].head(20).to_string(index=False))
        f.write("\n\n")

        # Overall statistics
        f.write("Overall Statistics (25°C):\n")
        f.write("-" * 80 + "\n")
        f.write(f"Number of pairs: {len(df_25)}\n")
        f.write(f"GE/RT range: {df_25['GE_RT_x05'].min():.4f} to {df_25['GE_RT_x05'].max():.4f}\n")
        f.write(f"Mean |GE/RT|: {df_25['GE_RT_x05'].abs().mean():.4f}\n")
        f.write(f"Pairs with positive GE/RT: {(df_25['GE_RT_x05'] > 0).sum()}\n")
        f.write(f"Pairs with negative GE/RT: {(df_25['GE_RT_x05'] < 0).sum()}\n")
        f.write("\n")

    print(f"  Saved: {output_path.name}")


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 70)
    print("OpenCOSMO-RS Results Visualization")
    print("=" * 70)
    print()

    # Load data
    all_file = RESULTS_DIR / "activity_coefficients_all.csv"
    summary_file = RESULTS_DIR / "summary_statistics.csv"

    if not all_file.exists():
        print(f"ERROR: Results file not found: {all_file}")
        print("Run calculate_binary_phase_diagrams.py first.")
        return

    print("Loading data...")
    df_all = pd.read_csv(all_file)
    df_summary = pd.read_csv(summary_file)

    print(f"  Total data points: {len(df_all)}")
    print(f"  Summary entries: {len(df_summary)}")
    print()

    # Load literature data
    lit_df = load_literature_data()
    if lit_df is not None:
        print(f"  Literature data points: {len(lit_df)}")

    # Create heatmaps
    print("\nCreating GE/RT heatmaps...")
    for T, label in zip(TEMPERATURES, TEMP_LABELS):
        create_ge_heatmap(df_summary, T, PLOTS_DIR / f"ge_heatmap_{label}.png")

    # Create phase diagrams for key systems
    print("\nCreating phase diagrams for key systems...")
    for s1, s2 in KEY_SYSTEMS:
        safe_name = f"{s1}_{s2}".replace(' ', '_').replace('-', '')
        plot_phase_diagram(df_all, s1, s2, PLOTS_DIR / f"phase_diagram_{safe_name}.png")

    # Create literature comparison
    print("\nCreating literature comparison...")
    create_literature_comparison(df_all, df_summary, lit_df, PLOTS_DIR / "literature_comparison.png")

    # Create statistics summary
    print("\nCreating statistics summary...")
    create_statistics_summary(df_summary, PLOTS_DIR / "summary_statistics.txt")

    print()
    print("=" * 70)
    print(f"All plots saved to: {PLOTS_DIR}")
    print("=" * 70)


if __name__ == '__main__':
    main()
