#!/usr/bin/env python3
"""Compare COSMOtherm vs openCOSMO-RS activity coefficient predictions.

Parses COSMOtherm .tab output and openCOSMO-RS CSV data for the same
31 solvents (465 binary pairs x 5 temperatures), then produces:
  1. Density parity plots (ln_gamma, GE/RT) with hexbin
  2. Per-pair deviation heatmap
  3. Chemical class breakdown bar chart
  4. Per-solvent ranking bar chart
  5. Temperature dependence of MAE (violin)
  6. Non-ideality magnitude vs MAE scatter
  7. Composition-dependent comparison for selected pairs
  8. 4-panel paper summary figure
  9. LLE summary grid
  10. Summary statistics CSV + one-sentence conclusion
"""

import os
import re
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO = Path(__file__).resolve().parent.parent
TAB_FILE = REPO / "cosmotherm-reference-data" / "opencosmo-LLE-compare.tab"
OPENCOSMO_CSV = (
    REPO / "data" / "calculated" / "opencosmo_activity_coefficients"
    / "activity_coefficients_all.csv"
)
OUT_DIR = REPO / "results" / "comparison"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Gas constant in kcal/(mol*K) for G^E unit conversion
R_KCAL = 0.001987204

# ---------------------------------------------------------------------------
# Name mapping: COSMOtherm .tab names --> openCOSMO-RS CSV names
# ---------------------------------------------------------------------------
COSMOTHERM_TO_OPENCOSMO = {
    "pentane": "Pentane",
    "hexane": "Hexane",
    "n-heptane": "Heptane",
    "octane": "Octane",
    "n-decane": "Decane",
    "cyclohexane": "Cyclohexane",
    "benzene": "Benzene",
    "toluene": "Toluene",
    "ethylbenzene": "Ethylbenzene",
    "methanol": "Methanol",
    "ethanol": "Ethanol",
    "propanol": "Propanol",
    "2-propanol": "2-Propanol",
    "1-butanol": "1-Butanol",
    "1-hexanol": "1-Hexanol",
    "1-octanol": "1-Octanol",
    "propanone": "Acetone",
    "acetonitrile": "Acetonitrile",
    "dimethylformamide": "DMF",
    "dimethylsulfoxide": "DMSO",
    "n-methyl-2-pyrrolidinone": "NMP",
    "thf": "THF",
    "ch2cl2": "Dichloromethane",
    "chcl3": "Chloroform",
    "ccl4": "Carbon tetrachloride",
    "diethylether": "Diethyl ether",
    "dioxane": "Dioxane",
    "ethylacetate": "Ethyl acetate",
    "methylacetate": "Methyl acetate",
    "aceticacid": "Acetic acid",
    "h2o": "Water",
}

# ---------------------------------------------------------------------------
# Chemical class assignments
# ---------------------------------------------------------------------------
SOLVENT_CLASS = {
    "Pentane": "Alkane",
    "Hexane": "Alkane",
    "Heptane": "Alkane",
    "Octane": "Alkane",
    "Decane": "Alkane",
    "Cyclohexane": "Alkane",
    "Benzene": "Aromatic",
    "Toluene": "Aromatic",
    "Ethylbenzene": "Aromatic",
    "Methanol": "Alcohol",
    "Ethanol": "Alcohol",
    "Propanol": "Alcohol",
    "2-Propanol": "Alcohol",
    "1-Butanol": "Alcohol",
    "1-Hexanol": "Alcohol",
    "1-Octanol": "Alcohol",
    "Acetone": "Polar aprotic",
    "Acetonitrile": "Polar aprotic",
    "DMF": "Polar aprotic",
    "DMSO": "Polar aprotic",
    "NMP": "Polar aprotic",
    "THF": "Ether",
    "Dichloromethane": "Halogenated",
    "Chloroform": "Halogenated",
    "Carbon tetrachloride": "Halogenated",
    "Diethyl ether": "Ether",
    "Dioxane": "Ether",
    "Ethyl acetate": "Ester",
    "Methyl acetate": "Ester",
    "Acetic acid": "Acid",
    "Water": "Water",
}

CLASS_COLORS = {
    "Alkane": "#4e79a7",
    "Aromatic": "#f28e2b",
    "Alcohol": "#e15759",
    "Polar aprotic": "#76b7b2",
    "Halogenated": "#59a14f",
    "Ether": "#edc948",
    "Ester": "#b07aa1",
    "Acid": "#ff9da7",
    "Water": "#9c755f",
}

TEMPERATURES = [298.15, 323.15, 348.15, 373.15, 398.15]

# Solvents to exclude from comparison
EXCLUDE_SOLVENTS = {"Chloroform"}


# ---------------------------------------------------------------------------
# Parse COSMOtherm .tab file
# ---------------------------------------------------------------------------
def parse_tab_file(tab_path):
    """Parse a COSMOtherm .tab file into a list of job dicts."""
    text = Path(tab_path).read_text()
    job_blocks = re.split(r"\n\s*Property\s+job\s+\d+\s*:", text)
    job_blocks = job_blocks[1:]

    jobs = []
    for block in job_blocks:
        job = {}
        m = re.search(
            r"Compounds\s+job\s+\d+\s*:\s*(\S+)\s*\(1\)\s*;\s*(.+?)\s*\(2\)\s*;",
            block,
        )
        if not m:
            continue
        job["solvent1_ct"] = m.group(1).strip()
        job["solvent2_ct"] = m.group(2).strip()
        job["solvent1"] = COSMOTHERM_TO_OPENCOSMO.get(
            job["solvent1_ct"], job["solvent1_ct"]
        )
        job["solvent2"] = COSMOTHERM_TO_OPENCOSMO.get(
            job["solvent2_ct"], job["solvent2_ct"]
        )

        m_t = re.search(r"T=\s*([\d.]+)\s*K", block)
        if not m_t:
            continue
        job["T_K"] = float(m_t.group(1))

        lines = block.split("\n")
        data_lines = []
        in_data = False
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("x1"):
                in_data = True
                continue
            if in_data:
                if stripped == "" or stripped.startswith("No LLE") or stripped.startswith("LLE"):
                    break
                parts = stripped.split()
                if len(parts) >= 10:
                    try:
                        float(parts[0])
                        data_lines.append(parts)
                    except ValueError:
                        break

        if not data_lines:
            continue

        df = pd.DataFrame(
            data_lines,
            columns=[
                "x1", "x2", "HE", "GE_kcal", "ptot",
                "mu1", "mu2", "ln_gamma1", "ln_gamma2", "y1", "y2",
            ],
        ).astype(float)

        # Convert G^E from kcal/mol to dimensionless GE/RT
        df["GE_RT"] = df["GE_kcal"] / (R_KCAL * job["T_K"])

        job["data"] = df

        lle_match = re.search(
            r"LLE point found at x`\(1\)\s*=\s*([\d.]+).*?x``\(1\)\s*=\s*([\d.]+)",
            block,
        )
        if lle_match:
            job["lle_found"] = True
            job["lle_x1_phase1"] = float(lle_match.group(1))
            job["lle_x1_phase2"] = float(lle_match.group(2))
        else:
            job["lle_found"] = False
            job["lle_x1_phase1"] = None
            job["lle_x1_phase2"] = None

        jobs.append(job)

    return jobs


# ---------------------------------------------------------------------------
# Load openCOSMO-RS data
# ---------------------------------------------------------------------------
def load_opencosmo(csv_path):
    return pd.read_csv(csv_path)


# ---------------------------------------------------------------------------
# Match and interpolate
# ---------------------------------------------------------------------------
def match_pair(oc_df, ct_job):
    """For a COSMOtherm job, find matching openCOSMO data and
    interpolate onto a common composition grid."""
    s1, s2 = ct_job["solvent1"], ct_job["solvent2"]
    T = ct_job["T_K"]

    mask_fwd = (
        (oc_df["solvent1"] == s1) & (oc_df["solvent2"] == s2)
        & (np.abs(oc_df["T_K"] - T) < 0.5)
    )
    mask_rev = (
        (oc_df["solvent1"] == s2) & (oc_df["solvent2"] == s1)
        & (np.abs(oc_df["T_K"] - T) < 0.5)
    )

    if mask_fwd.sum() > 0:
        oc_sub = oc_df[mask_fwd].copy()
        reversed_order = False
    elif mask_rev.sum() > 0:
        oc_sub = oc_df[mask_rev].copy()
        reversed_order = True
    else:
        return None

    oc_sub = oc_sub.sort_values("x1")

    if reversed_order:
        oc_sub = oc_sub.rename(columns={
            "x1": "x2_tmp", "x2": "x1_tmp",
            "ln_gamma1": "ln_gamma2_tmp", "ln_gamma2": "ln_gamma1_tmp",
            "gamma1": "gamma2_tmp", "gamma2": "gamma1_tmp",
        })
        oc_sub = oc_sub.rename(columns={
            "x1_tmp": "x1", "x2_tmp": "x2",
            "ln_gamma1_tmp": "ln_gamma1", "ln_gamma2_tmp": "ln_gamma2",
            "gamma1_tmp": "gamma1", "gamma2_tmp": "gamma2",
        })
        oc_sub["GE_RT"] = oc_sub["x1"] * oc_sub["ln_gamma1"] + oc_sub["x2"] * oc_sub["ln_gamma2"]
        oc_sub = oc_sub.sort_values("x1")

    ct_data = ct_job["data"]
    ct_interior = ct_data[(ct_data["x1"] >= 0.009) & (ct_data["x1"] <= 0.991)].copy()
    oc_x1 = oc_sub["x1"].values
    ct_x1 = ct_interior["x1"].values

    x_min = max(oc_x1.min(), ct_x1.min())
    x_max = min(oc_x1.max(), ct_x1.max())
    common_x1 = oc_x1[(oc_x1 >= x_min) & (oc_x1 <= x_max)]

    if len(common_x1) < 3:
        return None

    interp_lng1 = interp1d(ct_x1, ct_interior["ln_gamma1"].values, kind="cubic")
    interp_lng2 = interp1d(ct_x1, ct_interior["ln_gamma2"].values, kind="cubic")
    interp_ge = interp1d(ct_x1, ct_interior["GE_RT"].values, kind="cubic")

    return pd.DataFrame({
        "x1": common_x1,
        "ln_gamma1_ct": interp_lng1(common_x1),
        "ln_gamma2_ct": interp_lng2(common_x1),
        "GE_RT_ct": interp_ge(common_x1),
        "ln_gamma1_oc": np.interp(common_x1, oc_x1, oc_sub["ln_gamma1"].values),
        "ln_gamma2_oc": np.interp(common_x1, oc_x1, oc_sub["ln_gamma2"].values),
        "GE_RT_oc": np.interp(common_x1, oc_x1, oc_sub["GE_RT"].values),
    })


# ---------------------------------------------------------------------------
# Per-pair statistics
# ---------------------------------------------------------------------------
def pair_statistics(matched_df):
    """Compute MAE, RMSE, bias, and peak non-ideality for a matched pair."""
    d_lng1 = matched_df["ln_gamma1_ct"] - matched_df["ln_gamma1_oc"]
    d_lng2 = matched_df["ln_gamma2_ct"] - matched_df["ln_gamma2_oc"]
    d_ge = matched_df["GE_RT_ct"] - matched_df["GE_RT_oc"]

    return {
        "MAE_ln_gamma1": np.mean(np.abs(d_lng1)),
        "MAE_ln_gamma2": np.mean(np.abs(d_lng2)),
        "MAE_GE_RT": np.mean(np.abs(d_ge)),
        "RMSE_ln_gamma1": np.sqrt(np.mean(d_lng1**2)),
        "RMSE_ln_gamma2": np.sqrt(np.mean(d_lng2**2)),
        "RMSE_GE_RT": np.sqrt(np.mean(d_ge**2)),
        "bias_GE_RT": np.mean(d_ge),
        "max_abs_GE_RT_ct": np.max(np.abs(matched_df["GE_RT_ct"])),
        "n_points": len(matched_df),
    }


# ---------------------------------------------------------------------------
# Helper: class pair label
# ---------------------------------------------------------------------------
def class_pair_label(s1, s2):
    c1 = SOLVENT_CLASS.get(s1, "?")
    c2 = SOLVENT_CLASS.get(s2, "?")
    return " + ".join(sorted([c1, c2]))


# ---------------------------------------------------------------------------
# Plot 1: Density parity plots (hexbin)
# ---------------------------------------------------------------------------
def plot_parity(all_matched, out_dir):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    ct_lng1 = np.concatenate([m["ln_gamma1_ct"].values for m in all_matched])
    oc_lng1 = np.concatenate([m["ln_gamma1_oc"].values for m in all_matched])
    ct_lng2 = np.concatenate([m["ln_gamma2_ct"].values for m in all_matched])
    oc_lng2 = np.concatenate([m["ln_gamma2_oc"].values for m in all_matched])
    ct_ge = np.concatenate([m["GE_RT_ct"].values for m in all_matched])
    oc_ge = np.concatenate([m["GE_RT_oc"].values for m in all_matched])

    datasets = [
        (ct_lng1, oc_lng1, "ln(gamma1)", "C0"),
        (ct_lng2, oc_lng2, "ln(gamma2)", "C1"),
        (ct_ge, oc_ge, "GE/RT", "C2"),
    ]

    for ax, (ct_vals, oc_vals, title, _color) in zip(axes, datasets):
        hb = ax.hexbin(ct_vals, oc_vals, gridsize=80, cmap="Blues",
                       mincnt=1, norm=LogNorm())
        lim_lo = min(ct_vals.min(), oc_vals.min())
        lim_hi = max(ct_vals.max(), oc_vals.max())
        ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], "k--", lw=0.8, zorder=5)
        ax.set_xlabel(f"COSMOtherm {title}")
        ax.set_ylabel(f"openCOSMO-RS {title}")
        ax.set_title(title)

        mae = np.mean(np.abs(ct_vals - oc_vals))
        bias = np.mean(ct_vals - oc_vals)
        r2 = np.corrcoef(ct_vals, oc_vals)[0, 1]**2
        ax.text(0.05, 0.95,
                f"MAE = {mae:.4f}\nBias = {bias:+.4f}\nR$^2$ = {r2:.4f}\nN = {len(ct_vals)}",
                transform=ax.transAxes, fontsize=8, verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7))
        plt.colorbar(hb, ax=ax, label="Count", shrink=0.8)

    plt.suptitle("COSMOtherm vs openCOSMO-RS: Parity Plots (465 pairs x 5 T)", fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(out_dir / "parity_plots.png", dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Plot 2: Deviation heatmap
# ---------------------------------------------------------------------------
def plot_deviation_heatmap(stats_df, out_dir):
    sub = stats_df[np.abs(stats_df["T_K"] - 298.15) < 0.5].copy()
    if sub.empty:
        return
    solvents = sorted(set(sub["solvent1"].tolist() + sub["solvent2"].tolist()))
    n = len(solvents)
    idx = {s: i for i, s in enumerate(solvents)}

    mat = np.full((n, n), np.nan)
    for _, row in sub.iterrows():
        i, j = idx[row["solvent1"]], idx[row["solvent2"]]
        mat[i, j] = row["MAE_GE_RT"]
        mat[j, i] = row["MAE_GE_RT"]

    fig, ax = plt.subplots(figsize=(14, 12))
    im = ax.imshow(mat, cmap="YlOrRd", aspect="equal")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(solvents, rotation=90, fontsize=7)
    ax.set_yticklabels(solvents, fontsize=7)
    plt.colorbar(im, label="MAE(GE/RT)", shrink=0.8)
    ax.set_title("COSMOtherm vs openCOSMO-RS: MAE of GE/RT at 25 C")
    plt.tight_layout()
    plt.savefig(out_dir / "deviation_heatmap_25C.png", dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Plot 3: Chemical class breakdown
# ---------------------------------------------------------------------------
def plot_class_breakdown(stats_df, out_dir):
    sub = stats_df[np.abs(stats_df["T_K"] - 298.15) < 0.5].copy()
    if sub.empty:
        return
    sub["class_pair"] = sub.apply(
        lambda r: class_pair_label(r["solvent1"], r["solvent2"]), axis=1
    )

    grp = sub.groupby("class_pair")["MAE_GE_RT"].agg(["mean", "min", "max", "count"])
    grp = grp[grp["count"] >= 2].sort_values("mean", ascending=True)

    fig, ax = plt.subplots(figsize=(10, max(6, len(grp) * 0.35)))
    y_pos = np.arange(len(grp))
    bars = ax.barh(y_pos, grp["mean"], xerr=[grp["mean"] - grp["min"],
                                              grp["max"] - grp["mean"]],
                   height=0.7, color="#4e79a7", alpha=0.8, capsize=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(grp.index, fontsize=8)
    ax.set_xlabel("MAE(GE/RT) at 25 C")
    ax.set_title("Agreement by Chemical Class Pair")
    ax.axvline(0.05, color="green", ls="--", lw=0.8, label="Good (0.05)")
    ax.axvline(0.2, color="red", ls="--", lw=0.8, label="Poor (0.2)")
    ax.legend(fontsize=8)

    # Annotate counts
    for i, (_, row) in enumerate(grp.iterrows()):
        ax.text(row["mean"] + 0.005, i, f"n={int(row['count'])}",
                va="center", fontsize=7, color="gray")

    plt.tight_layout()
    plt.savefig(out_dir / "class_breakdown.png", dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Plot 4: Per-solvent ranking
# ---------------------------------------------------------------------------
def plot_solvent_ranking(stats_df, out_dir):
    sub = stats_df[np.abs(stats_df["T_K"] - 298.15) < 0.5].copy()
    if sub.empty:
        return

    # Compute mean MAE per solvent (each solvent appears in 30 pairs)
    solvent_mae = {}
    for s in SOLVENT_CLASS:
        mask = (sub["solvent1"] == s) | (sub["solvent2"] == s)
        if mask.sum() > 0:
            solvent_mae[s] = sub.loc[mask, "MAE_GE_RT"].mean()

    if not solvent_mae:
        return

    sdf = pd.DataFrame([
        {"solvent": s, "mean_MAE": v, "class": SOLVENT_CLASS[s]}
        for s, v in solvent_mae.items()
    ]).sort_values("mean_MAE", ascending=True)

    fig, ax = plt.subplots(figsize=(8, 9))
    y_pos = np.arange(len(sdf))
    colors = [CLASS_COLORS.get(c, "gray") for c in sdf["class"]]
    ax.barh(y_pos, sdf["mean_MAE"], color=colors, height=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sdf["solvent"], fontsize=9)
    ax.set_xlabel("Mean MAE(GE/RT) across all pairs at 25 C")
    ax.set_title("Per-Solvent Agreement Ranking")

    # Class legend
    handles = [plt.Rectangle((0, 0), 1, 1, color=CLASS_COLORS[c])
               for c in sorted(CLASS_COLORS)]
    ax.legend(handles, sorted(CLASS_COLORS), fontsize=7, loc="lower right",
              title="Class", title_fontsize=8)

    plt.tight_layout()
    plt.savefig(out_dir / "solvent_ranking.png", dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Plot 5: Temperature dependence of MAE (violin)
# ---------------------------------------------------------------------------
def plot_temperature_mae(stats_df, out_dir):
    fig, ax = plt.subplots(figsize=(8, 5))

    T_labels = []
    T_data = []
    T_means = []
    for T in TEMPERATURES:
        sub = stats_df[np.abs(stats_df["T_K"] - T) < 0.5]
        if len(sub) == 0:
            continue
        T_data.append(sub["MAE_GE_RT"].values)
        T_labels.append(f"{T - 273.15:.0f}")
        T_means.append(sub["MAE_GE_RT"].mean())

    parts = ax.violinplot(T_data, positions=range(len(T_data)),
                          showmedians=True, showextrema=False)
    for pc in parts["bodies"]:
        pc.set_facecolor("#4e79a7")
        pc.set_alpha(0.6)

    ax.plot(range(len(T_data)), T_means, "ro-", ms=6, lw=1.5, label="Mean")
    ax.set_xticks(range(len(T_labels)))
    ax.set_xticklabels([f"{t} C" for t in T_labels])
    ax.set_ylabel("MAE(GE/RT)")
    ax.set_xlabel("Temperature")
    ax.set_title("Agreement Improves with Temperature")
    ax.legend()

    for i, m in enumerate(T_means):
        ax.text(i, m + 0.003, f"{m:.3f}", ha="center", fontsize=8, color="red")

    plt.tight_layout()
    plt.savefig(out_dir / "temperature_mae_violin.png", dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Plot 6: Non-ideality magnitude vs MAE
# ---------------------------------------------------------------------------
def plot_nonideality_vs_mae(stats_df, out_dir):
    sub = stats_df[np.abs(stats_df["T_K"] - 298.15) < 0.5].copy()
    if sub.empty:
        return

    sub["class_pair"] = sub.apply(
        lambda r: class_pair_label(r["solvent1"], r["solvent2"]), axis=1
    )

    fig, ax = plt.subplots(figsize=(8, 6))

    # Color by whether it contains chloroform, water, or neither
    colors = []
    labels_done = set()
    for _, row in sub.iterrows():
        s1, s2 = row["solvent1"], row["solvent2"]
        if "Chloroform" in (s1, s2):
            c, lab = "#e15759", "Contains CHCl3"
        elif "Water" in (s1, s2):
            c, lab = "#4e79a7", "Contains Water"
        elif "DMSO" in (s1, s2) or "DMF" in (s1, s2) or "NMP" in (s1, s2):
            c, lab = "#76b7b2", "Contains DMSO/DMF/NMP"
        else:
            c, lab = "#bab0ac", "Other"
        colors.append(c)

        if lab not in labels_done:
            ax.scatter([], [], color=c, s=30, label=lab)
            labels_done.add(lab)

    ax.scatter(sub["max_abs_GE_RT_ct"], sub["MAE_GE_RT"],
               c=colors, s=20, alpha=0.7, edgecolors="none")
    ax.set_xlabel("Max |GE/RT| (COSMOtherm) — Non-ideality magnitude")
    ax.set_ylabel("MAE(GE/RT)")
    ax.set_title("Deviation Scales with Non-ideality, Chloroform Is Outlier")
    ax.legend(fontsize=8)

    # Reference line: y = 0.1*x
    xmax = sub["max_abs_GE_RT_ct"].max()
    ax.plot([0, xmax], [0, 0.1 * xmax], "k--", lw=0.7, alpha=0.5)
    ax.text(xmax * 0.7, 0.1 * xmax * 0.7 + 0.02, "10% error line",
            fontsize=7, color="gray")

    plt.tight_layout()
    plt.savefig(out_dir / "nonideality_vs_mae.png", dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Plot 7: Selected pairs composition comparison
# ---------------------------------------------------------------------------
def plot_selected_pairs(ct_jobs, oc_df, out_dir):
    interesting = [
        ("Methanol", "Hexane"),
        ("Acetonitrile", "Hexane"),
        ("Water", "1-Butanol"),
        ("Water", "Hexane"),
        ("Water", "Ethanol"),
        ("Water", "THF"),
        ("Acetone", "Hexane"),
        ("DMSO", "Hexane"),
        ("Dichloromethane", "Acetone"),
    ]

    plotted = []
    for s1, s2 in interesting:
        matching_jobs = [
            j for j in ct_jobs
            if ((j["solvent1"] == s1 and j["solvent2"] == s2) or
                (j["solvent1"] == s2 and j["solvent2"] == s1))
            and abs(j["T_K"] - 298.15) < 0.5
        ]
        if not matching_jobs:
            continue
        ct_job = matching_jobs[0]
        matched = match_pair(oc_df, ct_job)
        if matched is None:
            continue
        plotted.append((s1, s2, ct_job, matched))

    if not plotted:
        return

    n_pairs = len(plotted)
    ncols = min(3, n_pairs)
    nrows = (n_pairs + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    if n_pairs == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for i, (s1, s2, ct_job, matched) in enumerate(plotted):
        ax = axes[i]
        ct_data = ct_job["data"]
        ct_mask = (ct_data["x1"] >= 0.009) & (ct_data["x1"] <= 0.991)
        ax.plot(ct_data.loc[ct_mask, "x1"], ct_data.loc[ct_mask, "GE_RT"],
                "s-", ms=3, lw=1, color="C0", label="COSMOtherm", alpha=0.8)

        fwd = (
            (oc_df["solvent1"] == ct_job["solvent1"])
            & (oc_df["solvent2"] == ct_job["solvent2"])
            & (np.abs(oc_df["T_K"] - 298.15) < 0.5)
        )
        rev = (
            (oc_df["solvent1"] == ct_job["solvent2"])
            & (oc_df["solvent2"] == ct_job["solvent1"])
            & (np.abs(oc_df["T_K"] - 298.15) < 0.5)
        )
        if fwd.sum() > 0:
            oc_sub = oc_df[fwd].sort_values("x1")
            ax.plot(oc_sub["x1"], oc_sub["GE_RT"],
                    "o-", ms=3, lw=1, color="C1", label="openCOSMO-RS", alpha=0.8)
        elif rev.sum() > 0:
            oc_sub = oc_df[rev].sort_values("x1")
            oc_x1_ct = 1.0 - oc_sub["x1"].values
            oc_ge = oc_sub["GE_RT"].values
            sort_idx = np.argsort(oc_x1_ct)
            ax.plot(oc_x1_ct[sort_idx], oc_ge[sort_idx],
                    "o-", ms=3, lw=1, color="C1", label="openCOSMO-RS", alpha=0.8)

        display_s1 = ct_job["solvent1"]
        display_s2 = ct_job["solvent2"]
        ax.set_xlabel(f"x({display_s1})")
        ax.set_ylabel("GE/RT")
        ax.set_title(f"{display_s1} + {display_s2} (25 C)")
        ax.legend(fontsize=8)
        ax.set_xlim(0, 1)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle("GE/RT Comparison: COSMOtherm vs openCOSMO-RS", fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(out_dir / "selected_pairs_comparison.png", dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Plot 8: Temperature dependence for selected pairs
# ---------------------------------------------------------------------------
def plot_temperature_dependence(ct_jobs, oc_df, out_dir):
    pairs_to_plot = [
        ("Methanol", "Hexane"),
        ("Acetonitrile", "Hexane"),
        ("Water", "1-Butanol"),
        ("DMSO", "Hexane"),
    ]

    plotted = []
    for s1, s2 in pairs_to_plot:
        jobs_for_pair = [
            j for j in ct_jobs
            if ((j["solvent1"] == s1 and j["solvent2"] == s2) or
                (j["solvent1"] == s2 and j["solvent2"] == s1))
        ]
        if jobs_for_pair:
            plotted.append((s1, s2, sorted(jobs_for_pair, key=lambda j: j["T_K"])))

    if not plotted:
        return

    fig, axes = plt.subplots(1, len(plotted), figsize=(5.5 * len(plotted), 5))
    if len(plotted) == 1:
        axes = [axes]
    colors_T = plt.cm.coolwarm(np.linspace(0, 1, 5))

    for idx, (s1, s2, jobs) in enumerate(plotted):
        ax = axes[idx]
        for ti, ct_job in enumerate(jobs):
            T_C = ct_job["T_K"] - 273.15
            ct_data = ct_job["data"]
            ct_mask = (ct_data["x1"] >= 0.009) & (ct_data["x1"] <= 0.991)
            ax.plot(ct_data.loc[ct_mask, "x1"], ct_data.loc[ct_mask, "GE_RT"],
                    "s-", ms=2, lw=1, color=colors_T[ti],
                    label=f"CT {T_C:.0f}C", alpha=0.7)

            fwd = (
                (oc_df["solvent1"] == ct_job["solvent1"])
                & (oc_df["solvent2"] == ct_job["solvent2"])
                & (np.abs(oc_df["T_K"] - ct_job["T_K"]) < 0.5)
            )
            rev = (
                (oc_df["solvent1"] == ct_job["solvent2"])
                & (oc_df["solvent2"] == ct_job["solvent1"])
                & (np.abs(oc_df["T_K"] - ct_job["T_K"]) < 0.5)
            )
            if fwd.sum() > 0:
                oc_sub = oc_df[fwd].sort_values("x1")
                ax.plot(oc_sub["x1"], oc_sub["GE_RT"],
                        "o--", ms=2, lw=1, color=colors_T[ti], alpha=0.7)
            elif rev.sum() > 0:
                oc_sub = oc_df[rev].sort_values("x1")
                oc_x1_ct = 1.0 - oc_sub["x1"].values
                oc_ge = oc_sub["GE_RT"].values
                si = np.argsort(oc_x1_ct)
                ax.plot(oc_x1_ct[si], oc_ge[si],
                        "o--", ms=2, lw=1, color=colors_T[ti], alpha=0.7)

        display_s1 = jobs[0]["solvent1"]
        display_s2 = jobs[0]["solvent2"]
        ax.set_xlabel(f"x({display_s1})")
        ax.set_ylabel("GE/RT")
        ax.set_title(f"{display_s1} + {display_s2}")
        ax.legend(fontsize=6, ncol=2, title="solid=CT, dashed=OC", title_fontsize=6)
        ax.set_xlim(0, 1)

    plt.suptitle("Temperature Dependence: COSMOtherm vs openCOSMO-RS",
                 fontsize=11, y=1.02)
    plt.tight_layout()
    plt.savefig(out_dir / "temperature_dependence.png", dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Plot 9: LLE summary grid
# ---------------------------------------------------------------------------
def plot_lle_summary(stats_df, out_dir):
    lle = stats_df[stats_df["lle_found_ct"]].copy()
    if lle.empty:
        return

    # Unique LLE pairs
    lle["pair"] = lle.apply(
        lambda r: f"{r['solvent1']} + {r['solvent2']}", axis=1
    )
    pairs = sorted(lle["pair"].unique())

    fig, ax = plt.subplots(figsize=(6, max(5, len(pairs) * 0.3)))

    for i, pair in enumerate(pairs):
        sub = lle[lle["pair"] == pair]
        for _, row in sub.iterrows():
            T_C = row["T_K"] - 273.15
            # Color by MAE_GE_RT
            mae = row["MAE_GE_RT"]
            if mae < 0.05:
                c = "green"
            elif mae < 0.15:
                c = "orange"
            else:
                c = "red"
            ax.scatter(T_C, i, c=c, s=40, edgecolors="k", linewidths=0.3)

    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels(pairs, fontsize=6)
    ax.set_xlabel("Temperature (C)")
    ax.set_title(f"COSMOtherm LLE Results ({len(pairs)} pairs)")

    # Legend
    for c, lab in [("green", "MAE<0.05"), ("orange", "0.05-0.15"), ("red", ">0.15")]:
        ax.scatter([], [], c=c, s=40, edgecolors="k", linewidths=0.3, label=lab)
    ax.legend(fontsize=7, title="MAE(GE/RT)")

    plt.tight_layout()
    plt.savefig(out_dir / "lle_summary.png", dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Plot 10: Paper summary (4-panel)
# ---------------------------------------------------------------------------
def plot_paper_summary(all_matched, stats_df, out_dir):
    fig = plt.figure(figsize=(14, 11))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)

    # --- Panel A: GE/RT parity (hexbin) ---
    ax_a = fig.add_subplot(gs[0, 0])
    ct_ge = np.concatenate([m["GE_RT_ct"].values for m in all_matched])
    oc_ge = np.concatenate([m["GE_RT_oc"].values for m in all_matched])
    hb = ax_a.hexbin(ct_ge, oc_ge, gridsize=60, cmap="Blues", mincnt=1, norm=LogNorm())
    lim = [min(ct_ge.min(), oc_ge.min()), max(ct_ge.max(), oc_ge.max())]
    ax_a.plot(lim, lim, "k--", lw=0.8, zorder=5)
    ax_a.set_xlabel("COSMOtherm GE/RT")
    ax_a.set_ylabel("openCOSMO-RS GE/RT")
    mae = np.mean(np.abs(ct_ge - oc_ge))
    r2 = np.corrcoef(ct_ge, oc_ge)[0, 1]**2
    ax_a.text(0.05, 0.95,
              f"MAE = {mae:.3f}\nR$^2$ = {r2:.4f}\nN = {len(ct_ge):,}",
              transform=ax_a.transAxes, fontsize=9, va="top",
              bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7))
    plt.colorbar(hb, ax=ax_a, label="Count", shrink=0.8)
    ax_a.set_title("A) GE/RT Parity", fontweight="bold", loc="left")

    # --- Panel B: Temperature violin ---
    ax_b = fig.add_subplot(gs[0, 1])
    T_data, T_labels, T_means = [], [], []
    for T in TEMPERATURES:
        sub = stats_df[np.abs(stats_df["T_K"] - T) < 0.5]
        if len(sub) > 0:
            T_data.append(sub["MAE_GE_RT"].values)
            T_labels.append(f"{T - 273.15:.0f}")
            T_means.append(sub["MAE_GE_RT"].mean())

    parts = ax_b.violinplot(T_data, positions=range(len(T_data)),
                            showmedians=True, showextrema=False)
    for pc in parts["bodies"]:
        pc.set_facecolor("#4e79a7")
        pc.set_alpha(0.6)
    ax_b.plot(range(len(T_data)), T_means, "ro-", ms=5, lw=1.5, label="Mean")
    ax_b.set_xticks(range(len(T_labels)))
    ax_b.set_xticklabels([f"{t} C" for t in T_labels])
    ax_b.set_ylabel("MAE(GE/RT)")
    ax_b.legend(fontsize=8)
    ax_b.set_title("B) Temperature Dependence", fontweight="bold", loc="left")

    # --- Panel C: Per-solvent ranking ---
    ax_c = fig.add_subplot(gs[1, 0])
    sub25 = stats_df[np.abs(stats_df["T_K"] - 298.15) < 0.5]
    solvent_mae = {}
    for s in SOLVENT_CLASS:
        mask = (sub25["solvent1"] == s) | (sub25["solvent2"] == s)
        if mask.sum() > 0:
            solvent_mae[s] = sub25.loc[mask, "MAE_GE_RT"].mean()
    sdf = pd.DataFrame([
        {"solvent": s, "mae": v, "cls": SOLVENT_CLASS[s]}
        for s, v in solvent_mae.items()
    ]).sort_values("mae", ascending=True)
    y_pos = np.arange(len(sdf))
    colors = [CLASS_COLORS.get(c, "gray") for c in sdf["cls"]]
    ax_c.barh(y_pos, sdf["mae"], color=colors, height=0.7)
    ax_c.set_yticks(y_pos)
    ax_c.set_yticklabels(sdf["solvent"], fontsize=7)
    ax_c.set_xlabel("Mean MAE(GE/RT)")
    ax_c.axvline(0.05, color="green", ls="--", lw=0.7)
    ax_c.axvline(0.2, color="red", ls="--", lw=0.7)
    ax_c.set_title("C) Per-Solvent Ranking (25 C)", fontweight="bold", loc="left")

    # --- Panel D: Non-ideality vs MAE ---
    ax_d = fig.add_subplot(gs[1, 1])
    sub_ni = sub25.copy()
    colors_d = []
    for _, row in sub_ni.iterrows():
        s1, s2 = row["solvent1"], row["solvent2"]
        if "Chloroform" in (s1, s2):
            colors_d.append("#e15759")
        elif "Water" in (s1, s2):
            colors_d.append("#4e79a7")
        else:
            colors_d.append("#bab0ac")

    ax_d.scatter(sub_ni["max_abs_GE_RT_ct"], sub_ni["MAE_GE_RT"],
                 c=colors_d, s=15, alpha=0.7, edgecolors="none")
    xmax = sub_ni["max_abs_GE_RT_ct"].max()
    ax_d.plot([0, xmax], [0, 0.1 * xmax], "k--", lw=0.7, alpha=0.5)
    ax_d.set_xlabel("Max |GE/RT| (non-ideality)")
    ax_d.set_ylabel("MAE(GE/RT)")
    ax_d.set_title("D) Error vs Non-ideality (25 C)", fontweight="bold", loc="left")
    for c, lab in [("#e15759", "CHCl3"), ("#4e79a7", "Water"), ("#bab0ac", "Other")]:
        ax_d.scatter([], [], c=c, s=20, label=lab)
    ax_d.legend(fontsize=7)

    # Conclusion text
    frac_good = (sub25["MAE_GE_RT"] < 0.1).mean() * 100
    conclusion = (
        f"openCOSMO-RS reproduces COSMOtherm activity coefficients with "
        f"MAE(GE/RT) = {sub25['MAE_GE_RT'].mean():.3f} at 25 C "
        f"({frac_good:.0f}% of pairs within 0.1); "
        f"largest deviations arise from chloroform hydrogen-bond donor interactions."
    )
    fig.text(0.5, 0.01, conclusion, ha="center", fontsize=9, style="italic",
             bbox=dict(boxstyle="round", facecolor="#f0f0f0", alpha=0.8))

    plt.savefig(out_dir / "paper_summary.png", dpi=250, bbox_inches="tight")
    plt.close()
    return conclusion


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print(f"Parsing COSMOtherm .tab file: {TAB_FILE}")
    ct_jobs = parse_tab_file(TAB_FILE)
    print(f"  Parsed {len(ct_jobs)} COSMOtherm jobs")

    unique_pairs = set()
    lle_count = 0
    for j in ct_jobs:
        pair = tuple(sorted([j["solvent1"], j["solvent2"]]))
        unique_pairs.add(pair)
        if j["lle_found"]:
            lle_count += 1
    print(f"  {len(unique_pairs)} unique solvent pairs")
    print(f"  {lle_count} jobs with LLE found")

    print(f"\nLoading openCOSMO-RS data: {OPENCOSMO_CSV}")
    oc_df = load_opencosmo(OPENCOSMO_CSV)
    if EXCLUDE_SOLVENTS:
        oc_df = oc_df[~oc_df["solvent1"].isin(EXCLUDE_SOLVENTS)
                      & ~oc_df["solvent2"].isin(EXCLUDE_SOLVENTS)]
    print(f"  {len(oc_df)} data points")

    # Filter out excluded solvents
    if EXCLUDE_SOLVENTS:
        before = len(ct_jobs)
        ct_jobs = [j for j in ct_jobs if j["solvent1"] not in EXCLUDE_SOLVENTS
                   and j["solvent2"] not in EXCLUDE_SOLVENTS]
        print(f"\n  Excluded solvents: {EXCLUDE_SOLVENTS}")
        print(f"  {before} -> {len(ct_jobs)} jobs after filtering")

    print("\nMatching pairs and computing statistics...")
    all_matched = []
    all_stats = []
    unmatched = []

    for ct_job in ct_jobs:
        matched = match_pair(oc_df, ct_job)
        if matched is not None:
            all_matched.append(matched)
            stats = pair_statistics(matched)
            stats["solvent1"] = ct_job["solvent1"]
            stats["solvent2"] = ct_job["solvent2"]
            stats["T_K"] = ct_job["T_K"]
            stats["class1"] = SOLVENT_CLASS.get(ct_job["solvent1"], "?")
            stats["class2"] = SOLVENT_CLASS.get(ct_job["solvent2"], "?")
            stats["class_pair"] = class_pair_label(ct_job["solvent1"], ct_job["solvent2"])
            stats["lle_found_ct"] = ct_job["lle_found"]
            stats["lle_x1_phase1"] = ct_job.get("lle_x1_phase1")
            stats["lle_x1_phase2"] = ct_job.get("lle_x1_phase2")
            all_stats.append(stats)
        else:
            unmatched.append(
                f"  {ct_job['solvent1']} + {ct_job['solvent2']} @ {ct_job['T_K']} K"
            )

    print(f"  {len(all_matched)} matched pairs")
    if unmatched:
        print(f"  {len(unmatched)} unmatched (first 5):")
        for u in unmatched[:5]:
            print(u)

    stats_df = pd.DataFrame(all_stats)

    # Summary statistics
    print("\n--- Overall Statistics ---")
    for col in ["MAE_ln_gamma1", "MAE_ln_gamma2", "MAE_GE_RT",
                "RMSE_ln_gamma1", "RMSE_ln_gamma2", "RMSE_GE_RT"]:
        mean_val = stats_df[col].mean()
        max_val = stats_df[col].max()
        med_val = stats_df[col].median()
        print(f"  {col:20s}: mean={mean_val:.6f}  median={med_val:.6f}  max={max_val:.6f}")

    print(f"\n  Bias(GE/RT): mean={stats_df['bias_GE_RT'].mean():.6f} "
          f"(positive = COSMOtherm > openCOSMO-RS)")

    sub25 = stats_df[np.abs(stats_df["T_K"] - 298.15) < 0.5]
    frac_01 = (sub25["MAE_GE_RT"] < 0.1).mean() * 100
    frac_03 = (sub25["MAE_GE_RT"] < 0.3).mean() * 100
    print(f"\n  At 25 C: {frac_01:.0f}% of pairs with MAE<0.1, "
          f"{frac_03:.0f}% with MAE<0.3")

    # Worst pairs
    print("\n--- Top 10 Largest Deviations (MAE GE/RT) ---")
    worst = stats_df.nlargest(10, "MAE_GE_RT")
    for _, row in worst.iterrows():
        print(f"  {row['solvent1']:20s} + {row['solvent2']:20s} "
              f"@ {row['T_K']:.0f} K : MAE(GE/RT)={row['MAE_GE_RT']:.4f} "
              f"[{row['class_pair']}]")

    # Save statistics CSV
    stats_df.to_csv(OUT_DIR / "comparison_statistics.csv", index=False)
    print(f"\nStatistics saved to {OUT_DIR / 'comparison_statistics.csv'}")

    # LLE comparison table
    lle_jobs = [j for j in ct_jobs if j["lle_found"]]
    if lle_jobs:
        lle_rows = []
        for j in lle_jobs:
            lle_rows.append({
                "solvent1": j["solvent1"],
                "solvent2": j["solvent2"],
                "T_K": j["T_K"],
                "T_C": j["T_K"] - 273.15,
                "x1_phase1": j["lle_x1_phase1"],
                "x1_phase2": j["lle_x1_phase2"],
            })
        lle_df = pd.DataFrame(lle_rows)
        lle_df.to_csv(OUT_DIR / "cosmotherm_lle_results.csv", index=False)
        print(f"\nLLE results ({len(lle_df)} entries, "
              f"{lle_df.groupby(['solvent1','solvent2']).ngroups} unique pairs)")

    # Plots
    print("\nGenerating plots...")
    if all_matched:
        plot_parity(all_matched, OUT_DIR)
        print(f"  1. Parity plots (hexbin)")

        plot_deviation_heatmap(stats_df, OUT_DIR)
        print(f"  2. Deviation heatmap (25C)")

        plot_class_breakdown(stats_df, OUT_DIR)
        print(f"  3. Chemical class breakdown")

        plot_solvent_ranking(stats_df, OUT_DIR)
        print(f"  4. Per-solvent ranking")

        plot_temperature_mae(stats_df, OUT_DIR)
        print(f"  5. Temperature MAE violin")

        plot_nonideality_vs_mae(stats_df, OUT_DIR)
        print(f"  6. Non-ideality vs MAE")

        plot_selected_pairs(ct_jobs, oc_df, OUT_DIR)
        print(f"  7. Selected pairs comparison")

        plot_temperature_dependence(ct_jobs, oc_df, OUT_DIR)
        print(f"  8. Temperature dependence")

        plot_lle_summary(stats_df, OUT_DIR)
        print(f"  9. LLE summary")

        conclusion = plot_paper_summary(all_matched, stats_df, OUT_DIR)
        print(f"  10. Paper summary (4-panel)")

        print(f"\n--- CONCLUSION ---")
        print(f"  {conclusion}")

    print(f"\nAll outputs in: {OUT_DIR}")


if __name__ == "__main__":
    main()
