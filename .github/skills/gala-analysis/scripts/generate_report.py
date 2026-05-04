"""
GALAS Report Generator
Produces standard defect analysis figures and summary tables from GALAS pipeline output.

Usage:
    python generate_report.py --path ./data/ [--single_frame dump.0.txt] [--output ./data/figures/]
"""

import os
import os.path as op
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def ensure_dir(path):
    if not op.isdir(path):
        os.makedirs(path, exist_ok=True)


def load_all_component_data(data_path):
    csv_path = op.join(data_path, "all_component_data.csv")
    if not op.isfile(csv_path):
        return None
    return pd.read_csv(csv_path)


def load_component_csvs(data_path):
    """Load per-frame component CSVs into a dict keyed by frame name."""
    comp_csv_dir = op.join(data_path, "components", "csvs")
    if not op.isdir(comp_csv_dir):
        return {}
    results = {}
    for f in sorted(os.listdir(comp_csv_dir)):
        if f.endswith(".csv"):
            results[f] = pd.read_csv(op.join(comp_csv_dir, f))
    return results


# ---------- Defect classification helpers ----------

DEFECT_SIGNATURES_FCC = {
    "Mono-vacancy": (12, 24),
    "Di-vacancy": (18, 40),
    "Tri-vacancy (linear)": (24, 60),
    "Tri-vacancy (triangular)": (24, 54),
}

ADJACENT_MONOVAC_NODES = {
    "Adjacent mono-vac (square face, 4 shared)": 20,
    "Adjacent mono-vac (triangle face, 3 shared)": 21,
    "Adjacent mono-vac (face edge, 2 shared)": 22,
    "Adjacent mono-vac (face corner, 1 shared)": 23,
}


def classify_components(comp_df):
    """Classify components by matching (Nodes, Edges) signatures.
    Returns a dict of {defect_type: count} and a list of per-component records."""
    counts = {k: 0 for k in DEFECT_SIGNATURES_FCC}
    counts.update({k: 0 for k in ADJACENT_MONOVAC_NODES})
    counts["Unclassified"] = 0

    records = []
    for _, row in comp_df.iterrows():
        n, e = int(row["Nodes"]), int(row["Edges"])
        classified = False
        for name, (sig_n, sig_e) in DEFECT_SIGNATURES_FCC.items():
            if n == sig_n and e == sig_e:
                counts[name] += 1
                records.append({"Component": int(row["Component"]), "Nodes": n, "Edges": e, "Type": name})
                classified = True
                break
        if not classified:
            for name, sig_n in ADJACENT_MONOVAC_NODES.items():
                if n == sig_n:
                    counts[name] += 1
                    records.append({"Component": int(row["Component"]), "Nodes": n, "Edges": e, "Type": name})
                    classified = True
                    break
        if not classified:
            counts["Unclassified"] += 1
            records.append({"Component": int(row["Component"]), "Nodes": n, "Edges": e, "Type": "Unclassified"})
    return counts, records


# ---------- Plotting functions ----------

def plot_component_size_distribution(comp_df, output_dir, frame_label=""):
    """Histogram of component sizes (node counts), excluding the largest (GBS)."""
    sizes = comp_df["Nodes"].values
    if len(sizes) < 2:
        return
    sizes = sizes[1:]  # skip GBS (index 0, largest)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(sizes, bins=range(1, max(sizes) + 2), edgecolor="black", alpha=0.7)
    ax.set_xlabel("Component Size (number of atoms)")
    ax.set_ylabel("Count")
    ax.set_yscale("log")
    ax.set_title(f"In-Grain Defect Component Size Distribution{' — ' + frame_label if frame_label else ''}")
    plt.tight_layout()
    fig.savefig(op.join(output_dir, f"component_size_dist{'_' + frame_label if frame_label else ''}.png"), dpi=150)
    plt.close(fig)


def plot_defect_type_bar(counts, output_dir, frame_label=""):
    """Bar chart of defect type counts."""
    labels = [k for k, v in counts.items() if v > 0]
    values = [counts[k] for k in labels]
    if not labels:
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.barh(labels, values, color="steelblue", edgecolor="black")
    ax.set_xlabel("Count")
    ax.set_title(f"Defect Type Breakdown{' — ' + frame_label if frame_label else ''}")
    plt.tight_layout()
    fig.savefig(op.join(output_dir, f"defect_type_bar{'_' + frame_label if frame_label else ''}.png"), dpi=150)
    plt.close(fig)


def plot_trajectory_gbs_vs_defects(all_comp_df, output_dir):
    """Delta-v plot: change in GBS and total in-grain defect vertices relative to first frame."""
    if all_comp_df is None or len(all_comp_df) < 2:
        return

    gbs_sizes = all_comp_df["largest_grain"].values
    total_defect_nodes = all_comp_df["nodes"].values - gbs_sizes
    frames = np.arange(len(all_comp_df))

    delta_gbs = gbs_sizes - gbs_sizes[0]
    delta_defects = total_defect_nodes - total_defect_nodes[0]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(frames, delta_gbs, label=r"$\Delta v$ GBS ($C_{GBS}$)", color="tab:blue")
    ax.plot(frames, delta_defects, label=r"$\Delta v$ In-grain defects ($\Sigma C_k$)", color="tab:orange")
    ax.set_xlabel("Frame")
    ax.set_ylabel(r"$\Delta$ Vertices")
    ax.set_title("Evolution of GBS and In-Grain Defect Sizes")
    ax.legend()
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.5)
    plt.tight_layout()
    fig.savefig(op.join(output_dir, "trajectory_delta_v.png"), dpi=150)
    plt.close(fig)


def plot_trajectory_defect_counts(comp_csvs, output_dir):
    """Line plot of each defect type count across frames."""
    if not comp_csvs:
        return

    frame_names = sorted(comp_csvs.keys())
    type_counts_over_time = {k: [] for k in list(DEFECT_SIGNATURES_FCC) + list(ADJACENT_MONOVAC_NODES) + ["Unclassified"]}

    for fname in frame_names:
        counts, _ = classify_components(comp_csvs[fname])
        for dtype in type_counts_over_time:
            type_counts_over_time[dtype].append(counts.get(dtype, 0))

    frames = np.arange(len(frame_names))
    fig, ax = plt.subplots(figsize=(10, 6))
    for dtype, vals in type_counts_over_time.items():
        if max(vals) > 0:
            ax.plot(frames, vals, label=dtype)
    ax.set_xlabel("Frame")
    ax.set_ylabel("Count")
    ax.set_title("Defect Type Counts Over Trajectory")
    ax.legend(fontsize=8, loc="upper right")
    plt.tight_layout()
    fig.savefig(op.join(output_dir, "trajectory_defect_counts.png"), dpi=150)
    plt.close(fig)


def plot_trajectory_total_components(all_comp_df, output_dir):
    """Line plot of total number of distinct defect components per frame."""
    if all_comp_df is None or len(all_comp_df) < 2:
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    frames = np.arange(len(all_comp_df))
    ax.plot(frames, all_comp_df["components"].values, color="tab:green")
    ax.set_xlabel("Frame")
    ax.set_ylabel("Number of Components")
    ax.set_title("Total Defect Components Over Trajectory")
    plt.tight_layout()
    fig.savefig(op.join(output_dir, "trajectory_total_components.png"), dpi=150)
    plt.close(fig)


# ---------- Summary table ----------

def print_single_frame_summary(comp_df, frame_label=""):
    """Print a summary table for a single frame."""
    counts, records = classify_components(comp_df)
    total_defect_atoms = comp_df["Nodes"].sum()
    gbs_atoms = comp_df.iloc[0]["Nodes"] if len(comp_df) > 0 else 0
    n_components = len(comp_df) - 1  # exclude GBS

    print(f"\n{'='*60}")
    print(f"  GALAS Defect Summary{' — ' + frame_label if frame_label else ''}")
    print(f"{'='*60}")
    print(f"  Total defect atoms:           {total_defect_atoms:>10,}")
    print(f"  GBS atoms (C_GBS):            {gbs_atoms:>10,}")
    print(f"  In-grain defect components:   {n_components:>10,}")
    print(f"{'─'*60}")
    print(f"  {'Defect Type':<45} {'Count':>8}")
    print(f"{'─'*60}")
    for dtype, count in counts.items():
        if count > 0:
            print(f"  {dtype:<45} {count:>8,}")
    print(f"{'='*60}\n")
    return counts, records


# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(description="Generate GALAS defect analysis report")
    parser.add_argument("--path", type=str, default="./data/", help="Path to data folder")
    parser.add_argument("--single_frame", type=str, default=None,
                        help="Analyze only this frame (filename in components/csvs/, e.g. dump.0.txt.csv)")
    parser.add_argument("--output", type=str, default=None, help="Output directory for figures")
    args = parser.parse_args()

    output_dir = args.output or op.join(args.path, "figures")
    ensure_dir(output_dir)

    if args.single_frame:
        # Single-frame analysis
        csv_path = op.join(args.path, "components", "csvs", args.single_frame)
        if not op.isfile(csv_path):
            print(f"ERROR: Component CSV not found: {csv_path}")
            return
        comp_df = pd.read_csv(csv_path)
        label = args.single_frame.replace(".csv", "")
        counts, records = print_single_frame_summary(comp_df, label)
        plot_component_size_distribution(comp_df, output_dir, label)
        plot_defect_type_bar(counts, output_dir, label)

        # Write individual defect log
        rec_df = pd.DataFrame(records)
        log_path = op.join(output_dir, f"defect_log_{label}.csv")
        rec_df.to_csv(log_path, index=False)
        print(f"  Individual defect log saved to: {log_path}")
    else:
        # Full trajectory analysis
        all_comp_df = load_all_component_data(args.path)
        comp_csvs = load_component_csvs(args.path)

        if all_comp_df is not None:
            print(f"\nTrajectory: {len(all_comp_df)} frames found\n")

        # Per-frame summaries
        for fname, cdf in comp_csvs.items():
            label = fname.replace(".csv", "")
            counts, records = print_single_frame_summary(cdf, label)
            plot_component_size_distribution(cdf, output_dir, label)
            plot_defect_type_bar(counts, output_dir, label)
            rec_df = pd.DataFrame(records)
            rec_df.to_csv(op.join(output_dir, f"defect_log_{label}.csv"), index=False)

        # Trajectory-level plots
        plot_trajectory_gbs_vs_defects(all_comp_df, output_dir)
        plot_trajectory_defect_counts(comp_csvs, output_dir)
        plot_trajectory_total_components(all_comp_df, output_dir)

        print(f"\n  All figures saved to: {output_dir}/")


if __name__ == "__main__":
    main()
