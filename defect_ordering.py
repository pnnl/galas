# This material was prepared as an account of work sponsored by an agency of the 
# United States Government. Neither the United States Government nor the United 
# States Department of Energy, nor Battelle, nor any of their employees, nor any 
# jurisdiction or organization that has cooperated in the development of these 
# materials, makes any warranty, express or implied, or assumes any legal 
# liability or responsibility for the accuracy, completeness, or usefulness or 
# any information, apparatus, product, software, or process disclosed, or 
# represents that its use would not infringe privately owned rights. Reference 
# herein to any specific commercial product, process, or service by trade name, 
# trademark, manufacturer, or otherwise does not necessarily constitute or imply 
# its endorsement, recommendation, or favoring by the United States Government 
# or any agency thereof, or Battelle Memorial Institute. The views and opinions 
# of authors expressed herein do not necessarily state or reflect those of the 
# United States Government or any agency thereof.
#                    PACIFIC NORTHWEST NATIONAL LABORATORY
#                               operated by
#                                BATTELLE
#                                for the
#                      UNITED STATES DEPARTMENT OF ENERGY
#                       under Contract DE-AC05-76RL01830

import os
import os.path as op
import numpy as np
import pandas as pd
import pickle
import logging
import argparse
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default='./data/',
                    help='Path to data folder')
parser.add_argument('--lattice_constant', type=float, default=4.0559,
                    help='Lattice constant for normalizing distances')
parser.add_argument('--rdf_rmax', type=float, default=60.0,
                    help='Maximum radius for RDF computation (Angstrom)')
parser.add_argument('--rdf_bins', type=int, default=300,
                    help='Number of bins for RDF histogram')
parser.add_argument('--sq_qmax', type=float, default=2.0,
                    help='Maximum q for structure factor (1/Angstrom)')
parser.add_argument('--sq_nq', type=int, default=500,
                    help='Number of q points for structure factor')
parser.add_argument('--sro_shells', type=int, default=5,
                    help='Number of neighbor shells for Warren-Cowley parameters')
parser.add_argument('--sro_shell_width', type=float, default=0.0,
                    help='Shell width in Angstrom (0 = auto from lattice constant)')
args = parser.parse_args()

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

def compute_rdf(positions, rmax, nbins, box_volume):
    """Compute radial distribution function for a set of point positions."""
    N = len(positions)
    if N < 2:
        return np.zeros(nbins), np.linspace(0, rmax, nbins)
    
    dists = pdist(positions)
    density = N / box_volume
    
    bin_edges = np.linspace(0, rmax, nbins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    dr = bin_edges[1] - bin_edges[0]
    
    hist, _ = np.histogram(dists, bins=bin_edges)
    
    # Normalize: g(r) = hist / (N * density * shell_volume)
    shell_volumes = (4.0 / 3.0) * np.pi * (bin_edges[1:]**3 - bin_edges[:-1]**3)
    # Each pair counted once in pdist; total pairs = N*(N-1)/2
    # g(r) normalization: hist / (0.5 * N * density * shell_volume)
    gr = hist / (0.5 * N * density * shell_volumes)
    
    return gr, bin_centers


def compute_cross_rdf(positions_a, positions_b, rmax, nbins, box_volume):
    """Compute cross-RDF between two sets of positions (type A around type B)."""
    Na = len(positions_a)
    Nb = len(positions_b)
    if Na == 0 or Nb == 0:
        return np.zeros(nbins), np.linspace(0, rmax, nbins)
    
    # Compute all distances between A and B
    tree_b = cKDTree(positions_b)
    dists = []
    for pos in positions_a:
        d = tree_b.query_ball_point(pos, rmax, return_length=False)
        for idx in d:
            dist = np.linalg.norm(pos - positions_b[idx])
            if dist > 1e-10:  # exclude self if same set
                dists.append(dist)
    dists = np.array(dists)
    
    density_b = Nb / box_volume
    bin_edges = np.linspace(0, rmax, nbins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    hist, _ = np.histogram(dists, bins=bin_edges)
    shell_volumes = (4.0 / 3.0) * np.pi * (bin_edges[1:]**3 - bin_edges[:-1]**3)
    
    # Normalize: g_AB(r) = hist / (Na * density_b * shell_volume)
    gr = hist / (Na * density_b * shell_volumes + 1e-30)
    
    return gr, bin_centers


def compute_nn_distribution(positions, n_neighbors=1):
    """Compute nearest-neighbor distance distribution."""
    if len(positions) < 2:
        return np.array([])
    tree = cKDTree(positions)
    dists, _ = tree.query(positions, k=n_neighbors + 1)  # +1 for self
    nn_dists = dists[:, 1:]  # exclude self-distance (0)
    return nn_dists[:, 0]  # first nearest neighbor


def compute_structure_factor(positions, qmax, nq):
    """Compute radially-averaged structure factor S(q) via Debye formula."""
    N = len(positions)
    if N < 2:
        return np.zeros(nq), np.linspace(0.01, qmax, nq)
    
    q_values = np.linspace(0.01, qmax, nq)
    dists = pdist(positions)
    
    S_q = np.ones(nq)
    for i, q in enumerate(q_values):
        qr = q * dists
        # Debye formula: S(q) = 1 + (2/N) * sum(sin(qr)/(qr))
        sinc_term = np.where(qr > 1e-10, np.sin(qr) / qr, 1.0)
        S_q[i] = 1.0 + (2.0 / N) * np.sum(sinc_term)
    
    return S_q, q_values


def compute_warren_cowley(centroids_by_type, type_names, n_shells, shell_width):
    """Compute Warren-Cowley SRO parameters between defect types."""
    # Combine all centroids with type labels
    all_positions = []
    all_types = []
    for t_idx, (name, positions) in enumerate(zip(type_names, centroids_by_type)):
        for pos in positions:
            all_positions.append(pos)
            all_types.append(t_idx)
    
    all_positions = np.array(all_positions)
    all_types = np.array(all_types)
    N = len(all_positions)
    n_types = len(type_names)
    
    if N < 2:
        return np.zeros((n_shells, n_types, n_types)), np.zeros(n_shells)
    
    # Global concentrations
    concentrations = np.array([np.sum(all_types == t) / N for t in range(n_types)])
    
    # Build KD-tree
    tree = cKDTree(all_positions)
    
    # Define shell boundaries
    shell_boundaries = np.array([(n * shell_width, (n + 1) * shell_width) for n in range(n_shells)])
    shell_centers = 0.5 * (shell_boundaries[:, 0] + shell_boundaries[:, 1])
    
    # Compute alpha for each shell and type pair
    alpha = np.zeros((n_shells, n_types, n_types))
    
    for shell_idx in range(n_shells):
        r_min, r_max = shell_boundaries[shell_idx]
        
        for i_type in range(n_types):
            i_mask = all_types == i_type
            i_indices = np.where(i_mask)[0]
            Ni = len(i_indices)
            if Ni == 0:
                continue
            
            for j_type in range(n_types):
                if concentrations[j_type] < 1e-10:
                    continue
                
                # Count j-type neighbors in this shell for all i-type atoms
                total_neighbors = 0
                j_count = 0
                
                for idx in i_indices:
                    neighbors_in_shell = tree.query_ball_point(
                        all_positions[idx], r_max)
                    for n_idx in neighbors_in_shell:
                        if n_idx == idx:
                            continue
                        d = np.linalg.norm(all_positions[idx] - all_positions[n_idx])
                        if d >= r_min and d < r_max:
                            total_neighbors += 1
                            if all_types[n_idx] == j_type:
                                j_count += 1
                
                if total_neighbors > 0:
                    p_ij = j_count / total_neighbors
                    alpha[shell_idx, i_type, j_type] = 1.0 - p_ij / concentrations[j_type]
    
    return alpha, shell_centers


def hertz_distribution(r, density):
    """Analytical nearest-neighbor distribution for a random Poisson process."""
    return 4.0 * np.pi * density * r**2 * np.exp(-(4.0 / 3.0) * np.pi * density * r**3)


# ------------------------------------------------------------------
# Main analysis
# ------------------------------------------------------------------

# Load component info to identify defect types
data_path = args.path
components_csv_dir = op.join(data_path, 'components', 'csvs')
graphs_csv_dir = op.join(data_path, 'graphs', 'csvs')

# Find all component CSVs
comp_csvs = sorted([f for f in os.listdir(components_csv_dir) if f.endswith('.csv')])
if not comp_csvs:
    logging.error('No component CSVs found')
    raise SystemExit(1)

# Use the first frame (or only frame)
comp_csv = comp_csvs[0]
frame_name = comp_csv

logging.info(f'Analyzing defect ordering for frame: {frame_name}')

# Load component info and graph CSV
info = pd.read_csv(op.join(components_csv_dir, comp_csv))
graph_csv_name = comp_csv  # same name in graphs/csvs/
df = pd.read_csv(op.join(graphs_csv_dir, graph_csv_name))

# Load defect metadata to get atom indices per component
meta_files = [f for f in os.listdir(op.join(data_path, 'components')) if f.endswith('.meta.npz')]
if not meta_files:
    logging.error('No .meta.npz files found in components/')
    raise SystemExit(1)

meta = np.load(op.join(data_path, 'components', meta_files[0]))
defect_indices = meta['defect_indices']
labels = meta['labels']

# Estimate box volume from coordinate ranges
x_range = df['x'].max() - df['x'].min()
y_range = df['y'].max() - df['y'].min()
z_range = df['z'].max() - df['z'].min()
box_volume = x_range * y_range * z_range

# Define defect type signatures (FCC)
defect_signatures = {
    'Mono-vacancy': (12, 24),
    'Di-vacancy': (18, 40),
    'Tri-vacancy (linear)': (24, 60),
    'Tri-vacancy (triangular)': (24, 54),
}

# Classify components and compute centroids
defect_centroids = {}  # type_name -> list of centroid positions
unclassified_centroids = []

for _, row in info.iterrows():
    comp_id = int(row['Component'])
    n_nodes = int(row['Nodes'])
    n_edges = int(row['Edges'])
    
    # Skip GBS (component 0, the largest)
    if comp_id == 0:
        continue
    
    # Get atom indices for this component
    local_mask = labels == comp_id
    original_indices = defect_indices[local_mask]
    
    if len(original_indices) == 0:
        continue
    
    # Compute centroid
    centroid = df.iloc[original_indices][['x', 'y', 'z']].mean().values
    
    # Classify
    classified = False
    for defect_name, (sig_nodes, sig_edges) in defect_signatures.items():
        if n_nodes == sig_nodes and n_edges == sig_edges:
            if defect_name not in defect_centroids:
                defect_centroids[defect_name] = []
            defect_centroids[defect_name].append(centroid)
            classified = True
            break
    
    if not classified:
        unclassified_centroids.append(centroid)

# Add unclassified as a category if present
if unclassified_centroids:
    defect_centroids['Unclassified'] = unclassified_centroids

# Report defect counts
print("=" * 70)
print("DEFECT ORDERING ANALYSIS")
print("=" * 70)
print(f"\nFrame: {frame_name}")
print(f"Box volume: {box_volume:.1f} Å³")
print(f"\nDefect counts by type:")
all_centroids = []
for name, centroids in defect_centroids.items():
    print(f"  {name}: {len(centroids)}")
    all_centroids.extend(centroids)
all_centroids = np.array(all_centroids)
print(f"  Total (excl. GBS): {len(all_centroids)}")

# Create output directory
fig_dir = op.join(data_path, 'figures')
os.makedirs(fig_dir, exist_ok=True)
report_path = op.join(data_path, 'ordering_report.csv')

# ------------------------------------------------------------------
# 1. RDF of defect centers (grouped by type)
# ------------------------------------------------------------------
print("\n" + "-" * 70)
print("1. RADIAL DISTRIBUTION FUNCTION (RDF)")
print("-" * 70)

fig, ax = plt.subplots(figsize=(10, 6))

rdf_results = {}
for name, centroids in defect_centroids.items():
    positions = np.array(centroids)
    if len(positions) < 2:
        print(f"  {name}: too few defects for RDF (n={len(positions)})")
        continue
    
    gr, r = compute_rdf(positions, args.rdf_rmax, args.rdf_bins, box_volume)
    rdf_results[name] = (gr, r)
    ax.plot(r, gr, label=f'{name} (n={len(positions)})', linewidth=1.5)
    
    # Report first peak
    peak_idx = np.argmax(gr[5:]) + 5  # skip first few noisy bins
    print(f"  {name}: first peak at r = {r[peak_idx]:.2f} Å, g(r) = {gr[peak_idx]:.2f}")

# Also compute all-defect RDF
if len(all_centroids) >= 2:
    gr_all, r_all = compute_rdf(all_centroids, args.rdf_rmax, args.rdf_bins, box_volume)
    ax.plot(r_all, gr_all, 'k--', label=f'All defects (n={len(all_centroids)})', linewidth=2)
    peak_idx = np.argmax(gr_all[5:]) + 5
    print(f"  All defects: first peak at r = {r_all[peak_idx]:.2f} Å, g(r) = {gr_all[peak_idx]:.2f}")

ax.axhline(1.0, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('r (Å)')
ax.set_ylabel('g(r)')
ax.set_title('Radial Distribution Function of Defect Centers')
ax.legend()
ax.set_xlim(0, args.rdf_rmax)
plt.tight_layout()
plt.savefig(op.join(fig_dir, 'defect_rdf.png'), dpi=150)
plt.close()
print(f"\n  Figure saved: figures/defect_rdf.png")

# ------------------------------------------------------------------
# 2. Nearest-Neighbor Distance Distribution
# ------------------------------------------------------------------
print("\n" + "-" * 70)
print("2. NEAREST-NEIGHBOR DISTANCE DISTRIBUTION")
print("-" * 70)

fig, ax = plt.subplots(figsize=(10, 6))

nn_results = {}
for name, centroids in defect_centroids.items():
    positions = np.array(centroids)
    if len(positions) < 2:
        continue
    
    nn_dists = compute_nn_distribution(positions)
    nn_results[name] = nn_dists
    
    ax.hist(nn_dists, bins=80, alpha=0.5, density=True,
            label=f'{name} (mean={nn_dists.mean():.1f} Å)')
    
    print(f"  {name}: mean NN dist = {nn_dists.mean():.2f} Å, "
          f"std = {nn_dists.std():.2f} Å, min = {nn_dists.min():.2f} Å")

# Plot Hertz distribution for comparison (random expectation)
if len(all_centroids) >= 2:
    density_all = len(all_centroids) / box_volume
    r_hertz = np.linspace(0.1, args.rdf_rmax * 0.5, 200)
    p_hertz = hertz_distribution(r_hertz, density_all)
    ax.plot(r_hertz, p_hertz, 'k--', linewidth=2, label='Random (Hertz)')
    
    nn_all = compute_nn_distribution(all_centroids)
    ax.hist(nn_all, bins=80, alpha=0.3, density=True, color='gray',
            label=f'All defects (mean={nn_all.mean():.1f} Å)')
    print(f"  All defects: mean NN dist = {nn_all.mean():.2f} Å, "
          f"std = {nn_all.std():.2f} Å")
    print(f"  Random expectation (Hertz): mean = {(0.5541 / density_all**(1/3)):.2f} Å")

ax.set_xlabel('Nearest-Neighbor Distance (Å)')
ax.set_ylabel('Probability Density')
ax.set_title('Nearest-Neighbor Distance Distribution of Defect Centers')
ax.legend()
plt.tight_layout()
plt.savefig(op.join(fig_dir, 'defect_nn_distribution.png'), dpi=150)
plt.close()
print(f"\n  Figure saved: figures/defect_nn_distribution.png")

# ------------------------------------------------------------------
# 3. Warren-Cowley SRO Parameters
# ------------------------------------------------------------------
print("\n" + "-" * 70)
print("3. WARREN-COWLEY SHORT-RANGE ORDER PARAMETERS")
print("-" * 70)

# Determine shell width
shell_width = args.sro_shell_width if args.sro_shell_width > 0 else args.lattice_constant * 1.5
n_shells = args.sro_shells

type_names = list(defect_centroids.keys())
centroids_by_type = [np.array(defect_centroids[name]) for name in type_names]

# Only compute if multiple types exist
if len(type_names) >= 2 and all(len(c) >= 1 for c in centroids_by_type):
    alpha, shell_centers = compute_warren_cowley(
        centroids_by_type, type_names, n_shells, shell_width)
    
    print(f"\n  Shell width: {shell_width:.2f} Å")
    print(f"  Number of shells: {n_shells}")
    print(f"\n  α_ij values (shell 1, r = {shell_centers[0]:.1f} Å):")
    print(f"  {'':20s}", end='')
    for name in type_names:
        print(f"{name[:12]:>14s}", end='')
    print()
    
    for i, name_i in enumerate(type_names):
        print(f"  {name_i:20s}", end='')
        for j, name_j in enumerate(type_names):
            print(f"{alpha[0, i, j]:14.4f}", end='')
        print()
    
    print(f"\n  Interpretation:")
    print(f"    α < 0 → unlike defects attract (ordering)")
    print(f"    α = 0 → random distribution")
    print(f"    α > 0 → like defects cluster (segregation)")
    
    # Plot alpha vs shell for each pair
    fig, ax = plt.subplots(figsize=(10, 6))
    for i in range(len(type_names)):
        for j in range(i, len(type_names)):
            label = f'{type_names[i][:8]}-{type_names[j][:8]}'
            ax.plot(shell_centers, alpha[:, i, j], 'o-', label=label)
    
    ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Shell distance (Å)')
    ax.set_ylabel('Warren-Cowley α')
    ax.set_title('Warren-Cowley Short-Range Order Parameters')
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(op.join(fig_dir, 'defect_warren_cowley.png'), dpi=150)
    plt.close()
    print(f"\n  Figure saved: figures/defect_warren_cowley.png")
else:
    print("  Skipped: need at least 2 defect types for Warren-Cowley analysis")
    alpha = None

# ------------------------------------------------------------------
# 4. Structure Factor S(q) / Fourier Analysis
# ------------------------------------------------------------------
print("\n" + "-" * 70)
print("4. STRUCTURE FACTOR S(q)")
print("-" * 70)

fig, ax = plt.subplots(figsize=(10, 6))

sq_results = {}
for name, centroids in defect_centroids.items():
    positions = np.array(centroids)
    if len(positions) < 10:
        print(f"  {name}: too few defects for S(q) (n={len(positions)})")
        continue
    
    # For large sets, subsample to keep computation tractable
    if len(positions) > 5000:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(positions), 5000, replace=False)
        positions_sub = positions[idx]
        print(f"  {name}: subsampled to 5000 for S(q) computation")
    else:
        positions_sub = positions
    
    Sq, q = compute_structure_factor(positions_sub, args.sq_qmax, args.sq_nq)
    sq_results[name] = (Sq, q)
    ax.plot(q, Sq, label=f'{name} (n={len(positions)})', linewidth=1.5)

# All defects combined
if len(all_centroids) >= 10:
    if len(all_centroids) > 5000:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(all_centroids), 5000, replace=False)
        all_sub = all_centroids[idx]
    else:
        all_sub = all_centroids
    Sq_all, q_all = compute_structure_factor(all_sub, args.sq_qmax, args.sq_nq)
    ax.plot(q_all, Sq_all, 'k--', label='All defects', linewidth=2)
    
    # Find peaks in S(q)
    from scipy.signal import find_peaks
    peaks, properties = find_peaks(Sq_all, height=1.5, distance=20)
    if len(peaks) > 0:
        print(f"\n  Peaks in S(q) for all defects:")
        for p in peaks[:5]:
            real_space = 2 * np.pi / q_all[p]
            print(f"    q = {q_all[p]:.4f} Å⁻¹  →  d = {real_space:.1f} Å  "
                  f"(S = {Sq_all[p]:.2f})")
        print(f"\n  Interpretation: peaks indicate preferred periodic spacing")
        print(f"    Strong sharp peaks → long-range defect superlattice")
        print(f"    Broad peaks → short-range correlations only")
    else:
        print(f"\n  No significant peaks in S(q) → defects lack long-range periodicity")

ax.axhline(1.0, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('q (Å⁻¹)')
ax.set_ylabel('S(q)')
ax.set_title('Structure Factor of Defect Centers')
ax.legend()
ax.set_xlim(0, args.sq_qmax)
plt.tight_layout()
plt.savefig(op.join(fig_dir, 'defect_structure_factor.png'), dpi=150)
plt.close()
print(f"\n  Figure saved: figures/defect_structure_factor.png")

# ------------------------------------------------------------------
# Save numerical results
# ------------------------------------------------------------------
results_dir = op.join(data_path, 'ordering')
os.makedirs(results_dir, exist_ok=True)

# Save RDF data
for name, (gr, r) in rdf_results.items():
    safe_name = name.replace(' ', '_').replace('(', '').replace(')', '')
    pd.DataFrame({'r_angstrom': r, 'g_r': gr}).to_csv(
        op.join(results_dir, f'rdf_{safe_name}.csv'), index=False)

# Save NN distributions
for name, nn_dists in nn_results.items():
    safe_name = name.replace(' ', '_').replace('(', '').replace(')', '')
    pd.DataFrame({'nn_distance_angstrom': nn_dists}).to_csv(
        op.join(results_dir, f'nn_dist_{safe_name}.csv'), index=False)

# Save S(q) data
for name, (Sq, q) in sq_results.items():
    safe_name = name.replace(' ', '_').replace('(', '').replace(')', '')
    pd.DataFrame({'q_inv_angstrom': q, 'S_q': Sq}).to_csv(
        op.join(results_dir, f'sq_{safe_name}.csv'), index=False)

# Save Warren-Cowley parameters
if alpha is not None:
    rows = []
    for shell_idx in range(n_shells):
        for i, name_i in enumerate(type_names):
            for j, name_j in enumerate(type_names):
                rows.append({
                    'shell': shell_idx + 1,
                    'shell_center_angstrom': shell_centers[shell_idx],
                    'type_i': name_i,
                    'type_j': name_j,
                    'alpha': alpha[shell_idx, i, j]
                })
    pd.DataFrame(rows).to_csv(op.join(results_dir, 'warren_cowley.csv'), index=False)

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
print(f"\nNumerical results saved to: {results_dir}/")
print(f"Figures saved to: {fig_dir}/")
