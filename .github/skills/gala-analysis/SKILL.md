---
name: gala-analysis
description: "Run GALAS (Graph Analytics for Large Atomistic Simulations) defect analysis on LAMMPS dump files. Use when: analyzing defects in atomistic simulation data, identifying vacancies or grain boundaries, tracking defect evolution over MD trajectories, generating defect summary reports with counts and types, producing plots of defect metrics over time."
argument-hint: "Provide the path to your data directory containing LAMMPS dumps, lattice type (e.g. FCC), and lattice constant."
---

# GALA Defect Analysis

Run the full GALAS pipeline on LAMMPS dump files to identify, classify, and track structural defects in atomistic simulations.

## When to Use

- Analyze a single atomic structure snapshot for defect counts and types
- Process a trajectory of snapshots to track defect evolution over time
- Identify mono-vacancies, di-vacancies, tri-vacancies, vacancy clusters, and grain boundary superstructure (GBS)
- Generate summary statistics and plots of defect distributions

## Prerequisites

The following Python packages must be available in the active environment:

- `ovito` (atomic structure I/O, PTM/CNA classification, neighbor finding)
- `networkx` (graph construction and component analysis)
- `numpy`, `pandas` (data manipulation)
- `matplotlib` (plotting)
- `tqdm` (progress bars)

Also required: LAMMPS dump files placed in `<data_path>/dumps/`.

## Lattice Defaults

| Lattice | Lattice Constant (Å) | Ideal Neighbors | Cutoff Formula | Example Material |
|---|---|---|---|---|
| **FCC** | 4.0559 | 12 | $(a + a/\sqrt{2})/2$ | Al at 300 K |
| **BCC** | 2.8665 | 8 | $a \cdot \sqrt{3}/2$ | Fe at 300 K |
| **HCP** | 3.2094 (a), 5.2105 (c) | 12 | $(a + a/\sqrt{2})/2$ | Ti at 300 K |

When the user specifies a lattice type, use the corresponding defaults unless they provide explicit overrides. For BCC, the ideal neighbor count is 8 (first shell) or 14 (first + second shell) — confirm with the user which convention to use.

## Key Concepts

| Term | Meaning |
|------|---------|
| **PTM / CNA / a-CNA** | Methods to classify atoms by local structure (FCC, BCC, HCP, etc.). Atoms not matching the ideal lattice are "defective." |
| **Material graph** | Graph where each vertex is an atom and edges connect nearest neighbors within a cutoff radius. |
| **Connected component** | A disjoint subgraph of defective atoms. Each component represents one defect region. |
| **C_GBS** | The largest connected component — the grain boundary superstructure. |
| **C_k** | All remaining components — individual in-grain defects (vacancies, clusters, dislocations, stacking faults). |
| **Isomorphism check** | Structural comparison against template defect graphs to classify component type. |

### Defect Signatures (FCC)

| Defect Type | Vertices | Edges | Notes |
|---|---|---|---|
| Mono-vacancy | 12 | 24 | Each vertex has degree 4 in subgraph |
| Di-vacancy | 18 | 40 | Two adjacent vacant sites |
| Tri-vacancy (linear) | 24 | 60 | Three voids in a line |
| Tri-vacancy (triangular) | 24 | 54 | Three voids in a triangle |
| Adjacent mono-vac (square face) | 20 | — | 4 shared atoms |
| Adjacent mono-vac (triangle face) | 21 | — | 3 shared atoms |
| Adjacent mono-vac (face edge) | 22 | — | 2 shared atoms |
| Adjacent mono-vac (face corner) | 23 | — | 1 shared atom |

## Procedure

### Step 0: Validate inputs and set up directories

1. Confirm the user has provided or you can identify:
   - `data_path`: path to the data directory (default: `./data/`)
   - `lattice_type`: the ideal lattice (default: `FCC`)
   - `lattice_constant`: lattice constant in Å (default: `4.0559` for FCC Al at 300 K)
   - `ideal_neighbors`: number of neighbors in the ideal lattice (default: `12` for FCC)
   - `atom_assignment_method`: `PTM`, `CNA`, or `aCNA` (default: `PTM`)
2. Verify `<data_path>/dumps/` exists and contains LAMMPS dump files.
3. Create output directories if they don't exist: `neighbors/`, `graphs/`, `graphs/csvs/`, `components/`, `components/csvs/`, `defects/`.

### Step 1: Collect neighbors

Run [collect_neighbors.py](./../../collect_neighbors.py):

```
python collect_neighbors.py --path <data_path> --lattice_constant <lattice_constant>
```

**What it does:**
- Reads each LAMMPS dump via OVITO
- Applies PTM to classify atom structure types
- Computes Voronoi analysis for atomic volumes
- Finds all neighbor pairs within the cutoff radius: $r_{\text{cut}} = \frac{a + \frac{a}{\sqrt{2}}}{2}$
- Writes neighbor pair files to `<data_path>/neighbors/`
- Writes per-atom CSV data (coordinates, structure type, volume) to `<data_path>/graphs/csvs/`

### Step 2: Build graphs

Run [make_graphs.py](./../../make_graphs.py):

```
python make_graphs.py --path <data_path>
```

**What it does:**
- Reads neighbor pair files and constructs a NetworkX graph per frame
- Edge weights are neighbor distances
- Computes per-atom graph features: `n_neighbors`, `summed_neighbor_distances`, `norm_distances`
- Optionally computes extra features with `--extra`: triangles, weight stats, degree stats
- Saves graphs as `.gpickle` in `<data_path>/graphs/`
- Updates CSVs in `<data_path>/graphs/csvs/`

### Step 3: Generate components

Run [generate_components.py](./../../generate_components.py):

```
python generate_components.py --path <data_path> --ideal_neighbors <ideal_neighbors>
```

**What it does:**
- Loads each graph and removes atoms with ideally-coordinated neighbor counts
- Extracts connected components from the remaining defect subgraph
- The largest component is the GBS ($C_{\text{GBS}}$); all others are in-grain defects ($C_k$)
- Writes component graphs to `<data_path>/components/`
- Writes per-component statistics CSV to `<data_path>/components/csvs/`
- Writes summary file `<data_path>/all_component_data.csv` with columns: `frame, nodes, edges, components, largest_grain`

### Step 4: Collect and classify defects

Run [collect_defect_atoms.py](./../../collect_defect_atoms.py):

```
python collect_defect_atoms.py --path <data_path> --start 0 --n_nodes 12 --n_edges 24
```

**What it does:**
- Filters components matching a specific (nodes, edges) defect signature
- Collects the atom indices belonging to each matching component
- Saves a defect dictionary pickle to `<data_path>/defects/`

Run this step multiple times with different `--n_nodes` / `--n_edges` for each defect type of interest (see Defect Signatures table above).

### Step 5: Track defect changes over trajectory

Run [collect_changed_defects.py](./../../collect_changed_defects.py):

```
python collect_changed_defects.py --path <data_path> --start 0 --n_nodes 12 --n_edges 24
```

**What it does:**
- For each defect identified in Step 4, tracks its constituent atoms across all frames
- Identifies atoms whose neighbor count changed over the simulation (indicating defect migration, recombination, or transformation)
- Writes per-defect CSV with full trajectory data to `<data_path>/defects/`

### Step 6: Alternative — single-structure analysis via MaterialGraph class

For a single structure (not a trajectory), use [material_graph_class.py](./../../material_graph_class.py):

```python
from material_graph_class import MaterialGraph

mg = MaterialGraph(
    input="path/to/dump.txt",
    output="path/to/output_dir",
    atom_assignment_method="PTM",    # or 'CNA', 'aCNA'
    nondefect_lattice="FCC",
    cutoff=3.462
)
mg.subgraph()
mg.monovacancy_check()
```

This produces:
- `mg.G`: full defect graph
- `mg.C_gbs`: GBS component (largest component)
- `mg.C_k`: list of in-grain defect components
- `mg.C_k_sizes`: array of component sizes
- `mg.confirmed_monovac`: indices of confirmed mono-vacancy components

## Reporting

After completing the pipeline, compile and present the following outputs.

### Single Structure Report

1. **Total defect atom count**: number of atoms not matching the ideal lattice
2. **GBS size**: number of atoms in $C_{\text{GBS}}$ (vertex count of the largest component)
3. **Number of in-grain defect components**: count of $C_k$
4. **Defect type breakdown** (table):
   - Mono-vacancies (12 nodes, 24 edges)
   - Di-vacancies (18 nodes, 40 edges)
   - Tri-vacancies — linear (24 nodes, 60 edges) and triangular (24 nodes, 54 edges)
   - Adjacent mono-vacancy pairs grouped by shared-atom count (20–23 nodes)
   - Unclassified components (grouped by node/edge count)
5. **Component size distribution**: histogram of $|C_k|$ sizes
6. **Individual defect log**: for each defect component, report its index, node count, edge count, classified type, and mean graph properties (degree, weighted degree)

### Trajectory Report (multi-frame)

All items from the single-structure report, plus:

7. **$\Delta v$ plot**: change in vertex count for $C_{\text{GBS}}$ and total $C_k$ relative to the first frame (cf. Fig. 3a in the paper)
8. **Defect count evolution**: line/bar plot of each defect type count vs. frame/time (cf. Fig. 5a)
9. **Atom flow between regions**: FCC ↔ GBS ↔ mono-vacancy ↔ other defects (cf. Fig. 3b)
10. **Energy–structure correlation** (if energy data available): energy vs. $|C_{\text{GBS}}|$ colored by simulation regime (cf. Fig. 4)
11. **Changed-defect summary**: for defects whose neighbor counts changed, report the defect index, frame range of change, and nature of change (migration, recombination, merging)
12. **Defect fate tracking**: for mono-vacancies from the initial frame, classify final state (persisted, absorbed into GBS, recombined, merged into cluster) (cf. Fig. 5b)

### Visualization — Bundled Plotting Script

Use the bundled report generator to produce all standard figures:

```
python .github/skills/gala-analysis/scripts/generate_report.py --path <data_path>
```

For a single frame only:

```
python .github/skills/gala-analysis/scripts/generate_report.py --path <data_path> --single_frame dump.0.txt.csv
```

The script ([generate_report.py](./scripts/generate_report.py)) produces:
- Component size distribution histogram (log-scale y-axis)
- Defect type breakdown bar chart
- Trajectory Δv plot (GBS vs. in-grain defects)
- Defect type count evolution over frames
- Total component count over frames
- Per-frame defect log CSVs

All figures are saved to `<data_path>/figures/` by default (override with `--output`).

## Output Checklist

Before presenting results to the user, verify:

- [ ] All pipeline steps completed without error
- [ ] `all_component_data.csv` exists and is non-empty
- [ ] Defect counts are reported by type
- [ ] Individual defect locations/indices are listed
- [ ] For trajectories: time-evolution plots are generated
- [ ] Summary statistics are printed in a clear table format
