# GALAS: Graph Analytics for Large Atomistic Simulations

The graph-theoretical concept of connected components is employed to extract the evolution of defect configurations in large atomistic simulations. Building upon standard nearest neighbor analysis, graph theory and associated tools are used to reduce multi-million-atom systems into discrete component subgraphs that represent distinct structural defects. This method allows the automated identification, characterization, and tracking of defective regions within large volumes of data representing atomic-scale processes. Such analysis elucidates relationships between external stimuli and defect distributions, which have a large influence on material properties. 

## Scripts
* ```material_graph_class.py```: improved material graph creation method using PTM, CNA, or a-CNA and vacancy isomorphism check
* ```collect_neighbors.py```: collect neighboring atom pairs from LAMMPS dumps
* ```make_graphs.py```: create full graphs from neighbor pairs 
* ```generate_components.py```: break full graphs into components representing defect regions
* ```collect_defect_atoms.py```: aggregate all defects of a certain type
* ```collect_changed_defects.py```: collect defects that changed neighbors over the course of the simulation
* ```defect_ordering.py```: compute spatial ordering metrics (RDF, nearest-neighbor distributions, Warren-Cowley SRO parameters, structure factor) for defect centers grouped by type

## Data Directory Structure

The data directory structure and corresponding files types are as follows. Prior to running any of the scripts, only the ```./data/dumps/``` directory need to be created and populated with the raw simulation files (here, LAMMPS dumps).

```
./data/  
|
└───vacancy.edgelist.txt.gz
│
└───/dumps/
│       dump.0.txt
│       ...
│ 
└───/neighbors/
│       dump.0.txt.neighbors.txt.gz
│       ...
│   
└───/graphs/
|   │   dump.0.txt.gpickle
|   │   ...
│   └───/csvs/
│           dump.0.txt.csv
│           ...
│   
└───/components/
    │   dump.0.txt.gpickle
    │   ...
    └───/csvs/
            dump.0.txt.csv
            ...

└───/ordering/
        rdf_Mono-vacancy.csv
        nn_dist_Mono-vacancy.csv
        sq_Mono-vacancy.csv
        warren_cowley.csv
        ...
```

## MD simulation of severe shear deformation in polycrystalline Al

The GALAS codebase was applied to atomic structures obtained from a molecular dynamics (MD) simulation of a polycrystalline aluminum structure containing ~8.3 million atoms subjected to high shear strain and annealing. The graph-component approach detailed in this codebase aided the identification and categorization of various defect structures within the large simulation. Data are provided as LAMMPS dump files obtained at every 10 ps during shear and 100 ps during annealing at 300 K: https://figshare.com/projects/MD_simulation_of_severe_shear_deformation_in_polycrystalline_Al/140534


## Citation

J. A. Bilbrey, N. Chen, S. Hu, P. V. Sushko, Graph-component approach to defect identification in large atomistic simulations, <em>Computational Materials Science</em>, 214, 2022.
DOI: [10.1016/j.commatsci.2022.111700](https://www.sciencedirect.com/science/article/pii/S0927025622004244?dgcid=author)


## Copilot Skill: Automated Defect Analysis

A GitHub Copilot skill is included that can run the full GALAS pipeline and produce a defect report automatically via the chat interface.

### What it does

The skill guides an AI agent through the complete analysis workflow:

1. **Neighbor collection** — reads LAMMPS dump files, classifies atoms by structure type (PTM/CNA), and finds nearest-neighbor pairs within a cutoff radius.
2. **Graph construction** — builds a NetworkX graph for each frame where atoms are nodes and neighbor pairs are edges weighted by distance.
3. **Component extraction** — removes ideally-coordinated atoms and decomposes the remaining defect graph into connected components. The largest component is the grain boundary superstructure (GBS); all others are individual in-grain defects.
4. **Defect classification** — matches each component's (node count, edge count) signature against known defect templates (mono-vacancies, di-vacancies, tri-vacancies, etc.).
5. **Trajectory tracking** — follows defect atoms across frames to detect migration, recombination, or merging.
6. **Reporting** — prints a summary table of defect counts by type, generates matplotlib plots (component size distributions, defect evolution over time, GBS vs. in-grain vertex changes), and logs each individual defect with its location and classification.

### How to use it

In VS Code with GitHub Copilot Chat, type `/gala-analysis` and provide your parameters. Examples:

```
/gala-analysis Analyze ./data/ for FCC Al with lattice constant 4.0559
```

```
/gala-analysis Run single-frame defect analysis on ./data/dumps/dump.0.txt using PTM
```

```
/gala-analysis Process the full trajectory in ./data/ and report defect changes over time
```

The skill supports FCC, BCC, and HCP lattices with built-in defaults and produces:
- A defect type breakdown table with total counts
- Individual defect log (component index, node/edge counts, classified type)
- For trajectories: plots of defect count evolution, GBS size changes, and component distributions

### Requirements

The skill requires the same dependencies as the core GALAS scripts: `ovito`, `networkx`, `numpy`, `pandas`, `matplotlib`, and `tqdm`.

