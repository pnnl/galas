# GALAS: Graph Analytics for Large Atomistic Simulations

The graph-theoretical concept of connected components is employed to extract the evolution of defect configurations in large atomistic simulations. Building upon standard nearest neighbor analysis, graph theory and associated tools are used to reduce multi-million-atom systems into discrete component subgraphs that represent distinct structural defects. This method allows the automated identification, characterization, and tracking of defective regions within large volumes of data representing atomic-scale processes. Such analysis elucidates relationships between external stimuli and defect distributions, which have a large influence on material properties. 

## Scripts
* ```material_graph_class.py```: improved material graph creation method using PTM, CNA, or a-CNA and vacancy isomorphism check
* ```collect_neighbors.py```: collect neighboring atom pairs from LAMMPS dumps
* ```make_graphs.py```: create full graphs from neighbor pairs 
* ```generate_components.py```: break full graphs into components representing defect regions
* ```collect_defect_atoms.py```: aggregate all defects of a certain type
* ```collect_changed_defects.py```: collect defects that changed neighbors over the course of the simulation

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
```

## MD simulation of severe shear deformation in polycrystalline Al

The GALAS codebase was applied to atomic structures obtained from a molecular dynamics (MD) simulation of a polycrystalline aluminum structure containing ~8.3 million atoms subjected to high shear strain and annealing. The graph-component approach detailed in this codebase aided the identification and categorization of various defect structures within the large simulation. Data are provided as LAMMPS dump files obtained at every 10 ps during shear and 100 ps during annealing at 300 K: https://figshare.com/projects/MD_simulation_of_severe_shear_deformation_in_polycrystalline_Al/140534


## Citation

J. A. Bilbrey, N. Chen, S. Hu, P. V. Sushko, Graph-component approach to defect identification in large atomistic simulations, <em>Computational Materials Science</em>, 214, 2022.
DOI: [10.1016/j.commatsci.2022.111700](https://www.sciencedirect.com/science/article/pii/S0927025622004244?dgcid=author)

