# GALAS: Graph Analytics for Large Atomistic Simulations

## Scripts
* collect_neighbors.py
* make_graphs.py
* generate_components.py
* collect_defect_atoms.py
* collect_changed_defects.py

## Data Directory Structure
```
./data/  
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
