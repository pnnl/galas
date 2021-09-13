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
import gzip
import numpy as np
import pandas as pd
import networkx as nx
import logging
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default='./data/', 
                    help='Path to data folder')
parser.add_argument('--extra', action='store_true',
                    help='Compute extra graph properties')
args = parser.parse_args()

# NOTE: code assumes neighbors folder is present and populated
#       and file structure is in place

# Collect all LAMMPS dumps
data_path = op.join(args.path, 'neighbors')
all_frames = [f.replace('.neighbors.txt.gz','') for f in os.listdir(data_path)]
logging.info(f'{len(all_frames)} to total frames in {data_path}')

# Read through neighbor files and create graphs
for load_file in all_frames:
    # Define file names
    neighbor_file = op.join('neighbors', load_file + '.neighbors.txt.gz')
    graph_file = op.join('graphs', load_file + '.gpickle')
    csv_file = op.join('graphs/csvs', load_file + '.csv')

    logging.info(f'Creating graph from {neighbor_file}')
    
    # Read in csv containing coordinate data
    try:
        df = pd.read_csv(op.join(args.path,csv_file))
    except:
        logging.warning('... coordinate data not found')
        df=pd.DataFrame()
       
    # Create graph and store as gpickle
    try:
        G=nx.read_edgelist(gzip.open(op.join(args.path,neighbor_file), "rb"),
                           nodetype=int, data=(("weight", float),))
        nx.write_gpickle(G, op.join(args.path, graph_file))
        logging.info(f'... gpickle stored at {graph_file}')
    except:
        logging.warning(f'... failed to read {neighbor_file}')
        continue

    # Write neighbor info to csv
    df['n_neighbors'] = [x[1] for x in sorted(list(G.degree()))] 
    df['summed_neighbor_distances'] = [x[1] for x in sorted(list(G.degree(weight='weight')))] 
    df['norm_distances'] = df['summed_neighbor_distances']/df['n_neighbors']
    
    # Optional computation of additional graph properties
    if extra:
        # Compute number of triangles for each atom
        df['triangles']=df['idx'].apply(lambda x: nx.triangles(G, x))

        # Compute stats for neighbor distance (graph edge weight)
        weights = [[x[1]['weight'] for x in G[comp].items()] 
                   for comp in df['idx'].sort_values('idx').tolist()]
        df['weight_mean']=df['idx'].apply(lambda x: np.array(weights[x]).mean()) 
        df['weight_std']=df['idx'].apply(lambda x: np.array(weights[x]).std())
        df['weight_skew']=df['idx'].apply(lambda x: skew(weights[x]))

        # Compute stats for number of neighbors (graph node degree)
        degrees = [[G.degree[node] for node in [x for x in G.neighbors(comp)]] 
                   for comp in df['idx'].sort_values('idx').tolist()]
        df['neighbor_degree_mean']=df['idx'].apply(lambda x: np.array(degrees[x]).mean())
        df['neighbor_degree_std']=df['idx'].apply(lambda x: np.array(degrees[x]).std())
        df['neighbor_degree_skew']=df['idx'].apply(lambda x: skew(degrees[x]))
    
    df.to_csv(op.join(args.path, csv_file), index=False)
    logging.info(f'... info written to {csv_file}')

    # Clear graph from memory
    G.clear()
