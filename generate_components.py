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
import networkx as nx
import logging
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default='./data/', 
                    help='Path to data folder')
parser.add_argument('--ideal_neighbors', type=int, default=12,
                    help='Ideal number of neighbors in unit cell -- default for FCC')
args = parser.parse_args()

# Collect all graphs
data_path = op.join(args.path, 'graphs')
all_frames = [f for f in os.listdir(data_path)]
logging.info(f'{len(all_frames)} to total frames in {data_path}')

all_component_path = op.join(args.path, 'all_component_data.csv')
with open(all_component_path, 'w') as f:
    f.write("frame,nodes,edges,components,largest_grain\n")
    
# Read through graph files and generate components
for frame in all_frames:
    # Load graph
    logging.info(f'Loading graph from {frame}')
    try:
        G = nx.read_gpickle(frame)
    except:
        logging.warning('... graph could not be read')
        continue
        
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()

    # Lad graph data
    info_path = frame.replace('gpickle','csv').replace('graphs/','graphs/csvs/')
    try:
        df=pd.read_csv(info_path)
    except:
        logging.warning(f'... graph info could not be read from {info_path}')
        continue

    # Remove nodes with ideal number of neighbors
    nodes_to_remove = list(df.loc[df['n_neighbors']==args.ideal_neighbors].index)
    G.remove_nodes_from(nodes_to_remove)

    # Save full component subgraph
    component_path = frame.replace('graphs/','components/')
    nx.write_gpickle(G, component_path)
    logging.info(f'... component graph written to {component_path}')
              
    # Collect info on each atom in component
    a,b,c,d,e,f,g=[],[],[],[],[],[],[]
    d_std,e_std,f_std,g_std=[],[],[],[]
    n_components = 0
    for idx,g in enumerate([G.subgraph(c).copy() for c in sorted(nx.connected_components(G), key=len, reverse=True)]):
        nodes_to_keep = list(g.nodes())
        a.append(idx)
        b.append(g.number_of_nodes())
        c.append(g.number_of_edges())
        d.append(np.array([x[1] for x in g.degree()]).mean())
        d_std.append(np.array([x[1] for x in g.degree()]).std())
        e.append(np.array([x[1] for x in g.degree(weight='weight')]).mean())
        e_std.append(np.array([x[1] for x in g.degree(weight='weight')]).std())
        f.append(df.iloc[nodes_to_keep][f'summed_neighbor_distances'].to_numpy().mean())
        f_std.append(df.iloc[nodes_to_keep][f'summed_neighbor_distances'].to_numpy().std())
        g.append(df.iloc[nodes_to_keep][f'n_neighbors'].to_numpy().mean())
        g_std.append(df.iloc[nodes_to_keep][f'n_neighbors'].to_numpy().std())
        n_components += 1
        
    # Write component info csv   
    component_info_path = frame.replace('gpickle','csv').replace('graphs/','components/csvs/')
    info=pd.DataFrame({'Component':a,'Nodes':b,'Edges':c,'MeanDegree_sub':d, 'MeanDegree_sub_std':d_std, 
                       'MeanWDegree_sub':e, 'MeanWDegree_sub_std':e_std, 'MeanDegree_full':f, 
                       'MeanDegree_full_std':f_std, 'MeanWDegree_full':g, 'MeanWDegree_full_std':g_std})
    info.to_csv(component_info_path, index=False)
    logging.info(f'... detailed component graph info written to {component_info_path}')

    # Write data for each component
    with open(all_component_path, 'a') as f:
        f.write(f"{frame},{n_nodes},{n_edges},{n_components},{max(b)}\n")
    logging.info(f'... component info written to {all_component_path}')
        
    # Clear graph to free memory
    G.clear()

