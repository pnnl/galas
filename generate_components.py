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
from scipy.sparse import load_npz, save_npz
from scipy.sparse.csgraph import connected_components
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
all_frames = [f for f in os.listdir(data_path) if f.endswith('.npz')]
logging.info(f'{len(all_frames)} to total frames in {data_path}')

all_component_path = op.join(args.path, 'all_component_data.csv')
with open(all_component_path, 'w') as f:
    f.write("frame,nodes,edges,components,largest_grain\n")
    
# Read through graph files and generate components
for frame in all_frames:
    # Load sparse adjacency matrix
    logging.info(f'Loading graph from {frame}')
    try:
        A = load_npz(op.join(data_path, frame))
    except:
        logging.warning('... graph could not be read')
        continue
        
    n_nodes = A.shape[0]
    n_edges = A.nnz // 2  # each edge stored in both directions

    # Load graph data
    info_path = op.join(data_path, 'csvs', frame.replace('.npz', '.csv'))
    try:
        df = pd.read_csv(info_path)
    except:
        logging.warning(f'... graph info could not be read from {info_path}')
        continue

    # Identify defect atoms (non-ideal neighbor count)
    defect_mask = df['n_neighbors'].values != args.ideal_neighbors
    defect_indices = np.where(defect_mask)[0]

    # Extract defect subgraph (sparse submatrix)
    A_defect = A[defect_indices][:, defect_indices]

    # Save defect subgraph as sparse matrix
    component_path = op.join(args.path, 'components', frame)
    save_npz(component_path, A_defect.tocsr())
    logging.info(f'... component graph written to {component_path}')

    # Find connected components using scipy (much faster than NetworkX)
    n_comp, labels = connected_components(A_defect, directed=False)

    # Sort components by size (largest first) and relabel
    component_sizes = np.bincount(labels, minlength=n_comp)
    sorted_comp_ids = np.argsort(-component_sizes)
    label_map = np.empty(n_comp, dtype=int)
    label_map[sorted_comp_ids] = np.arange(n_comp)
    sorted_labels = label_map[labels]

    # Save component metadata (defect indices and labels) for downstream use
    meta_path = op.join(args.path, 'components', frame.replace('.npz', '.meta.npz'))
    np.savez(meta_path, defect_indices=defect_indices, labels=sorted_labels)

    # Collect info on each component
    comp_ids, comp_nodes, comp_edges = [], [], []
    deg_mean, deg_std, wdeg_mean, wdeg_std = [], [], [], []
    full_deg_mean, full_deg_std, full_wdeg_mean, full_wdeg_std = [], [], [], []

    # Precompute degree arrays for the defect subgraph
    A_defect_csr = A_defect.tocsr()
    sub_degrees = np.diff(A_defect_csr.indptr)
    sub_weighted_degrees = np.asarray(A_defect_csr.sum(axis=1)).flatten()

    for comp_id in range(n_comp):
        # Get local indices within defect subgraph for this component
        local_mask = sorted_labels == comp_id
        local_indices = np.where(local_mask)[0]
        original_indices = defect_indices[local_indices]

        # Node and edge count for this component
        comp_sub = A_defect_csr[local_indices][:, local_indices]
        n_comp_nodes = len(local_indices)
        n_comp_edges = comp_sub.nnz // 2

        comp_ids.append(comp_id)
        comp_nodes.append(n_comp_nodes)
        comp_edges.append(n_comp_edges)

        # Degree stats within the defect subgraph
        comp_degrees = sub_degrees[local_indices]
        deg_mean.append(comp_degrees.mean())
        deg_std.append(comp_degrees.std())

        # Weighted degree stats within the defect subgraph
        comp_wdegrees = sub_weighted_degrees[local_indices]
        wdeg_mean.append(comp_wdegrees.mean())
        wdeg_std.append(comp_wdegrees.std())

        # Stats from the full graph (stored in CSV)
        full_deg_mean.append(df.iloc[original_indices]['summed_neighbor_distances'].mean())
        full_deg_std.append(df.iloc[original_indices]['summed_neighbor_distances'].std())
        full_wdeg_mean.append(df.iloc[original_indices]['n_neighbors'].mean())
        full_wdeg_std.append(df.iloc[original_indices]['n_neighbors'].std())

    # Write component info csv   
    component_info_path = op.join(args.path, 'components', 'csvs', frame.replace('.npz', '.csv'))
    info = pd.DataFrame({
        'Component': comp_ids, 'Nodes': comp_nodes, 'Edges': comp_edges,
        'MeanDegree_sub': deg_mean, 'MeanDegree_sub_std': deg_std,
        'MeanWDegree_sub': wdeg_mean, 'MeanWDegree_sub_std': wdeg_std,
        'MeanDegree_full': full_deg_mean, 'MeanDegree_full_std': full_deg_std,
        'MeanWDegree_full': full_wdeg_mean, 'MeanWDegree_full_std': full_wdeg_std
    })
    info.to_csv(component_info_path, index=False)
    logging.info(f'... detailed component graph info written to {component_info_path}')

    # Write data for each component
    largest = comp_nodes[0] if comp_nodes else 0
    with open(all_component_path, 'a') as f:
        f.write(f"{frame},{n_nodes},{n_edges},{n_comp},{largest}\n")
    logging.info(f'... component info written to {all_component_path}')

    # Free memory
    del A, A_defect, A_defect_csr

