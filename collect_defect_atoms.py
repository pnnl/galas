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

import numpy as np
import pandas as pd
import os
import os.path as op
import pickle
import logging
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default='./data/', 
                    help='Path to data folder')
parser.add_argument('--start', type=int, default=0,
                    help='Index of structure to start collection')
parser.add_argument('--n_nodes', type=int, default=12,
                    help='Number of nodes ideal defect')
parser.add_argument('--n_edges', type=int, default=24,
                    help='Number of edges in ideal defect')
args = parser.parse_args()

if not op.isdir(op.join(args.path, 'defects')):
    os.mkdir(op.join(args.path, 'defects'))

data_path = op.join(args.path, 'graphs')
all_frames=sorted([f for f in os.listdir(data_path) if f.endswith('.npz')])

# load dataframe from t0
df=pd.read_csv(op.join(data_path, 'csvs', all_frames[args.start].replace('.npz','.csv')))

# load component metadata (defect indices and sorted labels)
meta = np.load(op.join(args.path, 'components', all_frames[args.start].replace('.npz', '.meta.npz')))
defect_indices = meta['defect_indices']
labels = meta['labels']

# Reconstruct components as lists of original atom indices, sorted by size (largest first)
n_comp = labels.max() + 1
components = []
for comp_id in range(n_comp):
    original_indices = defect_indices[labels == comp_id]
    components.append(original_indices.tolist())

# gather component info
info = pd.read_csv(op.join(args.path, 'components', 'csvs', all_frames[args.start].replace('.npz','.csv')))

# get all single vacancies
defect_index = info.loc[(info.Nodes==args.n_nodes)&(info.Edges==args.n_edges)].Component.tolist()

defect_dict={}
for d_index in defect_index:
    defect_dict[d_index] = components[d_index]
with open(op.join(args.path, 'defects', f'defect_dict_{args.n_nodes}nodes_{args.n_edges}edges_startat{args.start}.pickle'), 'wb') as handle:
    pickle.dump(defect_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
