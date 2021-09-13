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


with open(op.join(args.path, 'defects', f'defect_dict_{args.n_nodes}nodes_{args.n_edges}edges_startat{args.start}.pickle'), 'rb') as handle:
    defect_dict=pickle.load(handle)

component = list(defect_dict.keys())


for c in component:
    df_all=pd.DataFrame()

    data_path = op.join(args.path, 'graphs')
    all_frames=sorted([f for f in os.listdir(data_path) if op.isfile(op.join(data_path, f))])

    for structure in sorted([int(x.split('.')[n]) for x in all_frames]):
        # load dataframe from next time step
        df=pd.read_csv(op.join(data_path, 'csvs', f'dump.{structure}.txt.csv'))
        df['norm_distances'] = df['summed_neighbor_distances']/df['n_neighbors']
        df['step']=structure

        # merge with full dataframe    
        df_all=pd.concat([df_all, df.iloc[vacancy_dict[c]]])

    # check if any neighbor counts in the vacancy changed
    atoms = df_all.loc[(df_all.step==0)].idx.tolist()
    changed_atoms=[int(atom) for atom in atoms if df_all.loc[df_all.idx==atom].n_neighbors.std()>0]

    if len(changed_atoms)>0:
        df_all.to_csv(op.join(args.path,'defects', f'defect_info_{args.n_nodes}nodes_{args.n_edges}edges_component{c}.csv'), index=False)
