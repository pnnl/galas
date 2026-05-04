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
from ovito.io import import_file
from ovito.data import CutoffNeighborFinder
from ovito.modifiers import PolyhedralTemplateMatchingModifier, VoronoiAnalysisModifier
import numpy as np
import pandas as pd
import logging
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default='./data/', 
                    help='Path to data folder')
parser.add_argument('--lattice_constant', type=float, default=4.0559, 
                    help='Ideal lattice constant -- default for FCC Al at 300 K')
args = parser.parse_args()

# Set up file structure
# NOTE: code assumes dumps folder already present and populated
if not op.isdir(op.join(args.path, 'neighbors')):
    os.mkdir(op.join(args.path, 'neighbors'))
    
if not op.isdir(op.join(args.path, 'graphs')):
    os.mkdir(op.join(args.path, 'graphs'))

if not op.isdir(op.join(args.path, 'graphs/csvs')):
    os.mkdir(op.join(args.path, 'graphs/csvs'))
    
if not op.isdir(op.join(args.path, 'components')):
    os.mkdir(op.join(args.path, 'components'))
    
if not op.isdir(op.join(args.path, 'components/csvs')):
    os.mkdir(op.join(args.path, 'components/csvs'))

# Collect all LAMMPS dumps
data_path = op.join(args.path, 'dumps')
all_frames = [f for f in os.listdir(data_path)]
logging.info(f'{len(all_frames)} to total frames in {data_path}')

# Read through LAMMPS files and collect neighbors
for load_file in all_frames:
    # Define file names
    neighbor_file = op.join('neighbors', load_file + '.neighbors.txt.gz')
    csv_file = op.join('graphs/csvs', load_file + '.csv')

    # Read in LAMMPS dump
    logging.info(f'Collecing neighbors from {load_file}')
    try:
        pipeline = import_file(op.join(data_path, load_file))
        pipeline.modifiers.append(PolyhedralTemplateMatchingModifier())
        
        # Set up the Voronoi analysis modifier to compute atomic volumes
        voro = VoronoiAnalysisModifier(
            compute_indices = True,
            use_radii = True,
            edge_threshold = 0.1
        )
        pipeline.modifiers.append(voro)
        
        data = pipeline.compute()
    except:
        logging.warning(f'... file could not be read by Ovito')
        continue

    # Determine neighbor cutoff based on unit cell
    # NOTE: cutoff equation for FCC only!
    cutoff = np.sqrt(2)*args.lattice_constant/2 
    cutoff = (args.lattice_constant+cutoff)/2
    logging.info(f'... cutoff for neighbor distances is {cutoff:0.3f} A')

    # Initialize neighbor finder object
    finder = CutoffNeighborFinder(cutoff, data)

    # Find neighbors and accumulate into arrays
    row, col, dist = [], [], []
    for index in range(data.particles.count):
        for neigh in finder.find(index):
            row.append(index)
            col.append(neigh.index)
            dist.append(neigh.distance)

    # Bulk-write neighbor pairs to file
    edges = np.column_stack((row, col, dist))
    with gzip.open(op.join(args.path, neighbor_file), "wb") as f:
        np.savetxt(f, edges, fmt='%d %d %.6f')
    logging.info(f'... neighbor pairs written as {neighbor_file}') 

    # Collect coordinate data for graph/csvs  
    pos = np.asarray(data.particles['Position'])
    d={'idx': np.asarray(data.particles['Particle Identifier']), 
       'atom_type': np.asarray(data.particles['Particle Type']),
       'structure_type': np.asarray(data.particles['Structure Type']),
       'atomic_volume': np.asarray(data.particles['Atomic Volume']),
       'x': pos[:, 0],
       'y': pos[:, 1],
       'z': pos[:, 2]}
    df = pd.DataFrame(d)
    df.to_csv(op.join(args.path, csv_file), index=False)
    logging.info(f'... coordinate data written to {csv_file}')

