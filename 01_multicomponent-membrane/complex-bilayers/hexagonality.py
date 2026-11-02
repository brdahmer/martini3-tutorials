#!/usr/bin/env python3
# coding: utf-8

###Import required modules
import MDAnalysis as mda
from MDAnalysis.transformations import PositionAverager
import numpy as np
from MDAnalysis.lib.distances import capped_distance
from tqdm import tqdm
import argparse
import networkx

# ##### Lipid Hexagonality Determination

parser = argparse.ArgumentParser()
parser.add_argument('-f', dest="traj", default='traj_comp.xtc',
                    help = 'The trajectory to analyze.')
parser.add_argument('-s', dest="topol",
        help = '.tpr, .gro, or .pdb file with the same atom numbering as the trajectory.')
#parser.add_argument('-skip', dest="skip", type=int, default=1,
#        help = 'How many frames to skip while iterating.')
opts = parser.parse_args()
u = mda.Universe(opts.topol, opts.traj)
n_frame_avg = 5  # Averaging highlights the hexagonal packing
u.trajectory.add_transformations(PositionAverager(n_frame_avg))

### select the beads of interest
tail_beads = u.select_atoms('name C1A C1B')

### initialize the storage containers and some variables
#variables
neighbor_cutoff = 6. # Angstrom

fraction_hexagon = []
dist = np.empty((len(tail_beads), len(tail_beads)), dtype=np.float64)

### Iterate over the trajectory
for frame in tqdm(u.trajectory):
    if frame.frame < n_frame_avg:
        continue # Accumulating averaging frames

    cdists = capped_distance(tail_beads.positions,
                             tail_beads.positions,
                             min_cutoff=1.,
                             max_cutoff=neighbor_cutoff,
                             box=u.dimensions,
                             return_distances=False)
    graph = networkx.graph.Graph()
    graph.add_edges_from(cdists)
    graph.remove_nodes_from([n for n,d in graph.degree if d<2])
    hex_nodes = set()
    for n,d in graph.degree:
        if d == 6:
            hex_nodes.add(n)
            hex_nodes.update(graph[n].keys())
    fraction_hexagon.append((frame.time, len(hex_nodes)/len(tail_beads)))

fraction_hexagon = np.array(fraction_hexagon)
np.savetxt('hexagonality.xvg', fraction_hexagon, fmt= "%.4f")
