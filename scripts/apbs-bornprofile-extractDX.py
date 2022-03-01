import numpy as np
import pandas as pd
from gridData import Grid

import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    '--input', help='input filename')
parser.add_argument(
    '--output', help='output filename. .dx is auto added to name')
parser.add_argument(
    'grid', type=float, help='grid interval')
parser.add_argument(
    'gridL', type=float, help='loose grid interval')

args = parser.parse_args()

filename = args.input
outname = args.output

# Grid interval
grid_int= args.grid

# Edge offset for Grid
offset = grid_int/2.0

# For losse edge
grid_intL = args.gridL
offsetL = grid_intL/2.0

df = pd.read_csv(filename,delim_whitespace=True)
df = df.to_numpy()

data = df[:,0:4]
samples = data[:, 0:3]

coord = data[:, 0:3].T
weight = data[:,3]
weight = weight

# get dimensions
xmin = np.min(data[:,0])
ymin = np.min(data[:,1])
zmin = np.min(data[:,2])
xmax = np.max(data[:,0])
ymax = np.max(data[:,1])
zmax = np.max(data[:,2])

# get grid dimensions
xlen = int((xmax-xmin)/grid_int) + 1
ylen = int((ymax-ymin)/grid_int) + 1
zlen = int((zmax-zmin)/grid_int) + 1

# get edge.
xedge = np.linspace(xmin-offset, xmax+offset, xlen+1)
yedge = np.linspace(ymin-offset, ymax+offset, ylen+1)
zedge = np.linspace(zmin-offset, zmax+offset, zlen+1)
edges = [xedge, yedge, zedge]

# get loose edge.
xedge_L = np.arange(xmin-offsetL, xmax+offsetL, grid_intL)
yedge_L = np.arange(ymin-offsetL, ymax+offsetL, grid_intL)
zedge_L = np.arange(zmin-offsetL, zmax+offsetL, grid_intL)
edges_L = [xedge_L, yedge_L, zedge_L]

# Filling empty points with high value
filling_value = 2* weight.max()
grid = np.full([xlen, ylen, zlen], filling_value)

def get_ind(x, xmin, xint):
    ind = int((x-xmin)/xint)
    return ind

for i, c in enumerate(coord.T):
    xind = get_ind(c[0], xmin, grid_int)
    yind = get_ind(c[1], ymin, grid_int)
    zind = get_ind(c[2], zmin, grid_int)
    grid[xind, yind, zind] = weight[i]

# resample
g = Grid(grid, edges=edges)

# Needs GridDataFormat 0.6 or higher
g.interpolation_spline_order = 3
g_L = g.resample(edges_L)

# export
g.export(outname+'.dx')
g_L.export(outname+'-loose.dx')
