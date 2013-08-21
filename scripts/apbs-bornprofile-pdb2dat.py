#!/usr/bin/env python
# Code for merging any number of PDB files into a single .dat file takes the
# file names and --title= and writes the coordinates in rows into a new file
# with the specified title.
import numpy as np
import argparse
from MDAnalysis import *
parser=argparse.ArgumentParser()
parser.add_argument('pdbs', nargs='+', type=str)
parser.add_argument('--title')
args = parser.parse_args()
pdbs = args.pdbs
title = args.title
def pdb2dat(pdblist, name):
    xyz=[]
    for pdb in pdblist:
        u = Universe(pdb)
        positions= u.atoms.positions        
        for atom in positions:
            xyz.append(atom)
    dat = open('{x}.dat'.format(x = name),'w')
    for coord in xyz:
        dat.write("{x} {y} {z} \n".format(x = coord[0], y = coord[1], z = coord[2]))
    dat.close()
pdb2dat(pdbs,title)
