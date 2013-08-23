#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""
:Author: Lennard van der Feltz
:Year: 2013
:Licence: GPL 3
:Copyright: (c) 2013 Lennard van der Feltz
"""

usage = """%prog input pdbs --title=output dat file
 Code for merging any number of PDB files into a single .dat file takes the
 file names and --title= and writes the coordinates in rows into a new file
 with the specified title."""
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
