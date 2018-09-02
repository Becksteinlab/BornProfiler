#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#

"""
Code for merging any number of PDB files into a single .dat file takes
the file names and --title= and writes the coordinates in rows into a
new file with the specified title. If title not specified names output
'pdbs.dat'.
"""

import MDAnalysis as mda

def pdb2dat(pdblist, name):
    xyz = []
    for pdb in pdblist:
        u = mda.Universe(pdb)
        positions = u.atoms.positions
        for atom in positions:
            xyz.append(atom)
    with open('{x}.dat'.format(x = name),'w') as dat:
        coordinate_string = ""
        for coord in xyz:
            coordinate_string += "{x} {y} {z} \n".format(x=coord[0], y=coord[1], z=coord[2])
            dat.write(coordinate_string)


if __name__ == "__main__":
    import argparse

    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pdbs', nargs='+', type=str)
    parser.add_argument('--title', default="pdbs")
    args = parser.parse_args()

    pdb2dat(args.pdbs, args.title)
