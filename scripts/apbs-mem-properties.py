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
usage = """%prog PDBID --InputPDB

 Code takes PDBID of protein, finds the coordinate file on the opm database, &
 makes use of the dummy atoms to return the membrane thickness and position
 relative to the geometrical center of the protein. Also finds the maximal
 protein radius at membrane edges. In general membrane will need to be given
 a smaller radius than this. 
 Alternatively, if InputPDB is specified, still accesses opm database for membrane information, but finds radii from input PDB and gives outputs in terms of the input. This will only work if using a protomer as otherwise the membrane top and bottom will be displaced, and thus the radii will also be faulty."""
import numpy
import sys
import os
import urllib2
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('PDBID')
parser.add_argument('--InputPDB')
args = parser.parse_args()
def get_width(positions, center):
# Finds the maximum distance in the xy plane from center for a series of 
# positions. This is accomplished with the 0:2 portion of the index which
# returns only the x and y.
    disparray = positions - center
    distarray =(disparray * disparray)[:, 0:2].sum(axis=1)
    return (distarray.max()**0.5)

def memplacer(PDBID,InputPDB):
    pdburl = 'http://opm.phar.umich.edu/pdb/{pdb}.pdb'.format(pdb=PDBID)
    dl = urllib2.urlopen(pdburl)
    pdb = open('{pdb}.pdb'.format(pdb=PDBID), 'w')
    pdb.write(dl.read())
    pdb.close()
    try:
        import MDAnalysis
        import MDAnalysis.analysis
    except ImportError:
       print("MDAnalysis required for this script. Available from https://code.google.com/p/mdanalysis/")
    U = MDAnalysis.Universe('{pdb}.pdb'.format(pdb=PDBID))
    leaflet0 = MDAnalysis.analysis.leaflet.LeafletFinder(U, "resname DUM and prop z < 0").groups(0)
    leaflet1 = MDAnalysis.analysis.leaflet.LeafletFinder(U, "resname DUM and prop z > 0").groups(0)
    protein = U.selectAtoms("protein")
    centroid = protein.centroid()
    zbot = leaflet0.positions[0][2]
    ztop = leaflet1.positions[0][2]
    dist = centroid[2] - zbot
    thickness = ztop - zbot
    ztoptop = ztop + .5
    ztopbottom = ztop - .5
    zbottop = zbot + .5
    zbotbottom= zbot - .5
    if InputPDB is None:
# With the above tolerances in place, find the protein atoms from that section 
# and the maximum distance(in the xy plane) to gain an idea of membrane radius
        bottomlayer = protein.selectAtoms("prop z > {zbb} and prop z < {zbt}".format(zbb=zbotbottom, zbt = zbottop))
        toplayer = protein.selectAtoms("prop z > {ztb} and prop z < {ztt}".format(ztb=ztopbottom, ztt = ztoptop))
        botradius = get_width(bottomlayer.positions, leaflet0.centroid())
        topradius = get_width(toplayer.positions, leaflet1.centroid())
        print("Membrane thickness = {thicky}, beginning at {disty} below protein centroid \n max radius of bottom membrane exclusion = {bot} \n max radius of top membrane exclusion = {top} ".format(thicky=thickness, disty=dist, bot = botradius, top = topradius))
    else:
        print("Input PDB specified. Following values only valid if input PDB is protomer of PDBID")
        U2 = MDAnalysis.Universe(InputPDB)
        protein = U2.selectAtoms('protein')
        centroid = protein.centroid()
        zbot = centroid[2] - dist
        ztop = zbot + thickness
        ztoptop = ztop + .5
        ztopbottom = ztop - .5
        zbottop = zbot + .5
        zbotbottom= zbot - .5
# With the above tolerances in place, find the protein atoms from that section 
# and the maximum distance(in the xy plane) to gain an idea of membrane radius
        bottomlayer = protein.selectAtoms("prop z > {zbb} and prop z < {zbt}".format(zbb=zbotbottom, zbt = zbottop))
        toplayer = protein.selectAtoms("prop z > {ztb} and prop z < {ztt}".format(ztb=ztopbottom, ztt = ztoptop))
        botradius = get_width(bottomlayer.positions, bottomlayer.centroid())
        topradius = get_width(toplayer.positions, toplayer.centroid())
        print("Membrane thickness = {thicky}, beginning at z={zbottom} in provided input pdb \n max radius of bottom membrane exclusion = {bot} \n max radius of top membrane exclusion = {top} ".format(thicky=thickness, zbottom=zbot, bot = botradius, top = topradius))


memplacer(args.PDBID,args.InputPDB)   
