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
usage = """%prog PDBID

 Code takes PDBID of protein, finds the coordinate file on the opm database, &
 makes use of the dummy atoms to return the membrane thickness and position
 relative to the geometrical center of the protein. Also finds the maximal
 protein radius at membrane edges. In general membrane will need to be given
 a smaller radius than this."""
import numpy
import MDAnalysis
import MDAnalysis.analysis
import sys
import os
import urllib2
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('PDBID')

args = parser.parse_args()
def get_width(positions, center):
    disparray = positions - center
    distarray =(disparray * disparray)[:, 0:2].sum(axis=1)
    return (distarray.max()**0.5)

def memplacer(PDBID):
    pdburl = 'http://opm.phar.umich.edu/pdb/{pdb}.pdb'.format(pdb=PDBID)
    dl = urllib2.urlopen(pdburl)
    pdb = open('{pdb}.pdb'.format(pdb=PDBID), 'w')
    pdb.write(dl.read())
    pdb.close()
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
    bottomlayer = protein.selectAtoms("prop z > {zbb} and prop z < {zbt}".format(zbb=zbotbottom, zbt = zbottop))
    toplayer = protein.selectAtoms("prop z > {ztb} and prop z < {ztt}".format(ztb=ztopbottom, ztt = ztoptop))
    botradius = get_width(bottomlayer.positions, leaflet0.centroid())
    topradius = get_width(toplayer.positions, leaflet1.centroid())
    print("Membrane thickness = {thicky}, beginning at {disty} below protein centroid \n max radius of bottom membrane exclusion = {bot} \n max radius of top membrane exclusion = {top} ".format(thicky=thickness, disty=dist, bot = botradius, top = topradius))
memplacer(args.PDBID)   
