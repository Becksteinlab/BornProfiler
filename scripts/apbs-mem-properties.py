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
import bornprofiler
import urllib2
import sys
import traceback
import argparse
import logging
logger = logging.getLogger("bornprofiler")
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
    return numpy.array([numpy.mean(distarray)**0.5,numpy.max(distarray)**0.5])
def get_pdb_from_opm(PDBID):
    pdburl = 'http://opm.phar.umich.edu/pdb/{pdb}.pdb'.format(pdb=PDBID.lower())
    try:
        dl = urllib2.urlopen(pdburl)
        pdb = open('{pdb}.pdb'.format(pdb=PDBID), 'w')
        pdb.write(dl.read())
        pdb.close()
        return pdb
    except:
        traceback.print_exc()
        logger.fatal("File not found in opm database. Double check PDBID")
        sys.exit(1)
def memplacer(PDBID,InputPDB):
    bornprofiler.start_logging()
    try:
        import MDAnalysis
        import MDAnalysis.analysis
    except ImportError:
        traceback.print_exc()
        logger.fatal("MDAnalysis required for this script. Available from https://code.google.com/p/mdanalysis/")
        sys.exit(1)
    pdb = get_pdb_from_opm(PDBID)
    U = MDAnalysis.Universe('{pdb}.pdb'.format(pdb=PDBID))
    try: 
        leaflet0 = MDAnalysis.analysis.leaflet.LeafletFinder(U, "resname DUM and prop z < 0").groups(0)
    except IndexError:
        traceback.print_exc()
        logger.fatal("Membrane information not available for submitted protein")
        sys.exit(1)
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
# and the average and maximum distance(in the xy plane) to gain an idea of membrane radius
        bottomlayer = protein.selectAtoms("prop z > {zbb} and prop z < {zbt}".format(zbb=zbotbottom, zbt = zbottop))
        toplayer = protein.selectAtoms("prop z > {ztb} and prop z < {ztt}".format(ztb=ztopbottom, ztt = ztoptop))
        botradii = get_width(bottomlayer.positions, leaflet0.centroid())
        topradii = get_width(toplayer.positions, leaflet1.centroid())
        Meminfo = "Membrane thickness = {thicky}, beginning at {disty} below protein centroid \nminimum suggested radius of bottom membrane exclusion = {botsug}, maximum {botmax} \nminimum suggested radius of top membrane exclusion = {topsug}, maximum {topmax} ".format(thicky=thickness, disty=dist,botsug = botradii[0], botmax = botradii[1], topsug = topradii[0], topmax = topradii[1])

    else:
        logger.info("Input PDB specified. Resulting values only valid if input PDB is protomer of PDBID. See documentation of this script for further information.")
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
# and the average and maximum distance(in the xy plane) to gain an idea of membrane radius
        bottomlayer = protein.selectAtoms("prop z > {zbb} and prop z < {zbt}".format(zbb=zbotbottom, zbt = zbottop))
        toplayer = protein.selectAtoms("prop z > {ztb} and prop z < {ztt}".format(ztb=ztopbottom, ztt = ztoptop))
        botradii = get_width(bottomlayer.positions, bottomlayer.centroid())
        topradii = get_width(toplayer.positions, toplayer.centroid())
        Meminfo = "Membrane thickness = {thicky}, beginning at z={zbottom} in provided input pdb \nminimum suggested radius of bottom membrane exclusion = {botsug}, maximum {botmax} \nminimum suggested radius of top membrane exclusion = {topsug}, maximum {topmax} ".format(thicky=thickness, zbottom = zbot,botsug = botradii[0], botmax = botradii[1], topsug = topradii[0], topmax = topradii[1])
    Meminfo_filename = "Membrane_info_{pdbid}.txt".format(pdbid=PDBID)
    Meminfofile = open(Meminfo_filename, 'w')
    Meminfofile.write(Meminfo)
    Meminfofile.close()
    logger.info("Membrane information written to {filename}".format(filename=Meminfo_filename))
    bornprofiler.stop_logging()

memplacer(args.PDBID,args.InputPDB)   
