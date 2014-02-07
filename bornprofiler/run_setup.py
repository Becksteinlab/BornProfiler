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
:Year: 2014
:Licence: GPL 3
:Copyright: (c) 2014 Lennard van der Feltz
"""

"""Repository for various functions used repeatedly in run setup scripts. Many functions expect a logger under the title 'logger' to be in operation."""
import bornprofiler
import urllib2
import sys
import traceback
import argparse
import logging
import numpy

logger = logging.getLogger("bornprofiler")



def get_pdb_from_opm(PDBID):
    """Function for downloading pdb file from opm database."""
    pdburl = 'http://opm.phar.umich.edu/pdb/{pdb}.pdb'.format(pdb=PDBID.lower())
    try:
        dl = urllib2.urlopen(pdburl)
        pdb = open('{pdb}.pdb'.format(pdb=PDBID), 'w')
        pdb.write(dl.read())
        pdb.close()
        return pdb
    except:
        raise

def get_protein_pdb(PDBID):
    """Function for downloading pdb file and removing non-protein structures. Depends on get_pdb_from_opm."""
    try:
        import MDAnalysis
        import MDAnalysis.analysis
    except ImportError:
        traceback.print_exc()
        logger.fatal("MDAnalysis required for this script. Available from https://code.google.com/p/mdanalysis/")
        sys.exit(1)
    try:
        logger.info("Downloading pdb for {pdb} from opm database".format(pdb=PDBID))
        pdb = get_pdb_from_opm(PDBID)
    except:
        traceback.print_exc()
        logger.fatal("File not found in opm database. Double check PDBID {pdb}".format(pdb=PDBID))
        sys.exit(1)

    U = MDAnalysis.Universe('{pdb}.pdb'.format(pdb=PDBID))
    logger.info("Stripping waters, ligands, and other non-protein molecules")
    protein = U.selectAtoms("protein")
    logger.info("Writing protein info to {pdbid}_protein.pdb".format(pdbid=PDBID))
    protein.write("{pdbid}_protein.pdb".format(pdbid=PDBID))
    logger.info("Removing header section")
    f = open("{pdbid}_protein.pdb".format(pdbid=PDBID), 'r')
    outstring = ""
    for line in f:
        if line.split()[0] == 'ATOM':
            outstring+= line
        else:
            pass
    f.close()
    u = open("{pdbid}_protein.pdb".format(pdbid=PDBID), 'w')
    u.write(outstring)
    u.close()

def get_width(positions, center):
    """ Finds the maximum distance in the xy plane from center for a series of 
    positions. This is accomplished with the 0:2 portion of the index which
    returns only the x and y."""
    disparray = positions - center
    distarray =(disparray * disparray)[:, 0:2].sum(axis=1)
    return numpy.array([numpy.mean(distarray)**0.5,numpy.max(distarray)**0.5])

def memplacer(PDBID,InputPDB):
    try:
        import MDAnalysis
        import MDAnalysis.analysis
    except ImportError:
        traceback.print_exc()
        logger.fatal("MDAnalysis required for this script. Available from https://code.google.com/p/mdanalysis/")
        sys.exit(1)
    try:
        pdb = get_pdb_from_opm(PDBID)
    except:
        traceback.print_exc()
        logger.fatal("File not found in opm database. Double check PDBID")
        sys.exit(1)

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
    return [botradii[0],topradii[0],thickness,zbot]

