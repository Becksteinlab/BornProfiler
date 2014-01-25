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
usage = """%prog PDBID
Downloads PDB of protein with user input PDBID from opm database, then writes a new file with only the atom information of the protein itself. This removes all ligands, waters, and dummy atoms, as well as any header lines from the pdb."""

import bornprofiler
import urllib2
import sys
import traceback
import argparse
import logging
logger = logging.getLogger("bornprofiler")
parser = argparse.ArgumentParser()
parser.add_argument('PDBID')
args = parser.parse_args()


def get_pdb_from_opm(PDBID):
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

bornprofiler.start_logging()
get_protein_pdb(args.PDBID)
bornprofiler.stop_logging()
