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
from bornprofiler import run_setup
logger = logging.getLogger("bornprofiler")
parser = argparse.ArgumentParser()
parser.add_argument('PDBID')
args = parser.parse_args()

bornprofiler.start_logging()
run_setup.get_protein_pdb(args.PDBID)
bornprofiler.stop_logging()
