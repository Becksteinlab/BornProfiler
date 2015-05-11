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
from bornprofiler import run_setup
logger = logging.getLogger("bornprofiler")
parser = argparse.ArgumentParser()
parser.add_argument('PDBID')
parser.add_argument('--InputPDB')
args = parser.parse_args() 

bornprofiler.start_logging()
run_setup.memplacer(args.PDBID,args.InputPDB)
bornprofiler.stop_logging()   
