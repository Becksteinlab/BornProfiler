#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""
Code takes PDBID of protein, finds the coordinate file on the opm
database, & makes use of the dummy atoms to return the membrane
thickness and position relative to the geometrical center of the
protein. Also finds the maximal protein radius at membrane edges. In
general membrane will need to be given a smaller radius than this.
Alternatively, if InputPDB is specified, still accesses opm database
for membrane information, but finds radii from input PDB and gives
outputs in terms of the input. This will only work if using a protomer
as otherwise the membrane top and bottom will be displaced, and thus
the radii will also be faulty."""

import logging

import bornprofiler
from bornprofiler import run_setup

logger = logging.getLogger("bornprofiler")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('PDBID', help="Protein Databank ID of a membrane protein")
    parser.add_argument('--InputPDB', help="Use PQR from this file but download OPM orientation.")
    args = parser.parse_args()

    bornprofiler.start_logging()
    run_setup.memplacer(args.PDBID,args.InputPDB)
    bornprofiler.stop_logging()
