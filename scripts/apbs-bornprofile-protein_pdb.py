#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""Downloads PDB of protein with user input PDBID from OPM database
<https://opm.phar.umich.edu/>, then writes a new file with only the
atom information of the protein itself. This removes all ligands,
waters, and dummy atoms, as well as any header lines from the pdb.

"""

import logging

import bornprofiler
from bornprofiler import run_setup

logger = logging.getLogger("bornprofiler")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('PDBID', help="Protein Databank ID of a membrane protein.")
    args = parser.parse_args()

    bornprofiler.start_logging()
    run_setup.get_protein_pdb(args.PDBID)
    bornprofiler.stop_logging()
