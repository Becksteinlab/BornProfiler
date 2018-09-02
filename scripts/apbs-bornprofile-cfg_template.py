#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""Generates template cfg file while omitting sections based on flagged options"""

import bornprofiler
import logging
from bornprofiler.bpio import RunParameters
logger = logging.getLogger("bornprofiler")

def cfg_template(title, nomembrane, noplotting):
    RunParameters(title,nomembrane,noplotting)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('--title',default="default.cfg")
    parser.add_argument('-nomembrane', action='store_true')
    parser.add_argument('-noplotting', action='store_true')
    args = parser.parse_args()

    bornprofiler.start_logging()
    cfg_template(args.title, args.nomembrane, args.noplotting)
    bornprofiler.stop_logging()
