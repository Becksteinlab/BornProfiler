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
usage = """%prog [options]
Generates template cfg file while omitting sections based on flagged options"""
import argparse
import bornprofiler
import logging
from bornprofiler.bpio import RunParameters
logger = logging.getLogger("bornprofiler")
parser = argparse.ArgumentParser()
parser.add_argument('--title',default="default.cfg")
parser.add_argument('-nomembrane', action='store_true')
parser.add_argument('-noplotting', action='store_true')
args = parser.parse_args()
title = args.title
nomembrane = args.nomembrane
noplotting = args.noplotting
def cfg_template(title,nomembrane,noplotting):
    RunParameters(title,nomembrane,noplotting)

bornprofiler.start_logging()
cfg_template(title,nomembrane,noplotting)
bornprofiler.stop_logging()
