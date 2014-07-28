#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""
:Authors:  Oliver Beckstein & Lennard van der Feltz
:Year: 2010 2014
:License: GPL3
:Copyright: (c) 2010 Oliver Beckstein
:Copyright: (c) 2013 Oliver Beckstein
:Copyright: (c) 2014 Lennard van der Feltz
"""
from __future__ import with_statement
usage = """%prog [options] parameter-file

Setup Born profile calculation with or without a membrane. Default is with 
membrane, if option --nomembrane is specified, setup will proceed without 
membrane considerations. Parameters are read from the parameter file.

This script creates directories for ion positions, required input files to
APBS, and scripts that can be run locally or through a queuing system to
perform all calculations.

All parameters for the calculation must be set in the parameter file. In this
way, one parameter file completely describes the calculation. In particular,
each calculation should have a different job name (section [job] name =
JobName). For example, the ion is set in the section ::

  [bornprofile]
  ion = <name>

where `<name>` is one of the ions for which the package knows the Born
radius (see below). The input file also specifies the input PQR file
and the list of sample points at which the ion is placed
(``[bornprofile] points``). The sample points file should be formatted
as one white-space separated xyz coordinate per line, or a PDB file or
a HOLE sph file.

The "Born radii" for ions were taken from Table III in

  Alexander A. Rashin, Barry Honig (1985) J. Phys. Chem. 89(26):5588-5593
  http://dx.doi.org/10.1021/j100272a006

This paper suggests using the corrected covalent radius (Born radius)
and not the Pauling radius.

Born radii for H3O+, OH- (and H+... for testing) have been derived
from the solvation free energies in

  J.R. Pliego and J.M. Riveros. Chemical Physics Letters, 332(5-6): 597--602,
  2000.  http://dx.doi.org/10.1016/S0009-2614(00)01305-1.

directly via the Born equation. USE AT YOUR OWN RISK!!
"""

import bornprofiler

import os
import logging
logger = logging.getLogger('bornprofiler')

if __name__ == "__main__":
  import sys
  from optparse import OptionParser

  bornprofiler.start_logging()

  parser = OptionParser(usage=__doc__)
  parser.add_option("--run", dest="run", action="store_true",
                    help="immediately run apbs and draw_membrane2a to produce "
                    "all input files (can take a while); the default is to do "
                    "this as part of the individual jobs")
  parser.add_option("--nomembrane", dest = "no_membrane",action = "store_true",help="skip membrane steps")
  opts,args = parser.parse_args()

  try:
    filename = args[0]
  except:
    logger.fatal("Provide the parameter filename. See --help.")
    sys.exit(1)

  logger.info("run config = %(filename)r", vars())
  if opts.no_membrane:
     P = bornprofiler.core.Placeion(filename)
     P.generate()
  else:
     P = bornprofiler.core.MPlaceion(filename)
     P.generate(run=opts.run)

  bornprofiler.stop_logging()
