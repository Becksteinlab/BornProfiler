#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""
Setup Born profile calculation with or without a membrane. Parameters
are read from the parameter file.

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
from __future__ import with_statement

import bornprofiler

import os
import logging
logger = logging.getLogger('bornprofiler')

if __name__ == "__main__":
  import sys
  import argparse

  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("filename",
                      help="parameter file")
  parser.add_argument("--run", dest="run", action="store_true",
                      help="immediately run apbs and draw_membrane2a to produce "
                      "all input files (can take a while); the default is to do "
                      "this as part of the individual jobs")
  parser.add_argument("--nomembrane", dest="no_membrane", action="store_true",
                      help="do not use a membrane (skip membrane setup steps)")
  args = parser.parse_args()

  bornprofiler.start_logging()

  if not args.filename:
    logger.fatal("Provide the parameter filename. See --help.")
    sys.exit(1)

  logger.info("run config = %(filename)r", args.filename)
  if not args.no_membrane:
    P = bornprofiler.core.MPlaceion(args.filename)
  else:
    raise NotImplementedError("Omitting of membrane is not yet working.")
    # P = bornprofiler.core.Placeion(args.filename)
    # P.generate() ## can probably be omitted and use P.generate(run=args..run)

  P.generate(run=args.run)

  bornprofiler.stop_logging()
