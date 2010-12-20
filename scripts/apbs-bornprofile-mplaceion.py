#!/usr/bin/env python
# :Author: Oliver Beckstein <oliver.beckstein@bioch.ox.ac.uk>
# :Year: 2010
# :Licence: GNU Public Licence, version 3
#
# Copyright (c) 2010 Oliver Beckstein <oliver.beckstein@bioch.ox.ac.uk>

"""%prog [options] parameter-file

Setup Born profile calculation with a membrane. Parameters are read from the
parameter file. A new parameter-file can be generated with the --template
option; in this case only the file is written and no further actions are
performed.

All parameters for the calculation must be set in the parameter file. In this
way, one parameter file completely describes the calculation. In particular,
each calculation should have a different job name (section [job] name =
JobName).

This script creates directories for ion positions, required input files to
APBS, and scripts that can be run locally or through a queuing system to
perform all calculations.
"""

from __future__ import with_statement

import bornprofiler

import os
import logging
logger = logging.getLogger('bornprofiler') 

if __name__ == "__main__":
  import sys
  from optparse import OptionParser

  bornprofiler.start_logging()

  parser = OptionParser(usage=__doc__)
  parser.add_option("--template", dest="write_template", action="store_true",
                    help="write template parameter file and exit")
  parser.add_option("--run", dest="run", action="store_true",
                    help="immediately run apbs and draw_membrane2a to produce "
                    "all input files (can take a while); the default is to do "
                    "this as part of the individual jobs")
  opts,args = parser.parse_args()
  
  try:
    filename = args[0]
  except:
    logger.fatal("Provide the parameter filename. See --help.")
    sys.exit(1)


  if opts.write_template:
      bornprofiler.write_parameters(filename)
      sys.exit(0)

  logger.info("run config = %(filename)r", vars())

  P = bornprofiler.core.MPlaceion(filename)
  P.generate(run=opts.run)

  bornprofiler.stop_logging()
