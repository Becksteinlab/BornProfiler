#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#

import bornprofiler

import logging
logger = logging.getLogger('bornprofiler')

usage = """%prog [options]

Set up the BornProfiler configuration directories. This only has to be
done once (but it will not cause any damage to run this script again).
"""


if __name__ == "__main__":
  from optparse import OptionParser

  parser = OptionParser(usage=usage)
  opts, args = parser.parse_args()

  bornprofiler.start_logging()
  bornprofiler.config.setup()
  if bornprofiler.config.check_setup():
    logger.info("Init successful: you now have the template directories under %r.",
                bornprofiler.config.configdir)
    logger.info("The package can also be customized by editing %r.",
                bornprofiler.config.CONFIGNAME)
    logger.info("Questions and feedback: Oliver Beckstein <obeckste@asu.edu>")
  else:
    logger.error("Something is wrong: Failed to setup the template directories.")
    logger.warn("You can proceed but problems migh arise and you will not be able "
                "to easily customize generation of runs and submission scripts.")
  bornprofiler.stop_logging()
