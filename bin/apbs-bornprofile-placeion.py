#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
:Author: Kaihsu Tai, Oliver Beckstein
:Year: 2008, 2010
:Licence: GPL
:Copyright: (c) 2008 Kaihsu Tai
:Copyright: (c) 2010 Oliver Beckstein
:URL: http://en.wikiversity.org/wiki/Talk:Poisson%E2%80%93Boltzmann_profile_for_an_ion_channel

I wrote some Python code to automate this process. The job submission
requires a queuing system called Grid Engine. Copyright © 2008 Kaihsu
Tai. Moral rights asserted (why?). Hereby licensed under either GFDL
or GNU General Public License at your option.
"""
from __future__ import with_statement

import os
import logging
logger = logging.getLogger('bornprofile') 

usage = """%prog [options] pqr-file samplepoints-file

This script sets up input files for Born energy calculations for
APBS. It expects to find ion positions in the file samplepoints, one
xyz coordinate per line.

The "Born radii" for ions were taken from Table III in

  Alexander A. Rashin, Barry Honig (1985) J. Phys. Chem. 89(26):5588-5593
  http://dx.doi.org/10.1021/j100272a006

This paper suggests using the corrected covalent radius (Born radius)
and not the Pauling radius.
"""

from bornprofiler import IONS, JOBSCRIPTS, Placeion

   
if __name__ == "__main__":
  import sys
  from optparse import OptionParser

  logging.basicConfig()

  parser = OptionParser(usage=usage)
  parser.add_option("--name", dest="jobName",
                    metavar="STRING",
                    help="name for the job submission script [%default]")
  parser.add_option("--ion", dest="ionName", type="choice", choices=IONS.keys(),
                    metavar="NAME",
                    help="name of the ion to be sampled. Available values: %r. "
                    "Radii were taken from Table III in Rashin & Honig  1985. "
                    "The default ion is '%%default'." % (IONS.keys(),))
  parser.add_option("--dime", dest="dime", nargs=3, type="int",
                    metavar="NX NY NZ",
                    help="dimensions of the computational grid [97 97 193]")
  parser.add_option("--ionic-strength", dest="ionicStrength", type="float",
                    metavar="CONC",
                    help="set ionic strength of NaCl bath to the given concentration "
                    "CONC in mol/l [%default]")
  parser.add_option("--script", dest="script",
                    metavar="NAME",
                    help="name of a stored script template %r or (advanced usage!) a "
                    "filename that contains appropriate place holders [%%default]" % JOBSCRIPTS.keys())
  parser.set_defaults(ionicStrength=0.15, jobName="bornprofile", 
                      ionName="Na", dime=[97,97,193], script="local")

  opts,args = parser.parse_args()
  
  try:
    pqrfile, pointsfile = args
  except:
    logger.fatal("Needs PQR file and sample points. See --help.")
    sys.exit(1)

  try:
    script_template = JOBSCRIPTS[opts.script]
  except KeyError:
    try:
      script_template = open(opts.script).readlines()
    except:
      logger.fatal("--scripts=NAME must be either a filename or one of %r; see the " 
                   "source for the correct format of the file." % JOBSCRIPTS.keys())
      sys.exit(2)

  P = Placeion(pqrfile, pointsfile, ionName=opts.ionName, ionicStrength=opts.ionicStrength,
               dime=opts.dime, jobName=opts.jobName, script=script_template)
  P.generate()
