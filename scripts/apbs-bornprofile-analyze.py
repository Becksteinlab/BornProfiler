#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""
:Author: Kaihsu Tai, Oliver Beckstein
:Year: 2008, 2010
:Licence: GPL
:Copyright: (c) 2008 Kaihsu Tai
:Copyright: (c) 2010-2013 Oliver Beckstein
:URL: http://en.wikiversity.org/wiki/Talk:Poisson%E2%80%93Boltzmann_profile_for_an_ion_channel

I wrote some Python code to automate this process. The job submission
requires a queuing system called Grid Engine. Copyright © 2008 Kaihsu
Tai. Moral rights asserted (why?). Hereby licensed under either GFDL
or GNU General Public License at your option.

"""

import logging
logger = logging.getLogger('bornprofiler') 

usage = """%prog [options] parameter-file
       %prog [options] samplepoints-file *.out       

Extract the electrostatic free energy from the numbered APBS output files
(produced via the apbs-bornprofil-placeion.py script) and associate each energy
with the position of the ion. The prefered usage is to supply the run
configuration file and possibly the directory where the datafiles are stored
(--basedir).

.. Note:: The same samplepoints-file must be provided that was used for setting
   up the APBS calculations.

.. SeeAlso:: apbs-bornprofile-analyze3d.py can also be used and has a number of
   different options. (Eventually, the two programs will be merged.)
"""

import bornprofiler
from bornprofiler.analysis import AnalyzeElec
 
if __name__ == "__main__":
  import sys
  import os
  import glob
  from optparse import OptionParser

  bornprofiler.start_logging()

  parser = OptionParser(usage=usage)
  parser.add_option("--name", dest="jobName",
                    metavar="STRING",
                    help="name used in output files (welec_STRING.{dat,pdf})")
  parser.add_option("--ion", dest="ionName",
                    metavar="STRING",
                    help="Set the ion name for --pdb if not running from config file [%default]")
  parser.add_option("--pdb", "-p", dest="pdbfilename", 
                    metavar="FILE",
                    help="Export points to a PDB file with the energy in the B-factor. "
                    "If set to 'auto' then a filename is chosen.")
  parser.add_option("--basedir", "-B", dest="basedir",
                    metavar="DIR",
                    help="when using a run parameter file, the job output is found under "
                    "'DIR/<job.name>/w[0-9][0-9][0-9][0-9]/job*.out', with job.name taken from "
                    "the run parameter file [%default]")
  parser.add_option("--read", dest="create", action="store_false",
                    help="If set, read positions and enrgies from a previously created dat file.")
  parser.set_defaults(basedir=os.path.curdir, ionName="ION", create=True)


  opts,args = parser.parse_args()

  if len(args) == 0:
    logger.fatal("Needs samplepoints file and at least one APBS output file. See --help.")
    sys.exit(1)  
  elif len(args) == 1:
    # run parameter file
    f = bornprofiler.analysis.get_files(args[0], basedir=opts.basedir)
    opts.jobName = f['jobName']
    opts.ionName = f['ionName']
    args = [f['samplepoints']] + f['datafiles']
  elif len(args) == 2:
    # maybe the shell did not expand globs or we run in ipython?
    samplepoints,fileglob = args
    args = [samplepoints] + glob.glob(fileglob)

  if opts.jobName is None:
    opts.jobName = "bornprofile"

  kwargs = {'jobName': opts.jobName, 'create': opts.create}
  A = AnalyzeElec(*args, **kwargs)
  A.write()
  A.plot()

  if opts.pdbfilename:
    if opts.pdbfilename == "auto":
      opts.pdbfilename = None
    A.export(filename=opts.pdbfilename, format="pdb", ion=opts.ionName)    

  bornprofiler.stop_logging()
