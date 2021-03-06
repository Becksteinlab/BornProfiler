#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""
:Author: Oliver Beckstein
:Year: 2010
:Licence: GPL
:Copyright: (c) 2010-2013 Oliver Beckstein
"""

usage = """%prog [options] samplepoints-file *.out
       %prog [options] parameter-file

Extract the electrostatic free energy from the numbered APBS output files
(produced via the placeion.py script) and associate each energy with the position of the ion. 

A 3D electrostatic free energy landscape is constructed from the
values at the sample points.

.. Note:: The same samplepoints-file must be provided that was used for setting
   up the APBS calculations.
"""

import shutil
import os

import bornprofiler
from bornprofiler.analysis import AnalyzeElec3D

import logging
logger = logging.getLogger('bornprofiler') 
 
if __name__ == "__main__":
  import sys
  import os
  import glob
  from optparse import OptionParser


  parser = OptionParser(usage=usage)
  parser.add_option("--name", dest="jobName",
                    metavar="STRING",
                    help="name used in output files (welec_STRING.{dat,pdf})")
  parser.add_option("--delta", "-d", dest="delta", type="float",
                    metavar="FLOAT",
                    help="When exporting to a DX grid, sample points on bins with "
                    "size FLOAT Angstroem. Should be the value used for building the "
                    "sample points (e.g. Hollow's grid_spacing).")
  parser.add_option("--dx", "-x", dest="dxfilename", 
                    metavar="FILE",
                    help="Export grid to dx file FILE (see --delta). If set to 'auto' "
                    "then a filename is chosen.")
  parser.add_option("--pdb", "-p", dest="pdbfilename", 
                    metavar="FILE",
                    help="Export points to a PDB file with the energy in the B-factor. "
                    "If set to 'auto' then a filename is chosen.")
  parser.add_option("--copy-pqr", "-q", dest="copy_pqr", action="store_true",
                    help="copy the input pqr file to this directory for later visualisation")
  parser.add_option("--ion", dest="ionName",
                    metavar="STRING",
                    help="Set the ion name for --pdb if not running from config file [%default]")
  parser.add_option("--basedir", "-B", dest="basedir",
                    metavar="DIR",
                    help="when using a run parameter file, the job output is found under "
                    "'DIR/<job.name>/w[0-9][0-9][0-9][0-9].../job*.out', with job.name taken from "
                    "the run parameter file [%default]")
  parser.add_option("--read", dest="create", action="store_false",
                    help="If set, read positions and enrgies from a previously created dat file.")
  parser.set_defaults(basedir=os.path.curdir, ionName="ION", create=True, copy_pqr=False)

  opts,args = parser.parse_args()

  bornprofiler.start_logging()

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
    if opts.copy_pqr:
      logger.warn("--copy-pqr: No pqr information when running without cfg file. Ignored.")

  if opts.jobName is None:
    opts.jobName = "bornprofile"

  kwargs = {'jobName': opts.jobName, 'create': opts.create}
  A = AnalyzeElec3D(*args, **kwargs)
  A.write()

  if opts.dxfilename:
    if opts.delta:
      if opts.dxfilename == "auto":
        opts.dxfilename = None
      A.export(filename=opts.dxfilename, format="dx", delta=opts.delta)
    else:
      logger.warn("The --dx option was set but no --delta SPACING provided. No dx file will be written.")

  if opts.pdbfilename:
    if opts.pdbfilename == "auto":
      opts.pdbfilename = None
    A.export(filename=opts.pdbfilename, format="pdb", ion=opts.ionName)    

  if opts.copy_pqr:
    shutil.copy(f['pqr'], os.curdir)
    logger.info("cp %r %s", f['pqr'], os.curdir)

  bornprofiler.stop_logging()
