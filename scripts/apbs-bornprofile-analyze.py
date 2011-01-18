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
  parser.add_option("--plotter", dest="plotter", type="choice",
                    choices=('matplotlib','Gnuplot'),
                    help="plotting backend [%default]")
  parser.add_option("--basedir", dest="basedir",
                    metavar="DIR",
                    help="when using a run parameter file, the job output is found under "
                    "'DIR/<job.name>/w[0-9][0-9][0-9][0-9]/job*.out', with job.name taken from "
                    "the run parameter file [%default]")
  parser.set_defaults(plotter='matplotlib', basedir=os.path.curdir)


  opts,args = parser.parse_args()

  if len(args) == 0:
    logger.fatal("Needs samplepoints file and at least one APBS output file. See --help.")
    sys.exit(1)  
  elif len(args) == 1:
    # run parameter file
    from bornprofiler.io import RunParameters
    try:
      p = RunParameters(args[0])
      samplepoints = p.get_bornprofile_kwargs('points')
      fileglob = os.path.join(opts.basedir, p.get_bornprofile_kwargs('name'), 
                              'w[0-9][0-9][0-9][0-9]', 'job*.out')
      opts.jobName = p.get_bornprofile_kwargs('name')
    except:
      logger.fatal("Cannot obtain information about the sample points and directory from the "
                   "run parameter file %r.", args[0])
      raise
    args = [samplepoints] + glob.glob(fileglob)
  elif len(args) == 2:
    # maybe the shell did not expand globs or we run in ipython?
    samplepoints,fileglob = args
    args = [samplepoints] + glob.glob(fileglob)

  if opts.jobName is None:
    opts.jobName = "bornprofile"

  kwargs = {'jobName': opts.jobName}
  A = AnalyzeElec(*args, **kwargs)
  A.plot(plotter=opts.plotter)

  bornprofiler.stop_logging()
