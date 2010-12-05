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
import numpy

import logging
logger = logging.getLogger('bornprofiler') 

usage = """%prog [options] samplepoints-file *.out
       %prog [options] parameter-file

Extract the electrostatic free energy from the numbered APBS output files
(produced via the placeion.py script) and associate each energy with the position of the ion. 

.. Note:: The same samplepoints-file must be provided that was used for setting
   up the APBS calculations.
"""

from bornprofiler.core import BPbase
 
class AnalyzeElec(BPbase):
  "analyze APBS energy profiling results"

  def __init__(self, *args, **kwargs):
    self.pointsName = args[0]
    self.datafiles = args[1:]
    self.jobName = kwargs.pop('jobName', 'bornprofile')
 
    self.readPoints()
    if self.points.shape[0] != len(self.datafiles):
      raise ValueError("Number of sampled points (%d) does not match the number "
                       "of data files (%d). They MUST correspond 1-to-1." % 
                       (self.points.shape[0], len(self.datafiles)))
    self.accumulate()
    self.write()

  def accumulate(self):
    zE = []
    for num, point in enumerate(self.points):
      z = point[2]
      outName = self.outfilename(num)
      lines = ""
      # find file name in list of input files
      outPath = None
      for path in self.datafiles:
        if path.endswith(outName):
          outPath = path
          break
      if outPath is None:
        logger.warn("Sample point %d %r: Could not find file %s." % 
                    (num, point, outName))
        continue
      with open(outPath) as outFile:
        for line in outFile:
          if (line[0:18] == "  Local net energy"):
            zE.append([z, float(line.split()[6])])
    self.zE = numpy.array(zE).T
 
  def write(self):
    outName = self.datafile("welec")
    with open(outName, "w") as outFile:
      outFile.write("# z/angstrom E/(kJ/mol)\n")
      for z,E in self.zE.T:
        outFile.write("%(z)8.3f %(E)8.3e\n" % vars())
    logger.info("Wrote Born profile to %(outName)r.", vars(self)) 

  def plot(self, filename=None, plotter='matplotlib', **kwargs):
    """Plot Born profile.

    plot([filename[,plotter[,kwargs ...]]])

    :Keywords:
      *filename*
         name of image file to save the plot in; with *plotter* = 'Gnuplot'
         only eps files are supported (I think...)
      *plotter* either 'matplotlib' 
         (default, requires matplotlib) or 'gnuplot' (requires Gnuplot package)
      *kwargs*
         other keyword arguments that are passed on to :func:`pylab.plot`; ignored
         for Gnuplot
    """
    plotters = {'matplotlib': self._plot_matplotlib,
                'Gnuplot': self._plot_Gnuplot,
                }
    kwargs['filename'] = filename
    plotName = plotters[plotter](**kwargs)
    logger.info("Plotted grap %(plotName)r.", vars())
    return plotName
 
  def _plot_matplotlib(self, filename=None, **kwargs):
    from pylab import plot, xlabel, ylabel, savefig
    if filename is None:
      plotName = self.datafile("welec",".pdf")
    else:
      plotName = filename
    kwargs.setdefault('color', 'black')
    kwargs.setdefault('linewidth', 2)
    z,E = self.zE
    plot(z,E, **kwargs)
    xlabel(r'$z$ in nm')
    ylabel(r'$W$ in kJ$\cdot$mol$^{-1}$')
    savefig(plotName)
    return plotName

  def _plot_Gnuplot(self, filename=None, **kwargs):
    import Gnuplot    
    outName = self.datafile("welec")
    if filename is None:
      plotName = self.datafile("welec",".eps")
    else:
      plotName = filename
 
    # initialize
    g = Gnuplot.GnuplotProcess()
    cmd  = "set terminal postscript eps colour\n"
    cmd += 'set output "' + plotName + '"\n'
    cmd += """set style data lines
set xlabel "z / nm"
set ylabel "energy / (kJ/mol)"
"""
    cmd += 'plot "%s" using ($1/10):($2) title "%s"\n' (outName,self.jobName)
 
    # do it
    g(cmd)
    return plotName
 
if __name__ == "__main__":
  import sys
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
