#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:Author: Oliver Beckstein
:Year: 2010
:Licence: GPL
:Copyright: (c) 2010 Oliver Beckstein
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

A 3D electrostatic free energy landscape is constructred from the
values at the sample points.

.. Note:: The same samplepoints-file must be provided that was used for setting
   up the APBS calculations.
"""

import bornprofiler
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

  def accumulate(self):
    """Read the energy from each datafile and store with the coordinates.

    Populates :attr:`AnalyzeElec.data`, a (4,N) array where N is the
    number of windows.
    """
    data = []
    for num, point in enumerate(self.points):
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
            data.append(numpy.concatenate((point, [float(line.split()[6])])))
            break
      print "[%5.1f%%] Read point %6d/%6d  %r \r" % (100. * float(num+1)/len(self.points), 
                                                     num, len(self.points)-1, outPath),
    self.data = numpy.array(data).T
    print
 
  def ranges(self, padding=0):
    """Returns the range of values in each dimension.

    :Returns: Array *r* of shape (2,3): r[0] contains the smallest and 
              r[1] the largest values.
    """
    def minmax(dim):
      return numpy.array([self.data[dim,:].min(), self.data[dim,:].max()])
    ranges = numpy.array([minmax(dim) for dim in (0,1,2)]).T
    ranges[0] -= padding
    ranges[1] += padding
    return ranges
    
  def histogramdd(self, delta, fillfac=5):
    """Histogram PMF on a regular grid.

    The spacing *delta* must be the same or larger than used in Hollow; to be
    on the safe side, just use the value of *grid_spacing* (e.g. 2.0 for 2 A).
    The histogram grid is chosen large enough to encompass all data points.

    If *delta* is bigger than *grid_spacing* or points are not on a
    regular grid the values are averaged over each bin.

    Points without data are filled with fillfac * max(histogram).

    :Returns: histogram, edges (same as :func:`numpy.histogramdd`)
    """
    ranges = self.ranges(padding=delta/2.)
    bins = ((ranges[1] - ranges[0])/float(delta)).astype(int)
    h,e = numpy.histogramdd(self.data[:3].T, range=ranges.T, bins=bins, 
                            weights=self.data[3])
    N,e = numpy.histogramdd(self.data[:3].T, range=ranges.T, bins=bins)    
    h[N>0] /= N[N>0]  # average PMF
    h[N==0] = fillfac*h.max()
    return h,e

  def Grid(self, delta, resample_factor=None, interpolation_spline_order=3):
    """Package the PMF as a :class:`gridData.Grid` object.

    *delta* should be the original spacing of the points in angstroem.

    With a *resample_factor*, the data are interpolated on a new grid
    with *resample_factor* times as many bins as in the original
    histogram (See :meth:`gridData.Grid.resample_factor`).

    *interpolation_order* sets the interpolation order for the
     resampling procedure. 

    .. Warning:: Interpolating can lead to artifacts in the 3D PMF. If
                 in doubt, **do not resample**.
    """
    from gridData import Grid
    g = Grid(*self.histogramdd(delta), interpolation_spline_order=interpolation_spline_order)
    if not resample_factor:
      return g
    return g.resample_factor(resample_factor)

  def write(self):
    outName = self.datafile("welec")
    with open(outName, "w") as outFile:
      outFile.write("# x/A y/A z/A  W/(kJ/mol)\n")
      for x,y,z,E in self.data.T:
        outFile.write("%(x)8.3f %(y)8.3f %(z)8.3f %(E)8.3e\n" % vars())
    logger.info("Wrote Born PMF to %(outName)r.", vars()) 

  def export(self, delta, filename=None, **kwargs):
    """Write data to dx file *filename*.

    export(delta[,filename[,kwargs]])

    Requires grid spacing delta of original data points. *filename* is
    automatically chosen if not supplied. Other *kwargs* are passed on
    to :meth:`AnalElec.Grid`.
    """
    if filename is None:
      filename = self.datafile("welec", ext=".dx")
    g = self.Grid(delta, **kwargs)
    g.export(filename, format="dx")
    logger.info("Wrote Born PMF to dx file  %(filename)r with spacing %(delta)g A.",
                vars())

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

  def __repr__(self):
    return "<%s jobName=%r points=%r, %d data files>" % \
        (self.__class__.__name__, self.jobName, self.pointsName, len(self.datafiles))
 
if __name__ == "__main__":
  import sys
  import glob
  from optparse import OptionParser

  bornprofiler.start_logging()

  parser = OptionParser(usage=usage)
  parser.add_option("--name", dest="jobName",
                    metavar="STRING",
                    help="name used in output files (welec_STRING.{dat,pdf})")
  parser.add_option("--delta", "-d", dest="delta", type="float",
                    metavar="FLOAT",
                    help="When exporting to a DX grid, sample points on bins with "
                    "size FLOAT Angstroem. Should be the value used for building the "
                    "sample points (e.g. Hollow's grid_spacing).")
  parser.add_option("--export", "-x", dest="dxfilename", 
                    metavar="FILE",
                    help="Export grid to dx file FILE (see --delta). If set to 'auto' "
                    "then a filename is chosen.")
  parser.add_option("--plotter", dest="plotter", type="choice",
                    choices=('matplotlib',),
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
  A.write()

  if opts.delta and opts.dxfilename:
    if opts.dxfilename == "auto":
      opts.dxfilename = None
    A.export(opts.delta, filename=opts.dxfilename)

  bornprofiler.stop_logging()
