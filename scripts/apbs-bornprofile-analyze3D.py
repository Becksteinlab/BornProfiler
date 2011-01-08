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

A 3D electrostatic free energy landscape is constructed from the
values at the sample points.

.. Note:: The same samplepoints-file must be provided that was used for setting
   up the APBS calculations.
"""

import bornprofiler
from bornprofiler.core import BPbase
 
class AnalyzeElec(BPbase):
  """analyze APBS energy profiling results

  With create=True (default), reads position data from samplepoints
  file and energies from APBS output files.

  With create=False, reads positions and energies from a previous
  output file.
  """

  def __init__(self, *args, **kwargs):
    self.pointsName = args[0]
    self.datafiles = args[1:]
    self.jobName = kwargs.pop('jobName', 'bornprofile')

    create = kwargs.pop("create", True)
    if create:
      self.readPoints()
      if self.points.shape[0] != len(self.datafiles):
        raise ValueError("Number of sampled points (%d) does not match the number "
                         "of data files (%d). They MUST correspond 1-to-1." % 
                         (self.points.shape[0], len(self.datafiles)))
      self.accumulate()
    else:
      self.read()

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

  def read(self, filename=None):
    """Read datafile *filename* (format: x,y,z,W)"""
    if filename is None:
      filename = self.datafile("welec")
    self.data = numpy.loadtxt(filename).T
    logger.info("Read Born PMF from %(filename)r.", vars())

  def write(self):
    outName = self.datafile("welec")
    with open(outName, "w") as outFile:
      outFile.write("# x/A y/A z/A  W/(kJ/mol)\n")
      for x,y,z,E in self.data.T:
        outFile.write("%(x)8.3f %(y)8.3f %(z)8.3f %(E)8.3e\n" % vars())
    logger.info("Wrote Born PMF to %(outName)r.", vars()) 

  def export(self, filename=None, format="dx", **kwargs):
    """Export data to different file format.

    The format is deduced from the filename suffix or
    *format*. *kwargs* are set according to the exporter.

    dx
       histogram on a grid; must provide kwarg *delta* for the 
       grid spacing.
    pdb
       write PDB file with an ion for each sample point and the 
       energy as the B-factor. The kwarg *ion* can be used to 
       set the name/resName in the file (ION is the default).
    """
    formats = {'dx': {'ext': '.dx', 
                      'exporter': self._export_dx},
               'pdb': {'ext': '.pdb',
                       'exporter': self._export_pdb},
               }
    if filename is None:
      filename = self.datafile("welec", ext="."+format)
    else:
      format = os.path.splitext(filename)[1][1:]
    if not format in formats:
      raise ValueError("%r format is unsupported, only %r work" % (format, formats.keys()))    
    return formats[format]['exporter'](filename, **kwargs)

  def _export_dx(self, filename, delta=2.0, **kwargs):
    """Write data to dx file *filename*.

    export(delta[,filename[,kwargs]])

    Requires grid spacing delta of original data points. *filename* is
    automatically chosen if not supplied. Other *kwargs* are passed on
    to :meth:`AnalElec.Grid`.
    """
    g = self.Grid(delta, **kwargs)
    g.export(filename, format="dx")
    logger.info("Wrote Born PMF to dx file  %(filename)r with spacing %(delta)g A.",
                vars())

  def _export_pdb(self, filename, **kwargs):
    """Write datapoints to pdb file with energy in the B-factor field.

    :Keywords:
      *ion*
         name of the ion (name/resSeq in pdb)
    """
    
    # http://www.wwpdb.org/documentation/format32/sect9.html
    fmt = {'ATOM':   "ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f\n",
           'REMARK': "REMARK     %s\n",
           'TITLE':  "TITLE    %s\n",
           'CRYST1': "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
           }
    name = resName = kwargs.pop('ion', 'ION') or 'ION'  # None --> ION
    chainID = "Z"
    iCode = altLoc = " "
    occupancy = 1.0

    with open(filename, 'w') as pdb:
      pdb.write(fmt['TITLE'] % "Poisson-Boltzmann electrostatic free energy")
      pdb.write(fmt['REMARK'] % "W(x,y,z) in kJ/mol")
      pdb.write(fmt['REMARK'] % "Written by apbs-bornprofile-analyze3D.py")
      for serial,(x,y,z,W) in enumerate(self.data.T):
        serial += 1
        resSeq = serial
        tempFactor = W
        pdb.write(fmt['ATOM'] % vars())
    logger.info("Wrote ion positions for %(name)s to pdb file %(filename)r", vars())
                  
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
  parser.add_option("--dx", "-x", dest="dxfilename", 
                    metavar="FILE",
                    help="Export grid to dx file FILE (see --delta). If set to 'auto' "
                    "then a filename is chosen.")
  parser.add_option("--pdb", "-p", dest="pdbfilename", 
                    metavar="FILE",
                    help="Export points to a PDB file with the energy in the B-factor. "
                    "If set to 'auto' then a filename is chosen.")
  parser.add_option("--ion", dest="ionName",
                    metavar="STRING",
                    help="Set the ion name for --pdb if not running from config file [%default]")
  parser.add_option("--basedir", dest="basedir",
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
    from bornprofiler.io import RunParameters
    try:
      p = RunParameters(args[0])
      samplepoints = p.get_bornprofile_kwargs('points')
      fileglob = os.path.join(opts.basedir, p.get_bornprofile_kwargs('name'), 
                              'w[0-9][0-9][0-9][0-9]', 'job*.out')
      opts.jobName = p.get_bornprofile_kwargs('name')
      opts.ionName = p.get_bornprofile_kwargs('ion')
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

  kwargs = {'jobName': opts.jobName, 'create': opts.create}
  A = AnalyzeElec(*args, **kwargs)
  A.write()

  if opts.delta and opts.dxfilename:
    if opts.dxfilename == "auto":
      opts.dxfilename = None
    A.export(filename=opts.dxfilename, format="dx", delta=opts.delta)
  if opts.pdbfilename:
    if opts.pdbfilename == "auto":
      opts.pdbfilename = None
    A.export(filename=opts.pdbfilename, format="pdb", ion=opts.ionName)    

  bornprofiler.stop_logging()
