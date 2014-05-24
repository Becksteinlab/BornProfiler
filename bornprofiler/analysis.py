# BornProfiler -- analysis classes, typically used in scripts
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
# Copyright (c) 2008 Kaihsu Tai
# Copyright (c) 2011, 2012, 2013 Oliver Beckstein
from __future__ import with_statement

import os.path
import glob

import numpy
import bornprofiler
from bornprofiler.core import BPbase

import logging
logger = logging.getLogger('bornprofiler.analysis')

def get_files(runfile, basedir=os.path.curdir):
    """Read *runfile* and return dict with input files and names."""
    # convenience function --- one could also just read the cfg file
    from bornprofiler.bpio import RunParameters
    try:
        p = RunParameters(runfile,False,False)
        samplepoints = p.get_bornprofile_kwargs('points')
        fileglob = os.path.join(basedir, p.get_bornprofile_kwargs('name'),
                                'w[0-9][0-9][0-9][0-9]', 'job*.out')
        jobName = p.get_bornprofile_kwargs('name')
        ionName = p.get_bornprofile_kwargs('ion')
        pqr = p.get_bornprofile_kwargs('pqr')
    except:
        logger.fatal("Cannot obtain information about the sample points and "
                     "directory from the run parameter file %r.", runfile)
        raise
    return {'samplepoints': samplepoints,
            'datafiles': glob.glob(fileglob),  # gets all individual windows!
            'jobName': jobName,
            'ionName': ionName,
            'pqr': pqr,
            }

class Analyzer(BPbase):
    """Base class for analysis

    1. read points
    2. read energy for each point
    3. write data file(s)
    4. optional: plot

    Developer note: :attr:`exporters` contains the logic for exporting
    to various other formats besides plain column/text such as
    PDB. For additional exporters, add the entries in the local
    __init__ *after* calling ``super(name,
    self).__init__(*args,**kwargs)``.

    The Born energy W at each point (x,y,z) is stored in :attr:`data`,
    a (4,N) numpy array::

        X,Y,Z,W = data
    """
    def __init__(self, *args, **kwargs):
        self.pointsName = args[0]
        self.datafiles = args[1:]
        self.jobName = kwargs.pop('jobName', 'bornprofile')

        create = kwargs.pop("create", True)
        if create:
            self.readPoints()
            if self.points.shape[0] != len(self.datafiles):
                msg = "Number of sampled points (%d) does not match the number "\
                    "of data files (%d). Some windows are MISSING and will be skipped." % \
                    (self.points.shape[0], len(self.datafiles))
                logger.warn(msg)
                missing = self.find_missing_windows()
                for num,x,y,z,path in missing:
                    logger.warn("missing: %4d   (%8.3f,%8.3f,%8.3f)  taskid=%d",
                                num, x, y, z, self.get_taskid(num))
            self.accumulate()
        else:
            self.read()

        # other classes can add their exporters in __init__
        self.exporters = {'pdb': {'ext': '.pdb',
                                  'exporter': self._export_pdb},
                          }
        super(Analyzer, self).__init__(**kwargs)

    def read(self, filename=None):
        """Read datafile *filename* (format: x,y,z,W or z,W).

        The data are stored in :attr:`data`, a (4,N) array. If only z
        coordinates are provided then the x and y columns(data[0] and
        data[1]) are set to 0.
        """
        if filename is None:
            filename = self.datafile("welec")
        data = numpy.loadtxt(filename).T
        if data.shape[0] == 2:
            logger.warn("Reading old-style data file (z,W) %(filename)r.", vars())
            self.data = numpy.zeros((4,data.shape[-1]), dtype=data.dtype)
            self.data[[-2,-1]] = data
        else:
            self.data = data
        logger.info("Read Born PMF from %(filename)r.", vars())

    def write(self):
        outName = self.datafile("welec")
        with open(outName, "w") as outFile:
            outFile.write("# x/A y/A z/A  W/(kJ/mol)\n")
            for x,y,z,E in self.data.T:
                outFile.write("%(x)8.3f %(y)8.3f %(z)8.3f %(E)8.3e\n" % vars())
        logger.info("Wrote Born PMF to %(outName)r.", vars())

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
                    # maybe use regex o make this more robust...
                    if (line[0:18] == "  Local net energy"):
                        data.append(numpy.concatenate((point, [float(line.split()[6])])))
                        break
            print "[%5.1f%%] Read point %6d/%6d  %r \r" % \
                (100. * float(num+1)/len(self.points),
                 num, len(self.points)-1, outPath),
        print
        self.data = numpy.array(data).T

    def find_missing_windows(self):
        """Return a list of job numbers that have no data corresponding to a point."""
        missing = []
        for num, point in enumerate(self.points):
            outName = self.outfilename(num)
            # find file name in list of input files
            outPath = None
            for path in self.datafiles:   # could make self.datafiles a hash for faster search
                if path.endswith(outName):
                    outPath = path
                    break                 # good, found it
            if outPath is None:
                missing.append((num, point[0], point[1], point[2], outName))
        return numpy.rec.fromrecords(missing, names="num,x,y,z,path")

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
        if filename is None:
            filename = self.datafile("welec", ext="."+format)
        else:
            format = os.path.splitext(filename)[1][1:]
        if not format in self.exporters:
            raise ValueError("%r format is unsupported, only %r work" % (format, self.exporters.keys()))
        return self.exporters[format]['exporter'](filename, **kwargs)

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


class AnalyzeElec(Analyzer):
    "analyze APBS energy profiling results"

    def plot(self, filename=None, plotter='matplotlib', **kwargs):
        """Plot Born profile.

        plot([filename[,plotter[,kwargs ...]]])

        :Keywords:
          *filename*
             name of image file to save the plot in; with *plotter* = 'Gnuplot'
             only eps files are supported (I think...)
          *kwargs*
             other keyword arguments that are passed on to :func:`pylab.plot`
        """
        from pylab import plot, xlabel, ylabel, savefig

        if filename is None:
            plotName = self.datafile("welec",".pdf")
        else:
            plotName = filename
        kwargs.setdefault('color', 'black')
        kwargs.setdefault('linewidth', 2)
        z,W = self.data[-2], self.data[-1]  # should work for (4,N) and (2,N) data
        plot(z,W, **kwargs)
        xlabel(r'$z$ in nm')
        ylabel(r'$W$ in kJ$\cdot$mol$^{-1}$')
        savefig(plotName)
        logger.info("Plotted graph W(z) %(plotName)r.", vars())
        return plotName

    def __repr__(self):
        return "<%s jobName=%r points=%r, %d data files>" % \
            (self.__class__.__name__, self.jobName, self.pointsName, len(self.datafiles))

class MultiPlot1D(object):
    pass

class AnalyzeElec3D(Analyzer):
    """analyze APBS energy profiling results

    With create=True (default), reads position data from samplepoints
    file and energies from APBS output files.

    With create=False, reads positions and energies from a previous
    output file.
    """
    def __init__(self, *args, **kwargs):
        super(AnalyzeElec3D, self).__init__(*args, **kwargs)
        self.exporters['dx']=  {'ext': '.dx',
                                'exporter': self._export_dx}

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

    def Grid(self, delta, **kwargs):
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
        resample_factor = kwargs.pop('resample_factor', None)
        kwargs.setdefault('interpolation_spline_order', 3)
        g = Grid(*self.histogramdd(delta), **kwargs)
        if not resample_factor:
            return g
        return g.resample_factor(resample_factor)

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

