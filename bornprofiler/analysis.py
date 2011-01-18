# BornProfiler -- analysis classes, typically used in scripts
# Copyright (c) 2008 Kaihsu Tai
# Copyright (c) 2011 Oliver Beckstein
from __future__ import with_statement

import numpy
import bornprofiler
from bornprofiler.core import BPbase
 
import logging
logger = logging.getLogger('bornprofiler.analysis') 

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
        logger.info("Wrote Born profile to %(outName)r.", vars()) 

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

    def __repr__(self):
        return "<%s jobName=%r points=%r, %d data files>" % \
            (self.__class__.__name__, self.jobName, self.pointsName, len(self.datafiles))


class AnalyzeElec3D(BPbase):
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
        print
        self.data = numpy.array(data).T

 
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
 

