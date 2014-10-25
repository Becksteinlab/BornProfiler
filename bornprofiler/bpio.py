# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
# Copyright (c) 2005-2008 Kaihsu Tai, Oliver Beckstein
# Copyright (c) 2010-2013 Oliver Beckstein
"""
Input/Output for the BornProfiler modules
=========================================

The module contains functions and classes to handle common
functionality to read and write files.
"""

from __future__ import with_statement
import os.path, errno
from ConfigParser import SafeConfigParser
import numpy

import config

import logging
logger = logging.getLogger("bornprofiler.io")

def read_dat_to_array(datfile):
    """Reads the dat file into a (4,N) numpy array"""
    infolist = []
    dat_file = open(datfile)
    dat_file.readline()
    for line in dat_file:
       floatline = [float(x) for x in line.split()]
       infolist.append(floatline)
    return (numpy.array(infolist))

def write_dat_from_array(outName, datinfo):
    with open(outName, "w") as outFile:
        outFile.write("# x/A y/A z/A  W/(kJ/mol)\n")
        for x,y,z,E in datinfo:
            outFile.write("{x:>8,.3f} {y:>8,.3f} {z:>8,.3f} {E:>8,.3f}\n".format(x=x,y=y,z=z,E=E))
    logger.info("Wrote Born PMF to {outName}.".format(outName=outName))


def path(s):
    """Return *s* with user expansion."""
    return os.path.expanduser(s)

def float_or_None(s):
    """Return *s* as float or None when x == "None"."""
    if s == "None":
        return None
    return float(s)

class RunParameters(object):
    """All parameters for a BornProfiler or APBSmem run are stored in a INI-style file.

    This class accesses these parameters and can create a template with the
    default values.

    The class *guarantees* that all parameters exist; if needed they are
    populated with default values. Hence it is always possible to access a
    parameter without having to check if it is there.

   :Parameters:
        - zmem : membrane centre (A)
        - lmem : membrane thickness (A)
        - Vmem : cytosolic potential in kT/e (untested)
        - pdie : protein dielectric
        - sdie:  solvent dielectric
        - mdie:  membrane dielectric
        - Rtop : exclusion cylinder top
        - Rbot : exclusion cylinder bottom
        - x0_R : exclusion zone centre in X, ``None`` selects the default
        - y0_R : exclusion zone centre in Y, ``None`` selects the default
        - dx_R : shift centre of the exclusion zone in X
        - dy_R : shift centre of the exclusion zone in Y
        - cdie : dielectric in the channel (e.g. SDIE)
        - headgroup_l : thicknes of headgroup region
        - headgroup_die : dielectric for headgroup region
        - temperature : temperature
        - conc : ionic strength in mol/l
        - ... and more

    .. SeeAlso:; The example run input configurations file
       :doc:`examples/example_runinput.cfg` is the full specification and
       contains annotated values for *all* parameters.
    """
    # For each task (cf the apbs-* scripts!) we define the variables that we
    # want to pull from the run input config file in the parameter_selections
    # dict. Keys are tasks, values are dicts that correspond to sections in the
    # cfg file, together with parameter names *and* types (for the conversion
    # to python types).

    # When new parameter are added to the cfg file then they need to be added
    # in various places:
    #
    # 1. examples/example_runinput.cfg (documentation/specification! :-p )
    # 2. RunParameters.parameter_selections (possibly multiple times!)
    # 3. RunParameters._populate_default()
    # 4. membrane.BaseMem: set as attribute
    # 5. membrane.APBSnomem.vars: name of the attribute if needed for a task
    # 6. membrane.APBSmem.vars: name of the attribute if needed for a task
    # 7. membrane.BornAPBSmem.vars: name of the attribute if needed for a task

    # The following should all be in a template file, together with the
    # defaults and types...

    # Note: The parameter keys are converted to lowercase when accessing
    # the runparameter file but need to be properly cased in the kwargs
    # list because this is how they are being used in the downstream code.
    # section -> key, convertor_function
    #     float, int, str: simple values
    #     path:            expand ~ etc
    #     eval:            python expressions such as lists or tuples (NOT SAFE!!)
    parameter_selections = {
        'bornprofile':
            {'environment': [('temperature', float), ('conc', float), ('pdie', float),
                             ('sdie', float), ('pqr', path), ('runtype', str)],
             'membrane':    [('Rtop', float), ('Rbot', float), ('cdie', float),
                             ('x0_R', float_or_None), ('y0_R', float_or_None),
                             ('dx_R', float), ('dy_R', float),
                             ('headgroup_die', float), ('headgroup_l', float), ('Vmem', float),
                             ('lmem', float), ('zmem', float), ('mdie', float)],
             'bornprofile': [('ion', str), ('dime', eval), ('glen', eval), ('fglen', eval),
                             ('points', path)],
             'job': [('name', str), ('script', path), ('arrayscript', path)],
             'plotting': [('xcolumn',int),('ycolumn',int),('title',str),
                          ('xlabel',str), ('ylabel',str),('plot_label',str),
                          ('color',str),('protein_bottom',float),
                          ('protein_length', float)]
            },
        'bornprofilenomem':
            {'environment': [('temperature', float), ('conc', float), ('pdie', float),
                             ('sdie', float), ('pqr', path), ('runtype', str)],
             'bornprofile': [('ion', str), ('dime', eval), ('glen', eval), ('fglen', eval),
                             ('points', path)],
             'job': [('name', str), ('script', path), ('arrayscript', path)],
             'plotting': [('xcolumn',int),('ycolumn',int),('title',str),
                          ('xlabel',str), ('ylabel',str),('plot_label',str),
                          ('color',str),('protein_bottom',float),
                          ('protein_length', float)]
            },
        'apbsmem':
            {'environment': [('temperature', float), ('conc', float), ('pdie', float),
                             ('sdie', float), ('pqr', path),],
             'membrane':    [('Rtop', float), ('Rbot', float), ('cdie', float),
                             ('x0_R', float_or_None), ('y0_R', float_or_None),
                             ('dx_R', float), ('dy_R', float),
                             ('headgroup_die', float), ('headgroup_l', float), ('Vmem', float),
                             ('lmem', float), ('zmem', float), ('mdie', float)],
             'potential':   [('dime', eval), ('glen', eval),],
             },
        'apbsnomem':
            {'environment': [('temperature', float), ('conc', float), ('pdie', float),
                             ('sdie', float), ('pqr', path),],
             'potential':   [('dime', eval), ('glen', eval),],
             },
        }

    def __init__(self, filename,*omissions, **defaults):
        """Reads and parses the configuration file for a job."""
        self.filename = filename
        self.parser = SafeConfigParser(defaults)
        # set the defaults and override, omitting omissions...
        self._populate_default(*omissions)

        # start with user configuration:
        # (I could add the default values in a template file here... but still
        # need a way to determine the type for entries)
        self.parser.read(config.CONFIGNAME)
        if not os.path.exists(filename):
            self.write(filename)
            logger.info("Created new run input configuration %(filename)r", vars())
        else:
            # and then override with values from local file
            self.parser.readfp(open(filename))

        logger.info("Read run input configuration from %(filename)r", vars())
        if self.parser.get('bornprofile','ion') == 'H':
            logger.info("WARNING: submitted config file has H as ion. This will result in use of the mostly meaningless H+ proton. Try H30 for a more meaningful calculation")

    def _populate_default(self, nomembrane, **kwargs):
        # NOTE: - the parser turns all keys into *lowercase*
        #       - values must be strings
        #       - hack: python types are defined via external dicts
        #         (see self.bornprofile_parameters and
        #         get_bornprofile_kwargs())
        parser = kwargs.get('parser', None)
        
        if parser is None:
            parser = self.parser
        # can use %(basedir)s in other entries
        parser.set('DEFAULT', 'basedir', os.path.realpath(os.path.curdir))
        parser.set('DEFAULT', 'solvent_dielectric', '80')
        if nomembrane:
            pass
        else:
            parser.add_section('membrane')
            parser.set('membrane', 'Rtop', '0')
            parser.set('membrane', 'Rbot', '0')
            parser.set('membrane', 'x0_R', 'None')
            parser.set('membrane', 'y0_R', 'None')
            parser.set('membrane', 'dx_R', '0')
            parser.set('membrane', 'dy_R', '0')
            parser.set('membrane', 'cdie', '%(solvent_dielectric)s')
            parser.set('membrane', 'headgroup_die', '20')
            parser.set('membrane', 'headgroup_l', '0')
            parser.set('membrane', 'mdie', '2')
            parser.set('membrane', 'Vmem', '0')
            parser.set('membrane', 'lmem', '0')
            parser.set('membrane', 'zmem','0')
        parser.add_section('environment')
        parser.set('environment', 'pqr', 'protein.pqr')
        parser.set('environment', 'temperature', '298.15')
        parser.set('environment', 'conc', '0.1')
        parser.set('environment', 'pdie', '10')
        parser.set('environment', 'sdie', '%(solvent_dielectric)s')
        parser.set('environment', 'runtype', 'with_protein')  # alternative: mem_only
        parser.add_section('bornprofile')
        parser.set('bornprofile', 'ion', 'Na')
        parser.set('bornprofile', 'dime', '[(129,129,129),(129,129,129),(129,129,129)]')
        parser.set('bornprofile', 'glen', '[(250,250,250),(100,100,100),(50,50,50)]')
        parser.set('bornprofile', 'fglen', '(40,40,40)')
        parser.set('bornprofile', 'points', 'points.dat')
        parser.add_section('potential')   # for membrane.APBSmem
        parser.set('potential', 'dime', '(97,97,97)')
        parser.set('potential', 'glen', '(200,200,200)')
        parser.add_section('executables')
        #should come from global file
        #parser.set('executables', 'drawmembrane', 'draw_membrane2a')
        #parser.set('executables', 'apbs', 'apbs')
        parser.add_section('job')
        parser.set('job', 'name', 'mbornprofile')
        parser.set('job', 'script', 'q_local.sh')
        parser.set('job', 'arrayscript', 'q_array.sge')
        parser.add_section('plotting')
        parser.set('plotting','xcolumn','2')
        parser.set('plotting','ycolumn','3')
        parser.set('plotting','title','BP')
        parser.set('plotting','xlabel','z')
        parser.set('plotting','ylabel','W_elec')
        parser.set('plotting','plot_label','Ion')
        parser.set('plotting','color','black')
        parser.set('plotting','protein_bottom','-20')
        parser.set('plotting','protein_length','40')

    def _get_kwargs(self, *args, **kwargs):
        """Prepare kwargs for a specified task."""
        task = args[0]
        args = args[1:]
        parameter_selection = self.parameter_selections[task]

        kw = {}
        for section,parameters in parameter_selection.items():
            for option,convertor in parameters:
                try:
                    kw[option] = convertor(self.parser.get(section,option,vars=kwargs))
                except:
                    logger.error("Problem obtaining required parameter "
                                 "[%(section)s] %(option)s from "+str(self.filename), vars())
                    raise
        if len(args) == 1:
            return kw[args[0]]
        elif len(args) > 0:
            return [kw[k] for k in args]
        return kw

    def get_bornprofile_kwargs(self, *args, **kwargs):
        """Return a dict with kwargs appropriate for :class:`bornprofile`.

        Default values can be supplied in *kwargs*. This method picks unique
        parameter keys from the relevant sections of the run parameters file
        (i.e. *bornprofile*, *membrane*, and *environment*).

        If args are provided, then either a single value corresponding
        to the key or a list of such values is returned instead.
        """
        _args = ('bornprofile',) + args
        return self._get_kwargs(*_args, **kwargs)

    def get_bornprofilenomem_kwargs(self, *args, **kwargs):
        """Return a dict with kwargs appropriate for :class:`bornprofilenomem`.

        Default values can be supplied in *kwargs*. This method picks unique
        parameter keys from the relevant sections of the run parameters file
        (i.e. *bornprofile*, *membrane*, and *environment*).

        If args are provided, then either a single value corresponding
        to the key or a list of such values is returned instead.
        """
        _args = ('bornprofilenomem',) + args
        return self._get_kwargs(*_args,**kwargs)

    def get_apbsnomem_kwargs(self, *args, **kwargs):
        """Return a dict with kwargs appropriate for :class:`membrane.APBSnomem`.

        Default values can be supplied in *kwargs*. This method picks unique
        parameter keys from the relevant sections of the run parameters file
        (i.e. *bornprofile* and *environment*).

        If args are provided, then either a single value corresponding
        to the key or a list of such values is returned instead.
        """
        _args = ('apbsnomem',) + args
        return self._get_kwargs(*_args, **kwargs)

    def get_apbsmem_kwargs(self, *args, **kwargs):
        """Return a dict with kwargs appropriate for :class:`membrane.APBSmem`.

        Default values can be supplied in *kwargs*. This method picks unique
        parameter keys from the relevant sections of the run parameters file
        (i.e. *bornprofile*, *membrane*, and *environment*).

        If args are provided, then either a single value corresponding
        to the key or a list of such values is returned instead.
        """
        _args = ('apbsmem',) + args
        return self._get_kwargs(*_args, **kwargs)

    def write(self, filename=None):
        """Write the current parameters to *filename*."""
        if filename is None:
            filename = self.filename
        with open(filename, 'w') as f:
            self.parser.write(f)
        return filename

def readPoints(filename):
    """Read coordinates from either data file or pdb.

    Returns (N,3) array.
    """
    try:
        points = readPointsDat(filename)
    except ValueError:
        points = readPointsPDB(filename)
    logger.info("Read points from %(filename)r.", vars())
    return points

def readPointsDat(filename):
    """Read points from a simple data file.

    Example::
       # comment
       x y z
       x y z
       ...
    """
    # or could just use numpy.loadtxt(filename) ...
    points = []
    with open(filename) as pointsFile:
        for line in pointsFile:
            line = line.strip()
            if line.startswith('#') or len(line) == 0:
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError("%(filename)r must contain at least 3 entries x y z per line. Offending line was %(line)r." % vars())
            points.append(map(float, fields[0:3]))
    return numpy.array(points)

def readPointsPDB(filename):
    """Read points form a PDB formatted file.

    Takes x,y,z from any ATOM or HETATM record.
    """
    points = []
    with open(filename) as pointsFile:
        for line in pointsFile:
            line = line.strip()
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            points.append((x,y,z))
    return numpy.array(points)

def read_template(filename):
  """Return *filename* as one string.

  *filename* can be one of the template files.
  """
  # XXX: add a cache?? Would have to be careful with user supplied files.
  fn = config.get_template(filename)
  logger.debug("Reading file %(filename)r from %(fn)r.", vars())
  return "".join(file(fn).readlines())


class PQRReader(object):
    """Naive implementation of a PQR reader."""
    def __init__(self, filename, **kwargs):
        self.pqrName = filename
        self.pqrLines = []  # verbatim copy of all ATOM lines, some functions just use that
        with open(self.pqrName, "r") as pqrFile:
            for line in pqrFile:
                if (line[0:4] == "ATOM"):
                    self.pqrLines.append(line)
        _coords = []
        for line in self.pqrLines:
            fields = line.split()      # This depends on correct spacing in the PQR file.
            try:
                _coords.append(map(float, fields[5:8]))
            except:
                logger.fatal("Problem with PQR file format of file %(pqrName)r", vars(self))
                logger.fatal("Offending line: %s", line)
                raise

        self.coords = numpy.array(_coords)
        self.centroid = self.coords.mean(axis=0)

        logger.info("PQRReader: Read %(pqrName)r with centroid = %(centroid)r", vars(self))

    @property
    def coordinates(self):
        return self.coords

    @property
    def centerOfGeometry(self):
        return self.centroid
