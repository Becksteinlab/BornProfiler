# -*- coding: utf-8 -*-
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

def path(s):
    return os.path.expanduser(s)

class RunParameters(object):
    """All parameters for a BornProfiler or APBSmem run are stored in a INI-style file.

    This class accesses these parameters and can create a template with the
    default values.

   :Parameters:
        - zmem : membrane centre (A)
        - lmem : membrane thickness (A)
        - Vmem : cytosolic potential in kT/e (untested)
        - pdie : protein dielectric
        - sdie:  solvent dielectric
        - mdie:  membrane dielectric
        - Rtop : exclusion cylinder top
        - Rbot : exclusion cylinder bottom
        - cdie : dielectric in the channel (e.g. SDIE)
        - headgroup_l : thicknes of headgroup region
        - headgroup_die : dielectric for headgroup region
        - temperature : temperature
        - conc : ionic strength in mol/l  

        ... and others: see default file!
    """

    # this should all be in a template file, together with the defaults...

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
                             ('headgroup_die', float), ('headgroup_l', float), ('Vmem', float),
                             ('lmem', float), ('zmem', float), ('mdie', float)],
             'bornprofile': [('ion', str), ('dime', eval), ('glen', eval), ('fglen', eval),
                             ('points', path)],
             'job': [('name', str), ('script', path), ('arrayscript', path)],
             },
        'apbsmem':
            {'environment': [('temperature', float), ('conc', float), ('pdie', float),
                             ('sdie', float), ('pqr', path),],
             'membrane':    [('Rtop', float), ('Rbot', float), ('cdie', float),
                             ('headgroup_die', float), ('headgroup_l', float), ('Vmem', float),
                             ('lmem', float), ('zmem', float), ('mdie', float)],
             'potential':   [('dime', eval), ('glen', eval),],
             }
        }

    def __init__(self, filename, **defaults):
        """Reads and parses the configuration file for a job."""
        self.filename = filename
        self.parser = SafeConfigParser(defaults)

        # set the defaults and override ...
        self._populate_default()

        # start with user configuration:
        # (I could add the default values in a template file here... but still
        # need a way to determine the type for entries)
        self.parser.read(config.CONFIGNAME)

        if not os.path.exists(filename):
            self.write(filename)
        else:
            # and then override with values from local file
            self.parser.readfp(open(filename))    

        logger.info("Read run configuration from %(filename)r", vars())

    def _populate_default(self, parser=None):
        # NOTE: - the parser turns all keys into *lowercase*
        #       - values must be strings
        #       - hack: python types are defined via external dicts 
        #         (see self.bornprofile_parameters and
        #         get_bornprofile_kwargs())
        if parser is None:
            parser = self.parser
        # can use %(basedir)s in other entries
        parser.set('DEFAULT', 'basedir', os.path.realpath(os.path.curdir))
        parser.set('DEFAULT', 'solvent_dielectric', '80')
        parser.add_section('membrane')
        parser.set('membrane', 'Rtop', '0')
        parser.set('membrane', 'Rbot', '0')
        parser.set('membrane', 'cdie', '%(solvent_dielectric)s')
        parser.set('membrane', 'headgroup_die', '20')
        parser.set('membrane', 'headgroup_l', '0')
        parser.set('membrane', 'mdie', '2')
        parser.set('membrane', 'Vmem', '0')
        parser.set('membrane', 'lmem', '40')
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
        parser.set('executables', 'drawmembrane', 'draw_membrane2a')
        parser.set('executables', 'apbs', 'apbs')
        parser.add_section('job')
        parser.set('job', 'name', 'mbornprofile')
        parser.set('job', 'script', 'q_local.sh')
        parser.set('job', 'arrayscript', 'q_array.sge')
        
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
        """Return a dict with kwargs appropriate for :class:`membrane.BornAPBSmem`.

        Default values can be supplied in *kwargs*. This method picks unique
        parameter keys from the relevant sections of the run parameters file
        (i.e. *bornprofile*, *membrane*, and *environment*).

        If args are provided, then either a single value corresponding
        to the key or a list of such values is returned instead.
        """
        _args = ('bornprofile',) + args
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
