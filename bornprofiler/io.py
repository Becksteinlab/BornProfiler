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

import config

import logging
logger = logging.getLogger("bornprofiler.io")

def path(s):
    return os.path.expanduser(s)

class RunParameters(object):
    """All parameters for a BornProfiler or APBSmem run are stored in a INI-style file.

    This class accesses these parameters and can create a template with the
    default values.

    .. Note:: sdie=80 and mdie=2 are fixed in draw_membrane2.c at the moment

   :Parameters:
        - zmem : membrane centre (A)
        - lmem : membrane thickness (A)
        - Vmem (untested)
        - pdie : protein dielectric
        - Rtop : exclusion cylinder top
        - Rbot : exclusion cylinder bottom
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
    #     eval:            python expressions such as lists or tuples (NOT SAFE!!)
    bornprofile_parameters = \
        {'environment': [('temperature', float), ('conc', float), ('pdie', float),
                         ('sdie', float), ('pqr', path)],
         'membrane':    [('Rtop', float), ('Rbot', float), ('headgroup_die', float),
                         ('headgroup_l', float), ('mdie', float), ('Vmem', float),
                         ('lmem', float), ('zmem', float)],
         'bornprofile': [('ion', str), ('dime', eval), ('glen', eval), ('points', path)],
         'job': [('name', str), ('script', path), ('arrayscript', path)],
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
        parser.add_section('membrane')
        parser.set('membrane', 'Rtop', '0.0')
        parser.set('membrane', 'Rbot', '0.0')
        parser.set('membrane', 'headgroup_die', '20.0')
        parser.set('membrane', 'headgroup_l', '0.0')
        parser.set('membrane', 'mdie', '2.0')
        parser.set('membrane', 'Vmem', '0.0')
        parser.set('membrane', 'lmem', '40.0')
        parser.set('membrane', 'zmem','0.0')
        parser.add_section('environment')
        parser.set('environment', 'pqr', 'protein.pqr')
        parser.set('environment', 'temperature', '298.15')
        parser.set('environment', 'conc', '0.1')
        parser.set('environment', 'pdie', '10.0') # hardcoded!
        parser.set('environment', 'sdie', '80.0') # hardcoded in drawmembrane2a at the moment
        parser.add_section('bornprofile')
        parser.set('bornprofile', 'ion', 'Na')
        parser.set('bornprofile', 'dime', '[(129,129,129),(129,129,129),(129,129,129)]')
        parser.set('bornprofile', 'glen', '[(250,250,250),(100,100,100),(50,50,50)]')
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
        

    def get_bornprofile_kwargs(self, *args, **kwargs):
        """Return a dict with kwargs appropriate for :class:`membrane.BornAPBSmem`.

        Default values can be supplied in *kwargs*. This method picks unique
        parameter keys from the relevant sections of the run parameters file
        (i.e. *bornprofile*, *membrane*, and *environment*).

        If args are provided, then either a single value corresponding
        to the key or a list of such values is returned instead.
        """
        kw = {}
        for section,parameters  in self.bornprofile_parameters.items():
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

    def write(self, filename=None):
        """Write the current parameters to *filename*."""
        if filename is None:
            filename = self.filename
        with open(filename, 'w') as f:
            self.parser.write(f)
        return filename
