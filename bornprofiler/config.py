# APBS BornProfiler
# -*- coding: utf-8 -*-
"""
Configuration file handling and setup
=====================================

The user can set important paths in ~/.bornprofiler.cfg

  [executables]
  apbs = apbs
  drawmembrane = /path/to/drawmembrane2

Variables are accessed in the dict :data:`configuration`.

"""
from __future__ import with_statement

import os.path
from ConfigParser import SafeConfigParser

#: Default name for the configuration file.
CONFIGNAME = os.path.expanduser(os.path.join("~",".bornprofiler.cfg"))

#: Instance of :class:`ConfigParser.SafeConfigParser`.
cfg = SafeConfigParser()

def get_configuration(filename=CONFIGNAME):
    """Reads and parses the configuration file."""
    if not os.path.exists(filename):
        cfg.add_section('executables')
        cfg.set('executables', 'drawmembrane', 'drawmembrane2')
        cfg.set('executables', 'apbs', 'apbs')
        with open(filename, 'w') as configfile:
            cfg.write(configfile)  # write the default file
    else:
        cfg.readfp(open(filename))

    return {'apbs': cfg.get('executables', 'apbs'),
            'drawmembrane': cfg.get('executables', 'drawmembrane'),
            }
    
#: Dict containing important configuration variables, populated by 
#: :func:`get_configuration`.
configuration = get_configuration()    
