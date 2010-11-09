# APBS BornProfiler
# -*- coding: utf-8 -*-
"""
Configuration file handling and setup
=====================================

The user can set important paths in ``~/.bornprofiler.cfg``. The file
is automatically created with default values if it does not exist.

Edit the file with a text editor.

In order to restore the default values, simply delete the config file.


Format
------

The configuration file is in "INI" format: Section titles and variable
assignments::

  [executables]
  apbs = apbs
  drawmembrane = /path/to/drawmembrane2

  [membrane]
  class = APBSmem

Meaning of variables
--------------------

executables
   paths to binaries or just the name if they are found via
   :envvar:`PATH`
membrane
   configuration variables for apbs-mem-setup and friends
     - *class*: name of a Python class that is derived from
       :class:`bornprofiler.membrane.APBSmem`; needs to be in a file named
       ``custom.py`` in the current directory or in :mod:`bornprofiler.custom`
       (This is a bit of a hack...)


Accessing the configuration
---------------------------

Important variables are stored in the dict :data:`configuration`. Any
variable can be accessed via the getter method of the
:class:`ConfigParser.SafeConfigParser` instance, ``cfg``::

  from bornprofiler.config import cfg
  varname = cfg.get(section, varname)

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
        cfg.add_section('membrane')
        cfg.set('membrane', 'class','APBSmem')
        cfg.add_section('executables')
        cfg.set('executables', 'drawmembrane', 'drawmembrane2')
        cfg.set('executables', 'apbs', 'apbs')
        with open(filename, 'w') as configfile:
            cfg.write(configfile)  # write the default file
    else:
        cfg.readfp(open(filename))

    return {'apbs': cfg.get('executables', 'apbs'),
            'drawmembrane': cfg.get('executables', 'drawmembrane'),
            'configfilename': filename,
            }
    
#: Dict containing important configuration variables, populated by 
#: :func:`get_configuration`.
configuration = get_configuration()    
