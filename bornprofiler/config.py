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

from pkg_resources import resource_filename, resource_listdir

import utilities

import logging
logger = logging.getLogger("bornprofiler.config")

# Processing of the configuration file
# ------------------------------------ 

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

# User-accessible configuration
# -----------------------------
#
# TODO: move into config file!!

#: Directory to store user templates and rc files.
#: The default value is ``~/.bornprofiler``.
configdir = os.path.expanduser(os.path.join("~",".bornprofiler"))

#: Directory to store user supplied queuing system scripts.
#: The default value is ``~/.bornprofiler/qscripts``.
qscriptdir = os.path.join(configdir, 'qscripts')

#: Directory to store user supplied template files such as mdp files.
#: The default value is ``~/.bornprofiler/templates``.
templatesdir = os.path.join(configdir, 'templates')

#: Directory to store user python modules.
libdir = os.path.join(configdir, 'lib')

def setup():
    """Create the directories in which the user can store template and config files.
    
    This function can be run repeatedly without harm.
    """
    # setup() must be separate and NOT run automatically when config
    # is loaded so that easy_install installations work
    # (otherwise we get a sandbox violation)
    utilities.mkdir_p(configdir)
    utilities.mkdir_p(qscriptdir)
    utilities.mkdir_p(templatesdir)

def check_setup():
    """Check if templates directories are setup and issue a warning and help."""
    missing = [d for d in (configdir, qscriptdir, templatesdir)
               if not os.path.exists(d)]
    if len(missing) > 0:
         print "NOTE: Some configuration directories are not set up yet"
         print "      %r" % missing
         print "      You can create them with the command"
         print "      >>> bornprofiler.config.setup()"
    return len(missing) == 0
check_setup()


#: Search path for user queuing scripts and templates. The internal package-supplied
#: templates are always searched last via :func:`bornprofiler.config.get_templates`. 
#: Modify :data:`bornprofiler.config.path` directly in order to customize the template 
#: and qscript searching. By default it has the value ``['.', qscriptdir,
#: templatesdir]``. 
#: (Note that it is not a good idea to have template files and qscripts with the
#: same name as they are both searched on the same path.)
path = [os.path.curdir, qscriptdir, templatesdir]


# Location of template files
# --------------------------

def _generate_template_dict(dirname):
   """Generate a list of included files *and* extract them to a temp space.

   Templates have to be extracted from the egg because they are used
   by external code. All template filenames are stored in
   :data:`config.templates`.
   """
   return dict((resource_basename(fn), resource_filename(__name__, dirname+'/'+fn))
               for fn in resource_listdir(__name__, dirname)
               if not fn.endswith('~'))

def resource_basename(resource):
   """Last component of a resource (which always uses '/' as sep)."""
   if resource.endswith('/'):
        resource = resource[:-1]
   parts = resource.split('/')
   return parts[-1]


#: Registry of all template files that come with the package.
templates = _generate_template_dict('templates')

# Functions to access configuration data
# --------------------------------------

def get_template(t):
   """Find template file *t* and return its real path.

   *t* can be a single string or a list of strings. A string
   should be one of

   1. a relative or absolute path,
   2. a file in one of the directories listed in :data:`bornprofiler.config.path`,
   3. a filename in the package template directory (defined in the template dictionary
      :data:`bornprofiler.config.templates`) or
   4. a key into :data:`~bornprofiler.config.templates`.

   The first match (in this order) is returned. If the argument is a
   single string then a single string is returned, otherwise a list
   of strings.

   :Arguments: *t* : template file or key (string or list of strings)
   :Returns:   os.path.realpath(*t*) (or a list thereof)
   :Raises:    :exc:`ValueError` if no file can be located.

   """
   templates = [_get_template(s) for s in utilities.asiterable(t)]
   if len(templates) == 1:
        return templates[0]
   return templates

def get_templates(t):
   """Find template file(s) *t* and return their real paths.

   *t* can be a single string or a list of strings. A string should
   be one of

   1. a relative or absolute path,
   2. a file in one of the directories listed in :data:`bornprofiler.config.path`,
   3. a filename in the package template directory (defined in the template dictionary
      :data:`bornprofiler.config.templates`) or
   4. a key into :data:`~bornprofiler.config.templates`.

   The first match (in this order) is returned for each input argument.

   :Arguments: *t* : template file or key (string or list of strings)
   :Returns:   list of os.path.realpath(*t*) 
   :Raises:    :exc:`ValueError` if no file can be located.

   """
   return [_get_template(s) for s in utilities.asiterable(t)]

def _get_template(t):
   """Return a single template *t*."""
   if os.path.exists(t):           # 1) Is it an accessible file?
        pass
   else:                         
        _t = t
        _t_found = False
        for d in path:              # 2) search config.path
             p = os.path.join(d, _t)
             if os.path.exists(p):
                  t = p
                  _t_found = True
                  break
        _t = os.path.basename(t)
        if not _t_found:            # 3) try template dirs
             for p in templates.values():
                  if _t == os.path.basename(p):
                       t = p
                       _t_found = True     # NOTE: in principle this could match multiple
                       break               #       times if more than one template dir existed.
        if not _t_found:            # 4) try it as a key into templates
             try:
                  t = templates[t]
             except KeyError:
                  pass
             else:
                  _t_found = True
        if not _t_found:            # 5) nothing else to try...
            errmsg = "Failed to locate the template file %(t)r." % vars()
            logger.fatal(errmsg)
            raise ValueError(errmsg)
   return os.path.realpath(t)

def read_template(filename):
  """Return *filename* as one string.

  *filename* can be one of the template files.
  """
  import config
  fn = config.get_template(filename)
  return "".join(file(fn).readlines())


