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
  drawmembrane = /path/to/drawmembrane2a

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

import os.path, errno
from ConfigParser import SafeConfigParser

from pkg_resources import resource_filename, resource_listdir

import utilities

import logging
logger = logging.getLogger("bornprofiler.config")


APBS_MINIMUM_VERSION = 1,3   # want to be able to read gzipped files
DRAWMEMBRANE_REQUIRED_NAME = "draw_membrane2a.c"
DRAWMEMBRANE_MINIMUM_VERSION = 12,1,10  # date in MM/DD/YY (!)

# User-accessible configuration
# -----------------------------
#
# TODO: move into config file!!
# These are the default values; eventually they should be re-read from config

defaults = {}

#: Directory to store user templates and rc files.
#: The default value is ``~/.bornprofiler``.
defaults['configdir'] = os.path.expanduser(os.path.join("~",".bornprofiler"))

#: Directory to store user supplied queuing system scripts.
#: The default value is ``~/.bornprofiler/qscripts``.
defaults['qscriptdir'] = os.path.join(defaults['configdir'], 'qscripts')

#: Directory to store user supplied template files such as mdp files.
#: The default value is ``~/.bornprofiler/templates``.
defaults['templatesdir'] = os.path.join(defaults['configdir'], 'templates')


# Processing of the configuration file
# ------------------------------------ 

#: Default name for the configuration file.
CONFIGNAME = os.path.expanduser(os.path.join("~",".bornprofiler.cfg"))

#: Instance of :class:`ConfigParser.SafeConfigParser`.
class BPConfigParser(SafeConfigParser):
    def get_path(self, section, option):
        """Return option as an expanded path."""
        return os.path.expanduser(os.path.expandvars(self.get(section, option)))

#: :data:`cfg` is the instance of :class:`BPConfigParser` that makes all
#: global configuration data accessible
cfg = BPConfigParser()

def get_configuration(filename=CONFIGNAME):
    """Reads and parses the configuration file."""

    # defaults
    cfg.set('DEFAULT', 'configdir', defaults['configdir'])
    cfg.set('DEFAULT', 'qscriptdir', 
            os.path.join("%(configdir)s", os.path.basename(defaults['qscriptdir'])))
    cfg.set('DEFAULT', 'templatesdir', 
            os.path.join("%(configdir)s", os.path.basename(defaults['templatesdir'])))
    cfg.add_section('membrane')
    cfg.add_section('executables')
    cfg.set('executables', 'drawmembrane', 'draw_membrane2a')
    cfg.set('executables', 'apbs', 'apbs')

    if not os.path.exists(filename):
        with open(filename, 'w') as configfile:
            cfg.write(configfile)  # write the default file so that user can edit
        print(("NOTE: BornProfiler created the configuration file \n\t%r\n"+
               "      for you. Edit the file to customize the package.") % filename)
    else:
        # overwrite defaults
        cfg.readfp(open(filename))

    return {'apbs': cfg.get_path('executables', 'apbs'),
            'drawmembrane': cfg.get_path('executables', 'drawmembrane'),
            'configfilename': filename,
            }

#: Dict containing important configuration variables, populated by 
#: :func:`get_configuration` (mainly a shortcut; use :data:`cfg` in most cases)
configuration = get_configuration()    # also initializes cfg...

configdir = cfg.get_path('DEFAULT', 'configdir')
qscriptdir = cfg.get_path('DEFAULT', 'qscriptdir')
templatesdir = cfg.get_path('DEFAULT', 'templatesdir')


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
  fn = get_template(filename)
  logger.debug("Reading file %(filename)r from %(fn)r.", vars())
  return "".join(file(fn).readlines())


import subprocess
import re

def check_APBS(name=None):
    """Return ABPS version if apbs can be run and has the minimum required version.
    
    :Raises: error if it cannot be found (OSError ENOENT) or wrong version (EnvironmentError).
    """    
    APBS = name or cfg.get('executables', 'apbs')
    try:
        p = subprocess.Popen([APBS, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = p.communicate()
    except OSError:
        logger.critical("No APBS binary found --- set it in the config file %(configfilename)s", configuration)
        raise
    # do not check p.returncode because --verbose sets it to 13 (?) but no idea if this is a feature
    m = re.match('.*APBS\s*(?P<major>\d+)\.(?P<minor>\d+)', err)
    if m is None:
        raise EnvironmentError(errno.EIO, APBS,
                               "Cannot obtain APBS version string from %r." % err)        
    major,minor = int(m.group('major')), int(m.group('minor'))
    if not ((major,minor) >= APBS_MINIMUM_VERSION):
        raise EnvironmentError(errno.EIO, APBS, "APBS version %d.%d is too old, need at least %d.%d." % 
                               ((major,minor)+APBS_MINIMUM_VERSION))        
    return major,minor

def check_drawmembrane(name=None):
    """Return drawmembrane version or raise :exc:`EnvironmentError` if incompatible version of drawmembrane"""
    drawmembrane = name or cfg.get('executables', 'drawmembrane')
    try:
        p = subprocess.Popen([drawmembrane], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = p.communicate()
    except OSError:
        logger.critical("No drawmembrane binary found --- set it in the config file %(configfilename)s", configuration)
        raise
    m = re.search(r"\*\s*(?P<name>[^\s]+)\s+(?P<date>[/\d]+)", out)
    if m is None:
        raise EnvironmentError(errno.EIO, drawmembrane, "Cannot obtain version string from %r." % out)        
    if not m.group('name') == DRAWMEMBRANE_REQUIRED_NAME:
        raise EnvironmentError(errno.EIO, drawmembrane, "drawmembrane version %r does not work here, "
                               "need exactly %r (compile it from the src/drawmembrane directory)." % 
                               (m.group('name'), DRAWMEMBRANE_REQUIRED_NAME))
    month,day,year = map(int, m.group('date').split('/'))
    if not ((month,day,year) >= DRAWMEMBRANE_MINIMUM_VERSION):
        raise EnvironmentError(errno.EIO, drawmembrane, 
                               "version %d/%d/%d is too old, need at least %d/%d/%d" % 
                               ((month,day,year)+DRAWMEMBRANE_MINIMUM_VERSION))    
    return m.group('name'), (month,day,year) 

