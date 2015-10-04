# APBS BornProfiler python package
# Copyright (c) 2005-2008 Kaihsu Tai, Oliver Beckstein
# Copyright (c) 2010-2011 Oliver Beckstein
# Published under the GNU Public Licence, version 3

"""
Version information --- :mod:`bornprofiler.version`
===================================================

This module records the version of the software and also makes
functions available to access the version.
"""


__all__ = ['VERSION', 'get_version', 'get_version_tuple']

#: Is this a release? If not, the :data:`VERSION` will have
#: "-dev" appended.
RELEASE = False
#: Package version; this is the only place where it is set.
VERSION = 0,10,0

if not RELEASE:
    VERSION = VERSION + ('dev',)

def get_version():
    """Return current package version as a string."""
    s = ".".join(map(str, VERSION[:3]))
    if not RELEASE:
        s = s + "-" + VERSION[3]
    return s

def get_version_tuple():
    """Return current package version as a tuple (*MAJOR*, *MINOR*, *PATCHLEVEL*)."""
    return tuple(VERSION)

