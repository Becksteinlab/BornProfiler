# APBS BornProfiler python package
# Copyright (c) 2005-2008 Kaihsu Tai, Oliver Beckstein
# Copyright (c) 2010-2011 Oliver Beckstein
# Published under the GNU Public Licence, version 3

__all__ = ['VERSION', 'get_version', 'get_version_tuple']

#: Package version; this is the only place where it is set.
VERSION = 0,9,1

def get_version():
    """Return current package version as a string."""
    return ".".join(map(str,VERSION))

def get_version_tuple():
    """Return current package version as a tuple (*MAJOR*, *MINOR*, *PATCHLEVEL*)."""
    return tuple(VERSION)

