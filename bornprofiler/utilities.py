# GromacsWrapper: utilities.py
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
# Copyright (c) 2009-2013 Oliver Beckstein <orbeckst@gmail.com>
# See the file COPYING for details.

"""
:mod:`bornprofiler.utilities` -- Helper functions and classes
==============================================================

The module defines some convenience functions and classes that are
used in other modules.


Classes
-------

.. autoclass:: AttributeDict


Functions
---------

Some additional convenience functions that deal with files and
directories:

.. function:: openany(directory[,mode='r'])

   Context manager to open a compressed (bzip2, gzip) or plain file
   (uses :func:`anyopen`).

.. autofunction:: anyopen
.. autofunction:: realpath
.. function:: in_dir(directory[,create=True])

   Context manager to execute a code block in a directory.

   * The *directory* is created if it does not exist (unless
     *create* = ``False`` is set)   
   * At the end or after an exception code always returns to
     the directory that was the current directory before entering
     the block.

.. autofunction:: find_first
.. autofunction:: withextsep

Functions that improve list processing and which do *not* treat
strings as lists:

.. autofunction:: iterable
.. autofunction:: asiterable


Functions that help handling files:

.. autofunction:: unlink_f

"""

from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os
import errno
from contextlib import contextmanager
import bz2, gzip

import logging
logger = logging.getLogger('bornprofiler.utilities')

class AttributeDict(dict):
    """A dictionary with pythonic access to keys as attributes --- useful for interactive work."""
    def __getattribute__(self,x):
        try:
            return super(AttributeDict,self).__getattribute__(x)
        except AttributeError:
            return self[x]
    def __setattr__(self,name,value):
        try:
            super(AttributeDict,self).__setitem__(name, value)
        except KeyError:
            super(AttributeDict,self).__setattr__(name, value)

    def __getstate__(self):
        return self

    def __setstate__(self, state):
        self.update(state)

@contextmanager
def openany(datasource, mode='r'):
    """Open the datasource and close it when the context exits."""
    stream, filename = anyopen(datasource, mode=mode)
    try:
        yield stream
    finally:
        stream.close()

def anyopen(datasource, mode='r'):
    """Open datasource (gzipped, bzipped, uncompressed) and return a stream.

    :Arguments:
    - *datasource*: a file or a stream
    - *mode*: 'r' or 'w'
    """
    handlers = {'bz2': bz2.BZ2File, 'gz': gzip.open, '': file}
    
    if mode.startswith('r'):
        if hasattr(datasource,'next') or hasattr(datasource,'readline'):
            stream = datasource
            filename = '(%s)' % stream.name  # maybe that does not always work?
        else:
            stream = None
            filename = datasource
            for ext in ('bz2', 'gz', ''):   # file == '' should be last
                openfunc = handlers[ext]
                stream = _get_stream(datasource, openfunc, mode=mode)
                if not stream is None:
                    break
            if stream is None:
                raise IOError("Cannot open %(filename)r in mode=%(mode)r." % vars())
    elif mode.startswith('w'):
        if hasattr(datasource, 'write'):
            stream = datasource
            filename = '(%s)' % stream.name  # maybe that does not always work?
        else:
            stream = None
            filename = datasource
            name, ext = os.path.splitext(filename)
            if ext.startswith('.'):
                ext = ext[1:]
            if not ext in ('bz2', 'gz'):
                ext = ''   # anything else but bz2 or gz is just a normal file
            openfunc = handlers[ext]
            stream = openfunc(datasource, mode=mode)
            if stream is None:
                raise IOError("Cannot open %(filename)r in mode=%(mode)r with type %(ext)r." % vars())
    else:
        raise NotImplementedError("Sorry, mode=%(mode)r is not implemented for %(datasource)r" % vars())

    return stream, filename

def _get_stream(filename, openfunction=file, mode='r'):
    try:
        stream = openfunction(filename, mode=mode)
    except IOError:
        return None

    try:
        stream.readline()
        stream.close()
        stream = openfunction(filename,'r')
    except IOError:
        stream.close()
        stream = None
    return stream

@contextmanager
def in_dir(directory, create=True):
    """Context manager to execute a code block in a directory.

    * The directory is created if it does not exist (unless
      create=False is set)
    * At the end or after an exception code always returns to
      the directory that was the current directory before entering
      the block.
    """
    startdir = os.getcwd()
    try:
        try:
            os.chdir(directory)
            logger.info("Working in %(directory)r..." % vars())
        except OSError, err:
            if create and err.errno == errno.ENOENT:
                os.makedirs(directory)
                os.chdir(directory)
                logger.info("Working in %(directory)r (newly created)..." % vars())
            else:
                logger.exception("Failed to start working in %(directory)r." % vars())
                raise
        yield os.getcwd()
    finally:
        os.chdir(startdir)

def realpath(*args):
    """Join all args and return the real path, rooted at /.

    Expands '~', '~user', and environment variables such as $HOME.

    Returns ``None`` if any of the args is ``None``.
    """
    if None in args:
        return None
    return os.path.realpath(os.path.expanduser(os.path.expandvars(os.path.join(*args))))

def find_first(filename, suffices=None):
    """Find first *filename* with a suffix from *suffices*.

    :Arguments:
      *filename*
         base filename; this file name is checked first
      *suffices*
         list of suffices that are tried in turn on the root of *filename*; can contain the 
         ext separator (:data:`os.path.extsep`) or not

    :Returns: The first match or ``None``.
    """
    
    root,extension = os.path.splitext(filename)
    if suffices is None:
        suffices = []
    else:
        suffices = withextsep(suffices)
    extensions = [extension] + suffices  # native name is first
    for ext in extensions:
        fn = root + ext
        if os.path.exists(fn):
            return fn
    return None

def withextsep(extensions):
    """Return list in which each element is guaranteed to start with :data:`os.path.extsep`."""
    def dottify(x):
        if x.startswith(os.path.extsep):
            return x
        return os.path.extsep + x
    return [dottify(x) for x in asiterable(extensions)]


def iterable(obj):
    """Returns ``True`` if *obj* can be iterated over and is *not* a  string."""
    if type(obj) is str:
        return False    # avoid iterating over characters of a string

    if hasattr(obj, 'next'):
        return True    # any iterator will do 
    try: 
        len(obj)       # anything else that might work
    except TypeError: 
        return False
    return True

def asiterable(obj):
    """Returns obj so that it can be iterated over; a string is *not* treated as iterable"""
    if not iterable(obj):
        obj = [obj]
    return obj


# In utilities so that it can be safely used in tools, cbook, ...

def unlink_f(path):
    """Unlink path but do not complain if file does not exist."""
    try:
        os.unlink(path)
    except OSError, err:
        if err.errno != errno.ENOENT:
            raise

def mkdir_p(path):
    """Create a directory *path* with subdirs but do not complain if it exists.

    This is like GNU ``mkdir -p path``.
    """
    try:
        os.makedirs(path)
    except OSError, err:
        if err.errno != errno.EEXIST:
            raise
