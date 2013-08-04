# APBS BornProfiler python package
# Copyright (c) 2005-2008 Kaihsu Tai, Oliver Beckstein
# Copyright (c) 2010-2011 Oliver Beckstein
# Published under the GNU Public Licence, version 3

"""
APBS BornProfiler python package
================================

A collection of scripts to set up Poisson-Boltzmann electrostatic calculations
for the APBS_ package, in particular calculations of the electrostatic
solvation free energy of an ion along a pathway in a membrane protein (the
so-called *Born profile*).

.. _APBS: http://www.poissonboltzmann.org/apbs

"""
from version import get_version, get_version_tuple

import config, utilities

import logging
# see the advice on logging and libraries in
# http://docs.python.org/library/logging.html?#configuring-logging-for-a-library
class NullHandler(logging.Handler):
    def emit(self, record):
        pass
h = NullHandler()
logging.getLogger("bornprofiler").addHandler(h)
del h

def start_logging(logfile="bornprofiler.log"):
    """Start logging of messages to file and console."""
    import log
    log.create("bornprofiler", logfile=logfile)
    logging.getLogger("bornprofiler").info("BornProfiler STARTED logging to %r", logfile)

def stop_logging():
    """Stop logging to logfile."""
    import log
    logger = logging.getLogger("bornprofiler")
    logger.info("BornProfiler STOPPED logging")
    log.clear_handlers(logger)  # this _should_ do the job...

import core, io

def write_parameters(filename, **defaults):
    """Write a default parameter file to *filename*.

    key-value pairs in *kwargs* are written to the DEFAULT section of
    the parameter file. Edit the file with a text editor and then use
    it as input for :class:`core.MPlaceion`.

    .. Note:: If the file already exists then it will be updated but
              not reset to the default values of all sections.
    """
    p = io.RunParameters(filename, **defaults)
    return p.write()
