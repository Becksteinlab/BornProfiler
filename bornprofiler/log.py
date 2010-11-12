# log.py
# logging for bornprofiler module
# Copyright (c) 2010 Oliver Beckstein
# Published under the GNU Public Licence v3

"""
Setting up logging --- :mod:`bornprofiler.log`
====================================================

Configure logging for the BornProfiler. Import this module if logging is
desired in application code.

Logging to a file and the console.

See http://docs.python.org/library/logging.html?#logging-to-multiple-destinations

The top level logger of the library is named 'bornprofiler'.  Note that we
are configuring this logger with console output. If the root logger
also does this then we will get two output lines to the console. We'll
live with this because this is a simple convenience library and most
people will not bother with a logger (I think...)

In modules that use loggers get a logger like so::
     import logging
     logger = logging.getLogger('bornprofiler.MODULENAME')
"""

import logging

def create(logger_name="bornprofiler", logfile='bornprofiler.log'):
    """Create a top level logger.

    - The file logger logs everything (including DEBUG).
    - The console logger only logs INFO and above.

    Logging to a file and the console.
    
    See http://docs.python.org/library/logging.html?#logging-to-multiple-destinations
    
    The top level logger of the library is named 'bornprofiler'.  Note that
    we are configuring this logger with console output. If the root
    logger also does this then we will get two output lines to the
    console. We'll live with this because this is a simple
    convenience library...
    """

    logger = logging.getLogger(logger_name)

    logger.setLevel(logging.DEBUG)

    logfile = logging.FileHandler(logfile)
    logfile_formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logfile.setFormatter(logfile_formatter)
    logger.addHandler(logfile)

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console.setFormatter(formatter)

    logger.addHandler(console)

    return logger

def clear_handlers(logger):
    """clean out handlers in the library top level logger 

    (only important for reload/debug cycles...)
    """
    for h in logger.handlers:
        logger.removeHandler(h)


class NullHandler(logging.Handler):
    """Silent Handler.

    Useful as a default::
      h = NullHandler()
      logging.getLogger("bornprofiler").addHandler(h)
      del h
    """
    def emit(self, record):
        pass


