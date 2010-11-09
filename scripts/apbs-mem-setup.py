#!/usr/bin/env python
# (c) 2010 Oliver Beckstein
# Licensed under GPL
"""%%prog PQR-file

Runs customized apbs calculation of PQR file in a membrane.

A custom APBSmem class must be available in a python file ``custom.py``. 

Currently only the *CFTRmem* class, customized for Yohei's CFTR
structures (from the 30 ns MD) is available; see the source for
details.

Uses draw_membrane2 to add a low-dielectric (eps=2) region and sets
protein dielectric to eps=10.

Paths to draw_membrane2 and apbs are set in the configuration file
%(configfilename)r:

 - apbs = %(apbs)r
 - draw_membrane2 = %(drawmembrane)r

Commandline version of apbsmem, following
http://www.poissonboltzmann.org/apbs/examples/potentials-of-mean-force/the-polar-solvation-potential-of-mean-force-for-a-helix-in-a-dielectric-slab-membrane
"""

import os.path
from bornprofiler.config import cfg, configuration
import logging
logger = logging.getLogger('bornprofiler')

if __name__ == "__main__":
    import sys
    import errno
    from optparse import OptionParser

    logging.basicConfig()

    parser = OptionParser(usage=__doc__ % configuration)
    parser.add_option("-C","--class-name", dest="clsname",
                      help="Python name of a class derived from "
                      "bornprofiler.membrane.APBSmem [%default]")
    parser.add_option("-s", "--suffix", dest="suffix",
                      help="suffix for all generated files [%default]")    
    parser.set_defaults(clsname=cfg.get('membrane','class'), 
                        suffix="S")
    
    options, args = parser.parse_args()
    
    try:
        pqr = args[0]
    except IndexError:
        raise ValueError("Need PQR file as input")

    if not os.path.exists(pqr):
        raise IOError(errno.ENOENT, "PQR file not found", pqr)
        
    suffix = 'S'

    cls = None
    for modname in 'custom', 'bornprofiler.custom':
        try:
            mod = __import__(modname, fromlist=[opts.clsname])
            cls = mod.__getattribute__(opts.clsname)
        except (ImportError, AttributeError):
            pass
    try:
        C = cls(pqr, options.suffix)
    except TypeError:
        ermsg = "No setup class %r found in custom.py or bornprofiler.custom" % opts.clsname
        logger.fatal(errmsg)
        raise ValueError(errmsg)

    C.setup()
    C.run_apbs('solvation')
