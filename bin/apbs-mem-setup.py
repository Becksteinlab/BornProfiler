#!/usr/bin/env python
# (c) 2010 Oliver Beckstein
# Licensed under GPL
"""%%prog PQR-file

Runs customized apbs calculation of PQR file in a membrane. Currently
customized for Yohei's CFTR structures (from the 30 ns MD); see the
source for details.

Uses draw_membrane2 to add a low-dielectric (eps=2) region and sets
protein dielectric to eps=10.

.. Note:: Paths to draw_membrane2 and apbs are hard coded!
          apbs = %(APBS)r
          draw_membrane2 = %(DRAWMEMBRANE)r

Commandline version of apbsmem, following
http://www.poissonboltzmann.org/apbs/examples/potentials-of-mean-force/the-polar-solvation-potential-of-mean-force-for-a-helix-in-a-dielectric-slab-membrane
"""

import os.path
BASEDIR = os.path.expanduser("~/Projects/Channels/CFTR")

from membrane import APBSmem, APBS, DRAWMEMBRANE

class CFTRmem(APBSmem):
    """APBSmem with custom defaults"""
    def __init__(self, pqr, suffix, zmem=-2, lmem=40, mdie=2, sdie=80, pdie=10, 
                 headgroup_die=20, headgroup_l=0, Vmem=0,
                 Rtop=16, Rbot=10,  **kwargs):
        super(CFTRmem, self).__init__(pqr, suffix, zmem=-2, lmem=40, mdie=2, sdie=80, pdie=10, 
                 headgroup_die=20, headgroup_l=0, Vmem=0,
                 Rtop=16, Rbot=10,  **kwargs)

if __name__ == "__main__":
    import sys
    import errno
    from optparse import OptionParser

    parser = OptionParser(usage=__doc__ % vars())
    parser.add_option("-s", "--suffix", dest="suffix",
                      help="suffix for all generated files [%default]")
    parser.set_defaults(suffix="S")
    
    options, args = parser.parse_args()
    
    try:
        pqr = args[0]
    except IndexError:
        raise ValueError("Need PQR file as input")

    if not os.path.exists(pqr):
        raise IOError(errno.ENOENT, "PQR file not found", pqr)
        
    suffix = 'S'
    
    C = CFTRmem(pqr, options.suffix)
    C.setup()
    C.run_apbs('solvation')
