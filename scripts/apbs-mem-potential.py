#!/usr/bin/env python
# (c) 2010 Oliver Beckstein
# Licensed under GPL
"""%%prog [options] parameter-file [pqr]

Runs customized apbs calculation of a protein in a membrane. Because
the number of options is pretty large, everything must be specified in
a parameter file. The one exception is the pqr file: if provided as a
second argument, it override the setting of ``environment: pqr`` in
the configuration file.

A new parameter-file can be generated with the --template
option; in this case only the file is written and no further actions are
performed.

We use :program:`draw_membrane2a` to add a low-dielectric region
representing the membrane and optional a medium-dielectric region for
the headgroups. It is also possible to specify a channel exclusion
zone for a pore through a channel.

Paths to draw_membrane2 and apbs are set in the configuration file
%(configfilename)r or in the parameter file::

 - apbs = %(apbs)r
 - draw_membrane2 = %(drawmembrane)r

.. SeeAlso: apbsmem_ and the `APBS PMF tutorial`_.

.. _APBS PMF tutorial:
   http://www.poissonboltzmann.org/apbs/examples/potentials-of-mean-force/the-polar-solvation-potential-of-mean-force-for-a-helix-in-a-dielectric-slab-membrane
.. _apbsmem: http://mgrabe1.bio.pitt.edu/apbsmem/

"""

import os.path
import bornprofiler, bornprofiler.membrane
import logging
logger = logging.getLogger('bornprofiler')

if __name__ == "__main__":
    import sys
    import errno
    from optparse import OptionParser

    bornprofiler.start_logging()

    parser = OptionParser(usage=__doc__ % bornprofiler.config.configuration)
    parser.add_option("--template", dest="write_template", action="store_true",
                      help="write template parameter file and exit")
    parser.add_option("-s", "--suffix", dest="suffix",
                      help="suffix for all generated files [%default]")
    parser.add_option("--no-run", dest="run", action="store_false",
                      help="do not immediately run apbs for the Poisson-Boltzmann calculation "
                      "and only produce all required input files. apbs is still run in order "
                      "to obtain the dielectric, charge, and kappa maps needed for draw_membrane.")
    parser.set_defaults(suffix="S", run=True)
    
    opts,args = parser.parse_args()
    
    try:
        filename = args[0]
    except:
        logger.fatal("Provide the parameter filename. See --help.")
        sys.exit(1)
        
    if opts.write_template:
        bornprofiler.write_parameters(filename)
        sys.exit(0)
        
    params = bornprofiler.io.RunParameters(args[0])
    kw = params.get_apbsmem_kwargs()

    try:
        pqr = args[1]
        del kw['pqr']
    except IndexError:
        pqr = kw.pop('pqr')    
    
    A = bornprofiler.membrane.APBSmem(pqr, opts.suffix, **kw)
    A.generate()
    if opts.run:
        if A.unpack_dxgz:
            A.ungzip_dx()
        A.run_apbs('solvation')
        if A.unpack_dxgz:
            A.gzip_dx()

    bornprofiler.stop_logging()
