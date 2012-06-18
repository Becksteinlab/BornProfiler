#!/usr/bin/env python
# quick and dirty way to get the samplepoints.dat from HOLE .sph 
# (sphere) output

# $Id: mksample.py,v 1.10 2008/03/13 17:21:51 kaihsu Exp $

"""%prog [options] [FILE]

Create samplepoints for the Born profile scripts from a pdb file or a
HOLE sphere file, which includes points on the axis of a channel
pore. The coordinates of the pore axis are output as a simple list of
one point per line (x y z). 

By default, a sph file is read from stdin and output to stdout but
this can be changed by commandline options.

A sphere file is simply a pdb file and we only parse ATOM or HETATM
records. A sph file also contains the pore radius in the occupancy
filed (to 2 decimals), which can be written to the radius profile
R(z). If more accurate values are required one needs to grep HOLE's
output.
"""

from __future__ import with_statement

import string
import sys
import os

import logging
logger = logging.getLogger('bornprofiler') 

if __name__ == "__main__":
    from optparse import OptionParser
    
    logging.basicConfig()

    parser = OptionParser(usage=__doc__)
    parser.add_option("--outfile", "-o", dest="outfile",
                      metavar="FILE",
                      help="write sample points to FILE [%default]")
    parser.add_option("--pdb", "-s", dest="pdb",
                      metavar="FILE",
                      help="write path to pdb-formatted FILE [%default]")
    parser.add_option("--profile", "-p", dest="radiusfile",
                      metavar="FILE",
                      help="write radius profile R(z) to FILE (only HOLE sph "
                      "files; for any other we take what is in the occupancy field or write -1")
    parser.add_option("--skip", "-X", dest="skip", type="int",
                      metavar="N",
                      help="skip over N steps for the sample points [%default]")
    parser.add_option("--every", dest="every", type="int",
                      metavar="N",
                      help="use every Nth data point (same as --skip=N-1)")
    parser.set_defaults(outfile="/dev/stdout", pdb="path.pdb", skip=0)

    opts,args = parser.parse_args()

    try:
        infile = args[0]
    except IndexError:
        infile  = "/dev/stdin"

    if opts.every:
        skipN = opts.every - 1
    else:
        skipN = opts.skip

    try: 
        sphere = open(infile, "r")
    except IOError, details:
        sys.stderr.write("Error opening sph file '%s'. %s\n" % (infile,details))
        sys.exit(2)        
    logger.info("Opened sphere file '%s'." % infile)

    linecount = 0
    spheredata = []
    for line in sphere:
        record = line[0:6].strip()
        if not record in ('ATOM','HETATM'):
            continue
        if (string.find(line, "S 888") >= 0) or \
           (string.find(line, "S-888") >= 0) or \
           (string.find(line, "LAST-REC-END") >= 0):
            continue   # special HOLE sphere file records
        # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM ::
        # 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
        # HETATM   28  O   HOH X   0      -4.801   1.393 -31.205  1.00  1.00           O  
        x,y,z = map(float, [line[30:38], line[38:46], line[46:54]])
        occupancy = line[54:60]
        try:
            r = float(occupancy)
        except ValueError:   # file that does not use the occupancy field for radius?
            r = -1
        spheredata.append((x,y,z,r))
        linecount += 1
    sphere.close()
    logger.debug("total lines read:    %d", linecount)

    # sort on z
    # sort(key=...) requires python 2.4; otherwise rewrite with
    # sort(cmp=lambda x,y:int(x[2] - y[2])) or decorate-sort-undecorate
    spheredata.sort(key=lambda x:x[2])

    # write
    with open(opts.outfile, "w") as sample:
        logger.info("Opened output file '%s'.", opts.outfile)

        if opts.radiusfile:
            try: 
                radius = open(opts.radiusfile, "w")
            except IOError, details:
                sys.stderr.write("Error opening output radius profile file '%s'. %s\n"
                                 % (opts.radiusfile,details))
                sys.exit(2)
            radius.write('# HOLE radius along z-axis\n')
            radius.write('# sphere file: %s\n' % infile)
            logger.info("Opened radius file '%s'.", opts.radiusfile)

        if opts.pdb:
            try:
                pdb = open(opts.pdb, "w")
            except IOError, details:
                sys.stderr.write("Error opening output PDB file '%s'. %s\n"
                                 % (opts.pdb,details))
                sys.exit(2)

        skipN += 1
        printcount = linecount = 0
        for (x,y,z,r) in spheredata:
            if (0 == linecount % skipN):              
                sample.write("%(x)8.3f %(y)8.3f %(z)8.3f\n" % locals())
                printcount += 1
                if opts.pdb:
                    pdb.write("ATOM%(printcount)7i  C   XXX X   1    %(x)8.3f%(y)8.3f%(z)8.3f\n" % locals())
            if opts.radiusfile:
                radius.write("%(z)8.3f %(r)6.2f\n" % locals())
            linecount += 1

    logger.info("total lines written: %d", printcount)

    if opts.radiusfile:
        radius.close()
    if opts.pdb:
        pdb.close()
