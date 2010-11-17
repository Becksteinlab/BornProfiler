#!/usr/bin/env python
# quick and dirty way to get the samplepoints.dat from HOLE .sph 
# (sphere) output

# $Id: mksample.py,v 1.10 2008/03/13 17:21:51 kaihsu Exp $

"""
    Create samplepoints for the Born profile scripts from a
    HOLE sphere file, which includes points on the axis of
    a channel pore.
"""

import string
import getopt
import sys
import os

# global variables
verbose = 0

def usage(rc):
    str = """mksample

This script takes a .sph file from HOLE as input and outputs the coordinates
of the pore axis as a simple list of one point per line (x y z). By default,
a sph file is read from stdin and output to stdout but this can be changed by
commandline options.

Usage: mksample.py [options]
-h,--help             this help
-v,--verbose          verbose messages
--version             show the version
--debug=N             set verbosity level (0 default, verbose == 3)
-s,--sphere=FILE      sph file from HOLE (instead of stdin)
-o,--output=FILE      x y z list is written to FILE instead of stdout.
-p,--profile=FILE     write radius profile R(z) to FILE
-X,--skip=N           skip over N steps for the sample points 
                      (not the pore profile!); defaults to 0
--every=M             make sample point every Mth step; same as --skip=(M-1)

Note that the radii in the pore profile R(z) are only recorded with
two decimals in the sphere file. If more accurate values are required
one needs to grep HOLE's output.
"""
    sys.stderr.write(str)
    sys.exit(rc)

def msg(level,m):
    """
       Print a message if the verbosity level  >= the level of the message.
       The verbosity level 'verbose' must be declared globally
    """
    if (verbose >= level):
        if (verbose > 3):
            """If in debug mode then show debug level with message"""
            m = "[dbglvl%3d] %s" % (level,m)
        sys.stderr.write("%s\n" % m)
    

def mainCommand():
    """
        Main driver for running program from the command line.
    """
    global verbose
    
    shortOptlist = "hvs:o:p:X:D:"
    longOptlist = ["help","verbose","sphere=","output=","profile=","skip=","every=","debug=","version"]

    try: opts, args = getopt.getopt(sys.argv[1:], shortOptlist, longOptlist)
    except getopt.GetoptError, details:
        sys.stderr.write("GetoptError:  %s\n" % details)
        usage(2)

    infile  = "/dev/stdin"
    outfile = "/dev/stdout"
    radiusfile = None
    skipN = 0
    for o,a in opts:
        if o in ("-v","--verbose"):
            verbose = 3
        elif o in ("-D","--debug"):
            try: verbose = int(a)
            except ValueError:
                sys.stderr.write("Option error: Debug level must be a number\n")
                sys.exit(2)
        elif o in ("--version"):
            print "$Id: mksample.py,v 1.10 2008/03/13 17:21:51 kaihsu Exp $"
            sys.exit()
        elif o in ("-h","--help"):
            usage(0)
            sys.exit()
        elif o in ("-s","--sphere"):
            infile = a 
        elif o in ("-o","--output"):
            outfile = a
        elif o in ("-p","--profile"):
            radiusfile = a
        elif o in ("-X","--skip"):
            skipN = int(a)
        elif o in ("--every"):
            skipN = int(a) - 1

    try: sphere = open(infile, "r")
    except IOError, details:
        sys.stderr.write("Error opening sph file '%s'. %s\n" % (infile,details))
        sys.exit(2)        
    msg(3,"Opened sphere file '%s'." % infile)


    linecount = 0
    spheredata = []
    for line in sphere:
        if (string.find(line, "S 888") == -1) and \
           (string.find(line, "S-888") == -1) and \
           (string.find(line, "LAST-REC-END") == -1):
            spheredata.append(map(float, line[30:60].split()))  # split: not pdb standard:
                                  # I hope for spaces ... or rewrite with fields
            linecount += 1
    sphere.close()
    msg(3,"total lines read:    %d" % linecount)

    # sort on z
    # sort(key=...) requires python 2.4; otherwise rewrite with
    # sort(cmp=lambda x,y:int(x[2] - y[2])) or decorate-sort-undecorate
    spheredata.sort(key=lambda x:x[2])

    # write
    try: sample = open(outfile, "w")
    except IOError, details:
        sys.stderr.write("Error opening output sample points file '%s'. %s\n" % (outfile,details))
        sys.exit(2)        
    msg(3,"Opened output file '%s'." % outfile)

    if radiusfile:
        try: radius = open(radiusfile, "w")
        except IOError, details:
            sys.stderr.write("Error opening output radius profile file '%s'. %s\n"
                             % (radiusfile,details))
            sys.exit(2)
        radius.write('# HOLE radius along z-axis\n')
        radius.write('# sphere file: %s\n' % infile)
        msg(3,"Opened radius file '%s'." % radiusfile)

    skipN += 1
    printcount = linecount = 0
    for (x,y,z,r) in spheredata:
        if (0 == linecount % skipN):              
            sample.write("%(x)f %(y)f %(z)f\n" % locals())
            printcount += 1            
        if radiusfile:
            radius.write("%(z)f %(r).2f\n" % locals())
        linecount += 1

    msg(3,"total lines written: %d" % printcount)
    sample.close()

    if radiusfile:
        radius.close()
    sys.exit(0)


if __name__ == "__main__":
    mainCommand()
  

# $Log: mksample.py,v $
# Revision 1.10  2008/03/13 17:21:51  kaihsu
# refine
#
# Revision 1.9  2008/03/13 17:20:40  kaihsu
# oops; correct arithmetic
#
# Revision 1.8  2008/03/13 17:17:55  kaihsu
# refine usage message
#
# Revision 1.7  2008/03/13 17:16:58  kaihsu
# add --every option
#
# Revision 1.6  2008/03/13 17:13:30  kaihsu
# add --version option
#
# Revision 1.5  2007/11/02 22:02:55  oliver
# Added Sarah's feature request:
# * --profile option writes the HOLE profile R(z) in addition to the sample
#  points. The precision is limited to two decimals (B-factor filed in the
#  sph file)
# * reorganized code a bit: now we can properly sort the samplepoints on z
#   (otherwise the profiles would have nasty jumps, and it should also
#   provide for a ordered samplepoints list)
# * some cleanup (removed useless stuff)
#
# Revision 1.4  2004/12/16 11:53:07  oliver
# - --skip option to reduce number of sample points
# - debug level/messages works now (needed global verbose statement)
# - error checks on file opening
# - fixed option processing (after RTFM...)
#
# Revision 1.3  2004/12/13 13:46:01  oliver
# - help function
# - can read files from stdin or command line
# - write files to stdout or file given on command line
# - verbose is useless at the moment and so is the CGI stub
# (used pdb2pqr.py as a template)
#
# Revision 1.2  2004/07/27 18:14:51  kaihsu
# seems to work for nAChR
#
# Revision 1.1  2004/07/26 13:12:50  kaihsu
# quick-and-dirty samplepoints.dat maker
#
