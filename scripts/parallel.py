#!/usr/bin/env python
# $Id: parallel.py 833 2009-12-17 21:52:39Z oliver $
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""
:Author:  Oliver Beckstein
:Year: 2009
:License: GPL3
:Copyright: (c) 2009 Oliver Beckstein
:Copyright: (c) 2013 Oliver Beckstein
"""

usage ="""%prog NPROC command [args ...] --- [input ...]

Simple process-based parallelizer. Launches NPROC processes of
``command`` with the input list distributed equally over all
processes.

NOTE:

 - NO SANITY CHECKS!
 - It is not clear if there a restrictions on the length of command
   lines (in python's subprocess); if this is the case then not all
   input files are processed.
"""

import sys
import subprocess
import time
import numpy

WAIT_TIME = 1 # seconds to wait between polling processes

args = sys.argv[1:]
nproc = int(args.pop(0))

try:
    separator = args.index('---')
except ValueError:
    raise ValueError("EE Separate command with arguments from files with '---'.")

cmd = args[:separator]
files = args[separator+1:]

nfiles = len(files)

print "-- nproc:      %(nproc)r" % vars()
#print "-- files = %(files)r" % vars()
print "-- nfiles:     %(nfiles)d" % vars()
print "-- files/proc: %g" % (float(nfiles)/nproc, )

# one batch per process
batches = [[] for i in xrange(nproc)]    # need different [] !!
processes = numpy.array(nproc * [None])  # array of objects

# round robin distribution
i = -1
while len(files) > 0:
    i += 1
    i %= nproc
    x = files.pop(0)
    batches[i].append(x)

print "-- batch sizes: %r" % map(len, batches)

for i,batch in enumerate(batches):
    _cmd = cmd + batch
    #print "[%(i)2d] Launching process %(_cmd)r" % vars()
    print "-- [%2d] Launching process %r with %d inputs" % (i,cmd,len(batch))
    processes[i] = subprocess.Popen(_cmd)

def proc_running():
    """Returns True if at least one process is still running."""
    return numpy.any([(p.poll() is None) for p in processes])

def returncodes():
    return numpy.array([p.poll() for p in processes])

# poll status
while proc_running():
    time.sleep(WAIT_TIME)

rc = returncodes()

print "II all processes done"
print "-- return codes: %r" % rc

if numpy.all(rc == 0):
    exit_code = 0
else:
    print "WW Warning: some processes finished with non-zero exit codes."
    exit_code = 1

sys.exit(exit_code)
