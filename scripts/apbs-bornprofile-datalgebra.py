#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""
:Author: Lennard van der Feltz
:Year: 2014
:Licence: GPL 3
:Copyright: (c) 2013 Lennard van der Feltz
"""

import numpy 
import scipy.interpolate
import argparse
import logging
import bornprofiler.core

logger=logging.getLogger("bornprofiler")
parser = argparse.ArgumentParser()
parser.add_argument('--zdif',nargs=2,help="take the difference of the first dat file energy values from the second dat file energy values by interpolating")
parser.add_argument('--spacing',default=1,help="spacing between interpolation points.")
args = parser.parse_args()
zdif = args.zdif
spacing = args.spacing

def writeout(outName, datinfo):
    with open(outName, "w") as outFile:
        outFile.write("# x/A y/A z/A  W/(kJ/mol)\n")
        for x,y,z,E in datinfo:
            outFile.write("{x:>8,.3f} {y:>8,.3f} {z:>8,.3f} {E:>8,.3f}\n".format(x=x,y=y,z=z,E=E))
    logger.info("Wrote Born PMF to {outName}.".format(outName=outName))

def zdifference(dats,spacing):
    """Writes a dat file containing the difference between the energy of the two
    input dat files. This difference is done between the interpolated energy
    values at regular intervals with {spacing} distance between sample points.
    Does not assume the same path between dat files and output dat simply fills
    x and y values with zeroes. """
    datob1 = bornprofiler.core.datinfo(dats[0])
    datob2 = bornprofiler.core.datinfo(dats[1])
#Obtain lowest and highest z across dat files
    minz = min([datob1.minz,datob2.minz])
    maxz = max([datob2.maxz,datob1.maxz])
#Produce evenly spaced points between the min and max
    interarray = numpy.arange(minz,maxz,spacing)
    numpoints = interarray.shape[0]
    enerdif = datob1.z_interpolator()(interarray) - datob2.z_interpolator()(interarray)
#Convert energy difference array from shape (N) to shape (N,1)
    enerdif = numpy.reshape(enerdif,(numpoints,1))
#Convert point array from shape (N) to (N,3) with 0s to fill empty columns before the z values from interarray.
    interarray = numpy.reshape(interarray,(numpoints,1))
    interpoints = numpy.hstack((numpy.zeros((numpoints,2)),interarray))
#Stack the spacial and energy arrays together
    newdat = numpy.hstack((interpoints,enerdif))
    writeout("{dat0}_minus_{dat1}".format(dat0=dats[0].split('/')[-1].split('.')[0],dat1=dats[1].split('/')[-1]),newdat)

bornprofiler.start_logging()
zdifference(zdif,spacing)
bornprofiler.stop_logging()
