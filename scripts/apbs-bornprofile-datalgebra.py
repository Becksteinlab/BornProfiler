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
import bornprofiler.bpio as bpio
import bornprofiler.array_alg as array_alg

logger=logging.getLogger("bornprofiler")
parser = argparse.ArgumentParser()
parser.add_argument('--zdif',nargs=2,help="take the difference of the first dat file energy values from the second dat file energy values by interpolating along z")
parser.add_argument('--spacing',default=1.0,type = float,help="spacing between interpolation points.")
parser.add_argument('--zavg',nargs = '+',default = None, help ="Average the energy values of all dat files by interpolating along the z axis")
parser.add_argument('--title',help = "Name to give output file")
parser.add_argument('--zdiv',nargs = 2, help="take the quotient of the first dat file energy values divided by the second dat file energy values by interpolating along z")
args = parser.parse_args()
zdif = args.zdif
spacing = args.spacing
zavg = args.zavg
title = args.title
zdiv = args.zdiv

def zdifference(dats,spacing,title):
    """Writes a dat file containing the difference between the energy of the two
    input dat files. This difference is done between the interpolated energy
    values at regular intervals with {spacing} distance between sample points.
    Does not assume the same path between dat files and output dat simply fills
    x and y values with zeroes. """
    datarray1 = bpio.read_dat_to_array(dats[0])
    datarray2 = bpio.read_dat_to_array(dats[1])
    datob1 = array_alg.curve(datarray1,spacing=spacing)
    datob2 = array_alg.curve(datarray2,spacing=spacing)
    newdat = datob1 - datob2
    if title == None:
        title = "{dat0}_minus_{dat1}".format(dat0=dats[0].split('/')[-1],dat1=dats[1].split('/')[-1])
    else:
        title = title + ".dat"
    bpio.write_dat_from_array(title,newdat.datarray)

def zaverage(dats,spacing,title):
    """Writes a dat file containing the averaged energy value across all
    dat files. Interpolates energy values for every dat over the z values where
    all dat files exist. Uses {spacing} as distance between points for
    interpolation. Fills x and y of output file with zeroes."""
    datarrays = [bpio.read_dat_to_array(x) for x in dats]
    datobs = numpy.array([array_alg.curve(x,spacing=spacing) for x in datarrays])
    average = numpy.sum(datobs)/datobs.shape[0]
    if title == None:
       title = "Average.dat"
    else:
       title = title + ".dat"
    bpio.write_dat_from_array(title,average.datarray)

def zdivision(dats,spacing,title):
    """Writes a dat file containing the quotient of the energy of the two
    input dat files. This division is done between the interpolated energy
    values at regular intervals with {spacing} distance between sample points.
    Does not assume the same path between dat files and output dat simply fills
    x and y values with zeroes. """
    datarray1 = bpio.read_dat_to_array(dats[0])
    datarray2 = bpio.read_dat_to_array(dats[1])
    datob1 = array_alg.curve(datarray1,spacing=spacing)
    datob2 = array_alg.curve(datarray2,spacing=spacing)
    newdat = datob1/datob2
    if title == None:
        title = "{dat0}_divided_by_{dat1}".format(dat0=dats[0].split('/')[-1],dat1=dats[1].split('/')[-1])
    else:
        title = title + ".dat"
    bpio.write_dat_from_array(title,newdat.datarray)

bornprofiler.start_logging()
if zdif != None:
    zdifference(zdif,spacing,title)
elif zavg != None:
    zaverage(zavg,spacing,title)
elif zdiv != None:
    zdivision(zdiv,spacing,title)
bornprofiler.stop_logging()
