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
:Year: 2013
:Licence: GPL 3
:Copyright: (c) 2013 Lennard van der Feltz
"""
usage = """%prog pointx pointy pointz vecx vecy vecz length steplen

 Code for creating a straight line path of a given length starting at a point
 in the direction of a vec with points every steplen"""
import argparse
import numpy
parser = argparse.ArgumentParser()
parser.add_argument('point',type=float, nargs = 3)
parser.add_argument('vec',type=float, nargs = 3)
parser.add_argument('length',type=float)
parser.add_argument('steplen',type=float)
args = parser.parse_args()
vec = numpy.array(args.vec)
point = numpy.array(args.point)
length = args.length
steplen = args.steplen
def straightline(point, vec, length, steplen ):
    nsteps= length/steplen
    step = 0
    newpoint = point
    line = open("line.dat", 'w')
    line.write("{zero} {one} {two} \n".format(zero=point[0],one=point[1],two=point[2]))
    while step < nsteps:
        newpoint = newpoint+vec*steplen
        line.write("{zero} {one} {two} \n".format(zero=newpoint[0],one=newpoint[1],two=newpoint[2]))
        step += 1
    line.close()

straightline(point, vec, length, steplen )
