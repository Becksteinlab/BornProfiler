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

usage = """%prog pointx pointy pointz vecx vecy vecz length steplen --title

 Code for creating a straight line path pdb of a given length starting
 at a point in the direction of a vec with points every steplen.

"""

import argparse
import numpy

parser = argparse.ArgumentParser()
parser.add_argument('point',type=float, nargs = 3)
parser.add_argument('vec',type=float, nargs = 3)
parser.add_argument('length',type=float)
parser.add_argument('steplen',type=float)
parser.add_argument('--title', default="Line")

args = parser.parse_args()
vec = numpy.array(args.vec)
point = numpy.array(args.point)
length = args.length
steplen = args.steplen
title = args.title

def straightline(point, vec, length, steplen,Name ):
    nsteps= length/steplen
    step = 0
    newpoint = point
    line = open("{name}.pdb".format(name=Name), 'w')
    line.write("ATOM {atomnum:>6}   X  XXX X{atomnum:>4}{x:>12,.2f}{y:>8,.2f}{z:>8,.2f} \n".format(atomnum = step + 1,x=point[0],y=point[1],z=point[2]))
    coordinate_string = ""
    while step < nsteps:
        newpoint += vec*steplen
        coordinate_string +="ATOM {atomnum:>6}   X  XXX X{atomnum:>4}{x:>12,.2f}{y:>8,.2f}{z:>8,.2f} \n".format(atomnum = step + 2,x=newpoint[0],y=newpoint[1],z=newpoint[2])
        step += 1
    line.write(coordinate_string)
    line.close()

straightline(point,vec,length,steplen,title)

