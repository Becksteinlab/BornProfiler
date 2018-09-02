#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""Code for creating a straight line path pdb of a given length
starting at a point in the direction of a vec with points every
steplen.
"""

from __future__ import with_statement, division

import numpy

def straightline(point, vec, length, steplen, Name):
    point = numpy.asarray(point)
    vec = numpy.asarray(vec)
    nsteps = length//steplen
    step = 0
    newpoint = point
    with open("{name}.pdb".format(name=Name), 'w') as line:
        line.write("ATOM {atomnum:>6}   X  XXX X{atomnum:>4}{x:>12,.2f}{y:>8,.2f}{z:>8,.2f} \n".format(
            atomnum=step + 1, x=point[0], y=point[1], z=point[2]))
        coordinate_string = ""
        while step < nsteps:
            newpoint += vec*steplen
            coordinate_string += "ATOM {atomnum:>6}   X  XXX X{atomnum:>4}{x:>12,.2f}{y:>8,.2f}{z:>8,.2f} \n".format(
                atomnum=step + 2, x=newpoint[0], y=newpoint[1], z=newpoint[2])
            step += 1
        line.write(coordinate_string)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('point', type=float, nargs=3,
                        help="X Y Z coordinates of point on line")
    parser.add_argument('vec', type=float, nargs=3,
                        help="X Y Z coordinates of direction")
    parser.add_argument('length', type=float, help="length of "
                        "line. NOTE: N = length//steplen segments "
                        "will be generated, hence the final length can "
                        "differ from the requested length.")
    parser.add_argument('steplen', type=float,
                        help="length of a line segment")
    parser.add_argument('--title', default="Line")

    args = parser.parse_args()

    straightline(args.point, args.vec, args.length, args.steplen, args.title)

