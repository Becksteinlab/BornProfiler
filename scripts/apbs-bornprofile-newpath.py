#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
# Copyright (c) 2005-2008 Kaihsu Tai, Oliver Beckstein
# Copyright (c) 2010-2013 Oliver Beckstein
usage ="""%prog [options] path-datafile

Quick hack in order to generate a nice path from a manually changed
one. The path is assumed to go roughly along the z direction. The
script re-orders datapoints by connection so that the output file is
the sequence of data points (the putative pathway), even if the path
reverses direction.

The input path must consist of xyz coordinates, one per line, or a PDB
file. The points with highest and lowest z-value are considered the end
points.

Using the --mindist option, one can prune data points. The pruned path
consists of those closest points that have the prescribed minimum
distance.

In the code we are building a graph of the path and in this way we can
look at connectivity instead of having to rely on, e.g. the
z-coordinates.
"""

# lazy import... this is a quick hack
from numpy import *

try:
     # our own distance_array is at least as fast as cdist()
     from MDAnalysis.analysis.distances import distance_array
except ImportError:
     # see  http://stackoverflow.com/questions/1871536/euclidean-distance-between-points-in-two-different-numpy-arrays-not-within
     # for recipes
     try:
          from scipy.spatial.distance import cdist
          def distance_array(xy1,xy2):
               return cdist(xy1, xy2, 'euclidean')
     except ImportError:
          # Alex Martelli's solution (still good but slower that above)
          import numpy
          def distance_array(xy1,xy2):
               d0 = numpy.subtract.outer(xy1[:,0], xy2[:,0])
               d1 = numpy.subtract.outer(xy1[:,1], xy2[:,1])
               return numpy.hypot(d0, d1)

try:
     import networkx as NX
except ImportError:
     raise ImportError("NetworkX http://networkx.lanl.gov/ is required.")

import bornprofiler.bpio

def modify_path(infile, outfile="newpath.dat", mindist=1.0):
     """try to get a inter-point distance of  mindist or more if needed"""

     p = bornprofiler.bpio.readPoints(infile)

     # inplace sort by z (col=2):
     # http://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
     p[:] = p[p[:,2].argsort()]

     d = distance_array(p.astype(float32),p.astype(float32))
     G = NX.Graph(d < 20)   # cutoff 20 is sufficient to connect points
     for  (n,nbrdict) in G.adjacency_iter():
          for k,eattr in nbrdict.items():
               eattr['weight'] = d[n,k]**3   # favours nearest neighbours (mostly...)

     dij3 = NX.dijkstra_path(G, 0, G.nodes()[-1])
     newp = p[dij3]
     savetxt(outfile, newp, "%8.3f %8.3f %8.3f")

     # newp is ordered in traversal order;
     # get the distances between adjacent nodes
     deltas = map(linalg.norm, newp[1:] - newp[:-1])

     # build new linear graph with distances between nodes as edge attribute
     G2 = NX.Graph()
     G2.add_path(arange(0,len(newp)), weight=0)
     for  ((n,nbrdict),delta) in zip(G2.adjacency_iter(), deltas):
          for k,eattr in nbrdict.items():
               eattr['weight'] = delta

     # remove points so that distance is the shortest distance larger than mindist
     pruned = NX.Graph()

     n = m = k = 0
     while (n < len(G2)-1):
          dist = 0
          while (k < len(G2) and dist < mindist):
               m = k  # move along the chain
               k += 1
               dist += G2.edge[m][k]['weight']
               #print "segment: ", (m,k,dist)
          pruned.add_edge(n,k,weight=dist)
          #print "edge:", (n,k,dist)
          n = k

     # sort nodes so that the output pdb is also in linear order
     # (nodes() returns node numbers in arbirary order but by construction
     # we know that the linear graph goes from 0 -> last)
     pruned_coords = newp[sort(pruned.nodes())]

     import os.path
     root,ext = os.path.splitext(outfile)
     new_outfile = (root+"_dq%.1f"+ext) % mindist
     new_pdb = (root+"_dq%.1f.pdb") % mindist

     savetxt(new_outfile, pruned_coords, "%8.3f %8.3f %8.3f")
     print "Wrote pruned path to %(new_outfile)r" % locals()

     # write pdb
     with open(new_pdb, "w") as pdb:
          for i,(x,y,z) in enumerate(pruned_coords):
               atomnr = i+1
               pdb.write("ATOM%(atomnr)7i  C   XXX X   1    %(x)8.3f%(y)8.3f%(z)8.3f\n" % locals())
     print "Wrote pruned path to %(new_pdb)r" % locals()

     return new_outfile, new_pdb


if __name__ == "__main__":
     from optparse import OptionParser

     parser = OptionParser(usage=__doc__)
     parser.add_option("--outfile", "-o", dest="outfile",
                       metavar="FILE",
                       help="write sample points to FILE and FILE.pdb [%default]")
     parser.add_option("--mindist", "-d", dest="mindist", type="float",
                       metavar="FLOAT",
                       help="minimum distance between two points [%default]")
     parser.set_defaults(outfile="newpath.dat", mindist=1.0)

     opts,args = parser.parse_args()

     try:
          infile = args[0]
     except IndexError:
          infile  = "/dev/stdin"

     modify_path(infile, opts.outfile, mindist=opts.mindist)

