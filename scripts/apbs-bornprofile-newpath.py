#!/usr/bin/env python
"""%prog [options] path-datafile

Quick hack in order to generate a nice path from a manually changed
one. The input path must consist of xyz coordinates, one per line.

At the moment, it is assumed that the points are supposed to describe
a path roughly parallel to the z-axis --- this is used to find the
first and last point of the path.

"""

import MDAnalysis.analysis.distances
import networkx as NX
from numpy import *

def modify_path(infile, outfile="newpath.dat", mindist=1.0):
     #mindist = 1.0   # try to get a inter-point distance of  mindist or more if needed

     p = loadtxt(infile)
     # inplace sort by z (col=2): 
     # http://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
     p[:] = p[p[:,2].argsort()]

     d = MDAnalysis.analysis.distances.distance_array(p.astype(float32),p.astype(float32))
     G = NX.Graph(d < 20)   # cutoff 20 is sufficient to connect points
     for  (n,nbrdict) in G.adjacency_iter():
          for k,eattr in nbrdict.items():
               eattr['weight'] = d[n,k]**3   # favours nearest neighbours (mostly...)

     dij3 = NX.dijkstra_path(G, 0, G.nodes()[-1])
     newp = p[dij3]
     savetxt(outfile, newp, "%8.3f %8.3f %8.3f")

     # newp is ordered in traversal order
     deltas = map(linalg.norm, newp[1:] - newp[:-1])	

     G2 = NX.Graph()
     G2.add_path(arange(0,len(newp)), weight=0)
     for  ((n,nbrdict),delta) in zip(G2.adjacency_iter(), deltas):
          for k,eattr in nbrdict.items():
               eattr['weight'] = delta


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

     pruned_coords = newp[pruned.nodes()]

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

# def grow(n0, dist, n,k):
#     try:
#         w = G2.edge[n][k]['weight']
#         dist += w
#         if dist < mindist:
#             grow(n0, dist, n+1, k+1)
#     except KeyError:
#         pass
#     return n0, k, dist

if __name__ == "__main__":
     from optparse import OptionParser
     
     parser = OptionParser(usage=__doc__)
     parser.add_option("--outfile", "-o", dest="outfile",
                       metavar="FILE",
                       help="write sample points to FILE [%default]")
     parser.add_option("--mindist", "-d", dest="mindist", type="float",
                       metavar="FLOAT",
                       help="minimum distance between two points [%default]")
     parser.set_defaults(outfile="/dev/stdout", mindist=1.0)

     opts,args = parser.parse_args()

     try:
          infile = args[0]
     except IndexError:
          infile  = "/dev/stdin"

     modify_path(infile, opts.outfile, mindist=opts.mindist)
          
