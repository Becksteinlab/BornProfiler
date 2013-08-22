#!/usr/bin/env python
# Code takes PDBID of protein, finds the coordinate file on the opm database, &
# makes use of the dummy atoms to return the membrane thickness and position,
# relative to the geometrical center of the protein. Also finds the maximal
# protein radius at membrane edges. In general membrane will need to be given
# a smaller radius than this.
import numpy
import MDAnalysis
import MDAnalysis.analysis
import sys
import os
import urllib2
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('PDBID')

args = parser.parse_args()
def get_width(positions, center):
  max_diffsq = 0
  for position in positions:
    diffsq = 0
    i = 0
    while i < 2:
        diffsq += (position[i] - center[i])**2
        i += 1
    if diffsq > max_diffsq:
      max_diffsq= diffsq
  return max_diffsq**(0.5)


def memplacer(PDBID):
    pdburl = 'http://opm.phar.umich.edu/pdb/{pdb}.pdb'.format(pdb=PDBID)
    dl = urllib2.urlopen(pdburl)
    pdb = open('{pdb}.pdb'.format(pdb=PDBID), 'w')
    pdb.write(dl.read())
    pdb.close()
    U = MDAnalysis.Universe('{pdb}.pdb'.format(pdb=PDBID))
    leaflet0 = MDAnalysis.analysis.leaflet.LeafletFinder(U, "resname DUM and prop z < 0").groups(0)
    leaflet1 = MDAnalysis.analysis.leaflet.LeafletFinder(U, "resname DUM and prop z > 0").groups(0)
    protein = U.selectAtoms("protein")
    centroid = protein.centroid()
    zbot = leaflet0.positions[0][2]
    ztop = leaflet1.positions[0][2]
    dist = centroid[2] - zbot
    thickness = ztop - zbot
    ztoptop = ztop + .5
    ztopbottom = ztop - .5
    zbottop = zbot + .5
    zbotbottom= zbot - .5
    bottomlayer = protein.selectAtoms("prop z > {zbb} and prop z < {zbt}".format(zbb=zbotbottom, zbt = zbottop))
    toplayer = protein.selectAtoms("prop z > {ztb} and prop z < {ztt}".format(ztb=ztopbottom, ztt = ztoptop))
    botradius = get_width(bottomlayer.positions, leaflet0.centroid())
    topradius = get_width(toplayer.positions, leaflet1.centroid())
    print("Membrane thickness = {thicky}, beginning at {disty} below protein centroid \n max radius of bottom membrane exclusion = {bot} \n max radius of top membrane exclusion = {top} ".format(thicky=thickness, disty=dist, bot = botradius, top = topradius))
memplacer(args.PDBID)   
