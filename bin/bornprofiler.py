# -*- coding: utf-8 -*-
"""Bornprofiler module"""

from __future__ import with_statement

import os, errno
import logging
logger = logging.getLogger('bornprofile') 

# ionic radii (Angstrom) (all from Table III in Rashin & Honig)
#
#id symbol atomname radius charge
TABLE_IONS = """
Li  Li+ LI+    1.316 +1
Na  Na+ NA+    1.680 +1
K   K+  K+     2.172 +1
Rb  Rb+ RB+    2.311 +1
Cs  Cs+ CS+    2.514 +1
F   F-  F-     1.423 -1
Cl  Cl- CL-    1.937 -1
Br  Br- BR-    2.087 -1
I   I-  I-     2.343 -1
Mg  Mg+2  MG2+ 1.455 +2
Spermidine  SPM SPM 2.130  +1
"""

class Ion(dict):
  def __init__(self, name, symbol, atomname, radius, charge):
    super(Ion, self).__init__(name=name, symbol=symbol, atomname=atomname,
                              radius=float(radius), charge=float(charge))
  def __getattribute__(self, name):
    try:
      return self[name]
    except KeyError:
      return super(Ion, self).__getattribute__(name)

def _parse_TABLE_IONS():
  ions = {}
  lines = TABLE_IONS.split('\n')
  for line in lines:
    values = line.split()
    if len(values) == 0:
      continue
    ions[values[0]] = Ion(*values)
  return ions

IONS = _parse_TABLE_IONS()

# A job script template finds input and output filename in %(infile)s
# and %(outfile)s; the string is interpolated by python.



class BPbase(object):

  # number files, do not use z in filename
  def jobpath(self, *args):
    return os.path.join(self.jobName, *args)

  def jobscriptname(self, num):
    return self.jobpath(self.filename('job', num, '.bash'))

  def infilename(self, num):
    return self.filename("job",num,".in")

  def outfilename(self, num):
    return self.filename("job",num,".out")

  def datafile(self, prefix, ext='.dat'):
      return "%s_%s%s" % (prefix, self.jobName, ext)

  def filename(self, prefix, num, ext):
    return '%s_%04d%s' % (prefix,num,ext)

  def getPqrLine(self, aID, aType, rID, rType, x, y, z, q, r):
    # example: "ATOM  50000  NHX SPM 10000      43.624  57.177  58.408  1.0000 2.1300"
    line =  "ATOM  "
    line += "%5d  " % aID
    line += "%-4s" % aType
    line += "%3s " % rType
    line += "%5d    " % rID
    line += "%8.3f%8.3f%8.3f " % (x, y, z)
    line += "%7.4f%7.4f\n" % (q, r)
    return line
 
  def readPoints(self):
    with open(self.pointsName) as pointsFile:
      lines = pointsFile.readlines()
 
    points = []
    for line in lines:
      tokens = line.split()
      parsed = [float(tokens[0]), float(tokens[1]), float(tokens[2])]
      points.append(parsed)
    self.points = points
 
  def get_XYZ_dict(self, name, vec):
    return {name.upper()+'_XYZ': " ".join(map(str, vec))}

  def get_XYZ_str(self, vec):
    return " ".join(map(str, vec))

