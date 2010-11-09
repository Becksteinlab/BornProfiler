# -*- coding: utf-8 -*-
"""Bornprofiler module"""

from __future__ import with_statement

import os, errno
import numpy

import logging
logger = logging.getLogger('bornprofile') 

TABLE_IONS = """
#----------------------------------------
# ionic Born radii (Angstrom) 
#----------------------------------------
# from Table III in Rashin & Honig 1986
#
# values are white-space separated
#
#id symbol atomname radius charge

Li  Li+ LI+    1.316 +1
Na  Na+ NA+    1.680 +1
K   K+  K+     2.172 +1
Rb  Rb+ RB+    2.311 +1
Cs  Cs+ CS+    2.514 +1
F   F-  F-     1.423 -1
Cl  Cl- CL-    1.937 -1
Br  Br- BR-    2.087 -1
I   I-  I-     2.343 -1

Cu1 Cu+ Cu+    1.252 +1
Ag  Ag+ Ag+    1.434 +1
Cu2 Cu+2  Cu2+ 1.242 +2
 
Mg  Mg+2  MG2+ 1.455 +2
Ca  Ca+2  Ca2+ 1.862 +2
Sr  Sr+2  Sr2+ 2.054 +2
Ba  Ba+2  Ba2+ 2.119 +2
Zn  Zn+2  Zn2+ 1.338 +2
Cd  Cd+2  Cd2+ 1.509 +2
Hg  Hg+2  Hg2+ 1.541 +2

Al  Al+3  Al3+ 1.338 +3
Sc  Sc+3  Sc3+ 1.541 +3
Y   Y+3   Y3+  1.733 +3
La  La+3  La3+ 1.808 +3
Ce3 Ce+3  Ce3+ 1.761 +3
Ga  Ga+3  Ga3+ 1.338 +3 
In  In+3  In3+ 1.605 +3

Ce4 Ce+4  Ce4+ 1.761 +4
NH4 NH+4  NH4+ 2.130 +4

OH  OH-   OH-  1.498 -1
SH  SH-   SH-  1.969 -1
S   S-2   S2-  1.969 -2
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
    if len(values) == 0 or line.startswith("#"):
      continue
    ions[values[0]] = Ion(*values)
  return ions

IONS = _parse_TABLE_IONS()

def apbs_component(name, **kwargs):
  return APBS_SCRIPT_COMPONENTS[name] % kwargs

class BPbase(object):
  """Provide basic infra structure methods for bornprofiler classes.

  Defines the file name API.
  """
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
    points = []
    with open(self.pointsName) as pointsFile:
      for line in pointsFile:
        fields = line.split()
        points.append(map(float, fields[0:3]))
    self.points = numpy.array(points).T
 
  def get_XYZ_dict(self, name, vec):
    return {name.upper()+'_XYZ': " ".join(map(str, vec))}

  def get_XYZ_str(self, vec):
    return " ".join(map(str, vec))

# A job script template finds input and output filename in %(infile)s
# and %(outfile)s; the string is interpolated by python.

from config import read_template

JOBSCRIPTS = {
  'local': read_template('q_local.sh'),
  'SBCB':  read_template('q_SBCB.sh'),
  'darthtater': read_template('q_darthtater.sge'),
}

APBS_SCRIPT_COMPONENTS = {
  'header': read_template("000_placeion_header.in"),
  'read':   read_template("001_placeion_read.in"),
  # stage must be one of [cpx,ion,pro]
  'elec':   read_template("002_placeion_elec.in"),
  'printEnergy': read_template("003_placeion_printEnergy.in"),
}

class Placeion(BPbase):
  "preparing job for APBS energy profiling by placing ions"

  padding_xy = 40.0
  padding_z  = 80.0
 
  def __init__(self, pqrfile, pointsfile, ionName='Na', ionicStrength=0.15, 
               jobName='bornprofile', script=None, dime=None, temperature=300.0):
    self.pqrName = pqrfile
    self.pointsName = pointsfile
    self.jobName = jobName
    self.ion = IONS[ionName]
    self.ionicStrength = ionicStrength
    self.dime = dime or [97, 97, 193]
    self.temperature = temperature
    self.script = script

    self.cglen = [0, 0, 0]
    self.pqrLines = []

    self.readPQR()
    self.readPoints()

  def generate(self):
    """Generate all input files."""
    self.writePQRs()
    self.writeIn()
    return self.writeJob()

  def readPQR(self):
    with open(self.pqrName, "r") as pqrFile:
      lines = pqrFile.readlines()
      for line in lines:
        if (line[0:4] == "ATOM"):
          self.pqrLines.append(line)
 
    # finding the cglen (coarse grid lengths) automatically
    # This depends on correct spacing in the PQR file.
    minDim = [+100000, +100000, +100000]
    maxDim = [-100000, -100000, -100000]
    for line in lines:
      if (line[0:4] == "ATOM"):
        tokens = line.split()
        for i in range(0, 3):
          if float(tokens[i+5]) < minDim[i]:
            minDim[i] = float(tokens[i+5])
          if float(tokens[i+5]) > maxDim[i]:
            maxDim[i] = float(tokens[i+5])
    for i in [0, 1]:
      self.cglen[i] = (maxDim[i] - minDim[i]) + self.padding_xy
    self.cglen[2] = (maxDim[2] - minDim[2]) + self.padding_z
 
  def writePQRs(self):
    # make directory
    if (not os.path.exists(self.jobName)):
      os.mkdir(self.jobName)
 
    # write protein
    with open(self.jobpath("pro.pqr"), "w") as proFile:
      for line in self.pqrLines:
        proFile.write(line)
 
    # write complex and ion
    for num,point in enumerate(self.points.T):
      x,y,z = point
      aID = 50000
      rID = 10000
      # born ion is always resname ION
      ionLines = self.getPqrLine(aID, self.ion.atomname, rID, 'ION', x, y, z, self.ion.charge, self.ion.radius)
      # write ion
      with open(self.jobpath(self.filename("ion",num,".pqr")), "w") as ionFile:
        ionFile.write(ionLines)
      # write complex
      with open(self.jobpath(self.filename("cpx",num,".pqr")), "w") as cpxFile:
        for line in self.pqrLines:
          cpxFile.write(line)
        cpxFile.write(ionLines)
    logger.info("[%s] Wrote %d pqr files for protein, ion=%s and complex", 
                self.jobName, num+1, self.ion.atomname)
  
  def writeIn(self):
    for num,point in enumerate(self.points.T):
      z = point[2]
      with open(self.jobpath(self.infilename(num)), "w") as inFile:
        inFile.write(apbs_component('header', z=z))
        inFile.write(apbs_component('read', protein_pqr="pro.pqr",
                                    ion_pqr=self.filename("ion",num,".pqr"),
                                    complex_pqr=self.filename("cpx",num,".pqr")))
        for count,what in enumerate(["pro", "ion", "cpx"]):
          inFile.write(apbs_component('elec', stage=what, molnum=count+1,
                                      DIME_XYZ=self.get_XYZ_str(self.dime),
                                      CGLEN_XYZ=self.get_XYZ_str(self.cglen),
                                      conc=self.ionicStrength, temperature=self.temperature))
        inFile.write(apbs_component('printEnergy', complex='cpx', ion='ion', protein='pro'))

  def writeJob(self):
    if self.script is None:
      logger.warn("No script template provided; no queuing system scripts are written.")
      return None

    for num,point in enumerate(self.points.T):
      z = point[2]
      scriptargs = {
        'jobname': "n%04dz%.3f%s"%(num,z,self.jobName),
        'infile': self.infilename(num),
        'outfile': self.outfilename(num),
        }
      with open(self.jobscriptname(num), "w") as jobFile:
        jobFile.write(self.script % scriptargs)
 
    # TODO: template job array script, too
    qsubName = "qsub_" + self.jobName + ".bash"
    self.numJobs = self.points.shape[-1]
    with open(qsubName, "w") as jobFile:
      jobFile.write("""#$ -N BP%(jobName)s
#$ -S /bin/bash
#$ -l mem_free=500M,mem_total=500M
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-%(numJobs)d

declare -a job
 
""" % vars(self))
      for num,point in enumerate(self.points.T):
        z = point[2]
        jobFile.write('job[%d]="%s"\n' % (num+1, self.jobscriptname(num)))
      jobFile.write("""
run_d=$(dirname ${job[${SGE_TASK_ID}]})
script=$(basename ${job[${SGE_TASK_ID}]})
 
cd ${run_d} || { echo "Failed to cd ${run_d}. Abort."; exit 1; }
. ./${script}
""")
    return num+1


import membrane

class MPlaceion(BPbase):
    def __init__(self, *args, **kwargs):
      """Setup Born profile with membrane.

      MPlaceion(pqr,points[,memclass])

      :Arguments:
        *pqr*
           PQR file
        *points* 
           data file with sample points
      :Keywords:
        *memclass*
           a class or type :class:`APBSMem`.
      """
      self.pqrName = args[0]
      self.pointsName = args[1]
      self.MemSetup = kwargs.pop('memclass', membrane.APBSMem)  # to customize

      memSetups = dict([(suffix, self.MemSetup(self.pqrName, suffix, **kwargs)) for
                        suffix in ('L', 'M', 'S')])
        
        
