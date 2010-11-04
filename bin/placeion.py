#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
:Author: Kaihsu Tai, Oliver Beckstein
:Year: 2008, 2010
:Licence: GPL
:Copyright: (c) 2008 Kaihsu Tai
:Copyright: (c) 2010 Oliver Beckstein
:URL: http://en.wikiversity.org/wiki/Talk:Poisson%E2%80%93Boltzmann_profile_for_an_ion_channel

I wrote some Python code to automate this process. The job submission
requires a queuing system called Grid Engine. Copyright © 2008 Kaihsu
Tai. Moral rights asserted (why?). Hereby licensed under either GFDL
or GNU General Public License at your option.
"""

import os
import logging
logger = logging.getLogger('bornprofile') 

usage = """%programe [options] PQR samplepoints

This script sets up input files for Born energy calculations for
APBS. It expects to find ion positions in the file samplepoints.

The "Born radii" for ions (not Pauling radii!) were taken from Table III in

  Alexander A. Rashin, Barry Honig (1985) J. Phys. Chem. 89(26):5588-5593
  http://dx.doi.org/10.1021/j100272a006

This paper suggests using the corrected covalent radius (Born radius).
"""

# radius of the test ion (Angstrom)
# ionic radii (Angstrom)
# (all from Table III in Rashin & Honig)

class Ion(dict):
  def __init__(self, name, symbol, atomname, radius, charge):
    super(Ion, self).__init__(name=name, symbol=symbol, atomname=atomname,
                              radius=float(radius), charge=float(charge))
  def __getattribute__(self, name):
    try:
      return self[name]
    except KeyError:
      return super(Ion, self).__getattribute__(name)

#id symbol atomname radius charge
# I use OPLS/Gromacs atom names
_ions = """
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

def _setup_ions():
  lines = _ions.split('\n')
  ions = {}
  for line in lines:
    values = line.split()
    if len(values) == 0:
      continue
    ions[values[0]] = Ion(*values)
  return ions

IONS = _setup_ions()

# A job script template finds input and output filename in %(infile)s
# and %(outfile)s; the string is interpolated by python.

JOBSCRIPTS = {
  'local':
"""#!/bin/sh
echo "**  %(jobname)s"
echo "** APBS Born profile job running on $HOSTNAME"
apbs %(infile)s > %(outfile)s

""",
  'SBCB':
"""#$ -N %(jobname)s
#$ -S /bin/bash
#$ -l mem_free=500M,mem_total=500M
#$ -cwd
#$ -j y
#$ -r y

echo "APBS Born profile job running on $HOSTNAME"
ulimit -c 64
module load abps/32 
apbs %(infile)s > %(outfile)s

""",
}


class Placeion(object):
  "preparing job for APBS energy profiling by placing ions"

  padding_xy = 40.0
  padding_z  = 80.0
 
  def __init__(self, pqrfile, pointsfile, ionName='Na', ionicStrength=0.15, 
               jobName='bornprofile', script='SBCB'):
    self.pqrName = pqrfile
    self.pointsName = pointsfile
    self.jobName = jobName
    self.ion = IONS[ionName]
    self.ionicStrength = ionicStrength
    self.script = script

    self.cglen = [0, 0, 0]
    self.pqrLines = []
 
  def readPQR(self):
    pqrFile = open(self.pqrName, "r")
    lines = pqrFile.readlines()
    for line in lines:
      if (line[0:4] == "ATOM"):
        self.pqrLines.append(line)
    pqrFile.close()
 
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
    proFile = open(self.jobName + "/pro" + ".pqr", "w")
    for line in self.pqrLines:
      proFile.write(line)
    proFile.close()
 
    # write complex and ion
    for point in self.points:
      x = point[0]
      y = point[1]
      z = point[2]
      aID = 50000
      rID = 10000
      ionLines = ""

      # born ion is always resname ION
      ion = self.ion
      ionLines += self.getPqrLine(aID, ion.atomname, rID, 'ION', x, y, z, ion.charge, ion.radius)
 
      # write ion
      ionFile = open(self.jobName + "/ion_" + str(z) + ".pqr", "w")
      ionFile.write(ionLines)
      ionFile.close()
 
      # write complex
      cpxFile = open(self.jobName + "/cpx_" + str(z) + ".pqr", "w")
      for line in self.pqrLines:
        cpxFile.write(line)
      cpxFile.write(ionLines)
      cpxFile.close()
 
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
    pointsFile = open(self.pointsName, "r")
    lines = pointsFile.readlines()
    pointsFile.close()
 
    points = []
    for line in lines:
      tokens = line.split()
      parsed = [float(tokens[0]), float(tokens[1]), float(tokens[2])]
      points.append(parsed)
    self.points = points
 
  def writeIn(self):
    for point in self.points:
      z = point[2]
      with open(self.jobName + "/job_" + str(z) + ".in", "w") as inFile:
        inFile.write("""# APBS input file generated by """ + __file__ + """
# for the complex and the ion at z = """ + str(z) + "\n") 
 
        inFile.write("""
  read
    mol pqr pro.pqr
    mol pqr ion_""" + str(z) + """.pqr
    mol pqr cpx_""" + str(z) + """.pqr
  end
""")
 
        count = 1
        for what in ["pro", "ion", "cpx"]:
          inFile.write("""
  elec name """ + what + """
    mg-auto
    mol """ + str(count) + """
    dime 97 97 193
    cglen """ + str(self.cglen[0]) + " " + str(self.cglen[1]) + " " + str(self.cglen[2]) + """
    fglen 20 20 20    
    cgcent mol 3
    fgcent mol 2
    # NaCl ionic strength in mol/L
    ion  1 """ + str(self.ionicStrength) + """ 0.95 # sodium ions
    ion -1 """ + str(self.ionicStrength) + """ 1.81 # chloride ions
 
    lpbe
    bcfl mdh
    pdie 2.0 # protein and faux-lipid
    sdie 78.5 # Eisenberg and Crothers Phys. Chem. book 1979
    srfm smol
    chgm spl2
    srad 1.4
    swin 0.3
    sdens 10.0
    temp 300
    # gamma 0.105 # Uncomment for old versions of APBS: deprecated for APBS 1.0.0
    calcenergy total
    calcforce no
  end
""")
          count += 1
 
        inFile.write("""
  print energy cpx - ion - pro end
  quit
""")

 
  def writeJob(self):
    for point in self.points:
      z = point[2]
      scriptargs = {
        'jobname': "z" + str(z) + self.jobName,
        'infile': "job_" + str(z) + ".in",
        'outfile': "job_" + str(z) + ".out",
        }
      with open(self.jobName + "/job_" + str(z) + ".bash", "w") as jobFile:
        jobFile.write(self.script % scriptargs)
 
    # TODO: template job array script, too
    qsubName = "qsub_" + self.jobName + ".bash"
    with open(qsubName, "w") as jobFile:
      jobFile.write("""#$ -N """ + self.jobName + """
#$ -S /bin/bash
#$ -l mem_free=500M,mem_total=500M
#$ -cwd
#$ -j y
#$ -r y
 
declare -a job
 
""")
      count = 0
      for point in self.points:
        z = point[2]
        jobFile.write("job[" + str(count+1) + ']="' + self.jobName + "/job_" + str(z) + ".bash" + '"\n')
        count += 1
        jobFile.write("""
run_d=$(dirname ${job[${SGE_TASK_ID}]})
script=$(basename ${job[${SGE_TASK_ID}]})
 
cd ${run_d} || { echo "Failed to cd ${run_d}. Abort."; exit 1; }
chmod a+x ${script}
 
exec ./${script}
""")
    return count
 
  def generate(self):
    """Generate all input files."""
    self.readPQR()
    self.readPoints()
    self.writePQRs()
    self.writeIn()
    count = self.writeJob()
    return count
  
if __name__ == "__main__":
  import sys
  from optparse import OptionParser

  logging.basicConfig()

  parser = OptionParser(usage=usage)
  parser.add_option("--ionic-strength", dest="ionicStrength",
                    metavar="CONC",
                    help="set ionic strength of Na/Cl bath to the given concentration "
                    "CONC in mol/l [%default]")
  parser.add_option("--name", dest="jobName",
                    metavar="STRING",
                    help="name for the job submission script [%default]")
  parser.add_option("--ion", dest="ionName",
                    metavar="NAME",
                    help="name of the ion to be sampled. Available values:\n%r\n"
                    "Radii were taken from Table III in Rashin & Honig  1985. "
                    "The default ion is '%%default'" % (IONS.keys(),))
  parser.add_option("--script", dest="script",
                    metavar="NAME",
                    help="name of a stored script template or (advanced usage!) a "
                    "filename which contains appropriate place holders [%default]")
  parser.set_defaults(ionicStrength=0.15, jobName="bornprofile", 
                 ionName="Na", script="local")

  opts,args = parser.parse_args()
  
  try:
    pqrfile, pointsfile = args
  except:
    logger.fatal("Needs PQR file and sample points. See --help.")
    sys.exit(1)

  try:
    script_template = JOBSCRIPTS[opts.script]
  except KeyError:
    try:
      script_template = open(opts.script).readlines()
    except:
      logger.fatal("--scripts=NAME must be either a filename or one of %r; see the " 
                   "source for the correct format of the file." % JOBSCRIPTS.keys())
      sys.exit(2)

  P = Placeion(pqrfile, pointsfile, ionName=opts.ionName, ionicStrength=opts.ionicStrength,
               jobName=opts.jobName, script=script_template)
  P.generate()
