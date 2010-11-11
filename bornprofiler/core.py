# -*- coding: utf-8 -*-
"""Bornprofiler module"""

from __future__ import with_statement

import os, errno
import numpy

from config import configuration, read_template
from utilities import in_dir, asiterable

import logging
logger = logging.getLogger('bornprofile') 

TEMPLATES = {'born': read_template('mplaceion.in'),
             'q_array.sge': read_template('q_array.sge'),
             }
TABLE_IONS = read_template('bornions.dat')

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

  # naming scheme!!
  def get_ion_name(self, num):
    return self.filename("ion",num,".pqr")
  def get_complex_name(self, num):
    return self.filename("cpx",num,".pqr")
  def get_protein_name(self, num=None):
    return "pro.pqr"  
  def get_windowdirname(self, num):
    return "w%04d" % num


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

# change this to a monolitic script and hardcode stages!
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
      for line in pqrFile:
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
    with open(self.jobpath(self.get_protein_name()), "w") as proFile:
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
      with open(self.jobpath(self.get_ion_name(num)), "w") as ionFile:
        ionFile.write(ionLines)
      # write complex
      with open(self.jobpath(self.get_complex_name(num)), "w") as cpxFile:
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

    bash_jobarray = []
    for num,point in enumerate(self.points.T):
      z = point[2]
      scriptargs = {
        'jobname': "n%04dz%.3f%s"%(num,z,self.jobName),
        'infile': self.infilename(num),
        'outfile': self.outfilename(num),
        }
      scriptname = self.jobscriptname(num)
      bash_jobarray.append('job[%d]="%s"' % (num+1, scriptname))
      with open(scriptname, "w") as jobFile:
        jobFile.write(self.script % scriptargs)
 
    qsubName = "qsub_" + self.jobName + ".bash"
    self.numJobs = self.points.shape[-1]
    with open(qsubName, "w") as jobFile:
      jobFile.write(TEMPLATES['q_array.sge'] % {
          'jobName': self.jobName,
          'numJobs': self.numJobs,
          'jobArray': "\n".join(bash_jobarray),
          })
    return num+1


import membrane
from utilities import AttributeDict

def ngridpoints(c, nlev=4):
  """The allowed number of grid points.

  For mg-manual calculations, the arguments are dependent on the choice of *nlev*
  by the formula

    n = c*2**(nlev + 1) + 1

  where n is the dime argument, c is a non-zero integer, lev is the nlev
  value. The most common values for grid dimensions are 65, 97, 129, and 161
  (they can be different in each direction); these are all compatible with a
  nlev value of 4. If you happen to pick a "bad" value for the dimensions
  (i.e., mismatch with nlev), the APBS code will adjust the specified dime
  downwards to more appropriate values. This means that "bad" values will
  typically result in lower resolution/accuracy calculations!

  http://www.poissonboltzmann.org/apbs/user-guide/running-apbs/input-files/elec-input-file-section/elec-keywords/dime
  """
  return c*2**(nlev+1) + 1


class MPlaceion(BPbase):
  #: Schedule is run from first to last: L -> M -> S
  #: choose dime compatible with nlev=4 (in input file)
  schedule = {'dime': [(129, 129, 129),(129, 129, 129),(129, 129, 129)],
              'glen': [(250,250,250),(100,100,100),(50,50,50)],
              }

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
    self.pqrName = os.path.realpath(args[0])
    self.pointsName = os.path.realpath(args[1])

    # copied & pasted from Placeion because inheritance would be messy :-p
    self.jobName = kwargs.pop('jobName', "mbornprofile")
    self.ion = IONS[kwargs.pop('ionName', 'Na')]
    self.ionicStrength = kwargs.pop('ionicStrength', 0.15)
    self.temperature = kwargs.pop('temperature', 300.0)  # or use 'temp' ??
    self.script = kwargs.pop('script', None)

    # hack...
    self.SetupClass = kwargs.pop('memclass', membrane.BornAPBSmem)  # to customize

    self.readPQR()
    self.readPoints()

  def readPQR(self):
    self.pqrLines = []
    with open(self.pqrName, "r") as pqrFile:
      for line in pqrFile:
        if (line[0:4] == "ATOM"):
          self.pqrLines.append(line)

  def writePQRs(self):
    with in_dir(self.jobName):
      # make directory if it does not exist and cd into it
      
      for num,point in enumerate(self.points.T):
        windowdirname = self.get_windowdirname(num)
        with in_dir(windowdirname):
          # write protein (make hard-link to safe space)
          # HARDCODED name of protein pqr
          try:
            os.link(self.pqrName, self.get_protein_name())
          except OSError, err:
            if err.errno != errno.EEXIST:
              raise

          # write complex and ion 
          x,y,z = point
          aID = 99999
          rID = 9999
          # born ion is always resname ION
          ionLines = self.getPqrLine(aID, self.ion.atomname, rID, 'ION', x, y, z, self.ion.charge, self.ion.radius)
          # write ion
          with open(self.get_ion_name(num), "w") as ionFile:
            ionFile.write(ionLines)
          # write complex
          with open(self.get_complex_name(num), "w") as cpxFile:
            for line in self.pqrLines:
              cpxFile.write(line)
            cpxFile.write(ionLines)
    logger.info("[%s] Wrote %d pqr files for protein, ion=%s and complex", 
                self.jobName, num+1, self.ion.atomname)

  def generateMem(self, windows=None):
    """Generate special diel/kappa/charge files and mem_placeion.in for each window."""

    # use hard coded schedule for the moment, make it a table later
    #self.memSetups = [cls(self.pqrName, s.suffix, **s.as_kwargs()) for s in self.schedule]

    # TODO:
    # center i>0 on the ion
    # --> actually, need to set up for every single window because the focus
    # region changes and needs a new map
    # Also not sure if draw_membrane is still happy if there's no membrane in the
    # region.
    # - need a flexible ABPSmem and/or generate the template files on the
    #   fly
    # - do one directory per window because file naming is a desaster
    
    # do one BP calculation by hand and then work from there...

    if windows is None:
      windows = numpy.arange(1, self.numPoints+1)
    else:
      windows = asiterable(windows)

    for num in windows:
      windowdirpath = self.jobpath(self.get_windowdirname(num))
      with in_dir(windowdirpath, create=False):
        protein = self.get_protein_name()
        ion = self.get_ion_name(num)
        cpx = self.get_complex_name(num)
        # should get draw_membrane parameters as well, do this once we use cfg input files
        # TODO: add position of ion to comments
        kw = {'conc':self.ionicStrength, 'temperature':self.temperature,'comment':'window %d, z=XXX'%num,
              }
        kw.update(self.schedule)  # set dime and glen !
        # using a custom SetupClass pre-populates the parameters for draw_membrane2
        # (hack!! -- should be moved into a cfg input file)
        membornsetup = self.SetupClass(protein, ion, cpx, **kw)
        membornsetup.generate()
