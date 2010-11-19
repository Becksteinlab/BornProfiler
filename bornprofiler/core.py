# -*- coding: utf-8 -*-
"""Bornprofiler module"""

from __future__ import with_statement

import os, errno
import numpy

from config import configuration, read_template
from utilities import in_dir, asiterable

import logging
logger = logging.getLogger('bornprofiler.core') 

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

#: Rashin&Honig ion data as a dict, read from :file:`templates/bornions.dat`.
IONS = _parse_TABLE_IONS()

#: A job script template finds input and output filename in %(infile)s
#: and %(outfile)s; the string is interpolated by python.
JOBSCRIPTS = {
  'local': read_template('q_local.sh'),
  'SBCB':  read_template('q_SBCB.sh'),
  'darthtater': read_template('q_darthtater.sge'),
}

class BPbase(object):
  """Provide basic infra structure methods for bornprofiler classes.

  Defines the file name API.

  * simply number files, do not use z or other data in filename;
  * assume standard python 0-based indexing and naming
  * filenames are generated and typically only take *num*, the window
    number, as an argument

  .. Note:: This class cannot be used on its own and it relies on
            other attributes being set by the parent class.
  """
  def jobpath(self, *args):
    return os.path.join(self.jobName, *args)

  def jobscriptname(self, num):
    return self.filename('job', num, '.bash')

  def infilename(self, num):
    return self.filename("job",num,".in")

  def outfilename(self, num):
    return self.filename("job",num,".out")

  def datafile(self, prefix, ext='.dat'):
    return "%s_%s%s" % (prefix, self.jobName, ext)

  def filename(self, prefix, num, ext):
    return '%s_%04d%s' % (prefix,num,ext)

  def window_jobname(self, num):
    return "w%04d_%s" % (num,self.jobName)

  def get_ion_name(self, num):
    return self.filename("ion",num,".pqr")

  def get_complex_name(self, num):
    return self.filename("cpx",num,".pqr")

  def get_protein_name(self, num=None):
    return "pro.pqr"  

  def get_apbs_script_name(self, num=None):
    if num is None:
      return "mem_placeion.in"  # see :class:`membrane.BornAPBSmem`
    return self.infilename(num) # could probably make one of the two funcs redundant    

  def get_windowdirname(self, num):
    return "w%04d" % num

  def _process_window_numbers(self, windows=None):
    """Window numbers are 0-based."""
    if windows is None:
      windows = numpy.arange(self.numPoints)
    else:
      windows = numpy.unique(asiterable(windows))
    if windows.min() < 0  or windows.max() >= self.numPoints:
      errmsg = "window numbers must be between 0 and %d inclusive." % (self.numPoints-1)
      logger.fatal(errmsg)
      raise ValueError(errmsg)
    return windows

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
    self.points = numpy.array(points)
    self.numPoints = self.points.shape[0]
 
  def get_XYZ_dict(self, name, vec):
    return {name.upper()+'_XYZ': " ".join(map(str, vec))}

  def get_XYZ_str(self, vec):
    return " ".join(map(str, vec))

  def __repr__(self):
    return "<%s from pqr=%r points=%r>" % (self.__class__.__name__, 
                                           os.path.basename(self.pqrName), 
                                           os.path.basename(self.pointsName))


# For Placeion:
# change this to a monolithic script and hardcode stages!
# (The partite approach still exists for historic reasons)
APBS_SCRIPT_COMPONENTS = {
  'header': read_template("000_placeion_header.in"),
  'read':   read_template("001_placeion_read.in"),
  # stage must be one of [cpx,ion,pro]
  'elec':   read_template("002_placeion_elec.in"),
  'printEnergy': read_template("003_placeion_printEnergy.in"),
}
def apbs_component(name, **kwargs):
  return APBS_SCRIPT_COMPONENTS[name] % kwargs

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
    for num,point in enumerate(self.points):
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
    for num,point in enumerate(self.points):
      z = point[2]
      with open(self.jobpath(self.infilename(num)), "w") as inFile:
        inFile.write(apbs_component('header', z=z))
        inFile.write(apbs_component('read', protein_pqr=self.get_protein_name(num),
                                    ion_pqr=self.get_ion_name(num),
                                    complex_pqr=self.get_complex_name(num)))
        # TODO: make this a single monolithic template and get rid of 
        # the baroque split-template thing
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
    for num,point in enumerate(self.points):
      z = point[2]
      scriptargs = {
        'jobname': self.window_jobname(num),
        'infile': self.infilename(num),
        'outfile': self.outfilename(num),
        }
      scriptname = self.jobscriptname(num)
      bash_jobarray.append('job[%d]="%s"' % (num+1, scriptname))  # job array ids are 1-based!
      with open(self.jobpath(scriptname), "w") as jobFile:
        jobFile.write(self.script % scriptargs)
 
    qsubName = "qsub_" + self.jobName + ".bash"
    self.numJobs = self.points.shape[0]
    with open(qsubName, "w") as jobFile:
      jobFile.write(TEMPLATES['q_array.sge'] % {
          'jobName': self.jobName,
          'numJobs': self.numJobs,
          'jobArray': "\n".join(bash_jobarray),
          })
    return num+1


class MPlaceion(BPbase):
  """Generate all input files for a Born profile calculation WITH a membrane

  
  """

  #: Schedule is run from first to last: L -> M -> S
  #: choose dime compatible with nlev=4 (in input file)
  schedule = {'dime': [(129, 129, 129),(129, 129, 129),(129, 129, 129)],
              'glen': [(250,250,250),(100,100,100),(50,50,50)],
              }

  def __init__(self, *args, **kwargs):
    """Setup Born profile with membrane.

      MPlaceion(pqr,points[,memclass,jobName,ionName,ionicStrength,temperature,script,arrayscript,basedir])

      :Arguments:
        *pqr*
           PQR file
        *points* 
           data file with sample points

      :Keywords:
        *memclass*
           a class or type :class:`APBSMem` to customize draw_membrane
        *jobName*
           name of the run, used as top directory name and as unique identifier
           [mbornprofile]
        *ionName*
           name of the Rashin&Honig ion to calculate; any ion in IONS works [Na]
        *ionicStrength*
           concentration of monovalent NaCl solution in mol/l [0.15]
        *temperature*
           temperature in K [300.0]
        *script*
           template for a submission script (contains 'abps %(infile)s > %(outfile)s'
           or similar; %(jobname)s is also replaced). Can be (a) a local file, (b) a
           file stored in he user template dir, (c) a bundled template file.
        *arrayscript*
           template for a queuing system array script; contains the the placeholders
           %(jobName)s and %(jobArray)s; jobs are stored in the bash array 'job' so
           there should be a line 'declare -a job'. Window numbers correspond to the
           task ids (SGE) ar array ids (PBS) in the job array. See ``templates/q_array.sge``
           as an example.
        *basedir*
           top directory of set up, defaults to [.]           
    """
    self.__cache_MemBornSetup = {}

    self.pqrName = os.path.realpath(args[0])
    self.pointsName = os.path.realpath(args[1])

    # copied & pasted from Placeion because inheritance would be messy :-p
    self.jobName = kwargs.pop('jobName', "mbornprofile")
    self.ion = IONS[kwargs.pop('ionName', 'Na')]
    self.ionicStrength = kwargs.pop('ionicStrength', 0.15)
    self.temperature = kwargs.pop('temperature', 300.0)
    scriptname = kwargs.pop('script', None)
    if not scriptname is None:
      self.script = read_template(scriptname)
    else:
      self.script = None
    self.arrayscript = read_template(kwargs.pop('arrayscript', 'q_array.sge'))

    # hack for quickly customizing draw_membrane (uses custom classes)
    self.SetupClass = kwargs.pop('memclass', None)
    if self.SetupClass is None:
      import membrane
      self.SetupClass = membrane.BornAPBSmem  # to customize
      
    # sanity check
    assert len(self.schedule) != len(self.SetupClass.suffices), \
        "schedule does not correspond to naming scheme --- make sure that memclass is set up properly"

    # do some initial processing...
    self.readPQR()
    self.readPoints()

  def readPQR(self):
    self.pqrLines = []
    with open(self.pqrName, "r") as pqrFile:
      for line in pqrFile:
        if (line[0:4] == "ATOM"):
          self.pqrLines.append(line)

  def writePQRs(self, windows=None):
    """Generate input pqr files for all windows and store them in separate directories."""
    # NOTE: This is the only place where window numbers reference the sample points
    #       at the moment: point = self.points[num]
    windows = self._process_window_numbers(windows)
    with in_dir(self.jobName):
      # make directory if it does not exist and cd into it      
      for num in windows:
        windowdirname = self.get_windowdirname(num)
        with in_dir(windowdirname):
          # write protein (make hard-link to safe space)
          try:
            os.link(self.pqrName, self.get_protein_name())
          except OSError, err:
            if err.errno != errno.EEXIST:
              raise
          # write complex and ion 
          x,y,z = self.points[num]  # index points with window id
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
                self.jobName, len(windows), self.ion.atomname)

  def writeJob(self, windows=None):
    """Write the job script.

    Can be done selectively for a subset of windows and then the
    global array script will also only contain this subset of windows.
    """
    if self.script is None:
      logger.warn("No script template provided; no queuing system scripts are written.")
      return None

    windows = self._process_window_numbers(windows)

    bash_jobarray = []
    for num in windows:
      m = self.get_MemBornSetup(num)   # to access anything on a per-window basis
      scriptargs = {
        'jobname': self.window_jobname(num),
        'infile': self.infilename(num), # XXX: could take this from m, too...
        'outfile': self.outfilename(num),
        'unpack_dxgz': m.unpack_dxgz,   # XXX: these two really should be attributes of
        'apbs_version': m.apbs_version, # of MPlaceion but too lazy right now...
        'drawmembrane_script': m.filenames['born_setup_script'],
        }      
      scriptname = self.jobscriptname(num)
      scriptpath = self.jobpath(self.get_windowdirname(num), scriptname)
      bash_jobarray.append('job[%d]="%s"' % (num+1, scriptpath))    # job array ids are 1-based!
      # write scriptfile into the window directory
      with open(scriptpath, "w") as jobFile:
        jobFile.write(self.script % scriptargs)
 
    qsubName = "qsub_" + self.jobName + ".bash"
    self.numJobs = self.points.shape[0]  # always record maximum number
    with open(qsubName, "w") as jobFile:
      jobFile.write(self.arrayscript % {
          'jobName': self.jobName,
          'numJobs': self.numJobs,
          'jobArray': "\n".join(bash_jobarray),
          })
    return len(windows)


  def generateMem(self, windows=None, run=True):
    """Generate special diel/kappa/charge files and mem_placeion.in for each window.

    The APBS calculation is set up for manual focusing in three stages:
      1. **L** is a coarse grid and centered on the protein
      2. **M** is a medium grod, centered on the ion, and using focusing
         (the boundary values come from the **L** calculation)
      3. **S** is the finest grid, also centered and focused on the ion

    The ion positions ("windows") are simply sequentially numbered, starting at
    1, as they appear in the input file. Each window is set up in its own
    directiroy (called "wNNNN" where NNNN is the 4-digit, zero-padded number).

    .. Warning:: At the moment it is not checked that the "inner" focusing region M
       are always contained in the outer region L; it's on the TODO list.

    :Keywords:
       *windows* : ``None``, number, or list
          window number or list of window numbers to generate. If set to
          ``None`` (the default) then all windows are generated. Window
          numbers start at 0 and end at numPoints-1.
       *run* : bool
          ``True``: immediately generate files (can take a while); ``False`` defer
          file generation and just write a script
    """

    windows = self._process_window_numbers(windows)
    for num in windows:
      windowdirpath = self.jobpath(self.get_windowdirname(num))
      with in_dir(windowdirpath, create=False):
        # will throw an error if the directory does not exist -- currently users
        # responsibility to set it up properly (i.e. generate diel maps etc)
        membornsetup = self.get_MemBornSetup(num)
        membornsetup.generate(run=run)

  def get_MemBornSetup(self, num):
    """Return the setup class for window *num* (cached)."""
    try:
      return self.__cache_MemBornSetup[num]
    except KeyError:
      protein = self.get_protein_name()
      ion = self.get_ion_name(num)
      cpx = self.get_complex_name(num)
      apbs_script_name = self.get_apbs_script_name(num)
      # should get draw_membrane parameters as well, do this once we use cfg input files
      # TODO: add position of ion to comments
      kw = {'conc':self.ionicStrength, 'temperature':self.temperature,'comment':'window %d, z=XXX'%num,
            'basedir': os.path.realpath(self.jobpath(self.get_windowdirname(num))),
            'apbs_script_name': apbs_script_name}
      kw.update(self.schedule)  # set dime and glen !
      # using a custom SetupClass pre-populates the parameters for draw_membrane2
      # (hack!! -- should be moved into a cfg input file)
      self.__cache_MemBornSetup[num] = self.SetupClass(protein, ion, cpx, **kw)
    return self.__cache_MemBornSetup[num]

  def generate(self, windows=None, run=True):
    """Set up all input files for Born calculations with a membrane.

    generate([windows[,run]])

    The optional parameter *windows* allows one to select a subset of
    windows instead of all the points in the sample points file; it
    can be a single number or a list.

    Setting *run* to ``True`` (the default) immediately generates setup files,
    in particular it runs :program:`apbs` and :program:`draw_membrane2a` in
    order to add the membrane to the system. Because this can take a long time
    for a large number of windows it is also possible to only generate a bas
    script for each window and defer the setup (*run* = ``False``).
    """
    windows = self._process_window_numbers(windows)
    self.writePQRs(windows=windows)
    self.generateMem(windows=windows, run=run)
    self.writeJob(windows=windows)


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
