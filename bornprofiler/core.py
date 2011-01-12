# -*- coding: utf-8 -*-
"""Bornprofiler module"""

from __future__ import with_statement

import os, errno
import numpy

import io
from config import configuration, read_template
from utilities import in_dir, asiterable

import logging
logger = logging.getLogger('bornprofiler.core') 

# TODO: get rid of this TEMPLATE dict and hard-code directly in WriteIn()?
TEMPLATES = {'placeion': read_template('placeion.in'), }

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
    # PQR is white space separated!
    # but try to be close to PDB... http://www.wwpdb.org/documentation/format32/sect9.html
    #         1         2         3         4         5         6         7
    #123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
    #ATOM  seria name res reSeq    x-------y-------z-------occup-tempFa
    fmt =  "ATOM  %(aID)5d %(aType)-4s %(rType)3s %(rID)5d    %(x)7.3f %(y)7.3f %(z)7.3f %(q)5.3f %(r)5.3f\n"
    return fmt % vars()
 
  def readPoints(self):
    """Read positions for Born ions from data file.

    Tries to be smart and autodetect standard x-y-z dat file or pdb.
    """
    self.points = io.readPoints(self.pointsName)
    self.numPoints = self.points.shape[0]

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
          # write protein (make hard-link to save space)
          try:
            os.link(self.pqrName, self.get_protein_name())
          except OSError, err:
            if err.errno == errno.EEXIST:
              pass
            elif err.errno == errno.EXDEV:
              import shutil
              shutil.copy(self.pqrName, self.get_protein_name())
            else:
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
 
  def get_XYZ_dict(self, name, vec):
    return {name.upper()+'_XYZ': " ".join(map(str, vec))}

  def get_XYZ_str(self, vec):
    return " ".join(map(str, vec))

  def __repr__(self):
    return "<%s from pqr=%r points=%r>" % (self.__class__.__name__, 
                                           os.path.basename(self.pqrName), 
                                           os.path.basename(self.pointsName))


class Placeion(BPbase):
  "preparing job for APBS energy profiling by placing ions"

  padding_xy = 40.0  # xy only used with use_cubic_boundaries = False
  padding_z  = 80.0
 
  def __init__(self, *args, **kwargs):
    if len(args) == 1:
      # new style
      import io
      params = io.RunParameters(args[0])
      self.bornprofile_kwargs = kw = params.get_bornprofile_kwargs()
      self.pqrName = os.path.realpath(kw.pop('pqr'))
      self.pointsName = os.path.realpath(kw.pop('points'))
      self.ion = IONS[kw.pop('ion', 'Na')]
      self.jobName = kw.pop('name', "bornprofile")
      self.ionicStrength = kw.pop('conc', 0.15)
      self.temperature = kw.pop('temperature', 300.0)
      self.sdie = kw.pop('sdie', 78.5)
      self.pdie = kw.pop('pdie', 10.0)
      self.arrayscript = read_template(kw.pop('arrayscript', 'q_array.sge'))
      self.script = read_template(kw.pop('script', 'q_local.sh'))
      dime = numpy.array(kw.pop('dime', [129,129,129]))
      if len(dime.shape) == 2:
        self.dime = dime[0]   # only take the first one in a triplet
      else:
        self.dime = dime
      self.fglen = numpy.array(kw.pop('fglen', [40,40,40]))
      # glen not needed, auto-set from pqr and padding
    else:
      import warnings
      warnings.warn("Using deprecated Placeion(pqr,points) call", DeprecationWarning)
      self.pqrName = os.path.realpath(args[0])
      self.pointsName = os.path.realpath(args[1])
      self.jobName = kwargs.pop('jobName', "bornprofile")
      self.ion = IONS[kwargs.pop('ionName', 'Na')]
      self.ionicStrength = kwargs.pop('ionicStrength', 0.15)
      self.temperature = kwargs.pop('temperature', 300.0)
      self.sdie = kwargs.pop('sdie', 78.5)
      self.pdie = kwargs.pop('pdie', 10.0)
      self.dime = numpy.array(kwargs.pop('dime', [129,129,129]))
      self.fglen = numpy.array(kwargs.pop('fglen', [40,40,40]))
      self.script = read_template(kwargs.pop('script', 'q_local.sh'))
      self.arrayscript = read_template(kwargs.pop('arrayscript', 'q_array.sge'))

    self.cglen = numpy.array([0, 0, 0])  # automatically set by readPQR() + padding!
    self.use_cubic_boundaries = True     # False uses old placeion.py algorithm: manually adjust fglen!!
    self.pqrLines = []

    logger.info("Placeion: pqr=%(pqrName)r", vars(self))
    logger.info("Placeion: points=%(pointsName)r", vars(self))
    logger.info("Placeion: ion=%(ion)r", vars(self))

    self.readPQR()
    self.readPoints()

  def generate(self):
    """Generate all input files."""
    self.writePQRs()
    self.writeIn()
    return self.writeJob()

  def readPQR(self):
    """Read PQR file and determine the bounding box."""
    with open(self.pqrName, "r") as pqrFile:
      for line in pqrFile:
        if (line[0:4] == "ATOM"):
          self.pqrLines.append(line)
 
    # finding the cglen (coarse grid lengths) automatically
    _coords = []
    for line in self.pqrLines:
      if (line[0:4] == "ATOM"):
        fields = line.split()      # This depends on correct spacing in the PQR file.
        _coords.append(map(float, fields[5:8]))
    coords = numpy.array(_coords)
    if self.use_cubic_boundaries:
      # use largest cubic box (note: fglen is cubic in templates/placeion.in)
      # but can be set in run parameters bornprofile.fglen.
      paddings = numpy.max([self.padding_xy, self.padding_z])
      L = (coords.max() - coords.min()) + paddings
      self.cglen = numpy.array([L,L,L])
    else:
      # original placeion.py algorithm
      # (note: for this to work better, the fglen in templates/placeion.in should
      # match the ratio of box lengths)
      paddings = numpy.array([self.padding_xy, self.padding_xy, self.padding_z])
      self.cglen = (coords.max(axis=0) - coords.min(axis=0)) + paddings
         
  def writeIn(self, windows=None):
    windows = self._process_window_numbers(windows)
    for num in windows:
      x,y,z = self.points[num]
      # old-style ... manual setting up :-p
      params = {'x': x, 'y': y, 'z': z, 'protein_pqr': self.get_protein_name(num),
                'ion_pqr': self.get_ion_name(num), 'complex_pqr': self.get_complex_name(num),
                'DIME_XYZ': self.get_XYZ_str(self.dime), 'CGLEN_XYZ': self.get_XYZ_str(self.cglen),
                'FGLEN_XYZ': self.get_XYZ_str(self.fglen),
                'conc': self.ionicStrength, 'temperature': self.temperature,
                'sdie': self.sdie, 'pdie': self.pdie}
      inPath = self.jobpath(self.get_windowdirname(num), self.infilename(num))
      with open(inPath, "w") as inFile:
        inFile.write(TEMPLATES['placeion'] % params)

  def writeJob(self, windows=None):
    # TODO: consolidate with MPlaceion.writeJob
    if self.script is None:
      logger.warn("No script template provided; no queuing system scripts are written.")
      return None

    windows = self._process_window_numbers(windows)

    bash_jobarray = []
    for num in windows:
      scriptargs = {
        'jobname': self.window_jobname(num),
        'infile': self.infilename(num),
        'outfile': self.outfilename(num),
        'drawmembrane_script': 'false',  # hack so that we can use th MPlaceion scripts...
        'unpack_dxgz': 'false',          # hack so that we can use th MPlaceion scripts...
        'apbs_version': '*',             # hack so that we can use th MPlaceion scripts...
        }
      scriptname = self.jobscriptname(num)
      scriptpath = self.jobpath(self.get_windowdirname(num), scriptname)
      bash_jobarray.append('job[%d]="%s"' % (num+1, scriptpath))  # job array ids are 1-based!
      with open(scriptpath, "w") as jobFile:
        jobFile.write(self.script % scriptargs)
      os.chmod(scriptpath, 0755)

    qsubName = "qsub_" + self.jobName + ".bash"
    self.numJobs = self.points.shape[0]
    with open(qsubName, "w") as jobFile:
      jobFile.write(self.arrayscript % {
          'jobName': self.jobName,
          'numJobs': self.numJobs,
          'jobArray': "\n".join(bash_jobarray),
          })
    return len(windows)


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

      MPlaceion(paramfile[,basedir])

      MPlaceion(pqr,points[,memclass,jobName,ionName,ionicStrength,temperature,script,arrayscript,basedir])

      :Arguments:
        *parameterfile* 
           ini-style file containing all run parameters
        *pqr*
           PQR file                        **DEPRECATED**
        *points* 
           data file with sample points    **DEPRECATED**

      :Keywords:        
        *memclass*
           a class or type :class:`APBSMem` to customize draw_membrane
           **DEPRECATED**
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

    if len(args) == 1:
      # new style
      import io
      params = io.RunParameters(args[0])
      self.bornprofile_kwargs = kw = params.get_bornprofile_kwargs()
      self.pqrName = os.path.realpath(kw.pop('pqr'))
      self.pointsName = os.path.realpath(kw.pop('points'))
      self.ion = IONS[kw.pop('ion', 'Na')]
      self.jobName = kw.pop('name', "mbornprofile")
      self.ionicStrength = kw.pop('conc', 0.15)
      self.temperature = kw.pop('temperature', 300.0)
      self.arrayscript = read_template(kw.pop('arrayscript', 'q_array.sge'))
      self.script = read_template(kw.pop('script', 'q_local.sh'))

      # filter stuff... probably should do this in a cleaner manner (e.g. different
      # sections in the parameter file?)
      kw.pop('fglen', None)

      import membrane
      self.SetupClass = membrane.BornAPBSmem  # use parameters to customize (see get_MemBornSetup())
    else:
      import warnings
      warnings.warn("Using deprecated MPlaceion(pqr,points) call", DeprecationWarning)
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

    logger.info("MPlaceion: pqr=%(pqrName)r", vars(self))
    logger.info("MPlaceion: points=%(pointsName)r", vars(self))
    logger.info("MPlaceion: ion=%(ion)r", vars(self))
      
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
      os.chmod(scriptpath, 0755)
 
    qsubName = "qsub_" + self.jobName + ".bash"
    self.numJobs = self.points.shape[0]  # always record maximum number
    with open(qsubName, "w") as jobFile:
      jobFile.write(self.arrayscript % {
          'jobName': self.jobName,
          'numJobs': self.numJobs,
          'jobArray': "\n".join(bash_jobarray),
          })
    return len(windows)


  def generateMem(self, windows=None, run=False):
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
          file generation and just write a script [``False``]
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
      # OLD-STYLE (still used ...):
      kw = {'conc':self.ionicStrength, 'temperature':self.temperature,'comment':'window %d, z=XXX'%num,
            'basedir': os.path.realpath(self.jobpath(self.get_windowdirname(num))),
            'apbs_script_name': apbs_script_name}
      kw.update(self.schedule)  # set dime and glen !
      # NEW-STYLE (overrides old-style)
      try:
        kw.update(self.bornprofile_kwargs)
      except AttributeError:
        pass
      # using a custom SetupClass pre-populates the parameters for draw_membrane2
      # (hack!! -- should be moved into a cfg input file)
      # Cfg file sets remaining kw args [2010-11-19] but still messy;
      # but no custom classes needed anymore, just membrane.BornAPBSmem)
      self.__cache_MemBornSetup[num] = self.SetupClass(protein, ion, cpx, **kw)
    return self.__cache_MemBornSetup[num]

  def generate(self, windows=None, run=False):
    """Set up all input files for Born calculations with a membrane.

    generate([windows[,run]])

    The optional parameter *windows* allows one to select a subset of windows
    instead of all the points in the sample points file; it can be a single
    number or a list.

    Setting *run* to ``True`` immediately generates setup files, in particular
    it runs :program:`apbs` and :program:`draw_membrane2a` in order to add the
    membrane to the system. Because this can take a long time for a large
    number of windows it is also possible to only generate a bash script for
    each window and defer the setup (*run* = ``False``, the default).
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
