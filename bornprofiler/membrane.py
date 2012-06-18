# BornProfiler -- integration with draw_membrane
# Copyright (c) 2010-2011 Oliver Beckstein
# Published under the GNU Public Licence, version 3

"""
APBS calculations: Membrane simulations
=======================================

Requires the :program:`draw_membrane2a` binary.

Uses :programe:`draw_membrane2a` to add a low-dielectric (eps=2) region and sets
protein dielectric to eps=10.

.. Note:: Paths to draw_membrane2a and apbs are set in the configuration file
   ``~/.bornprofiler.cfg``::
          apbs = %(APBS)r
          draw_membrane2 = %(DRAWMEMBRANE)r

Commandline version of apbsmem, following
http://www.poissonboltzmann.org/apbs/examples/potentials-of-mean-force/the-polar-solvation-potential-of-mean-force-for-a-helix-in-a-dielectric-slab-membrane
and modified (see ``src/drawmembrane/draw_membrane2a.c`` in the BornProfiler distribution).

"""
from __future__ import with_statement

import os, errno
import warnings
from subprocess import call
import numpy
from itertools import izip

import logging
logger = logging.getLogger("bornprofiler.membrane")

from io import read_template

import config
from config import configuration
# executables; must be found on PATH by the shell or full path in config file
# drawmembrane MUST be draw_membrane2a
DRAWMEMBRANE = configuration["drawmembrane"]
APBS = configuration["apbs"]

#: cache for template files
TEMPLATES = {'dummy': read_template('dummy.in'),
             'solvation': read_template('solvation.in'),
             'born_dummy': read_template('mdummy.in'),
             'with_protein': read_template('mplaceion.in'),
             'memonly': read_template('mplaceion_memonly.in'),
             'born_setup_script': read_template('drawmembrane2.bash'),
}

class BaseMem(object):
    drawmembrane = DRAWMEMBRANE
    apbs = APBS
    def __init__(self, *args, **kwargs):
        """draw_membrane and common APBS run  parameters

        :Keywords:
           - zmem : membrane centre (A)
           - lmem : membrane thickness (A)
           - Vmem : cytosolic membrane potential in kT/e (untested)
           - pdie : protein dielectric
           - sdie: solvent dielectric
           - mdie: membrane dielectric
           - Rtop : exclusion cylinder top
           - Rbot : exclusion cylinder bottom
           - x0_R : X exclusion center 
           - y0_R : Y exclusion center 
           - cdie : channel solvent dielectric (by default same as sdie)
           - headgroup_l : thickness of headgroup region (A)
           - headgroup_die : headgroup dielectric
           - temperature : temperature
           - conc : ionic strength in mol/l
           - basedir: full path to the directory from which files are read and written;
             by default this is realpath('.')

        """
        self.versions = {}
        # check draw_membrane: raises a stink if not the right one
        self.versions['drawmembrane'] = config.check_drawmembrane(self.drawmembrane)

        # dx file compression
        # http://www.poissonboltzmann.org/apbs/user-guide/running-apbs/input-files/elec-input-file-section/elec-keywords/write
        self.versions['APBS'] = config.check_APBS(self.apbs)
        self.unpack_dxgz = False     # hack for v1.3, see below and http://sourceforge.net/support/tracker.php?aid=3108761
        if self.versions['APBS'] >= (1,3):
            self.dxformat = "gz"     # format in an APBS read/write statement
            self.dxsuffix = "dx.gz"  # APBS automatically adds .dx.gz when writing
            if self.versions['APBS'] == (1,3) \
                    and not config.cfg.getboolean('executables', 'apbs_always_read_dxgz'):
                # note that 1.3 'READ gz' is broken (at least on Mac OS X) so we hack around
                # this by gunzipping in the queuing script and replacing the gz format in readstatements
                # with the dx one :-p
                # svn 1623+ are fixed and work so one can set the hack in the config file:
                # [executables]
                # apbs_always_read_dxgz = True
                self.unpack_dxgz = True
        else:
            self.dxformat = "dx"
            self.dxsuffix = "dx"
        self.apbs_version = ".".join(map(str, self.versions['APBS']))

        # all the variables needed for draw_membrane2a
        # (TODO: having to add ALL these variables as attributes sucks; this
        # should be done automagically... makes it hard adding new parameters)
        self.zmem = kwargs.pop('zmem', 0.0)
        self.lmem = kwargs.pop('lmem', 40.0)
        self.Vmem = kwargs.pop('Vmem', 0.0)  # potential (untested)
        self.mdie = kwargs.pop('mdie', 2.0)  # membrane
        self.sdie = kwargs.pop('sdie', 80.0)  # idie (solvent)
        #self.sdie = kwargs.pop('sdie', 78.5)  # idie (solvent), Eisenberg and Crothers Phys. Chem. book 1979
        self.pdie = kwargs.pop('pdie', 10.0)  # protein
        self.headgroup_die = kwargs.pop('headgroup_die', 20.0)
        self.headgroup_l = kwargs.pop('headgroup_l', 0.0)  # geo3
        self.Rtop = kwargs.pop('Rtop', 0)  # geo1
        self.Rbot = kwargs.pop('Rbot', 0)  # geo2
        self.x0_R = kwargs.pop('x0_R', 0)  # new in draw_membrane2a.c 04/26/11
        self.y0_R = kwargs.pop('y0_R', 0)  # new in draw_membrane2a.c 04/26/11
        self.cdie = kwargs.pop('cdie', self.sdie)
        self.temperature = kwargs.pop('temperature', 298.15)
        self.conc = kwargs.pop('conc', 0.1)   # monovalent salt at 0.1 M
        self.basedir = kwargs.pop('basedir', os.path.realpath(os.path.curdir))

        #logger.debug("BornProfiler: detected APBS version %(apbs_version)r", vars(self))

        super(BaseMem, self).__init__(*args, **kwargs)

    def get_var_dict(self, stage):
        """Load required values for stage from vars(self) into d."""
        keys = [k for k in self.vars[stage].split(',') if k.strip() in self.__dict__]
        d = dict((k,self.__dict__[k.strip()]) for k in keys)
        return d

    def vec2str(self, vec):
        return " ".join(map(str, vec))

    def write(self, stage, **kwargs):
        """General template writer.

        *stage* is a key into :data:`TEMPLATES` and :attr:`filenames`.
        """
        vardict = self.get_var_dict(stage)
        vardict.update(kwargs)
        with open(self.infile(stage), 'w') as f:
            f.write(TEMPLATES[stage] % vardict)
        return self.infile(stage)

    def outfile(self, stage):
        return "%s%s.out" % (stage, self.suffix)

    def infile(self, stage):
        return self.filenames[stage]

    def run_apbs(self, stage, **kwargs):
        outfile = self.outfile(stage)
        infile = self.infile(stage)
        logger.info("APBS starting: %s --output-file=%s %s", self.apbs, outfile, infile)
        rc = call([self.apbs, '--output-file=%s' % outfile, infile])
        if rc != 0:
            errmsg = "ABPS failed [returncode %d]. Investigate error messages in output and io.mc" % rc
            logger.fatal(errmsg)
            raise EnvironmentError(errno.ENODATA, self.apbs, errmsg)
        logger.info("APBS completed run successfully. The peasants rejoice.")
        return rc

    def run_drawmembrane(self, **kwargs):
        stage = kwargs.pop('stage', "drawmembrane2")
        v = self.get_var_dict(stage)
        # hardcoded names dielx<infix>.dx etc!
        infix = kwargs.pop('infix', self.suffix)
        v.update(kwargs)  # override with kwargs
        cmdline = [self.drawmembrane] + \
            map(str, ["-z", v['zmem'], "-d", v['lmem'], "-m", v['mdie'],   # membrane
                      "-a", v['headgroup_l'], "-i", v['headgroup_die'],    # headgrops
                      "-s", v["sdie"], "-p", v['pdie'],                    # environment
                      "-V", v['Vmem'], "-I", v['conc'],
                      "-R", v['Rtop'], "-r", v['Rbot'], "-c", v['cdie'],   # channel exclusion
                      "-X", v['x0_R'], "-Y", v['y0_R'],
                      ])
        if self.dxformat == "gz":
            cmdline.append('-Z')   # special version draw_membrane2a that can deal with gz
        cmdline.append(infix)
        logger.info("COMMAND: %s", " ".join(cmdline))
        rc = call(cmdline)
        if rc != 0:
            errmsg = "drawmembrane failed [returncode %d]" % rc
            logger.fatal(errmsg)
            raise EnvironmentError(errno.ENODATA, self.drawmembrane, errmsg)
        logger.info("Drawmembrane finished. Look for dx files with 'm' in their name.")
        return rc

    def write_infile(self, name):
        if self.unpack_dxgz:
            extra = {'dxformat': 'dx', 'dxsuffix': 'dx'}
            logger.debug("Working around bug in APBS 1.3: reading gunzipped files")
        else:
            extra = {}
        self.write(name, **extra)

    def _gzip_dx(self, name="gzip", options=None):
        from subprocess import call
        from glob import glob
        if options is None:
            options = []
        if name == "gunzip" or '-d' in options:
            dxfiles = glob("*.dx.gz")
        else:
            dxfiles = glob("*.dx")
        if len(dxfiles) == 0:
            logger.warn("No dx files for %(name)s.", vars())
            return 0
        cmdline = [name] + options + dxfiles
        logger.info("Beginning to %s %s all %d dx files... patience.",
                    name, " ".join(options), len(dxfiles))
        logger.debug("COMMAND: %r", cmdline)
        rc = call(cmdline)
        if rc == 0:
            logger.info("Completed %s operation.", name)
        else:
            logger.error("%s error: returncode %d", name, rc)
        return rc

    def gzip_dx(self):
        """Run external gzip on all dx files.

        Saves about 98% of space.
        """
        return self._gzip_dx()

    def ungzip_dx(self):
        """Run external ungzip on all dx files."""
        return self._gzip_dx('gunzip')


class APBSmem(BaseMem):
    """Represent the apbsmem tools.

    Run for S, M, and L grid (change suffix).

    .. Note:: see code for kwargs
    """
    # used by apbs-mem-potential.py

    def __init__(self, *args, **kwargs):
        """Set up calculation.

        APBSmem(pqr, suffix[, kwargs])

        :Arguments:
          - arguments for drawmembrane (see source)
          - temperature: temperature [298.15]
          - conc: ionic concentration of monovalent salt with radius 2 A
                  in mol/l [0.1]
          - dime: grid dimensions, as list [(97,97,97)]
          - glen: grid length, as list in  Angstroem [(250,250,250}]
        """
        self.pqr = args[0]
        self.suffix = args[1]
        self.dime = kwargs.pop('dime', (97,97,97))     # grid points -- check allowed values!!
        self.glen = kwargs.pop('glen', (250,250,250))  # Angstroem
        self.filenames = {'dummy': 'dummy%(suffix)s.in' % vars(self),
                          'solvation': 'solvation%(suffix)s.in' % vars(self),
                          }
        # names for diel, kappa, and charge maps hard coded; only suffix (=infix) varies
        # (see templates/dummy.in)

        #: "static" variables required for a  calculation
        self.vars = {'dummy':
                     "pqr,suffix,"
                     "pdie,sdie,conc,temperature,dxformat,dxsuffix,"
                     "DIME_XYZ,GLEN_XYZ",
                     'solvation':
                     "pqr,suffix,pdie,sdie,conc,temperature,dxformat,dxsuffix,"
                     "DIME_XYZ,GLEN_XYZ",
                     'drawmembrane2':
                     "zmem,lmem,pdie,sdie,mdie,Vmem,conc,"
                     "Rtop,Rbot,x0_R,y0_R,cdie,"
                     "headgroup_l,headgroup_die,suffix",
                     }
        # generate names
        d = {}
        d['DIME_XYZ'] = self.vec2str(self.dime)
        d['GLEN_XYZ'] = self.vec2str(self.glen)
        self.__dict__.update(d)

        # process the drawmembrane parameters
        super(APBSmem, self).__init__(*args[2:], **kwargs)

    def generate(self):
        """Setup solvation calculation.

        1. create exclusion maps (runs apbs)
        2. create membrane maps (drawmembrane)
        3. create apbs run input file

        """
        self.write('dummy')
        self.run_apbs('dummy')
        self.run_drawmembrane()
        self.write_infile('solvation')
        if self.dxformat == "dx":
            logger.info("Manually compressing all dx files (you should get APBS >= 1.3...)")
            self.gzip_dx()

class BornAPBSmem(BaseMem):
    """Class to prepare a single window in a manual focusing run."""

    #: Suffices of file names are hard coded in templates and should not
    #: be changed; see ``templates/mdummy.in`` and ``templates/mplaceion.in``.
    #: The order of the suffices corresponds to the sequence in the schedule.
    suffices = ('L','M','S')
    default_runtype = "with_protein"

    def __init__(self, *args, **kwargs):
        """Set up calculation.

        BornAPBSmem(protein_pqr, ion_pqr, complex_pqr[, kwargs])

        Additional arguments:
          - apbs_script_name: name on the xxx.in script [mem_placeion.in]
          - drawmembrane arguments (see source)
          - temperature: temperature [298.15]
          - conc: ionic concentration of monovalent salt with radius 2 A
                  in mol/l [0.1]
          - dime: grid dimensions, as list with 1 or 3 entries; if o1 entry
                  then the same dime is used for all focusing stages [(97,97,97)]
          - glen: grid length, as list in  Angstroem, must contain three triplets
                  [(250,250,250),(100,100,100),(50,50,50)]
          - runtype: for standard Born calculations use 'with_protein'; if one only wants to
                  look at geometrically defined dielectric regions and uncharged proteins
                  (see the ``example/Parsegian``) then use 'memonly' ['with_protein']
        """
        self.protein_pqr = args[0]
        self.ion_pqr = args[1]
        self.complex_pqr = args[2]
        # names for diel, kappa, and charge maps hard coded; only infix varies
        # (see templates/mdummy.in and templates/mplaceion.in)
        # processing requires draw_membrane2a with changes to the filename handling code
        self.infices = ('_prot_L','_prot_M','_prot_S', '_cpx_L','_cpx_M','_cpx_S')
        self.suffix = ""  # hack to keep parent class happy :-p[
        self.comment = kwargs.pop('comment', "BORNPROFILE")
        self.runtype = kwargs.pop('runtype', self.default_runtype)

        dime = numpy.array(kwargs.pop('dime', (97,97,97)))     # grid points -- check allowed values!!
        if not (dime.shape == (3,) or dime.shape == (3,3)):
            raise TypeError("dime must be a triplet of values or a triplet of such triplets")
        if len(dime.shape) == 1:
            self.dime = numpy.concatenate([dime, dime, dime]).reshape(3,3)
        else:
            self.dime = dime
        glen = numpy.array(kwargs.pop('glen', [(250,250,250),(100,100,100),(50,50,50)]))
        if not glen.shape == (3,3):
            raise TypeError("glen must be a triplet of triplets.")
        self.glen = glen  # Angstroem

        # names for diel, kappa, and charge maps hard coded; only infix varies
        # TODO: proper naming and/or directories
        self.filenames = {'born_dummy': 'mem_dummy.in',
                          self.runtype: kwargs.pop('apbs_script_name', 'mem_placeion.in'),
                          'born_setup_script': kwargs.pop('dummy_script_name', 'run_drawmembrane.bash'),
                          }

        #: "static" variables required for generating a file from a template;
        #: The variable names can be found in self.__dict__
        self.vars = {'born_dummy':
                     "protein_pqr,ion_pqr,complex_pqr,"
                     "pdie,sdie,conc,temperature,dxformat,dxsuffix,"
                     "DIME_XYZ_L,DIME_XYZ_M,DIME_XYZ_S,"
                     "GLEN_XYZ_L,GLEN_XYZ_M,GLEN_XYZ_S",
                     self.runtype:
                     "protein_pqr,ion_pqr,complex_pqr,"
                     "pdie,sdie,conc,temperature,comment,dxformat,dxsuffix,"
                     "DIME_XYZ_L,DIME_XYZ_M,DIME_XYZ_S,"
                     "GLEN_XYZ_L,GLEN_XYZ_M,GLEN_XYZ_S",
                     'drawmembrane2':
                     "zmem,lmem,pdie,sdie,mdie,Vmem,conc,"
                     "Rtop,Rbot,x0_R,y0_R,cdie,"
                     "headgroup_l,headgroup_die,suffix",
                     }
        self.vars['born_setup_script'] = \
            ",".join([self.vars['drawmembrane2'],
                      "dxformat","infices","born_dummy_in","born_dummy_out"])
        # generate names
        d = {}
        for suffix,dime,glen in izip(self.suffices, self.dime, self.glen):
            d['DIME_XYZ_%s' % suffix] = self.vec2str(dime)
            d['GLEN_XYZ_%s' % suffix] = self.vec2str(glen)
        self.__dict__.update(d)

        # process the drawmembrane parameters
        super(BornAPBSmem, self).__init__(*args[3:], **kwargs)

    def generate(self, run=False):
        """Setup solvation calculation.

        If *run* = ``True`` then runs :program:`apbs` and
        :program:`draw_membrane2` (which can take a few minutes); otherwise
        just generate scripts (default).

        *run* = ``True``
          1. create exclusion maps (runs :program:`apbs` through :meth:`run_apbs`)
          2. create membrane maps (:program:`draw_membrane2a` through
             :meth:`run_drawmembrane2`)
          3. create apbs run input file

        *run* = ``False``
          1. write a bash script that calls :program:`apbs` and :program:`draw_membrane2`
             and which can be integrated into a window run script for parallelization.
          2. create apbs run input file
        """
        if run:
            return self._generate_locally()
        else:
            return self._generate_scripted()

    def _generate_locally(self):
        self.write('born_dummy')
        self.run_apbs('born_dummy')
        for infix in self.infices:
            self.run_drawmembrane(infix=infix)
        self.write_infile(self.runtype)
        if self.dxformat == "dx":
            logger.info("Manually compressing all dx files (you should get APBS >= 1.3...)")
            self.gzip_dx()
        return None

    def _generate_scripted(self):
        """Write input files and bash script to run drawmembrane and APBS.

        :Returns: name of the bash script
        """
        self.write('born_dummy')
        # hack: MUST provide filenames as kwargs and infices as a "bash" list (i.e. just spaces);
        # produces script == self.infile('born_setup_script')
        script = self.write('born_setup_script', infices=self.vec2str(self.infices),
                            born_dummy_in=self.infile('born_dummy'),
                            born_dummy_out=self.outfile('born_dummy'))
        self.write_infile(self.runtype)

        os.chmod(script, 0755)
        return script

