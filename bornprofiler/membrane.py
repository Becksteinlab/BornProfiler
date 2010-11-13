"""
APBS calculations: Membrane simulations
=======================================

Requires the :program:`draw_membrane2` binary.

Uses :programe:`draw_membrane2` to add a low-dielectric (eps=2) region and sets
protein dielectric to eps=10.

.. Note:: Paths to draw_membrane2 and apbs are hard coded!
          apbs = %(APBS)r
          draw_membrane2 = %(DRAWMEMBRANE)r

Commandline version of apbsmem, following
http://www.poissonboltzmann.org/apbs/examples/potentials-of-mean-force/the-polar-solvation-potential-of-mean-force-for-a-helix-in-a-dielectric-slab-membrane

"""

import os, errno
import warnings
from subprocess import call
import numpy
from itertools import izip

import logging
logger = logging.getLogger("bornprofiler.membrane")

import config
from config import configuration, read_template
# executables; must be found on PATH by the shell or full path in config file
# drawmembrane MUST be draw_membrane2a
DRAWMEMBRANE = configuration["drawmembrane"]
APBS = configuration["apbs"]

TEMPLATES = {'dummy': read_template('dummy.in'),
             'solvation': read_template('solvation.in'),
             'born_dummy': read_template('mdummy.in'),
             'born_run': read_template('mplaceion.in'),
}

class BaseMem(object):
    #drawmembrane = "draw_membrane4"  # from the 1.04 tarball -- segfaults
    drawmembrane = DRAWMEMBRANE  # from the APBS website (no headgroups)
    apbs = APBS
    def __init__(self, *args, **kwargs):
        """draw_membrane and common APBS run  parameters

        .. Note:: sdie=80 and mdie=2 are fixed in draw_membrane2.c at the moment

        :Keywords:
           - zmem : membrane centre (A)
           - lmem : membrane thickness (A)
           - Vmem (untested)
           - pdie : protein dielectric
           - Rtop : exclusion cylinder top
           - Rbot : exclusion cylinder bottom
           - temperature : temperature
           - conc : ionic strength in mol/l  
           - basedir: full path to the directory from which files are read and written; 
             by default this is realpath('.')                     
                        
        """
        self.zmem = kwargs.pop('zmem', 0.0)
        self.lmem = kwargs.pop('lmem', 40.0)
        self.Vmem = kwargs.pop('Vmem', 0.0)  # potential (untested)
        self.mdie = kwargs.pop('mdie', 2.0)  # membrane -- HARDCODED in draw_membrane2.c !!
        self.sdie = kwargs.pop('sdie', 80.0)  # idie (solvent) -- HARDCODED in draw_membrane2.c !!
        #self.sdie = kwargs.pop('sdie', 78.5)  # idie (solvent), Eisenberg and Crothers Phys. Chem. book 1979
        self.pdie = kwargs.pop('pdie', 10.0)  # protein
        self.hdie = kwargs.pop('headgroup_die', 20.0)
        self.lhgp = kwargs.pop('headgroup_l', 0.0)  # geo3
        self.Rtop = kwargs.pop('Rtop', 0)  # geo1
        self.Rbot = kwargs.pop('Rbot', 0)  # geo2
        self.temperature = kwargs.pop('temperature', 298.15)
        self.conc = kwargs.pop('conc', 0.1)   # monovalent salt at 0.1 M
        self.basedir = kwargs.pop('basedir', os.path.realpath(os.path.curdir))

        # check draw_membrane: raises a stink if not the right one
        config.check_drawmembrane(self.drawmembrane)

        # dx file compression
        # http://www.poissonboltzmann.org/apbs/user-guide/running-apbs/input-files/elec-input-file-section/elec-keywords/write
        apbsversion = config.check_APBS(self.apbs)
        if apbsversion >= (1,3):
            self.dxformat = "gz"     # format in an APBS read/write statement
            self.dxsuffix = "dx.gz"  # APBS automatically adds .dx.gz when writing
        else:
            self.dxformat = "dx"
            self.dxsuffix = "dx"

        super(BaseMem, self).__init__(*args, **kwargs)

    def get_var_dict(self, stage):
        d = dict((k,self.__dict__[k.strip()]) for k in self.vars[stage].split(','))
        # modify?
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
        # hardcoded dielx!
        infix = kwargs.pop('infix', self.suffix)
        v.update(kwargs)  # override with kwargs
        cmdline = [self.drawmembrane, infix] + \
            map(str, [v['zmem'], v['lmem'], v['pdie'], v['Vmem'], v['conc'],
                      v['Rtop'], v['Rbot']])
        if self.dxformat == "gz":
            cmdline.append('gz')   # special version draw_membrane2a that can deal with gz
        logger.info("COMMAND: %s", " ".join(cmdline))
        rc = call(cmdline)
        if rc != 0:
            errmsg = "drawmembrane failed [returncode %d]" % rc
            logger.fatal(errmsg)
            raise EnvironmentError(errno.ENODATA, self.drawmembrane, errmsg)
        logger.info("Drawmembrane finished. Look for dx files with 'm' in their name.")
        return rc

    def _gzip_dx(self, name="gzip", options=None):
        from subprocess import call
        from glob import glob
        if options is None:
            options = []
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

    # XXX: probably broken at the moment, check
    # - dx vs dx.gz format
    # - make template monolithic
    # - use separate window dirs ... ie make it more like BornAPBSmem

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
        self.vars = {'dummy': "pqr,suffix,pdie,sdie,conc,temperature,DIME_XYZ,GLEN_XYZ",
                     'solvation': "pqr,suffix,pdie,sdie,conc,temperature,,DIME_XYZ,GLEN_XYZ",
                     'drawmembrane2': "zmem,lmem,pdie,Vmem,conc,Rtop,Rbot,suffix"}
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
        self.write('solvation')
        if gz:
            self.gzip_dx()        

class BornAPBSmem(BaseMem):
    """Class to prepare a single window in a manual focusing run."""    

    #: Suffices of file names are hard coded in templates and should not
    #: be changed; see ``templates/mdummy.in`` and ``templates/mplaceion.in``.
    #: The order of the suffices corresponds to the sequence in the schedule.
    suffices = ('L','M','S')

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
                          'born_run': kwargs.pop('apbs_script_name', 'mem_placeion.in'),
                          }

        #: "static" variables required for generating a file from a template;
        #: The variable names can be found in self.__dict__
        self.vars = {'born_dummy': 
                     "protein_pqr,ion_pqr,complex_pqr,"
                     "pdie,sdie,conc,temperature,dxformat,dxsuffix,"
                     "DIME_XYZ_L,DIME_XYZ_M,DIME_XYZ_S,"
                     "GLEN_XYZ_L,GLEN_XYZ_M,GLEN_XYZ_S",
                     'born_run': 
                     "protein_pqr,ion_pqr,complex_pqr,"
                     "pdie,sdie,conc,temperature,comment,dxformat,dxsuffix,"
                     "DIME_XYZ_L,DIME_XYZ_M,DIME_XYZ_S,"
                     "GLEN_XYZ_L,GLEN_XYZ_M,GLEN_XYZ_S",
                     'born_jobscript': "jobname,infile,outfile",
                     'drawmembrane2': 
                     "zmem,lmem,pdie,Vmem,conc,Rtop,Rbot",}

        # generate names
        d = {}
        for suffix,dime,glen in izip(self.suffices, self.dime, self.glen):
            d['DIME_XYZ_%s' % suffix] = self.vec2str(dime)
            d['GLEN_XYZ_%s' % suffix] = self.vec2str(glen)
        self.__dict__.update(d)

        # process the drawmembrane parameters
        super(BornAPBSmem, self).__init__(*args[3:], **kwargs)
        
    def generate(self):
        """Setup solvation calculation.

        1. create exclusion maps (runs apbs)
        2. create membrane maps (drawmembrane)
        3. create apbs run input file

        """
        self.write('born_dummy')
        self.run_apbs('born_dummy')
        for infix in self.infices:
            self.run_drawmembrane(infix=infix)
        self.write('born_run')
        if self.dxformat == "dx":
            logger.info("Manually compressing all dx files (you should get APBS >= 1.3...)")
            self.gzip_dx()

