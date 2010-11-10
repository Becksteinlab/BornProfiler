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
from subprocess import call
import numpy
from itertools import izip

import logging
logger = logging.getLogger("bornprofile.membrane")

from config import configuration, read_template
# executables; must be found on PATH by the shell or full path in config file
DRAWMEMBRANE = configuration["drawmembrane"]
APBS = configuration["apbs"]

TEMPLATES = {'dummy': read_template('dummy.in'),
             'solvation': read_template('solvation.in'),
             'born_dummy': read_template('mdummy.in'),
             'born_run': read_template('mplaceion.in'),
}

class APBSmem(object):
    """Represent the apbsmem tools.

    Run for S, M, and L grid (change suffix).

    .. Note:: see code for kwargs
    """
    #drawmembrane = "draw_membrane4"  # from the 1.04 tarball -- segfaults
    drawmembrane = DRAWMEMBRANE  # from the APBS website (no headgroups)
    apbs = APBS
    def __init__(self, pqr, suffix, zmem=0, lmem=40, mdie=2, sdie=80, pdie=10, 
                 headgroup_die=20, headgroup_l=0, Vmem=0,
                 Rtop=0, Rbot=0,  **kwargs):
        """Set up calculation.

        Additional arguments:
          - temp: temperature [298.15]
          - conc: ionic concentration of monovalent salt with radius 2 A
                  in mol/l [0.1]
          - dime: grid dimensions, as list [(97,97,97)]
          - glen: grid length, as list in  Angstroem [(250,250,250}]
        """
        if len(suffix) != 1:
            errmsg = "suffix can only be a single character due to limitations in %r" % self.drawmembrane
            logger.fatal(errmsg)            
            raise ValueError(errmsg)
        self.pqr = pqr
        self.suffix = suffix
        self.zmem = zmem
        self.lmem = lmem
        self.Vmem = Vmem  # potential (untested)
        self.mdie = mdie  # membrane
        self.sdie = sdie  # idie (solvent)
        self.pdie = pdie  # protein
        self.hdie = headgroup_die
        self.lhgp = headgroup_l  # geo3
        self.Rtop = Rtop  # geo1
        self.Rbot = Rbot  # geo2

        self.temperature = kwargs.pop('temp', 298.15)
        self.conc = kwargs.pop('conc', 0.1)   # monovalent salt at 0.1 M
        self.dime = kwargs.pop('dime', (97,97,97))     # grid points -- check allowed values!!
        self.glen = kwargs.pop('glen', (250,250,250))  # Angstroem
        self.filenames = {'dummy': 'dummy%(suffix)s.in' % vars(self),
                          'solvation': 'solvation%(suffix)s.in' % vars(self),
                          }
        # names for diel, kappa, and charge maps hard coded; only suffix varies

        #: "static" variables required for a  calculation
        self.vars = {'dummy': "pqr,suffix,pdie,sdie,conc,temperature",
                     'solvation': "pqr,suffix,pdie,sdie,conc,temperature",
                     'drawmembrane2': "zmem,lmem,pdie,Vmem,conc,Rtop,Rbot,suffix"}

    def get_var_dict(self, stage):
        d = dict((k,self.__dict__[k.strip()]) for k in self.vars[stage].split(','))
        # modify?
        return d

    def get_XYZ_dict(self, name, vec):
        return {name.upper()+'_XYZ': " ".join(map(str, vec))}

    def get_XYZ(self, attr):
        """Join items of *attr* with spaces."""
        return " ".join(map(str, self.__getattibute__(attr)))

    def vec2str(self, vec):
        return " ".join(map(str, vec))

    def write(self, stage, **kwargs):
        vardict = self.get_var_dict(stage)
        vardict.update(self.get_XYZ_dict('dime', self.dime))
        vardict.update(self.get_XYZ_dict('glen', self.glen))
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
        rc = call([self.apbs, '--output-file=%s' % outfile, infile])
        return rc

    def run_drawmembrane(self, **kwargs):
        stage = kwargs.pop('stage', "drawmembrane2")
        v = self.get_var_dict(stage)
        # hardcoded dielx!
        dielx = kwargs.pop('dielx', "dielx%(suffix)s.dx" % vars(self))
        v.update(kwargs)  # override with kwargs
        if not os.path.exists(dielx):
            errmsg = "File %(dielx)r missing." % vars()
            logger.fatal(errmsg)
            raise OSError(errno.ENOENT, dielx, errmsg)
        logger.info("Running drawmembrane...")
        rc = call([self.drawmembrane, dielx] + 
                   map(str, [v['zmem'], v['lmem'], v['pdie'], v['Vmem'], v['conc'],
                             v['Rtop'], v['Rbot']]))
        if rc != 0:
            errmsg = "%r failed with return code %d." % (self.drawmembrane, rc)
            logger.fatal(errmsg)
            raise EnvironmentError(rc, "execution failed", self.drawmembrane)
        logger.info("Drawmembrane finished. Look for dx files with 'm' in their name.")
        return rc

    def setup(self, **kwargs):
        """Setup solvation calculation.

        1. create exclusion maps (runs apbs)
        2. create membrane maps (drawmembrane)
        3. create apbs run input file

        """
        self.write('dummy')
        self.run_apbs('dummy')
        self.run_drawmembrane()
        self.write('solvation')
        

class BornAPBSmem(APBSmem):
    """Class to prepare a single window in a manual focusing run."""    
    def __init__(self, protein_pqr, ion_pqr, complex_pqr,
                 zmem=0, lmem=40, mdie=2, sdie=80, pdie=10, 
                 headgroup_die=20, headgroup_l=0, Vmem=0,
                 Rtop=0, Rbot=0,  **kwargs):
        """Set up calculation.

        Additional arguments:
          - temp: temperature [298.15]
          - conc: ionic concentration of monovalent salt with radius 2 A
                  in mol/l [0.1]
          - dime: grid dimensions, as list with 1 or 3 entries; if o1 entry
                  then the same dime is used for all focusing stages [(97,97,97)]
          - glen: grid length, as list in  Angstroem, must contain three triplets
                  [(250,250,250),(100,100,100),(50,50,50)]
        """
        self.protein_pqr = protein_pqr
        self.ion_pqr = ion_pqr
        self.complex_pqr = complex_pqr
        self.suffices = ('L','M','S')   # hard coded in templates
        self.suffix = ""  # hack to keep parent class happy :-p[
        self.comment = kwargs.pop('comment', "BORNPROFILE")
        self.zmem = zmem
        self.lmem = lmem
        self.Vmem = Vmem  # potential (untested)
        self.mdie = mdie  # membrane
        self.sdie = sdie  # idie (solvent)
        self.pdie = pdie  # protein
        self.hdie = headgroup_die
        self.lhgp = headgroup_l  # geo3
        self.Rtop = Rtop  # geo1
        self.Rbot = Rbot  # geo2

        self.temperature = kwargs.pop('temp', 298.15)
        self.conc = kwargs.pop('conc', 0.1)   # monovalent salt at 0.1 M
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

        # names for diel, kappa, and charge maps hard coded; only suffix varies
        # TODO: proper naming and/or directories
        self.filenames = {'born_dummy': 'mdummy.in' % vars(self),
                          'born_run': 'mplaceion.in' % vars(self),
                          }

        #: "static" variables required for a  calculation
        self.vars = {'born_dummy': "protein_pqr,ion_pqr,complex_pqr,"
                     "pdie,sdie,conc,temperature,"
                     "DIME_XYZ_L,DIME_XYZ_M,DIME_XYZ_S,"
                     "GLEN_XYZ_L,GLEN_XYZ_M,GLEN_XYZ_S",
                     'born_run': "protein_pqr,ion_pqr,complex_pqr,"
                     "pdie,sdie,conc,temperature,comment,"
                     "DIME_XYZ_L,DIME_XYZ_M,DIME_XYZ_S,"
                     "GLEN_XYZ_L,GLEN_XYZ_M,GLEN_XYZ_S",
                     'drawmembrane2': "zmem,lmem,pdie,Vmem,conc,Rtop,Rbot"}

        # generate names
        d = {}
        for suffix,dime,glen in izip(self.suffices, self.dime, self.glen):
            d['DIME_XYZ_%s' % suffix] = self.vec2str(dime)
            d['GLEN_XYZ_%s' % suffix] = self.vec2str(glen)
        self.__dict__.update(d)
        

    def setup(self, **kwargs):
        """Setup solvation calculation.

        1. create exclusion maps (runs apbs)
        2. create membrane maps (drawmembrane)
        3. create apbs run input file

        """
        self.write('born_dummy')
        self.run_apbs('born_dummy')
        for dielx in ('dielxL.dx','dielxM.dx','dielxS.dx',
                      'dieCxL.dx','dieCxM.dx','dieCxS.dx'):
            self.run_drawmembrane(dielx=dielx)
        self.write('born_run')
