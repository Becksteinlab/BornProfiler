#!/usr/bin/env python
# (c) 2010 Oliver Beckstein
# Licensed under GPL
"""%%prog PQR-file

Runs customized apbs calculation of PQR file in a membrane. Currently
customized for Yohei's CFTR structures (from the 30 ns MD); see the
source for details.

Uses draw_membrane2 to add a low-dielectric (eps=2) region and sets
protein dielectric to eps=10.

.. Note:: Paths to draw_membrane2 and apbs are hard coded!
          apbs = %(APBS)r
          draw_membrane2 = %(DRAWMEMBRANE)r

Commandline version of apbsmem, following
http://www.poissonboltzmann.org/apbs/examples/potentials-of-mean-force/the-polar-solvation-potential-of-mean-force-for-a-helix-in-a-dielectric-slab-membrane
"""

from subprocess import call

import os.path

BASEDIR = os.path.expanduser("~/Projects/Channels/CFTR")

# executables; must be found on PATH by the shell or full path
DRAWMEMBRANE = os.path.expanduser(os.path.join("~", "bin", "draw_membrane2"))
APBS = os.path.expanduser("apbs")


TEMPLATES = {'dummy': """
# APBSmem: prepare the coefficient maps
# (modified by OB)
read
    mol pqr "%(pqr)s"
end

elec name solv0
    mg-dummy
    dime %(DIME_XYZ)s
    glen %(GLEN_XYZ)s
    lpbe
    pdie %(pdie)f
    sdie %(sdie)f
    bcfl zero
    ion 1 %(conc)f 2
    ion -1 %(conc)f 2
    gcent 0 0 0
    mol 1
    chgm spl2
    srfm mol
    srad 1.4
    swin 0.3
    sdens 10
    temp %(temperature)f
    calcenergy no
    calcforce no
    write dielx dx dielx%(suffix)s
    write diely dx diely%(suffix)s
    write dielz dx dielz%(suffix)s
    write kappa dx kappa%(suffix)s
    write charge dx charge%(suffix)s
end

quit
""",
             'solvation': """
# APBSmem: run simple 1st solvation
# (modified by OB)

read
    mol pqr "%(pqr)s"
    # Read Maps
    diel dx dielx%(suffix)sm.dx diely%(suffix)sm.dx dielz%(suffix)sm.dx
    kappa dx kappa%(suffix)sm.dx
    charge dx charge%(suffix)sm.dx
end

elec name solv0
    mg-manual
    dime %(DIME_XYZ)s
    glen %(GLEN_XYZ)s
    lpbe
    pdie %(pdie)f
    sdie %(sdie)f
    bcfl zero
    ion charge  1  conc %(conc)f  radius 0.95  # sodium
    ion charge -1  conc %(conc)f  radius 1.81  # chloride
    gcent 0 0 0
    mol 1
    chgm spl2
    srfm mol
    srad 1.4
    swin 0.3
    sdens 10
    temp %(temperature)f
    calcenergy total
    calcforce no
    write pot dx pot_solv0%(suffix)s
end

elec name ref0
    mg-manual
    dime %(DIME_XYZ)s
    glen %(GLEN_XYZ)s
    lpbe
    pdie %(pdie)f
    sdie %(sdie)f
    bcfl zero
    ion 1 %(conc)f 2
    ion -1 %(conc)f 2
    gcent 0 0 0
    mol 1
    chgm spl2
    srfm mol
    srad 1.4
    swin 0.3
    sdens 10
    temp %(temperature)f
    calcenergy total
    calcforce no
    usemap diel 1
    usemap kappa 1
    usemap charge 1
    write pot dx pot_ref%(suffix)s
end

print elecEnergy ref0 - solv0
end

quit
""",
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

        self._vecnames = ["dime", "glen"]  # see _update_self_names
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

    def write(self, stage, **kwargs):
        vardict = self.get_var_dict(stage)
        vardict.update(self.get_XYZ_dict('dime', self.dime))
        vardict.update(self.get_XYZ_dict('glen', self.glen))
        f = open(self.infile(stage), 'w')
        try:
            f.write(TEMPLATES[stage] % vardict)
        finally:
	    f.close()
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
        rc = call([self.drawmembrane, dielx] + 
                   map(str, [v['zmem'], v['lmem'], v['pdie'], v['Vmem'], v['conc'],
                             v['Rtop'], v['Rbot']]))
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
        

class CFTRmem(APBSmem):
    """APBSmem with custom defaults"""
    def __init__(self, pqr, suffix, zmem=-2, lmem=40, mdie=2, sdie=80, pdie=10, 
                 headgroup_die=20, headgroup_l=0, Vmem=0,
                 Rtop=16, Rbot=10,  **kwargs):
        super(CFTRmem, self).__init__(pqr, suffix, zmem=-2, lmem=40, mdie=2, sdie=80, pdie=10, 
                 headgroup_die=20, headgroup_l=0, Vmem=0,
                 Rtop=16, Rbot=10,  **kwargs)

if __name__ == "__main__":
    import sys
    import errno
    from optparse import OptionParser

    parser = OptionParser(usage=__doc__ % vars())
    parser.add_option("-s", "--suffix", dest="suffix",
                      help="suffix for all generated files [%default]")
    parser.set_defaults(suffix="S")
    
    options, args = parser.parse_args()
    
    try:
        pqr = args[0]
    except IndexError:
        raise ValueError("Need PQR file as input")

    if not os.path.exists(pqr):
        raise IOError(errno.ENOENT, "PQR file not found", pqr)
        
    suffix = 'S'
    
    C = CFTRmem(pqr, options.suffix)
    C.setup()
    C.run_apbs('solvation')
