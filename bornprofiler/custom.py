from bornprofiler.membrane import APBSmem, BornAPBSmem, MemonlyAPBSmem

class _BornAPBSmemCustomizer(BornAPBSmem):
    """Derive and mixin a defaults class to provide the dict defaults."""
    def __init__(self, *args, **kwargs):
        assert hasattr(self, 'defaults') == True, "Derive a class together with a mixin default class."
        kw = self.defaults.copy()
        kw.update(kwargs)
        super(_BornAPBSmemCustomizer, self).__init__(*args, **kw)

class _MemonlyAPBSmemCustomizer(MemonlyAPBSmem):
    """Derive and mixin a defaults class to provide the dict defaults."""
    def __init__(self, *args, **kwargs):
        assert hasattr(self, 'defaults') == True, "Derive a class together with a mixin default class."
        kw = self.defaults.copy()
        kw.update(kwargs)
        super(_MemonlyAPBSmemCustomizer, self).__init__(*args, **kw)

class _APBSmemCustomizer(APBSmem):
    """Derive and mixin a defaults class to provide the dict defaults."""
    def __init__(self, *args, **kwargs):
        kw = self.defaults.copy()
        kw.update(kwargs)
        super(_APBSmemCustomizer, self).__init__(*args, **kw)

#------------------------------------------------------------
# CFTR (Ivetac model)
#------------------------------------------------------------
class CFTRdefaults(object):
    """Mixin class to set defaults for CFTR APBS calcs."""
    defaults = dict(zmem=-2, lmem=40, mdie=2, sdie=80, pdie=10, 
                    headgroup_die=20, headgroup_l=0, Vmem=0,
                    Rtop=16, Rbot=10)

class CFTRmem(_APBSmemCustomizer, CFTRdefaults):
    """APBSmem with custom defaults"""

class CFTRborn(_BornAPBSmemCustomizer, CFTRdefaults):
    """BornAPBSmem with custom defaults"""


#------------------------------------------------------------
# Parsegian test case
#------------------------------------------------------------
class ParsegianDefaults(object):
    """Defaults for the example/parsegian case.

    The membrane is a low dielectric of epsilon=2 with a thickness of
    40 A (in the paper he uses 70 A), centered at (0,0,0).

    .. Note:: The fake protein is required to trick the BornProfiler
       scripts. The charges are 0 and the pdie is the same as that of the
       membrane so it should be electrostatically invisible.
    """
    defaults = dict(zmem=-20, lmem=40, mdie=2, pdie=2, Rtop=0, Rbot=0,
                    headgroup_l=0, Vmem=0,)

class ParsegianSlab(_MemonlyAPBSmemCustomizer, ParsegianDefaults):
    """Set up for an Na ion through a low dielectric slab.

    Thickness is 40 A instead of the 70 A of the paper.

    (See A. Parsegian, Nature 221 (1969), 844)
    """

class ParsegianPore(_MemonlyAPBSmemCustomizer, ParsegianDefaults):
    """Set up for an Na ion through an aqueous pore through a low dielectric slab.

    Pore radius is 5 A, and pore dielectric is 80.

    (See A. Parsegian, Nature 221 (1969), 844)
    """
    def __init__(self, *args, **kwargs):
        kw = self.defaults.copy()
        kw.update({'Rbot': 5, 'Rtop': 5})
        kw.update(kwargs)
        super(ParsegianPore, self).__init__(*args, **kw)
