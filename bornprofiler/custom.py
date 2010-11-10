from bornprofiler.membrane import APBSmem, BornAPBSmem

class CFTRdefaults(object):
    """Mixin class to set defaults for CFTR APBS calcs."""
    defaults = dict(zmem=-2, lmem=40, mdie=2, sdie=80, pdie=10, 
                    headgroup_die=20, headgroup_l=0, Vmem=0,
                    Rtop=16, Rbot=10)

class CFTRmem(APBSmem, CFTRdefaults):
    """APBSmem with custom defaults"""
    def __init__(self, *args, **kwargs):
        kw = self.defaults.copy()
        kw.update(kwargs)
        super(CFTRmem, self).__init__(*args, **kw)


class CFTRborn(BornAPBSmem, CFTRdefaults):
    """BornAPBSmem with custom defaults"""
    def __init__(self, *args, **kwargs):
        kw = self.defaults.copy()
        kw.update(kwargs)
        super(CFTRborn, self).__init__(*args, **kw)


