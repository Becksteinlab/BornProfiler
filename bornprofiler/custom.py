from bornprofiler.membrane import APBSmem

class CFTRmem(APBSmem):
    """APBSmem with custom defaults"""
    def __init__(self, pqr, suffix, zmem=-2, lmem=40, mdie=2, sdie=80, pdie=10, 
                 headgroup_die=20, headgroup_l=0, Vmem=0,
                 Rtop=16, Rbot=10,  **kwargs):
        super(CFTRmem, self).__init__(pqr, suffix, zmem=-2, lmem=40, mdie=2, sdie=80, pdie=10, 
                 headgroup_die=20, headgroup_l=0, Vmem=0,
                 Rtop=16, Rbot=10,  **kwargs)
