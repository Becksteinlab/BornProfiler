============================
 README Parsegian Example
============================

Classic paper of ion through a membrane. ~1969

Fake-protein: square of uncharged H-atoms with length 40 A.

Generate a straight path::
  ../../scripts/samplepoints.pl --com=0 --com=0 --lengthunit=Angstrom path.txt > path.dat


Dielectric slab
===============

In ipython::

  M = bornprofiler.core.MPlaceion("fakeprotein.pqr", "path.dat", jobName="ParsegianSlab",memclass=bornprofiler.custom.ParsegianSlab, script='q_SBCB.sh')
  M.generate(run=False)

Then run locally (or submit through a queuing system). On my 8-core
Mac::

   ../../scripts/parallel.py 8 ../../scripts/fake_qsub ./qsub_ParsegianSlab.bash --- `jot 19`

and analyze when done::

   apbs-bornprofile-analyze.py path.dat ParsegianSlab/w00*/job*.out


Aqueous pore
============

Pore of radius 5 A through the slab::

   P = bornprofiler.core.MPlaceion("fakeprotein.pqr", "path.dat", jobName="ParsegianPore",memclass=bornprofiler.custom.ParsegianPore, script='q_SBCB.sh')
   P.generate(run=False)

and run::

  ../../scripts/parallel.py 8 ../../scripts/fake_qsub ./qsub_ParsegianPore.bash --- `jot 19`  



