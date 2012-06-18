============================
 README Parsegian Example
============================

Classic paper of ion through a membrane [A. Parsegian, Nature 221
(1969), 844].

Fake-protein: square of uncharged H-atoms with length 40 A.

Generate a straight path::
  ../../scripts/samplepoints.pl --com=0 --com=0 --lengthunit=Angstrom path.txt > path.dat

All input parameters are set in the input files **ParsegianSlab.cfg**
and  **ParsegianPore.cfg**

.. Note:: The "fake protein" is required to trick the BornProfiler
   scripts into accepting a system that only contains dielectric
   regions. The charges are 0 and the pdie is the same as that of the
   membrane so it should be electrostatically invisible.


Dielectric slab
===============

Simple dielectric slab.

The membrane is a low dielectric of epsilon=2 with a thickness of 40 A
(in the paper he uses 70 A), centered at (0,0,0)::

  apbs-bornprofile-mplaceion.py ParsegianSlab.cfg

Then run locally (or submit through a queuing system). On my 8-core
Mac (you need to figure out how to run all the individual window
calculations on your machine; maybe you can substitute the
:program:`seq` command for :program:`jot`)::

   ../../scripts/parallel.py 8 ../../scripts/fake_qsub ./qsub_ParsegianSlab.bash --- `jot 19`

and analyze when done::

   apbs-bornprofile-analyze.py --name=ParsegianSlab path.dat ParsegianSlab/w00*/job*.out


Aqueous pore
============

Set up for an Na ion through an aqueous pore through a low dielectric
slab. Pore radius is 5 A, and pore dielectric is 80::

  apbs-bornprofile-mplaceion.py ParsegianPore.cfg

and run (this is specific for my Mac --- you need to figure out how to
run all the individual window calculations on your machine; maybe you
can substitute the :program:`seq` command for :program:`jot`)::

  ../../scripts/parallel.py 8 ../../scripts/fake_qsub ./qsub_ParsegianPore.bash --- `jot 19`  

Analyze::
  apbs-bornprofile-analyze.py --name=ParsegianPore path.dat ParsegianPore/w00*/job*.out



