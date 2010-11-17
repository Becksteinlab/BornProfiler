============================
 README Parsegian Example
============================

Classic paper of ion through a membrane. ~1969

Fake-protein: square of uncharged H-atoms with length 40 A.

../../scripts/samplepoints.pl --com=0 --com=0 --lengthunit=Angstrom path.txt > path.dat


 M = bornprofiler.core.MPlaceion("fakeprotein.pqr", "path.dat", jobName="ParsegianSlab",memclass=bornprofiler.custom.ParsegianSlab, script='q_SBCB.sh')


* "protein" gives NAN because there are no charges in the system (and
  hence the energy is 0)

data for ion at centre of the low-dielectric slab::

 grep Total job_0010.out 
  Total electrostatic energy = NAN kJ/mol
  Total electrostatic energy = NAN kJ/mol
  Total electrostatic energy = NAN kJ/mol
  Total electrostatic energy = 6.984003508929E+01 kJ/mol
  Total electrostatic energy = 3.964936713784E+02 kJ/mol
  Total electrostatic energy = 0.000000000000E+00 kJ/mol
  Total electrostatic energy = 2.257804856385E+02 kJ/mol
  Total electrostatic energy = 5.836501527120E+02 kJ/mol
  Total electrostatic energy = 1.178935989250E+03 kJ/mol

::
   W = [6.984003508929E+01,
   3.964936713784E+02,
   0.000000000000E+00,
   2.257804856385E+02,
   5.836501527120E+02,
   1.178935989250E+03,
   ]
