==================================================
 Example: nicotinic acetylcholine receptor (2BG9)
==================================================

Generation of input files
=========================

From the example directory, prepare the input::

  apbs-bornprofile-placeion.py Na_nAChR.cfg
  apbs-bornprofile-placeion.py Cl_nAChR.cfg

These commands will generate input files in two directories
corresponding to the *job:name* values in the run configuration (cfg)
files. 


Running APBS calculations
=========================

You can run the individual "windows" (=ion positions) by executing the
``job_NNNN.sh`` script in each sub-directory (``<job:name>/wNNNN``)
or, if your Gridengine queuing system is set up accordingly, by
submitting a GE array job::

  qsub qsub_<job:name>.sge

You can also provide your own template for the qeuing system
submission script by changing the value of *job:arrayscript* to your
script. Your script can be placed in your personal library of
BornProfiler files in ``~/.bornprofiler/qscripts``.


Analyzing Results
=================

The total electrostatic energy is read from the output file of APBS
(the ``job_NNNN.out`` file in each window directory). A script makes
this easy (and plots the result if the Python package
:mod:`matplotlib` is installed)::

  apbs-bornprofile-analyze.py Na_nAChR.cfg

The output would be

`welec_Na_nAChR.dat``
   Simple data file of the Born profile with *z* in column 1 and the
   energy in kJ/mol in column 2.

`welec_Na_nAChR.pdf``
   The Born profile plot in PDF format.

