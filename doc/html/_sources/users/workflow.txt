==========
 Workflow
==========

TODO: describe the individual steps

Sample point generation
=======================

Path (1D)
---------

* str8path
* HOLE


Volume (3D)
-----------

* HOLLOW (custom version)


Parameter settings
==================

Radii and Charges
-----------------

Use pdb2pqr to generate the input `PQR`_ file.

.. _PQR: http://www.poissonboltzmann.org/file-formats/biomolecular-structurw/pqr

BornProfiler run input file
---------------------------

* ...
* membrane position
* exclusion zone


Queuing system script
---------------------

Discuss example, highlight what needs to be customized and how.


Generate window input files
===========================

::
  apbs-bornprofile-mplaceion


Run jobs
========

Manually or typically through a queuing system.


Analyze and visualize data
==========================

::
  apbs-bornprofile-analyze

Distinguish 1D/3D.

Talk about using Chimera_ to analyze 3D energy landscapes and add some
tips & tricks.

.. _Chimera: http://www.cgl.ucsf.edu/chimera/

