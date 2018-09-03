.. -*- coding: utf-8 -*-

==================
 Basic operations
==================

Potential surface with membrane
===============================

PQR file
--------

Have a PQR file for your membrane protein of interest
available. Getting it from  the *Orientations of Proteins in Membranes
(OPM_) database* is convenient because you can get a sense for the
likely position of the membrane. In the following we assume we use the
GLIC channel (a proton gated cation channel) with OPM/PDB id `3p4w
<https://opm.phar.umich.edu/proteins/831>`_ as an example. Download
the file from OPM and create the PQR file with :program:`pdb2pqr`:

.. code-block:: bash
		
   wget https://storage.googleapis.com/opm-assets/pdb/3p4w.pdb
   pdb2pqr --ff CHARMM --whitespace --drop-water 3p4w.pdb 3p4w.pqr


Membrane position
-----------------

The membrane will be simulated as a low-dielectric brick-shaped slab.

Obtain the boundaries of the membrane that OPM suggested: The file
contains "DUM" atoms that show the position of the membrane
surfaces. The protein is oriented such that the membrane is in the X-Y
plane so we extract the minimum and maximum Z component of the DUM
atoms using some Unix command line utilities:

.. code-block:: bash

   grep '^HETATM.*DUM'  3p4w.pdb | cut -b 47-54 | sort -n | uniq		

which gives

     -16.200
      16.200

i.e., the membrane is located between -16.2 Å and +16.2 Å with a
thickness of 32.4 Å. The center of the membrane is at ``z = 0`` and
the bottom is at ``z = -16.2``.

Pore dimensions
---------------

Ion channels have aqueous pathways. They should be represented by a
high dielectric environment that is also accessible to ions.

In order to add the pore regions we need to get approximate dimensions
of the pore, which is represented as a truncated cone with a center in
the X-Y plane and potentially differing radii at the top and bottom of
the membrane plane.

You can run :program:`HOLE` to get radii and positions.

For this tutorial we simply load the structure in a visualization tool
such a :program:`VMD`, :program:`pymol` or :program:`Chimera` and
estimate a diameter of about 20 Å at the center of the protein, which
is at ``x=0`` and ``y=0`` in the plane of the membrane.



Set up run
----------

Generate a template input file ``3p4w.cfg`` for the protein::

  apbs-mem-potential.py --template 3p4w.cfg

Edit the template in ``[environment]`` section and the set pqr file.

.. code-block:: inifile

  [environment]
  pqr = 3p4w.pqr

  
Edit the template in the ``[membrane]`` section to add data for the
membrane.

* ``lmem`` is the thickness, 32.4 Å.
* ``zmem`` is the bottom of the membrane, ``z=-16.2``
* ``vmem`` is the membrane potential (in ???); here we leave it at 0.
  


.. code-block:: inifile

   [membrane]
   rtop = 10
   rbot = 10
   x0_r = None
   y0_r = None
   dx_r = 0
   dy_r = 0
   cdie = %(solvent_dielectric)s
   headgroup_die = 20
   headgroup_l = 0
   mdie = 2
   vmem = 0
   lmem = 32.4
   zmem = -16.2

   
Run calculation
---------------

Once all information is collected in the cfg file, one runs

.. code-block:: bash

   apbs-mem-potential.py 3p4w.cfg		

This will create input files for :program:`apbs` and run
:program:`drawmembrane2a` when necessary.

The output consists of *dx* files of the potential (in kT/e).
	 
   
.. _OPM: https://opm.phar.umich.edu/

A simple Born profile
=====================

TODO: Outline the problem of ion permeation, discuss simple example
and show how this package can solve the problem. Choose something very
simple such as nAChR or GLIC.

