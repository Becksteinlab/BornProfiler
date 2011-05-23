================================
 README for BornProfiler [2010]
================================

Copyright (c) 2008 Kaihsu Tai
Copyright (c) 2010-2011 Oliver Beckstein

A small Python package to set up "Born" calculations of an ion in a
membrane protein and calculate the electrostatic free energy with
APBS_.

For more details on this *internal SBCB release* see the Wiki at

    https://sbcb.bioch.ox.ac.uk/wiki/index.php/BornProfiler


.. _APBS: http://www.poissonboltzmann.org/apbs


Features
========

The BornProfiler package helps setting up Poisson-Boltzmann
calculations of the electrostatic potential of mean force of an ion in
a pore or channel under the influence of a membrane. The membrane is
modelled as a dielectric slab of epsilon=2.

 * Provide a path (list of coordinates) and a PQR file of the protein as input.
 * A membrane can be defined with arbitrary thickness, z-position, and
   dielectric. A headgroup region can also be defined with a different
   dielectric constant. 
 * Define all input parameters in a compact parameter file so that
   there is always a record of the exact calculation setup available. 
 * Born radii for all ions from the Rashin & Honig paper [Rashin1985]_
   are included; just select the ion in the input file.
 * Born radii for H3O+, OH- (and H+... for testing) have been derived
   from the solvation free energies in [Pliego2000]_ directly via the
   Born equation. USE AT YOUR OWN RISK!!
 * Customize run scripts and queuing system submission scripts by
   providing your own templates.  

.. [Rashin1985] A.Rashin & B.Honig, J Phys Chem B 89 (1985), 5588
.. [Pliego2000] J.R. Pliego and J.M. Riveros. Chemical Physics
                Letters, 332(5-6): 597--602, 2000. 
		doi:10.1016/S0009-2614(00)01305-1.  



History and Contributions
=========================

Based on Kaihsu Tai's Python rewrite
(`Poisson-Boltzmann profile for an ion channel`_) of the original placeion.sh and
analyze.sh bash scripts by Kaihsu Tai and Oliver Beckstein.

Uses material from the APBS Wiki (`PMF of a helix in a membrane`_) and contains a modified version of
Michael Grabe's ``draw_membrane2`` from APBSmem_.

.. _Poisson-Boltzmann profile for an ion channel:
   http://en.wikiversity.org/wiki/Poisson%E2%80%93Boltzmann_profile_for_an_ion_channel

.. _PMF of a helix in a membrane:
   http://www.poissonboltzmann.org/apbs/examples/potentials-of-mean-force/the-polar-solvation-potential-of-mean-force-for-a-helix-in-a-dielectric-slab-membrane

.. _APBSmem: 
   http://mgrabe1.bio.pitt.edu/apbsmem/



