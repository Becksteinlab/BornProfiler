=======
 USAGE
=======

Notes on how to use the BornProfiler package.

The basic philosophy of the package is

    1. provide input
        * a pqr file of the protein
        * a list of coordinates at which to calculate the Born energy of an ion 

    2. generate a range of files that can be used to run apbs

    3. run the generated files manually, or better, through the queuing system of a cluster

    4. extract data from finished runs and provide as graphs and/or
       data files for further processing  

The user runs scripts that are named ``apbs-bornprofile-TASK.py``. The
input is provided by a run input file, which is also referred to as a
*run input configuration file*. It typically bears the suffix *.cfg*. All
information about a calculation is stored in such a cfg file with the
idea that this gives a good summary of the calculation, which is
important for repeating it or writing it up.

The input is documented in the example run input file
``examples/example_runinput.cfg``.


Simple dielectric regions
=========================

See ``examples/parsegian``.


Ion channel
===========

See ``examples/2BG9``.


Generic
=======

Prepare a pqr file of, say, an ion channel and create a list of points
(e.g. along the pore axis) on which you want to calculate the
electrostatic potential of mean force. 

The list of points should be either a simple data file with one set of
white-space separate coordinates (x, y, and z in angstroem) per line
(comment lines starting with "#" and empty lines are ignored) or a PDB
file (the x, y, and z coordinates of all ATOM and HETATM records are
used). For simple, straight pores one can compute the list of points
along the straight pore axis. For more complicated pores one can use a
program such as HOLE (``apbs-bornprofile-mksample.py`` can process a
range of points lists such as the ``sph`` output from HOLE.)

See the help functions of the ``apbs-bornprofiler-mplaceion.py`` and
other ``apbs-bornprofiler-*.py`` scripts.


With a membrane
---------------

First generate a blank parameter file::

  apbs-bornprofile-mplaceion.py --template Mine.cfg

Edit the parameter file (use your favorite text editor instead of
``vi`` )::

  vi Mine.cfg

For instance, set the *pqr* file and the position of the membrane. The
*qscript* and *arrayscript* can also be your own templates. For what's
available see ``BornProfiler/bornprofiler/templates/``.

The file format of a run input config file is explained in the comments of
``examples/example_runinput.cfg``. You can use it as a basis for your own
input, too. In any case, you should have a look at this file in order to
understand how to use the BornProfiler package.


Set-up the individual windows::

  apbs-bornprofile-mplaceion.py Mine.cfg

Run the individual windows. For instance, if you have a working SGE
system, run the whole job array with ::

  qsub qsub_mine.sge


Analyze::

  apbs-bornprofile-analyze.py Mine.cfg

In bulk solvent
---------------

If you do not want a membrane then use ``apbs-bornprofiler-placeion.py``. The
``[membrane]`` section of the run configuration file is ignored. 

A second major difference of ``apbs-bornprofiler-placeion.py`` compared to
``apbs-bornprofiler-mplaceion.py`` is that ``apbs-bornprofiler-placeion.py``
uses automatic focusing whereas ``apbs-bornprofiler-mplaceion.py`` must do
manual focusing. Therefore, ``apbs-bornprofiler-placeion.py`` looks at the run
configuration variables glen_ *and* fglen_ to determine the focusing
schedule. (See the APBS documentation for explanations of the variables.)


``apbs-bornprofiler-placeion.py`` behaves similarly to the old
``placeion.sh`` and ``placeion.py`` scripts.


Potential for visualization
===========================

Use APBSmem_ or edit the **[potential]** section in the parameter file to
run a single APBS calculation on a large and fairly coarse grid::

  apbs-mem-potential.py Mine.cfg

The main output consists of two electrostatics potential files:

pot_membrane.dx
    potential in kT/e of the protein in the membrane (low dielectric 2
    or the value of *mdie* in the parameter file)  
pot_bulksolvent.dx
    potential in kT/e of the protein in the solvent (typically, a
    dielectric of ~80 or whatever was set for *sdie* in the parameter
    file)  

Other files of interest:

dielxSm.dx.gz
    dielectric map *with* membrane ("m"). This file can be visualized as a
    density together with the input PQR or PDB structure in order to check that
    the membrane was placed correctly. The values of the "density" are *mdie*
    (e.g. 2) for the hydrophobic core of the membrane, *pdie* (10) for the
    protein, *headgroup_die* (20) for the lipid headgroup region, *sdie* (80)
    for the solvent and *cdie* (typically the same as *sdie*) for the channel
    region.

Note that you probably have to gunzip the resulting dx files so that you can
read it in Chimera, VMD or PyMOL.

The default values for the grid probably work for most proteins but
depending on size you might have to change the values in ::

  [potential]
  glen = (200,200,200)
  dime = (97,97,97)

glen_ is the size of the box in �; the box is centered on the
protein. dime_ are the grid dimensions. For more accurate calculations,
increase dime but be warned that the memory requirements rise
accordingly.


.. _APBSmem: 
   http://mgrabe1.bio.pitt.edu/apbsmem/
.. _glen:
   http://www.poissonboltzmann.org/apbs/user-guide/running-apbs/input-files/elec-input-file-section/elec-keywords/glen
.. _fglen:
   http://www.poissonboltzmann.org/apbs/user-guide/running-apbs/input-files/elec-input-file-section/elec-keywords/fglen
.. _dime:
   http://www.poissonboltzmann.org/apbs/user-guide/running-apbs/input-files/elec-input-file-section/elec-keywords/dime
