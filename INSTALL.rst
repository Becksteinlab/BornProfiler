======================
 INSTALL BornProfiler
======================

Required pre-requisites
=======================

* Python 2.5 or better
* NumPy_ 
* a C compiler such as GNU gcc
* APBS_

.. _NumPy: http://numpy.scipy.org
.. _APBS: http://www.poissonboltzmann.org 


Installation
============

Get the source distribution and unpack the tar ball.

Install the python module and scripts::

  python setup.py install --user

(``--user`` might only work for Python 2.6; look at the output of
``python setup.py install --help`` for guidance on what your options
are.)

Compile the customized (and improved) version of ``draw_membrane``::

  mkdir BUILD && cd BUILD
  cmake -D CMAKE_INSTALL_PREFIX=$HOME -D CMAKE_BUILD_TYPE=Release ../src/drawmembrane
  make
  make install

The ``make install`` step will install the executable
``draw_membrane2a`` under ``CMAKE_INSTALL_PREFIX/bin``; change
``CMAKE_INSTALL_PREFIX`` if you prefer another location.

(``cmake`` is not really needed; if you don't have it try the
following::

   gcc ./src/drawmembrane/draw_membrane2a.c -o draw_membrane2a -lm -lz

and install manually in a place where you or your shell can find it.)


Configuration
=============

Finalize your installation by running ::

  apbs-bornprofile-init.py

This should tell you that it set up a configuration file
``~/.bornprofiler.cfg`` and a number of directories.

The default ``~/.bornprofiler.cfg`` looks like this::

   [DEFAULT]
   configdir = ~/.bornprofiler
   templatesdir = %(configdir)s/templates
   qscriptdir = %(configdir)s/qscripts

   [executables]
   apbs = apbs
   drawmembrane = draw_membrane2a

The file can be edited in a text editor. For instance, one can add the
full path to the ``apbs`` and ``draw_membrane2a`` executable binaries.

Any other variables used in run configuration input files can also be
added here and will be used as defaults.

Advanced use: You can drop templates for run scripts into *qscriptdir*
and have the BornProfiler package pick them up automagically. 
