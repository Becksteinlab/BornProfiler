======================================
 Building and installing BornProfiler
======================================

**BornProfiler** consists of a Python package :mod:`bornprofiler` and a
a stand-alone executable ``draw_membrane2``. The Python package is
needed to set up the calculations and ``draw_membrane2`` is needed as
a helper tool for ``apbs``; hence it needs to be installed on the same
machine where ``apbs`` is going to run.


Required pre-requisites
=======================

* python >= 2.7 and < 3
* NumPy_
* a C compiler such as GNU gcc
* APBS_ >= 1.3

.. _NumPy:  http://numpy.scipy.org
.. _APBS: http://www.poissonboltzmann.org


Installation from Source
========================

Unpack the tar ball::

  tar zxvf BornProfiler-0.9.2.tar.gz

Install the python module and scripts::

  cd BornProfiler
  python setup.py install

(see ``python setup.py install --help`` for guidance on what your options
are.)

Compile the customized (and improved) version of ``draw_membrane``
named ``draw_membrane2``::

  mkdir BUILD && cd BUILD
  cmake -D CMAKE_INSTALL_PREFIX=$HOME -D CMAKE_BUILD_TYPE=Release ../src/drawmembrane
  make
  make install

The ``make install`` step will install the executable
``draw_membrane2a`` under ``CMAKE_INSTALL_PREFIX/bin``; change
``CMAKE_INSTALL_PREFIX`` if you prefer another location.

(cmake_ is not really needed; if you don't have it try the following::

   gcc ./src/drawmembrane/draw_membrane2a.c -o draw_membrane2a -lm -lz

and install manually in a place where you or your shell can find it.)

.. Note:: ``draw_membrane2a`` also needs to be installed on the machine
          where you want to run your BornProfiler jobs: it will run
          together with ``apbs``. If you are going to run you
          calculations on a cluster then ``draw_membrane2a`` (and
          ``apbs``) need to be *both* installed on the cluster.

.. _cmake: http://www.cmake.org/


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
