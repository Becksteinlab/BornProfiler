# APBS BornProfiler
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
# Copyright (c) 2005-2008 Kaihsu Tai <k@kauha.eu>
# Copyright (c) 2005-2013 Oliver Beckstein <orbeckst@gmail.com>
# Copyright (c) 2013-2015 Oliver Beckstein, Lennard van der Feltz
# See the file COPYING for details.
#
# setuptools installation of APBS BornProfiler

from setuptools import setup, find_packages
import versioneer


with open("README.rst") as summary:
    LONG_DESCRIPTION = summary.read()

setup(version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      name="APBS-BornProfiler",
      description="Setting up of Born profile calculations for APBS",
      long_description=LONG_DESCRIPTION,
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="GPLv3",
      url="https://github.com/Becksteinlab/BornProfiler",
      keywords="science",
      packages=find_packages(exclude=['scripts']),
      package_data = {'bornprofiler': ['templates/*.sh', 'templates/*.bash', 'templates/*.ge', 'templates/*.sge',
                                       'templates/*.in', 'templates/*.dat',
                                       'templates/*.cfg']},
      scripts = ["scripts/apbs-bornprofile-potential.py",
                 "scripts/apbs-bornprofile-init.py",
                 "scripts/apbs-bornprofile-analyze.py",
                 "scripts/apbs-bornprofile-analyze3D.py",
                 # "scripts/apbs-bornprofile-datalgebra.py",
                 "scripts/apbs-bornprofile-placeion.py",
                 "scripts/apbs-bornprofile-mksample.py",
                 "scripts/apbs-bornprofile-newpath.py",
                 # "scripts/apbs-bornprofile-NEBpath.py",
                 "scripts/apbs-bornprofile-protein_pdb.py",
                 "scripts/apbs-bornprofile-1Dgraph.py",
                 "scripts/parallel.py", "scripts/fake_qsub",
                 "scripts/apbs-mem-potential.py",
                 "scripts/apbs-mem-properties.py",
                 "scripts/apbs-bornprofile-pdb2dat.py",
                 "scripts/apbs-bornprofile-straightpath.py",
                 "scripts/apbs-bornprofile-BPauto.py",
                 "scripts/apbs-bornprofile-autoanalyze.py",
                 "scripts/apbs-bornprofile-cfg_template.py",
           ],
      install_requires=[
          'numpy>=1.0.3',
          'scipy',
          'networkx'
        ],
      extras_require = {
        'analysis': ['GridDataFormats>=0.2.2',
                     'MDAnalysis',
                     ],
        },
      zip_safe=True,
)
