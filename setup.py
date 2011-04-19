# APBS BornProfiler
# 
# Copyright (c) 2005-2011 Oliver Beckstein <orbeckst@gmail.com>
# Copyright (c) 2005-2008 Kaihsu Tai <k@kauha.eu>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
#
# setuptools installation of APBS BornProfiler

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

# Dynamically calculate the version based on bornprofiler.VERSION.
# (but requires that we can actually import the package BEFORE it is
# properly installed!)
version = __import__('bornprofiler').get_version()

setup(name="APBS-BornProfiler",
      version=version,
      description="Setting up of Born profile calculations for APBS",
      long_description="""
""",
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="GPLv3",
      url="http://sbcb.bioch.ox.ac.uk/oliver/software/",
      keywords="science",
      packages=find_packages(exclude=['scripts']),
      package_data = {'bornprofiler': ['templates/*.sh', 'templates/*.bash', 'templates/*.sge',
                                       'templates/*.in', 'templates/*.dat',
                                       'templates/*.cfg']},
      scripts = ["scripts/apbs-bornprofile-init.py",
                 "scripts/apbs-bornprofile-analyze.py",
                 "scripts/apbs-bornprofile-analyze3D.py",
                 "scripts/apbs-bornprofile-placeion.py",
                 "scripts/apbs-bornprofile-mplaceion.py",
                 "scripts/apbs-bornprofile-mksample.py",
                 "scripts/apbs-bornprofile-newpath.py",
                 "scripts/parallel.py", "scripts/fake_qsub",
                 "scripts/apbs-mem-potential.py",
                 ],
      install_requires=[ 
        # commented out because I get fed up with easy_install's broken(?) dependency tracking...
        #'numpy>=1.0.3', 
        ], 
      extras_require = {
        'analysis': ['gridData>=0.3',
                     ],
        },
      dependency_links = ["http://sbcb.bioch.ox.ac.uk/oliver/download/Python/"],
      zip_safe=True,
)
