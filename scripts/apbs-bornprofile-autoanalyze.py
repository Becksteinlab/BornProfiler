#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""Runs analysis script apbs-bornprofile-analyze.py on the files in
specified prot/pdb/ion directories. Intended to follow use of BPauto
script. (EXPERIMENTAL)

"""

import os
import subprocess
import shutil

def autoanalyze(protein, pdbids, ions):
    os.chdir(protein)
    for pdbid in pdbids:
        os.chdir(pdbid)
        pdbidstart = pdbid[0]
        for ion in ions:
            os.chdir(ion)
            subprocess.call(["apbs-bornprofile-analyze.py","{pdbid}_{ion}.cfg".format(pdbid=pdbid,ion=ion)])
            os.chdir("..")
        os.chdir("..")
    os.chdir("..")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-protein")
    parser.add_argument("-pdbids", nargs = '+')
    parser.add_argument("-ions", nargs = '+')
    args = parser.parse_args()

    autoanalyze(args.protein, args.pdbids, args.ions)
