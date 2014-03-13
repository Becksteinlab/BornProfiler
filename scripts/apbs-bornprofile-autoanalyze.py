import os
import subprocess
import argparse
import shutil
parser=argparse.ArgumentParser()
parser.add_argument("-protein")
parser.add_argument("-pdbids", nargs = '+')
parser.add_argument("-ions", nargs = '+')
args = parser.parse_args()
protein = args.protein
pdbids = args.pdbids
ions = args.ions

def autoanalyze(protein,pdbids,ions):
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

autoanalyze(protein,pdbids,ions)
