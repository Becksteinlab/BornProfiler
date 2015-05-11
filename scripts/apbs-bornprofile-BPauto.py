#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""
:Author: Lennard van der Feltz
:Year: 2014
:Licence: GPL 3
:Copyright: (c) 2014 Lennard van der Feltz
"""
usage = """%prog -protein -pdbids -ions --pqr_forcefield --Nomembrane"""
import os
import argparse
import logging
import bornprofiler
import sys
import subprocess
from bornprofiler.config import cfg
import bornprofiler.run_setup as run_setup

logger = logging.getLogger("bornprofiler")

parser = argparse.ArgumentParser()
parser.add_argument("-protein")
parser.add_argument("-pdbids",nargs = '+')
parser.add_argument("-ions",nargs = '+')
parser.add_argument("--pqr_forcefield",default = 'CHARMM')
parser.add_argument("--pH", default = None)
parser.add_argument("--Nomembrane",action='store_true')
parser.add_argument("--path", default = True)
parser.add_argument("--pathres", default = 1)
parser.add_argument("--script", default = 'q_ASU.sh')
parser.add_argument("--arrayscript",default = 'array_ASU_workstations.sge')
parser.add_argument("--submit_com",default = None)
args = parser.parse_args()
protein = args.protein
pdbids = args.pdbids
ions = args.ions
forcefield = args.pqr_forcefield
pH = args.pH
Nomembrane = args.Nomembrane
path = args.path
pathres= args.pathres
script = args.script
arrayscript = args.arrayscript
submit_com = args.submit_com

def generate_directory(directory):
    try:
        os.mkdir(directory)
    except OSError:
        logger.warning("Directory {drc} already exists. Did not delete".format(drc=directory))

def prepare_run(protein,pdbids,ions,forcefield,pH,Nomembrane,path,pathres,script,arrayscript,submit_com):
    try:
        import MDAnalysis
    except ImportError:
        logger.fatal("Unable to import MDAnalysis. MDAnalysis required for many steps in run preparation. Available from https://code.google.com/p/mdanalysis/")    
        sys.exit(1)
    logger.info("generating {prot} directory".format(prot=protein))
    generate_directory(protein)
    os.chdir(protein)
    logger.info("generating cfg sections")
    cfg.add_section('plotting')
    cfg.add_section('environment')
    cfg.add_section('bornprofile')
    cfg.add_section('job')
    cfg.set("job","script",script)
    cfg.set("job","arrayscript",arrayscript)
    for pdbid in pdbids:
        logger.info("generating {prot}/{pdb} directory".format(prot=protein,pdb=pdbid))
        generate_directory(pdbid)
        os.chdir(pdbid)
        run_setup.get_protein_pdb(pdbid)
        if Nomembrane:
            pass
        else:
            bot_rad,top_rad,thickness,zbot = run_setup.memplacer(pdbid,None)
        try:
            if pH == None:
                subprocess.call(["pdb2pqr.py","--whitespace","--ff={force}".format(force=forcefield),"{pdb}_protein.pdb".format(pdb=pdbid),"{pdb}.pqr".format(pdb=pdbid)])
            else:
                subprocess.call(["pdb2pqr.py","--whitespace","--ff={force}".format(force=forcefield),"--with-ph={pH}".format(pH=pH),"{pdb}_protein.pdb".format(pdb=pdbid),"{pdb}.pqr".format(pdb=pdbid)])
        except:
            logger.info("pqr generation for {pdb}_protein.pdb failed. This pqr must be generated before running apbs calculations. Make sure pdb2pqr is available.".format(pdb=pdbid))
        logger.info("Calculating protein size and location information")
        u = MDAnalysis.Universe("{pdb}_protein.pdb".format(pdb=pdbid))
        bmin,bmax = u.atoms.bbox()
        centroidz = u.atoms.centroid()[2]
        if path == True:
            logger.info("Creating path")
            pmin = bmin[2]-10
            pmax = bmax[2] + 10
            subprocess.call(["apbs-bornprofile-straightpath.py","0", "0", "{z}".format(z = pmin), "0", "0", "1", "{leng}".format(leng = pmax-pmin), "{stepleng}".format(stepleng = pathres), "--title","Centerline"])
        else:
            logger.info("Calculating path size information")
            P = MDAnalysis.Universe("../../{path}".format(path=path))
            topcorner,botcorner = P.atoms.bbox()
            pmin=botcorner[2]
            pmax=topcorner[2]           
        for ion in ions:
            logger.info("generating {prot}/{pdb}/{ion} directory".format(prot=protein,pdb=pdbid,ion=ion))
            generate_directory(ion)
            os.chdir(ion)
            logger.info("Setting cfg elements based on protein, ion, membrane, and path information")
            if Nomembrane:
                cfg.set('membrane','rtop','0')
                cfg.set('membrane','rbot','0')
                cfg.set('membrane','lmem','0')
                cfg.set('membrane','zmem','0')
            else: 
                cfg.set('membrane','rtop','{}'.format(top_rad))
                cfg.set('membrane','rbot','{}'.format(bot_rad))
                cfg.set('membrane','lmem','{}'.format(thickness))
                cfg.set('membrane','zmem','{}'.format(zbot))
            cfg.set('environment','pqr','../{pdb}.pqr'.format(pdb=pdbid))
            cfg.set('bornprofile','ion',ion)
            #section to ensure fine grids always contained within coarse grids
            if pmax - centroidz > 75:
                box = (abs(centroidz - pmax) + 100)*2
                if box/161 - 250/129 < box/129 -250/129:
                    dime = 161
            elif pmin - centroidz < -75:
                box = (abs(centroidz - pmin) + 100)*2
                if box/161 - 250/129 < box/129 -250/129:
                    dime = 161
            else:
                box = 250
                dime = 129
            cfg.set('bornprofile','glen','[({box},{box},{box}),(100,100,100),(50,50,50)]'.format(box=box))
            cfg.set('bornprofile','dime','[({dime},{dime},{dime}),({dime},{dime},{dime}),({dime},{dime},{dime})]'.format(dime=dime))
            if path==True:
                cfg.set('bornprofile','points','../Centerline.pdb')
            else:
                cfg.set('bornprofile','points','../../../{path}'.format(path=path))

            cfg.set('plotting','protein_bottom','{bot}'.format(bot = bmin[2]))
            cfg.set('plotting','protein_length','{leng}'.format(leng = (bmax-bmin)[2]))
            cfg.set('job','name','{pdbid}line{ion}'.format(pdbid=pdbid,ion=ion))
            logger.info("writing cfg file.")
            with open('{pdbid}_{ion}.cfg'.format(pdbid=pdbid,ion=ion), 'wb') as config_file:
                cfg.write(config_file)
#            if Nomembrane:
#                subprocess.call(['apbs-bornprofile-placeion.py','{pdbid}_{ion}.cfg'.format(pdbid=pdbid,ion=ion),'--nomembrane'])
#            else:
            subprocess.call(['apbs-bornprofile-placeion.py','{pdbid}_{ion}.cfg'.format(pdbid=pdbid,ion=ion)])
            if submit_com != None:
                subprocess.call(['{submit_com}'.format(submit_com=submit_com),'qsub_{pdbid}line{ion}.bash'.format(pdbid=pdbid,ion=ion)])
            os.chdir('..')
        os.chdir("..")
        

bornprofiler.start_logging()
prepare_run(protein,pdbids,ions,forcefield,pH,Nomembrane,path,pathres,script,arrayscript,submit_com)
bornprofiler.stop_logging()
