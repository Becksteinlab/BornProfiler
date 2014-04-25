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
usage = """%prog dat_files --xcolumn --ycolumn --title --xlabel --ylabel --plot_labels --colors
Produces a matplotlib.pyplot plot of the data of the xcolumn vs the data of the ycolumn from a column-orgranized dat file with the labels, colors, and title. Defaults set for apbs analysis of energy vs z coordinate."""
import bornprofiler
import argparse
import logging
import matplotlib
matplotlib.use('AGG')
import sys
from bornprofiler import plotting
from bornprofiler.config import cfg
logger = logging.getLogger("bornprofiler")
parser = argparse.ArgumentParser()
parser.add_argument('dat_files', nargs = '+')
parser.add_argument('--xcolumn',default = 2)
parser.add_argument('--cfgs', nargs = '+')
parser.add_argument('--ycolumn',default = 3)
parser.add_argument('--title', default = 'Born Energy')
parser.add_argument('--xlabel',default = r'z ($\AA$)')
parser.add_argument('--ylabel',default = r'$\mathcal{W}$$_\mathrm{elec}$ (kJ/mol)')
parser.add_argument('--plot_labels',default=None, nargs = '+')
parser.add_argument('--colors', default=None ,nargs = '+')
parser.add_argument('--file_title', default = 'Born_Energy')
args = parser.parse_args()
dat_files = args.dat_files
xcolumn = args.xcolumn
cfgs = args.cfgs
ycolumn = args.ycolumn
title = args.title
xlabel = args.xlabel
ylabel = args.ylabel
plot_labels = args.plot_labels
colors = args.colors
file_title = args.file_title
bornprofiler.start_logging()
# Checks to ensure matching numbers of labels, colors, and files.
if  plot_labels == None or len(dat_files) == len(plot_labels):
    if colors == None or len(dat_files) == len(colors):
        pass
    else:
        logger.fatal("Number of colors does not match number of data files.")
        sys.exit(1)
else:
    logger.fatal("Number of plot labels does not match number of data files.")
    sys.exit(1)


# If no labels are specified, generates them from the file names, removing directory and format information

if plot_labels == None:
    plot_labels = [dat_file.split('/')[-1].split('.')[0] for dat_file in dat_files]    
# Follows automatic pylab color scheme if none specified. Additionally determines whether to apply cfg methods or not.
if cfgs ==None:
    if colors == None:
       plotting.graph_mult_data(dat_files,xcolumn,ycolumn,plot_labels,xlabel,ylabel,title,file_title)
    else:
        plotting.graph_mult_data_colors(dat_files,xcolumn,ycolumn,plot_labels,colors,xlabel,ylabel,title,file_title) 
else:
    if colors == None:
        plotting.graph_mult_data_cfgs(dat_files,xcolumn,ycolumn,plot_labels,xlabel,ylabel,title,file_title,cfgs)
    else:
        plotting.graph_mult_data_colors_cfgs(dat_files,xcolumn,ycolumn,plot_labels,colors,xlabel,ylabel,title,file_title,cfgs)

bornprofiler.stop_logging()   
