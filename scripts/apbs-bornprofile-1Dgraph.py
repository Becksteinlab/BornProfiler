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
Produces a pylab plot of the data of the xcolumn vs the data of the ycolumn from a column-orgranized dat file with the labels, colors, and title. Defaults set for apbs analysis of energy vs z coordinate."""
import bornprofiler
import argparse
import logging
import pylab
import sys

logger = logging.getLogger("bornprofiler")
parser = argparse.ArgumentParser()
parser.add_argument('dat_files', nargs = '+')
parser.add_argument('--xcolumn',default = 2)
parser.add_argument('--ycolumn',default = 3)
parser.add_argument('--title', default = 'Energy')
parser.add_argument('--xlabel',default = 'z(angstrom)')
parser.add_argument('--ylabel',default = 'kJ/mol')
parser.add_argument('--plot_labels',default=None, nargs = '+')
parser.add_argument('--colors', default=None ,nargs = '+')
args = parser.parse_args()
dat_files = args.dat_files
xcolumn = args.xcolumn
ycolumn = args.ycolumn
title = args.title
xlabel = args.xlabel
ylabel = args.ylabel
plot_labels = args.plot_labels
colors = args.colors
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

def graph_mult_data_colors(file_list,xcolumn,ycolumnplot_labels,colors,x_label,ylabel,title):
    logger.info("unpacking data")
    datalist = [[pylab.loadtxt(filename),label,color] for filename,label,color in zip(file_list, plot_labels,colors)]
    logger.info("plotting")
    for data,label,color in datalist:
        pylab.plot( data[:,xcolumn],data[:,ycolumn],label=label,color=color)
    pylab.legend()
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    logger.info("saving figure to {title}.png".format(title=title))
    pylab.savefig(title)

def graph_mult_data(file_list,xcolumn,ycolumn,plot_labels,x_label,ylabel,title):
    logger.info("unpacking data")
    datalist = [[pylab.loadtxt(filename),label] for filename,label in zip(file_list, plot_labels)]
    logger.info("plotting")
    for data,label in datalist:
        pylab.plot( data[:,2],data[:,3],label=label)
    pylab.legend()
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    logger.info("saving figure to {title}.png".format(title=title))
    pylab.savefig(title)

# If no labels are specified, generates them from the file names, removing directory and format information

if plot_labels == None:
    plot_labels = [dat_file.split('/')[-1].split('.')[0] for dat_file in dat_files]    
# Follows automatic pylab color scheme if none specified
if colors == None:
    graph_mult_data(dat_files,xcolumn,ycolumn,plot_labels,xlabel,ylabel,title)
else:
    graph_mult_data_colors(dat_files,xcolumn,ycolumn,plot_labels,colors,xlabel,ylabel,title) 

bornprofiler.stop_logging()   
