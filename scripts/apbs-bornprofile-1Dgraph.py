#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# BornProfiler --- A package to calculate electrostatic free energies with APBS
# Written by Kaihsu Tai, Lennard van der Feltz, and Oliver Beckstein
# Released under the GNU Public Licence, version 3
#
"""Produces a matplotlib.pyplot plot of the data of the xcolumn vs the
data of the ycolumn from a column-orgranized dat file with the labels,
colors, and title. Defaults set for apbs analysis of energy vs z
coordinate. Saves to a png and pdf of the title provided in
file_title. better_labels compares submitted plot_labels to a list of
standard ions and replaces them with strings containing charge and
bonding information. seaborn turns on better plotting via the seaborn
module."""

import sys
import logging
import matplotlib
matplotlib.use('AGG')

import bornprofiler
from bornprofiler import plotting
from bornprofiler.config import cfg

logger = logging.getLogger("bornprofiler")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dat_files', nargs = '+')
    parser.add_argument('--xcolumn',default = 2)
    parser.add_argument('--cfgs', nargs = '+',default=None)
    parser.add_argument('--ycolumn',default = 3)
    parser.add_argument('--title', default = 'Born Energy')
    parser.add_argument('--xlabel',default = r'z ($\AA$)')
    parser.add_argument('--ylabel',default = r'$\mathcal{W}$$_\mathrm{elec}$ (kJ/mol)')
    parser.add_argument('--plot_labels',default=None, nargs = '+')
    parser.add_argument('--colors', default=None ,nargs = '+')
    parser.add_argument('--file_title', default = 'Born_Energy')
    parser.add_argument('--better_labels',action='store_true')
    parser.add_argument('--seaborn',action='store_true')
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
    better_labels = args.better_labels
    seaborn = args.seaborn

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
    #Converts known ion names to custom labels which include charge and other information
    better_labels_dict = {'Ca':r'Ca$^{2+}$','Na':r"Na$^{+}$",'Cl':r"Cl$^{-}$",'K':r"K$^{+}$",'H30':r"H$_{3}$O$^{+}$",'H':r"H$^{+}$",'Li':r'Li$^{+}$','Rb':r'Rb$^{+}$','Cs':r'Cs$^{+}$','F':r'F$^{-}$','Br':r'Br$^{-}$','I':r'I$^{-}$'}
    if better_labels:
        better_labels = []
        for item in plot_labels:
            if item in better_labels_dict:
                 better_labels.append(better_labels_dict[item])
            else:
                better_labels.append(item)
        plot_labels=better_labels
    # If no labels are specified, generates them from the file names, removing directory and format information

    if plot_labels == None:
        plot_labels = [dat_file.split('/')[-1].split('.')[0] for dat_file in dat_files]
    # Follows automatic pylab color scheme if none specified. Additionally determines whether to apply cfg methods or not.
    plotting.graph_mult_data(dat_files,xcolumn,ycolumn,plot_labels,xlabel,ylabel,title,file_title,colors=colors,cfgs=cfgs,seaborn=seaborn)

    bornprofiler.stop_logging()
