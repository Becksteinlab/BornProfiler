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

"""Repository for various functions used repeatedly in plotting scripts. Many functions expect a logger under the title 'logger' to be in operation."""
import bornprofiler
import sys
import traceback
import argparse
import logging
import numpy
import matplotlib.pyplot as plt
from bornprofiler.config import cfg
logger = logging.getLogger("bornprofiler")


def graph_mult_data(file_list,xcolumn,ycolumn,plot_labels,x_label,ylabel,title):
    logger.info("unpacking data")
    datalist = [[numpy.loadtxt(filename),label] for filename,label in zip(file_list, plot_labels)]
    logger.info("plotting")
    for data,label in datalist:
        plt.plot( data[:,2],data[:,3],label=label)
    plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    logger.info("saving figure to {title}.png".format(title=title))
    plt.savefig(title)

def graph_mult_data_colors(file_list,xcolumn,ycolumn,plot_labels,colors,x_label,ylabel,title):
    logger.info("unpacking data")
    datalist = [[numpy.loadtxt(filename),label,color] for filename,label,color in zip(file_list, plot_labels,colors)]
    logger.info("plotting")
    for data,label,color in datalist:
        plt.plot( data[:,xcolumn],data[:,ycolumn],label=label,color=color)
    plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    logger.info("saving figure to {title}.png".format(title=title))
    plt.savefig(title)

def graph_mult_data_cfgs(file_list,xcolumn,ycolumn,plot_labels,x_label,ylabel,title,cfgs):
    logger.info("unpacking data")
    datalist = [[numpy.loadtxt(filename),label] for filename,label in zip(file_list, plot_labels)]
    logger.info("plotting")
    for data,label in datalist:
        plt.plot( data[:,2],data[:,3],label=label)
    plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    logger.info("saving figure to {title}.png".format(title=title))
    plt.savefig(title)

def graph_mult_data_colors_cfgs(file_list,xcolumn,ycolumn,plot_labels,colors,x_label,ylabel,title,cfgs):
    logger.info("unpacking data")
    datalist = [[numpy.loadtxt(filename),label,color] for filename,label,color in zip(file_list, plot_labels,colors)]
    logger.info("plotting")
    for data,label,color in datalist:
        plt.plot( data[:,xcolumn],data[:,ycolumn],label=label,color=color)
    plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    logger.info("saving figure to {title}.png".format(title=title))
    plt.savefig(title)
