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


def graph_mult_data(file_list,xcolumn,ycolumn,plot_labels,x_label,y_label,title):
    logger.info("unpacking data")
    datalist = [[numpy.loadtxt(filename),label] for filename,label in zip(file_list, plot_labels)]
    logger.info("plotting")
    for data,label in datalist:
        plt.plot( data[:,2],data[:,3],label=label)
    plt.legend(loc='best')
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    logger.info("saving figure to {title}.png".format(title=title))
    plt.savefig(title)

def graph_mult_data_colors(file_list,xcolumn,ycolumn,plot_labels,colors,x_label,y_label,title):
    logger.info("unpacking data")
    datalist = [[numpy.loadtxt(filename),label,color] for filename,label,color in zip(file_list, plot_labels,colors)]
    logger.info("plotting")
    for data,label,color in datalist:
        plt.plot( data[:,xcolumn],data[:,ycolumn],label=label,color=color)
    plt.legend(loc='best')
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    logger.info("saving figure to {title}.png".format(title=title))
    plt.savefig(title)

def graph_mult_data_cfgs(file_list,xcolumn,ycolumn,plot_labels,x_label,y_label,title,cfgs):
    logger.info("unpacking data")
    datalist = [[numpy.loadtxt(filename),label] for filename,label in zip(file_list, plot_labels)]
    logger.info("reading cfg")
    logger.warning("Assuming membrane and protein locations identical in all cfgs")
    cfg.readfp(open(cfgs[0]))
    protein_bottom = float(cfg.get('plotting','protein_bottom'))
    protein_top = protein_bottom + float(cfg.get('plotting','protein_length'))
    membrane_bottom = float(cfg.get('membrane','zmem'))
    membrane_top = membrane_bottom + float(cfg.get('membrane','lmem'))
    logger.info("plotting")
    fig = plt.figure()
    axis = fig.add_subplot(111)
    for data,label in datalist:
        axis.plot( data[:,2],data[:,3],label=label)
    axis.legend(loc='best')
    axis.set_title(title)
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    plot_bot,plot_top = axis.get_ylim()
    protein_low,protein_high = [plot_bot*.6 + plot_top*.4,plot_bot*.4+plot_top*.6]
    axis.fill_between([membrane_bottom,membrane_top],[plot_top,plot_top],[plot_bot,plot_bot],facecolor='yellow',alpha=0.3)
    axis.fill_between([protein_bottom,protein_top],[protein_high,protein_high],[protein_low,protein_low],facecolor='red',alpha=0.3)
    logger.info("saving figure to {title}.png".format(title=title))
    fig.savefig(title)

def graph_mult_data_colors_cfgs(file_list,xcolumn,ycolumn,plot_labels,colors,x_label,y_label,title,cfgs):
    logger.info("unpacking data")
    datalist = [[numpy.loadtxt(filename),label,color] for filename,label,color in zip(file_list, plot_labels,colors)]
    logger.info("reading cfg")
    logger.warning("Assuming membrane and protein locations identical in all cfgs")
    cfg.readfp(open(cfgs[0]))
    protein_bottom = float(cfg.get('plotting','protein_bottom'))
    protein_top = protein_bottom + float(cfg.get('plotting','protein_length'))
    membrane_bottom = float(cfg.get('membrane','zmem'))
    membrane_top = membrane_bottom + float(cfg.get('membrane','lmem'))
    logger.info("plotting")
    fig = plt.figure()
    axis = fig.add_subplot(111)
    for data,label,color in datalist:
        axis.plot( data[:,xcolumn],data[:,ycolumn],label=label,color=color)
    axis.legend(loc='best')
    axis.set_title(title)
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    plot_bot,plot_top = axis.get_ylim()
    protein_low,protein_high = [plot_bot*.6 + plot_top*.4,plot_bot*.4+plot_top*.6]
    axis.fill_between([membrane_bottom,membrane_top],[plot_top,plot_top],[plot_bot,plot_bot],facecolor='yellow',alpha=0.3)
    axis.fill_between([protein_bottom,protein_top],[protein_high,protein_high],[protein_low,protein_low],facecolor='red',alpha=0.3)

    logger.info("saving figure to {title}.png".format(title=title))
    fig.savefig(title)
