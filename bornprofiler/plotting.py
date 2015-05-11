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
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from bornprofiler.config import cfg
logger = logging.getLogger("bornprofiler")

def test_seaborn():
    try:
        import seaborn as sns
        return True
    except ImportError:
        logger.info("seaborn import failed. Proceeding with uglier graphs")
        return False

def protein_membrane_plot(axis,cfgs):
    """function for plotting membrane and protein on subplot axis using information from the first cfg in cfgs""" 
    plot_bot,plot_top = axis.get_ylim()
    logger.info("reading cfg")
    logger.warning("Assuming membrane and protein locations identical in all cfgs")
    cfg.readfp(open(cfgs[0]))
    protein_bottom = float(cfg.get('plotting','protein_bottom'))
    protein_top = protein_bottom + float(cfg.get('plotting','protein_length'))
    membrane_bottom = float(cfg.get('membrane','zmem'))
    membrane_top = membrane_bottom + float(cfg.get('membrane','lmem'))
    protein_low,protein_high = [plot_bot*.6 + plot_top*.4,plot_bot*.4+plot_top*.6]
    axis.fill_between([membrane_bottom,membrane_top],[plot_top,plot_top],[plot_bot,plot_bot],facecolor='yellow',alpha=0.3)
    axis.fill_between([protein_bottom,protein_top],[protein_high,protein_high],[protein_low,protein_low],facecolor='red',alpha=0.3)

def plot_markups(fig,axis,title,x_label,y_label,file_title,cfgs=None,seaborn=False):
    """function for applying labels,titles,membrane boxes and saving a matplotlib figure fig with subplot axis, previously plotted"""
    axis.legend(loc='best')
    axis.set_title(title)
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    if cfgs==None:
        pass
    else:
        protein_membrane_plot(axis,cfgs)
    if seaborn:
        import seaborn as sns
        sns.despine(trim=True, fig=fig)
        fig.tight_layout()
    logger.info("saving figure to {file_title}.pdf and {file_title}.png".format(file_title=file_title))
    fig.savefig("{file_title}.png".format(file_title=file_title),format='png')
    fig.savefig("{file_title}.pdf".format(file_title=file_title),format='pdf')


def graph_mult_data(file_list,xcolumn,ycolumn,plot_labels,x_label,y_label,title,file_title,colors=None,cfgs=None,seaborn=False):
    if seaborn:
        seaborn = test_seaborn()
        import seaborn as sns
    logger.info("unpacking data")
    if colors==None:
        datalist = [[numpy.loadtxt(filename),label] for filename,label in zip(file_list, plot_labels)]
    else:
        datalist = [[numpy.loadtxt(filename),label,color] for filename,label,color in zip(file_list, plot_labels,colors)]
    if seaborn:
        sns.set_style("ticks", rc={'font.family': 'Helvetica'})
        sns.set_context("paper")
    logger.info("plotting")
    fig = plt.figure()
    axis = fig.add_subplot(111)
    if seaborn:
        sns.offset_spines(fig=fig)
    if colors==None:
        for data,label in datalist:
            axis.plot( data[:,xcolumn],data[:,ycolumn],label=label)
    else:
        for data,label,color in datalist:
            axis.plot( data[:,xcolumn],data[:,ycolumn],label=label,color=color)
    plot_markups(fig,axis,title,x_label,y_label,file_title,cfgs=cfgs,seaborn=seaborn)

