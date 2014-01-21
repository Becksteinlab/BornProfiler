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
:Year: 2013
:Licence: GPL 3
:Copyright: (c) 2013 Lennard van der Feltz
"""
usage = """ %prog dat_file path_pdb --spacing --savestep
Code for determining minimum energy path between a fixed beginning and endpoint.
Employs the nudged elastic band method. This method takes the starting and 
ending point of a given path and uses a given energy landscape to determine how the points in between should move. Forces perpendicular to the path are a result
of the free energy gradient while a simple elastic force between each point and 
its highest energy neighbor is used to determine the forces parallel to the path
An additional damping factor helps motion reach equilibrated state. The forces 
are then integrated using a verlet integrator. Once all points are stable the 
resulting path is returned, along with the path from each step divisible by savestep."""
import argparse
import numpy
import logging
import itertools
import traceback
import os
import sys
import pathmetrics_par
import bornprofiler
import MDAnalysis
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
logger=logging.getLogger("bornprofiler")
parser = argparse.ArgumentParser()
parser.add_argument('dat_file')
parser.add_argument('path_pdb')
parser.add_argument('--spacing',type=float, default=1.0)
parser.add_argument('--savestep',type=float, default=1.0)
parser.add_argument('--title',default = "NEB")
args = parser.parse_args()
path_pdb = args.path_pdb
dat_file = args.dat_file
spacing = args.spacing 
savestep = args.savestep
title = args.title

def delaunay_interpolation(point, Polyhedra, field, fieldmax):
    """Returns the interpolated value of the field.
Finds the delaunay polyhedron which contains point, then interpolates the 
value from the data values at the vertices based on proximity to each vertex
in comparison to the total distance from all vertices."""
    simplex = Polyhedra.find_simplex(point)
    if simplex == -1:
        interpol = fieldmax*4
    else:
        vertices = Polyhedra.vertices[simplex]
        vals = field[vertices]
        interpol = griddata(Polyhedra.points[vertices],vals,point,'linear')
    return interpol

def inbounds(vect, lowbounds,highbounds):
    """Returns True if vect is within bounds, and False otherwise."""
    if numpy.all(vect > lowbounds) and numpy.all(vect < highbounds):
        return True
    else:
        return False


def gradient_at_point(point, Interpolator, grid_spacing):
    """Returns the gradient of grid_object obtained by use of the two-point 
    central difference algorithm"""
    point = numpy.reshape(point,(1,3))    
    spacings = grid_spacing/4.0*numpy.identity(3)
    xgrad = (Interpolator.__call__(point + spacings[0]) - Interpolator(point - spacings[0]))/(grid_spacing/2.0)
    ygrad = (Interpolator(point + spacings[1]) - Interpolator(point - spacings[1]))/(grid_spacing/2.0)
    zgrad = (Interpolator(point + spacings[2]) - Interpolator(point - spacings[2]))/(grid_spacing/2.0)
    return numpy.append(xgrad,[ygrad,zgrad])

def gradient_array(array,Interpolator,grid_spacing):
    """Returns the array of the gradients of the field grid"""
    gradients = numpy.array(map(lambda point: gradient_at_point(point, Interpolator,grid_spacing), array))
    return gradients

def energy_at_point(point, points, energies,energymax):
    return griddata(points,energies,point,fill_value=4*energymax,method='linear')
        

def energies_along_array(array, points,energies,energymax):
    return  Interpolator.__call__(array)

def tangential_vectors(path,path_energies):
    """Returns an array of tangential vectors tangential to the path defined as 
    the vector that points to the highest energy neighbor. Skips endpoints and 
    simply uses zeroes for them."""
    n = 1
    tangent_vectors = numpy.array([0,0,0])
    while n < len(path)-1:
        if path_energies[n-1] > path_energies[n+1]:
            tangent_vector = path[n-1] - path[n]
            tangent_vectors = numpy.vstack([tangent_vectors,tangent_vector])
        else:
            tangent_vector = path[n+1] - path[n]
            tangent_vectors = numpy.vstack([tangent_vectors,tangent_vector])
        n += 1
    tangent_vectors = numpy.vstack([tangent_vectors, numpy.array([0,0,0])])
    return tangent_vectors

def verlet(point, old_point, acceleration, time_step):
    """verlet integrates to give new position after time step"""  
    new_point = 2*point - old_point + acceleration*(time_step**2)
    return new_point

def verlet_array(current_array, old_array, accelerations, time_step):
    """Returns positions resulting from application of verlet integration
    to all entries in array assuming they correspond to the members
    of old_array and accelerations at matching index"""
    new_path = map(verlet, current_array, old_array, accelerations,[time_step]*len(current_array))
    return numpy.array(new_path)

def elastic_forces_along_array(path, equil_dists, k):
    """returns elastic forces on each point in a path. Returns zeroes for 
    endpoint forces."""
    n = 1
    elastic_forces = numpy.array([0,0,0])
    while n < len(path) - 1:
        force1 = elastic_force(path[n],path[n-1],k,equil_dists[n-1])
        force2 = elastic_force(path[n],path[n+1],k,equil_dists[n])
        net_force = force1+force2
        elastic_forces = numpy.vstack([elastic_forces, net_force])
        n += 1
    elastic_forces = numpy.vstack([elastic_forces, numpy.array([0,0,0])])
    return elastic_forces

def elastic_force(point1,point2,k,equil_dist):
    """Returns elastic force from point 1 to point 2 given equilibrium distance"""
    displacement = (point2-point1)
    distance = (numpy.sum(displacement*displacement))**0.5
    force = (displacement/distance)*(distance-equil_dist)*k
    return force

def divide_vect_by_scalar(vect,scalar):
    if scalar == 0:
        return vect
    else:
        return vect/scalar

def normalize_vects(array_of_vects):
    """Returns a normalized version of the input array"""
    norm_array = numpy.sum(array_of_vects*array_of_vects,axis=1)
    normalized_vects = numpy.array(map(divide_vect_by_scalar,array_of_vects,norm_array))
    return normalized_vects


def mult_vect_by_scalar(vect,scalar):
    return vect*scalar

def neighboring_displacements(array):
    n = 1
    disps = array[n] - array[n-1]
    while n < len(array) - 1.0:
        disp = array[n+1] - array[n]
        disps = numpy.vstack([disps,disp])
        n += 1
    return disps

def NEB(initial_path, Interpolator, spacing,savestep,k,time_step,title):
    """Returns a path(numpy array of points) of the initial path relaxed
    upon the energy grid. NEB algorithm essentially places beads at each point
    along the path joined together by elastic connectors whose equilibrium 
    length is the initial spacing. See script documentation for additional info"""
# bootstrapping
    tolerance = spacing/100.0
    damping_factor = 0.05
    all_paths = numpy.array([initial_path])
    path_energies = Interpolator.__call__(initial_path)
    all_paths_energies = numpy.array([path_energies])
#initial displacements to neighboring beads(points)
    initial_disps = neighboring_displacements(initial_path)
    equil_dists = numpy.sum(initial_disps*initial_disps,axis=1)
    current_path = initial_path
    previous_path = initial_path
    n = 0
#This initializes the system before we begin checking for equilibration
#by running through the algorithm 50 times
    logger.info("initializing NEB")
    while n < 50:
        tangent_vectors = tangential_vectors(current_path,path_energies)
        normalized = normalize_vects(tangent_vectors)
        elastic_forces = elastic_forces_along_array(current_path,equil_dists,k)
        tangent_force_mags = numpy.sum(normalized*elastic_forces,axis=1)
        tangent_forces = normalized*tangent_force_mags[:, numpy.newaxis]
        gradient_along_path = gradient_array(current_path,Interpolator,spacing)
        tangent_mags_of_gradient = numpy.sum(normalized*gradient_along_path,axis=1)
        tangent_gradient_component = normalized * tangent_mags_of_gradient[:, numpy.newaxis]
        perpendicular_forces = gradient_along_path - tangent_gradient_component
        perpendicular_forces[0] = numpy.array([0,0,0])
        perpendicular_forces[len(perpendicular_forces)-1] = numpy.array([0,0,0])
        total_forces = tangent_forces + perpendicular_forces
        new_path = verlet_array(current_path, previous_path, total_forces, time_step)
        previous_path = current_path
        current_path = new_path
        path_energies = Interpolator.__call__(current_path)
        all_paths = numpy.vstack([all_paths,numpy.array([current_path])])
        all_paths_energies = numpy.vstack([all_paths_energies, numpy.array([path_energies])])
        n += 1
    hausdorff = pathmetrics_par.d_H((current_path, previous_path))
    logger.info("changing to directory for path pdbs")
    os.chdir("{title}_path_pdbs".format(title=title))
    logger.info("beginning equilibration check and damping")
#Runs the algorithm until hausdorf distances between current and previous 
#fall within tolerance.
    while hausdorff > tolerance:
        if n%savestep == 0:
            logger.info("Saving a pdb of the path")
            write_pdb_from_positions_array(current_path,"{nth}path".format(nth = n),path_energies)
        tangent_vectors = tangential_vectors(current_path,path_energies)
        normalized = normalize_vects(tangent_vectors)
        elastic_forces = elastic_forces_along_array(current_path,equil_dists,k)
        tangent_force_mags = numpy.sum(normalized*elastic_forces,axis=1)
        tangent_forces = normalized*tangent_force_mags[:, numpy.newaxis]
        gradient_along_path = gradient_array(current_path,Interpolator,spacing)
        tangent_mags_of_gradient = numpy.sum(normalized*gradient_along_path,axis=1)
        tangent_gradient_component = normalized * tangent_mags_of_gradient[:, numpy.newaxis]
        perpendicular_forces = gradient_along_path - tangent_gradient_component
        perpendicular_forces[0] = numpy.array([0,0,0])
        perpendicular_forces[len(perpendicular_forces)-1] = numpy.array([0,0,0])
        total_forces = tangent_forces + perpendicular_forces
        new_path = verlet_array(current_path, previous_path, total_forces, time_step)   
        hausdorff = pathmetrics_par.d_H((new_path, current_path))
#Damping to allow system to reach equilibration. Effectively reduces velocity.  
        previous_path = current_path + (new_path - current_path)*damping_factor
        current_path = new_path
        path_energies = Interpolator.__call__(current_path)
        all_paths = numpy.vstack([all_paths,numpy.array([current_path])])
        all_paths_energies = numpy.vstack([all_paths_energies, numpy.array([path_energies])])
        n += 1
        logger.info("Hausdorff distance = {haus}".format(haus=hausdorff))
    logger.info("NEB completed, number of steps = {num}, final hausdorff distance = {haus}".format(num = n,haus = hausdorff))
    logger.info("Writing final path pdb")
    write_pdb_from_positions_array(all_paths[-1],"final{path}".format(path = title),all_paths_energies[-1])
    os.chdir("..")
    return [all_paths,all_paths_energies]
        


def positions_array_from_pdb(pdb):
    """Returns numpy array of positions from PDB extracted via MDAnalysis"""    
    U = MDAnalysis.Universe(pdb)
    Points = U.selectAtoms("all")
    return Points.positions

def write_pdb_from_positions_array(positions, pdb_title, energies):
    """Writes a pdb of dummy atoms with given title""" 
    pdb_string = ""
# Initializing atomnum
    atomnum = 1
#establishing progress along path length
    displacements = neighboring_displacements(positions)
    distances = numpy.sum((displacements*displacements)**0.5,axis=1)
    cumulative_distance = 0
    for position in positions:
        pdb_string +="ATOM {atom:>6}   X  XXX X{atom:>4}{x:>12,.2f}{y:>8,.2f}{z:>8,.2f}{progress:>6,.2f}{energy:>6} \n".format(atom = atomnum,x=position[0],y=position[1],z=position[2],progress=cumulative_distance, energy = energies[atomnum-1])
        atomnum += 1
        try:
            cumulative_distance += distances[atomnum-1]
        except:
            pass
    pdb = open("{title}.pdb".format(title=pdb_title), 'w')
    pdb.write(pdb_string)
    pdb.close()

def read_energy_dat(dat_file):
    """Returns an array of floats.
    Reads a dat file of format x y z energy into an array of floats whose
    members correspond to the lines of the file."""
    energy_file = open(dat_file, "r")
#reading first line because it merely contains labels
    energy_file.readline()
#reading second line to initialize arrays
    line1 = energy_file.readline()
    points = numpy.array([line1.split()[0:3]])
    energies = numpy.array([line1.split()[3]])
    for line in energy_file:
        line_array = numpy.array(line.split())
        points = numpy.vstack((points,line_array[0:3]))
        energies = numpy.vstack([energies,line_array[3]])
    energy_floats = energies.astype(numpy.float)
    point_floats = points.astype(numpy.float)
    data_list = [point_floats,energy_floats]
    return data_list

def write_dat_of_paths_energies(paths,dat_title,energies_arrays):
    """writes dat files of the format x y z energy for each path in paths with
    corresponding energies"""
    n=0
    dat_string = "#{x:>7}{y:>10}{z:>10}{number:>13}{progress:>10}{energy:>10} \n".format(number = "atomnum",x="x",y="y",z="z",progress="distance", energy = "energy")
    for path in paths:
        energies=energies_arrays[n]
# Initializing atom number
        atomnum = 1
#establishing progress along path length
        displacements = neighboring_displacements(path)
        distances = numpy.sum((displacements*displacements)**0.5,axis=1)
        cumulative_distance = 0
        current_bead = 0
        for position in path:
            dat_string +="{x:>10,.2f}{y:>10,.2f}{z:>10,.2f}{number:>10,.2f}{progress:>10,.2f}{energy:>10,.2f} \n".format(number = atomnum,x=position[0],y=position[1],z=position[2],progress=cumulative_distance, energy = energies[atomnum-1][0])
            
            atomnum += 1
            try:
                cumulative_distance += distances[atomnum-1]
            except:
                pass
        n += 1
        dat_string += "&\n"
    dat = open("{title}.dat".format(title=dat_title), 'w')
    dat.write(dat_string)
    dat.close()

def init(path_pdb,dat_file,spacing,savestep,title):
    """Starts up the logger and imports the necessary libraries,
    Checks to ensure that the pdb directory does not already exist.
    loads the dat file and parses into a delaunay object of the points while
    producing an array of the corresponding energies. 
    Loads the path from a pdb into a numpy array,
    then passes them into the NEB process. Writes the final path as a pdb 
    with the title final{title}."""
    logger.info("creating directory for path pdbs")
    try:
        os.mkdir("{title}_path_pdbs".format(title=title))
    except OSError:
        traceback.print_exc()
        logger.fatal("directory {title}_path_pdbs already exists. Remove or submit a different title with the --title flag.".format(title=title))
        sys.exit(1)
#    try:
#        from gridData import Grid
#    except ImportError:
#        traceback.print_exc()
#        logger.fatal("gridData library required for this script. Available from https://github.com/orbeckst/GridDataFormats")
#        sys.exit(1)    
    try:
        import MDAnalysis
        import MDAnalysis.analysis
    except ImportError:
        traceback.print_exc()
        logger.fatal("MDAnalysis required for this script. Available from https://code.google.com/p/mdanalysis/")
        sys.exit(1)

    logger.info("Reading energy file")
    data_list = read_energy_dat(dat_file)
    points = data_list[0]
    energies = data_list[1]
    logger.info("Setting constants based on energy properties")
    energymax = numpy.amax(energies)
    energymin = numpy.amin(energies)
    smallest_possible_period = 2.0*numpy.pi/(((energymax-energymin)/(spacing**2))**0.5)
    time_step = smallest_possible_period/50.0
    k_to_overcome_steepest_descent = ((50.0 - energymin)/(spacing**2))
    k = k_to_overcome_steepest_descent
    logger.info("Initializing interpolator for energy data.")
    Interpolator = LinearNDInterpolator(points,energies,fill_value=energymax*4)
    logger.info("Reading path file")
    initial_path = positions_array_from_pdb(path_pdb)
    logger.info("Applying NEB algorithm")
    all_paths,all_paths_energies = NEB(initial_path, Interpolator, spacing,savestep,k,time_step,title)
#    third_path = NEB(numpy.delete(numpy.delete(second_path,0,0),-1,0),free_energy_grid,spacing,k,time_step)
#    fourth_path = NEB(numpy.delete(numpy.delete(third_path,0,0),-1,0),free_energy_grid,spacing,k,time_step)
#    final_path = numpy.delete(numpy.delete(third_path,0,0),-1,0)
    logger.info("Writing all paths and energies to {title}.dat".format(title=title))
    write_dat_of_paths_energies(all_paths,title,all_paths_energies)
    
bornprofiler.start_logging()
init(path_pdb,dat_file,spacing,savestep,title)
bornprofiler.stop_logging()
