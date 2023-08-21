#!/bin/env python
#Copyright (C) 2023, GITRpy contributors
#Authors: Abdou Diaw
"""
This module is part of GITRpy (Global Impurity TRansport in Python) code.
It generates a netCDF file containing all impurities data.
"""
import os
import netCDF4 as nc
import numpy as np
from periodictable import oxygen, hydrogen, tungsten, deuterium, helium, carbon

def createParticleSourceFile(filename, data_dict):
    """
    Create a netCDF file for containing all impurities data
    Inputs:
    filename -- name of the netCDF file
    data_dict -- a dictionary containing all impurities data
    Outputs:
    filename -- netCDF file
    """
    # Load data from data_dict
    nParticles = data_dict['nParticles']
    # Create a netCDF file for containing all impurities data
    if os.path.exists(filename):
        os.remove(filename)
    rootgrp = nc.Dataset(filename, "w", format="NETCDF4")
    npp = rootgrp.createDimension("nP", nParticles)
    rootgrp.createVariable("impurity_id","f8",("nP"))
    rootgrp.createVariable("x","f8",("nP"))
    rootgrp.createVariable("y","f8",("nP"))
    rootgrp.createVariable("z","f8",("nP"))
    rootgrp.createVariable("vx","f8",("nP"))
    rootgrp.createVariable("vy","f8",("nP"))
    rootgrp.createVariable("vz","f8",("nP"))
    rootgrp.createVariable("charge","f8",("nP"))
    rootgrp.createVariable("mass","f8",("nP"))
    rootgrp.createVariable("IonizationState","f8",("nP"))
    rootgrp.createVariable("materialName", "S10", ("nP"))
    # Populate variables with data
    rootgrp.variables["impurity_id"][:] = data_dict['impurity_id']
    rootgrp.variables["x"][:] = data_dict['x']
    rootgrp.variables["y"][:] = data_dict['y']
    rootgrp.variables["z"][:] = data_dict['z']
    rootgrp.variables["vx"][:] = data_dict['vx']
    rootgrp.variables["vy"][:] = data_dict['vy']
    rootgrp.variables["vz"][:] = data_dict['vz']
    rootgrp.variables["charge"][:] = data_dict['charge']
    rootgrp.variables["mass"][:] = data_dict['mass']
    rootgrp.variables["IonizationState"][:] = data_dict['IonizationState']
    rootgrp.variables["materialName"][:] = data_dict['materialName']
    rootgrp.close()
    return filename

def readParticleSourceFile(filename):
    """
    Read a netCDF file containing all impurities data
    Inputs:
    filename -- name of the netCDF file
    Outputs:
    data_dict -- a dictionary containing all impurities data
    """
    # Load data from data_dict
    data_dict = {}
    rootgrp = nc.Dataset(filename, "r", format="NETCDF4")
    data_dict['nParticles'] = len(rootgrp.dimensions['nP'])
    data_dict['charge'] = rootgrp.variables["charge"][:]
    data_dict['materialName'] = rootgrp.variables["materialName"][:]
    data_dict['mass'] = rootgrp.variables["mass"][:]
    data_dict['IonizationState'] = rootgrp.variables["IonizationState"][:]
    data_dict['impurity_id'] = rootgrp.variables["impurity_id"][:]
    data_dict['x'] = rootgrp.variables["x"][:]
    data_dict['y'] = rootgrp.variables["y"][:]
    data_dict['z'] = rootgrp.variables["z"][:]
    data_dict['vx'] = rootgrp.variables["vx"][:]
    data_dict['vy'] = rootgrp.variables["vy"][:]
    data_dict['vz'] = rootgrp.variables["vz"][:]
    rootgrp.close()
    return data_dict

#


if __name__ == "__main__":
    filename = "input/particleSource.nc"
    te = 20.0
    ti = 20.0
    material1  = carbon
    material2 = oxygen

    # Domain:
    xmin = 0.01
    xmax = 0.01
    ymin = 0.01
    ymax =  0.05
    zmin = 0.05
    zmax = 0.05
    Z1=material1.number

    Z2=material2.number
    A1=material1.mass
    A2=material2.mass
    vth1 = np.sqrt(1.602e-19*te/(1.67e-27*A1))
    vth2 = np.sqrt(1.602e-19*ti/(1.67e-27*A2))
    nParticles = int(100)
    data_dict = {}
    data_dict['nParticles'] = nParticles
    data_dict['charge'] = np.zeros(nParticles)
    data_dict['mass'] = np.zeros(nParticles)
    data_dict['IonizationState'] = np.zeros(nParticles)
    data_dict['impurity_id'] = np.zeros(nParticles)
    data_dict['x'] = np.zeros(nParticles)
    data_dict['y'] = np.zeros(nParticles)
    data_dict['z'] = np.zeros(nParticles)
    data_dict['vx'] = np.zeros(nParticles)
    data_dict['vy'] = np.zeros(nParticles)
    data_dict['vz'] = np.zeros(nParticles)
    data_dict['vz'] = np.zeros(nParticles)
    data_dict['materialName'] = np.empty(nParticles, dtype='S10')


    # Create to particles species with same number and different mass and charge
    for  i in range(0,nParticles):
#        if i < int(nParticles/2):
        if i < int(nParticles/2):
            data_dict['charge'][i] = 1 #material1.number
            data_dict['mass'][i] = material1.mass
            data_dict['IonizationState'][i] = int(material1.number)
            data_dict['materialName'][i] ='C'
   
            data_dict['impurity_id'][i] = i
            data_dict['x'][i] = np.random.uniform(xmin,xmax)
            data_dict['y'][i] = np.random.uniform(ymin,ymax)
            data_dict['z'][i] = np.random.uniform(zmin,zmax)
 
            vn = vth1  * np.sqrt(-np.log(np.random.uniform(0, 1)))
            vr = vth1 * np.sqrt(-np.log(np.random.uniform(0, 1)))
            theta1  = 2 * np.pi * np.random.uniform(0, 1)
            theta2 = 2 * np.pi * np.random.uniform(0, 1)

            data_dict['vx'][i] = np.cos(theta1)*vn
            data_dict['vy'][i] = np.cos(theta2)*vr  
            data_dict['vz'][i] = np.sin(theta2)*vr

        else:
            data_dict['charge'][i] = 1
            data_dict['mass'][i] = material2.mass
            data_dict['IonizationState'][i] = int(material2.number)
            data_dict['materialName'][i] = 'O'
            data_dict['impurity_id'][i] = i

            data_dict['x'][i] = np.random.uniform(xmin,xmax)
            data_dict['y'][i] = np.random.uniform(ymin,ymax)
            data_dict['z'][i] = np.random.uniform(zmin,zmax)

            vn = vth2  * np.sqrt(-np.log(np.random.uniform(0, 1)))
            vr = vth2 * np.sqrt(-np.log(np.random.uniform(0, 1)))
            theta1  = 2 * np.pi * np.random.uniform(0, 1)
            theta2 = 2 * np.pi * np.random.uniform(0, 1)

            data_dict['vx'][i] = np.cos(theta1)*vn
            data_dict['vy'][i] = np.cos(theta2)*vr  
            data_dict['vz'][i] = np.sin(theta2)*vr

    createParticleSourceFile(filename, data_dict)
#
#    # Create to particles species with same number and different mass and charge
#    for  i in range(0,nParticles):
#        data_dict['charge'][i] = 1
#        data_dict['mass'][i] = hydrogen.mass
#        data_dict['IonizationState'][i] = int(hydrogen.number)
#        data_dict['materialName'][i] = str('H') #'W' # material2.symbol
#        data_dict['impurity_id'][i] = i
#        data_dict['x'][i] = np.random.uniform(xmin,xmax)
#        data_dict['y'][i] = np.random.uniform(ymin,ymax)
#        data_dict['z'][i] = np.random.uniform(zmin,zmax)
#        data_dict['vx'][i] = np.random.normal(0,vth2)
#        data_dict['vy'][i] = np.random.normal(0,vth2)
#        data_dict['vz'][i] = np.random.normal(0,vth2)
#    createParticleSourceFile(filename, data_dict)
    
#
#data = nc.Dataset('OW.nc', "r")
#print(data['x'][:])
##print(data)
