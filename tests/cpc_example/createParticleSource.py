#!/bin/env python
#Copyright (C) 2023, GITR contributors
#Authors: Abdou Diaw
"""
This module is part of GITR (Global Impurity Transport) code.
It generates a netCDF file containing all impurities data.
"""
import os
import netCDF4 as nc
import numpy as np
from periodictable import oxygen, hydrogen, tungsten

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

# test createParticleSourceFile
if __name__ == "__main__":
    filename = "O.nc"
    te = 10.0
    ti = 20.0
    material1  = oxygen
#    material2 = tungsten

    # Domain:
    xmin, xmax =    -0.0149994, 0.0149994
    ymin, ymax = -0.0149994, 0.0149994
    zmin, zmax = -0.00865992, 0.0299997

#    if mat in material1:
#        names ='O'
#    else:
#        names ='W'
    Z1=material1.number

#    Z2=material2.number
    A1=material1.mass
#    A2=material2.mass
    vth1 = np.sqrt(2*1.602e-19*te/(1.67e-27*Z1))
    nParticles = int(1)
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
    data_dict['materialName'] = np.empty(nParticles, dtype='S10')


    # concentration of each species of Oygen: 1 to 8
    ntotal = 1
    c1 = [ntotal]
    N1 = c1[0] * nParticles
    listN = [N1]
    charges = [1]

    # Create to particles species with same number and different mass and charge
    for  i, charge, N in zip(range(8), charges, listN):
            data_dict['charge'][i] = 1 #charge
            data_dict['mass'][i] = 1 #material1.mass
            data_dict['IonizationState'][i] = 1 #int(material1.number)
            data_dict['materialName'][i] = 'O'
            data_dict['impurity_id'][i] = i
            data_dict['x'][i] = np.random.uniform(xmin,xmax) * 0
            data_dict['y'][i] = np.random.uniform(ymin,ymax) * 0
            data_dict['z'][i] = np.random.uniform(zmin,zmax) * 0
            data_dict['vx'][i] = 1 #np.random.normal(0,vth1)
            data_dict['vy'][i] = 1 #np.random.normal(0,vth1)
            data_dict['vz'][i] = 0 #np.random.normal(0,vth1)
    createParticleSourceFile(filename, data_dict)
#
#data = nc.Dataset('O.nc', "r")
#
#print(data["x"][:])
#

