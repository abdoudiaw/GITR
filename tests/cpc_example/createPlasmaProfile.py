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

    rootgrp.createVariable("ne","f8",("nP"))
    rootgrp.createVariable("te","f8",("nP"))
    rootgrp.createVariable("ni","f8",("nP"))
    rootgrp.createVariable("ti","f8",("nP"))
    rootgrp.createVariable("vx","f8",("nP"))
    rootgrp.createVariable("vy","f8",("nP"))
    rootgrp.createVariable("vz","f8",("nP"))
    
    rootgrp.variables["ne"][:] = data_dict['ne']
    rootgrp.variables["te"][:] = data_dict['te']
    rootgrp.variables["ni"][:] = data_dict['ni']
    rootgrp.variables["ti"][:] = data_dict['ti']
    rootgrp.variables["vx"][:] = data_dict['vx']
    rootgrp.variables["vy"][:] = data_dict['vy']
    rootgrp.variables["vz"][:] = data_dict['vz']
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
    data_dict['ne'] = rootgrp.variables["ne"][:]
    data_dict['te'] = rootgrp.variables["te"][:]
    data_dict['ti'] = rootgrp.variables["ti"][:]
    data_dict['ni'] = rootgrp.variables["ni"][:]
    data_dict['vx'] = rootgrp.variables["vx"][:]
    data_dict['vy'] = rootgrp.variables["vy"][:]
    data_dict['vz'] = rootgrp.variables["vz"][:]
    rootgrp.close()
    return data_dict

# test createParticleSourceFile
if __name__ == "__main__":
    filename = "plasmaProfile.nc"
    te = 1.0
    ti = 1.0
    material1  = oxygen
    material2 = tungsten

    Z1=material1.number

    Z2=material2.number
    A1=material1.mass
    A2=material2.mass
    vth1 = np.sqrt(2*1.602e-19*te/(1.67e-27*Z1))
    nParticles = int(1)
    data_dict = {}
    data_dict['nParticles'] = nParticles
    data_dict['ne'] = np.zeros(nParticles)
    data_dict['te'] = np.zeros(nParticles)
    data_dict['ni'] = np.zeros(nParticles)
    data_dict['ti'] = np.zeros(nParticles)
    data_dict['vx'] = np.zeros(nParticles)
    data_dict['vy'] = np.zeros(nParticles)
    data_dict['vz'] = np.zeros(nParticles)


    # Create to particles species with same number and different mass and charge
    for  i in range(0,nParticles):
        data_dict['ne'][i] = 1e19
        data_dict['te'][i] = te
        data_dict['ni'][i] = 1e19
        data_dict['ti'][i] = ti
        data_dict['vx'][i] = np.random.normal(0,vth1)
        data_dict['vy'][i] = np.random.normal(0,vth1)
        data_dict['vz'][i] = np.random.normal(0,vth1)

    createParticleSourceFile(filename, data_dict)

#data = nc.Dataset('OW.nc', "r")
# #print(data)
# print(data['IonizationState'][:])

# xmin = -0.0149994
# xmax = 0.0149994
# ymin = -0.0149994
# ymax = 0.0149994
# zmin = -0.00865992
# zmax = 0.0299997

# volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
# volume = abs(volume)
# print("volume ", volume)
# print("nParticles ", nParticles)
# print("density ", nParticles/volume)
