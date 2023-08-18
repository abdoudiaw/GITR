
#!/bin/env python
#Copyright (C) 2023, GITRpy contributors
#Authors: Abdou Diaw
"""
This module is part of GITRpy (Global Impurity TRansport in Python) code.
It generates a netCDF file containing surface reactions data.
"""
import os
import netCDF4 as nc
import numpy as np

def createSurfaceReactionsData(filename, data_dict):
    """
    Create a netCDF file containing material data
    Inputs:
    filename -- name of the netCDF file
    data_dict -- a dictionary containing all species data
    Outputs:
    filename -- netCDF file
    """
    # Create a netCDF file for containing all impurities data
    if os.path.exists(filename):
        os.remove(filename)
    rootgrp = nc.Dataset(filename, "w", format="NETCDF4")

    # Store the case name as an attribute of the root group
    
    # Load data from data_dict
    nA = len(data_dict['A'])
    nE = len(data_dict['E'])
    nEdistBins = data_dict['energyDist'].shape[0]
    nEdistBinsRef = data_dict['energyDistRef'].shape[0]
    nAdistBins = data_dict['cosXDist'].shape[0] 
    nAdistBinsRef = data_dict['cosXDistRef'].shape[0]

    # Create dimensions
    rootgrp.createDimension("nA", nA)
    rootgrp.createDimension("nE", nE)
    rootgrp.createDimension("nEdistBins", nEdistBins)
    rootgrp.createDimension("nEdistBinsRef", nEdistBinsRef)
    rootgrp.createDimension("nAdistBins", nAdistBins)
    rootgrp.createDimension("nAdistBinsRef", nAdistBinsRef)

    # Create variables
    A_var = rootgrp.createVariable("A", "f4", ("nA",))
    A_var[:] = data_dict['A']
    E_var = rootgrp.createVariable("E", "f4", ("nE",))
    E_var[:] = data_dict['E']
    spyld_var = rootgrp.createVariable("spyld", "f4", ("nA", "nE"))
    spyld_var[:] = data_dict['spyld']
    cosYDist_var = rootgrp.createVariable("cosXDist", "f4", ("nAdistBins", "nA", "nE"))
    cosYDist_var[:] = data_dict['cosYDist']
    cosXDist_var = rootgrp.createVariable("cosYDist", "f4", ("nAdistBins", "nA", "nE"))
    cosXDist_var[:] = data_dict['cosXDist']
    cosYDistRef_var = rootgrp.createVariable("cosXDistRef", "f4", ("nAdistBinsRef", "nA", "nE"))
    cosYDistRef_var[:] = data_dict['cosYDistRef']
    cosXDistRef_var = rootgrp.createVariable("cosYDistRef", "f4", ("nAdistBinsRef", "nA", "nE"))
    cosXDistRef_var[:] = data_dict['cosXDistRef']
    energyDist_var = rootgrp.createVariable("energyDist", "f4", ("nEdistBins", "nA", "nE"))
    energyDist_var[:] = data_dict['energyDist']
    energyDistRef_var = rootgrp.createVariable("energyDistRef", "f4", ("nEdistBinsRef", "nA", "nE"))
    energyDistRef_var[:] = data_dict['energyDistRef']
    eDistEgrid_var = rootgrp.createVariable("eDistEgrid", "f4", ("nEdistBins",))
    eDistEgrid_var[:] = data_dict['eDistEgrid']
    eDistEgridRef_var = rootgrp.createVariable("eDistEgridRef", "f4", ("nEdistBinsRef",))
    eDistEgridRef_var[:] = data_dict['eDistEgridRef']
    phiGrid_var = rootgrp.createVariable("phiGrid", "f4", ("nAdistBins",))
    phiGrid_var[:] = data_dict['phiGrid']
    thetaGrid_var = rootgrp.createVariable("thetaGrid", "f4", ("nAdistBins",))
    thetaGrid_var[:] = data_dict['thetaGrid']

    rootgrp.close()
    return filename


if __name__ == "__main__":
    # load simple case data
    data = nc.Dataset("simpleSurfaceModel8ev.nc", "r")
    E = data.variables['E'][:]
    A = data.variables['A'][:]
    spyld = data.variables['spyld'][:, :]
    rfyld = data.variables['rfyld'][:, :]
    cosXDist = data.variables['cosXDist'][:, :, :]
    cosYDist = data.variables['cosYDist'][:, :, :]
    cosXDistRef = data.variables['cosXDistRef'][:, :, :]
    cosYDistRef = data.variables['cosYDistRef'][:, :, :]
    energyDist = data.variables['energyDist'][:, :, :]
    energyDistRef = data.variables['energyDistRef'][:, :, :]
    eDistEgrid = data.variables['eDistEgrid'][:]
    eDistEgridRef = data.variables['eDistEgridRef'][:]
    phiGrid = data.variables['phiGrid'][:]
    thetaGrid = data.variables['thetaGrid'][:]

    data.close()
 # Example data dictionary and group names
    data_dict = {
            'E': E,
            'A': A,
            'spyld': spyld,
            'rfyld': rfyld,
            'cosXDist': cosXDist,
            'cosYDist': cosYDist,
            'cosXDistRef': cosXDistRef,
            'cosYDistRef': cosYDistRef,
            'energyDist': energyDist,
            'energyDistRef': energyDistRef,
            'eDistEgrid': eDistEgrid,
            'eDistEgridRef': eDistEgridRef,
            'phiGrid': phiGrid,
            'thetaGrid': thetaGrid,
            }
    
    # We are using the simple data case above for the 3 cases D_W, W_W, O_W
    case_names = ['D_on_W', 'W_on_W', 'O_on_W']
    for case_name in case_names:
        filename = "surfaceReactions_"+str(case_name)+".nc"
        createSurfaceReactionsData(filename, data_dict)
