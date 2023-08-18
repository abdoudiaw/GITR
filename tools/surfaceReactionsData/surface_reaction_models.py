#!/bin/env python
#Copyright (C) 2023, GITR contributors
#Authors: Abdou Diaw
"""
This module is part of GITR (Global Impurity Transport) code.
It contains sputtering and reflection models
"""
import math  
import numpy as np
import netCDF4 as nc
from scipy import constants
import materials
#Constants
Q = constants.physical_constants["elementary charge"][0]
PI = constants.pi
AMU = constants.physical_constants["unified atomic mass unit"][0]
ANGSTROM = constants.angstrom
MICRON = constants.micro
NM = constants.nano
CM = constants.centi
EPS0 = constants.epsilon_0
A0 = constants.physical_constants["Bohr radius"][0]
K = constants.physical_constants["atomic unit of permittivity"][0]
ME = constants.physical_constants["electron mass"][0]
SQRTPI = np.sqrt(PI)
SQRT2PI = np.sqrt(2 * PI)
C = constants.physical_constants["speed of light in vacuum"][0]

# This file comes from RustBCA repository
# species, material, kinetic_energy


def bohdansky_heavy_ion(ion, target, energy_eV):
    '''
    Bohdansky sputtering yield formula in the heavy ion (M1/M2 > 0.5) regime.
    Returns None if the target does not have a surface binding energy.

    https://doi.org/10.1063/1.327954
    https://doi.org/10.1016/0168-583X(84)90271-4

    Args:
        ion (dict): a dictionary with the fields Z (atomic number), m (mass)
        target (dict): a dictionary with the fields Z (atomic number), m (mass), Es (surface binding energy)
        energy_eV (float): energy in electron-volts

    Returns:
        Y (float): sputtering yield in atoms/ion
    '''

            # Find target material from materials list
    if hasattr(materials, ion):
        _ion = getattr(materials, ion)
    else:
        raise ValueError("Material not found")
    # FInd particle material from materials list
    if hasattr(materials, target):
        _target = getattr(materials, target)
    else:
        raise ValueError("Material not found")
    

    # Z1 = _ion['Z']
    # M1 =  _ion['m']
    # z2 = target_material['Z']
    # M2 = target_material['m']
    # Us = target_material['Es']
    # Q  = target_material['Q']
    # Z2 = target_material['Z']

    z1 = _ion['Z']
    z2 = _target['Z']
    m1 = _ion['m']
    m2 = _target['m']
    Us = _target['Es']



    # z1 = ion['Z']
    # z2 = target['Z']
    # m1 = ion['m']
    # m2 = target['m']
    # Us = target['Es']

    if Us == 0.: return None
    alpha = 0.3*(m2/m1)**(2./3.)

    reduced_mass_2 = m2/(m1 + m2)
    reduced_mass_1 = m1/(m1 + m2)

    #Following assumptions are for very light ions (m1/m2<0.5)
    K = 0.4
    R_Rp = K*m2/m1 + 1.

    Eth = (1.9 + 3.8*(m1/m2) + 0.134*(m2/m1)**1.24)*Us

    a0 = 0.529*ANGSTROM
    a = 0.885*a0*(z1**(2./3.) + z2**(2./3.))**(-1./2.)
    reduced_energy = 0.03255/(z1*z2*(z1**(2./3.) + z2**(2./3.))**(1./2.))*reduced_mass_2*energy_eV
    sn = 3.441*np.sqrt(reduced_energy)*np.log(reduced_energy + 2.718)/(1. + 6.355*np.sqrt(reduced_energy) + reduced_energy*(-1.708 + 6.882*np.sqrt(reduced_energy)))
    Sn = 8.478*z1*z2/(z1**(2./3.) + z2**(2./3.))**(1./2.)*reduced_mass_1*sn

    sputtering_yield = alpha*Sn*(1 - (Eth/energy_eV)**(2./3.))*(1. - Eth/energy_eV)**2

    return sputtering_yield

def bohdansky_light_ion(ion, target, energy_eV):
    '''
    Bohdansky sputtering yield formula in the light ion (M1/M2 < 0.5) limit.
    Returns None if the target does not have a surface binding energy.

    https://doi.org/10.1063/1.327954
    https://doi.org/10.1016/0168-583X(84)90271-4

    Args:
        ion (dict): a dictionary with the fields Z (atomic number), m (mass)
        target (dict): a dictionary with the fields Z (atomic number), m (mass), Es (surface binding energy)
        energy_eV (float): energy in electron-volts

    Returns:
        Y (float): sputtering yield in atoms/ion
    '''

        # Find target material from materials list
    if hasattr(materials, ion):
        _ion = getattr(materials, ion)
    else:
        raise ValueError("Material not found")
    # FInd particle material from materials list
    if hasattr(materials, target):
        _target = getattr(materials, target)
    else:
        raise ValueError("Material not found")
    
    z1 = _ion['Z']
    z2 = _target['Z']
    m1 = _ion['m']
    m2 = _target['m']
    Us = _target['Es']



    if Us == 0.: return None
    alpha = 0.2

    reduced_mass_2 = m2/(m1 + m2)
    reduced_mass_1 = m1/(m1 + m2)

    #Following assumptions are for very light ions (m1/m2<0.5)
    K = 0.4
    R_Rp = K*m2/m1 + 1.

    Eth = (1.9 + 3.8*(m1/m2) + 0.134*(m2/m1)**1.24)*Us
    
    a0 = 0.529*ANGSTROM
    a = 0.885*a0*(z1**(2./3.) + z2**(2./3.))**(-1./2.)
    reduced_energy = 0.03255/(z1*z2*(z1**(2./3.) + z2**(2./3.))**(1./2.))*reduced_mass_2*energy_eV
    sn = 3.441*np.sqrt(reduced_energy)*np.log(reduced_energy + 2.718)/(1. + 6.355*np.sqrt(reduced_energy) + reduced_energy*(-1.708 + 6.882*np.sqrt(reduced_energy)))
    Sn = 8.478*z1*z2/(z1**(2./3.) + z2**(2./3.))**(1./2.)*reduced_mass_1*sn

    sputtering_yield = 0.042/Us*(R_Rp)*alpha*Sn*(1-(Eth/energy_eV)**(2./3.))*(1-(Eth/energy_eV))**2

    if sputtering_yield > 0:
        return sputtering_yield
    else:
        return 0.

def sputteringYieldYamamura(particle, target, energy_eV):
    '''
    Yamamura sputtering yield formula for normal incidence.
    https://doi.org/10.1080/01422448208226913
    Args:
        particle (object): a particle object
        target (material): material at this target
        energy_eV (float): energy in electron-volts
    Returns:
        Y (float): sputtering yield in atoms/ion
    '''
    # Find the material from materials list
    if hasattr(materials, target):
        target_material = getattr(materials, target)
    else:
        raise ValueError("Material not found")

    # FInd particle material from materials list
    if hasattr(materials, particle):
        projectile_material = getattr(materials, particle)
    else:
        raise ValueError("Material not found")
    

    z1 = projectile_material['Z']
    m1 =  projectile_material['m']

    z2 = target_material['Z']
    m2 = target_material['m']
    Us = target_material['Es']
    Q  = target_material['Q']

    reduced_mass_2 = m2/(m1 + m2)
    reduced_mass_1 = m1/(m1 + m2)

    #Lindhard's reduced energy
    reduced_energy = 0.03255/(z1*z2*(z1**(2./3.) + z2**(2./3.))**(1./2.))*reduced_mass_2*energy_eV

    #Yamamura empirical constants
    K = 8.478*z1*z2/(z1**(2./3.) + z2**(2./3.))**(1./2.)*reduced_mass_1
    a_star = 0.08 + 0.164*(m2/m1)**0.4 + 0.0145*(m2/m1)**1.29
    #Sputtering threshold energy
    Eth = (1.9 + 3.8*(m1/m2) + 0.134*(m2/m1)**1.24)*Us
    #Lindhard-Scharff-Schiott nuclear cross section
    sn = 3.441*math.sqrt(reduced_energy)*math.log(reduced_energy + 2.718)/(1. + 6.355*math.sqrt(reduced_energy) + reduced_energy*(-1.708 + 6.882*math.sqrt(reduced_energy)))
    #Lindhard-Scharff electronic cross section
    k = 0.079*(m1 + m2)**(3./2.)/(m1**(3./2.)*m2**(1./2.))*z1**(2./3.)*z2**(1./2.)/(z1**(2./3.) + z2**(2./3.))**(3./4.)
    se = k*math.sqrt(reduced_energy)

    return 0.42*a_star*Q*K*sn/Us/(1. + 0.35*Us*se)*(1. - math.sqrt(Eth/energy_eV))**2.8

# This file comes from RustBCA repository
def reflectionWierzbickiBiersack(particle, target, energy_eV):
    '''
    Wierzbicki-Biersack empirical reflection coefficient (1994); not as widely
        applicable as Thomas et al.
    https://doi.org/10.1080/10420159408221042
    Args:
        ion (dict): a dictionary with the fields Z (atomic number), m (mass)
        target (dict): a dictionary with the fields Z (atomic number), m (mass)
        energy_eV (float): energy in electron-volts
    Returns:
        R (float): reflection coefficient of ion on target with energy_eV
    '''
    #Wierzbicki and Biersack empirical reflection coefficient (1994)

    # Find target material from materials list
    if hasattr(materials, target):
        target_material = getattr(materials, target)
    else:
        raise ValueError("Material not found")
    # FInd particle material from materials list
    if hasattr(materials, particle):
        projectile_material = getattr(materials, particle)
    else:
        raise ValueError("Material not found")
    

    Z1 = projectile_material['Z']
    M1 =  projectile_material['m']
    z2 = target_material['Z']
    M2 = target_material['m']
    Us = target_material['Es']
    Q  = target_material['Q']
    Z2 = target_material['Z']

    energy_keV = energy_eV/1E3

    #Thomas-Fermi reduced energy
    reduced_energy = 32.55*energy_keV*M2/((M1 + M2)*Z1*Z2*(Z1**0.23 + Z2**0.23))
    mu = M2/M1

    #Here are some empirical coefficients
    a1 = 5.9638
    b1 = 0.0646
    c1 = 52.211

    a2 = 1.26E-3
    b2 = -0.9305
    c2 = 1.235

    #Wierzbicki and Biersack found that you can separate the dependence on mu, e
    RN_mu = math.exp(a1*math.sqrt(1. - b1*(math.log(mu/c1))**2.))
    RN_e = a2*math.exp(b2*(math.log(reduced_energy + 1.))**c2)

    if not (1.03 < mu <= 240):
        print("Warning: Wierzbicki-Biersack may not be accurate for this ion-target pair")
        print(f'False: 1.03 < {mu} <= 240')

    if not (1 < reduced_energy < 10):
        print("Warning: Wierzbicki-Biersack may not be accurate at this energy")
        print(f'False: 1 < {reduced_energy} <= 10')

    return RN_mu*RN_e

def ThompsonAngularDistribution(theta):
    """Thompson angular distribution function

    Args:
        theta (float): azimuthal angle
    Returns:
        f(theta) (float): angular distribution function
    """
    b = 1.0
    f_theta = (1+b)*(np.cos(theta*np.pi/180) )**b * np.sin(theta*np.pi/180)
    f_theta /= np.trapz(f_theta, theta)
    return f_theta

def ThompsonEnergyDistribution(E, Emax, Eb):
    """Thomas energy distribution function
    
    Args:
        E (float): energy of the sputtered atom
        Emax (float): maximum energy of the sputtered atom
        Eb (float): binding energy of the sputtered atom
    
    Returns:
        f(E) (float): energy distribution function
    """
    Ec= Emax/ Eb
    C = 2 * Eb / (1 - (1 + 2* Ec) / (1 + Ec)**2)    
    f_E = C * E / ( Eb + E )**3
    f_E /= np.trapz(f_E, E)
    return np.where(E < Emax, f_E, 0)

def ThompsonEnergyAngularDistribution(E, Emax, theta, Eb):
    """Thomas energy angular distribution function

    Args:
        E (float): energy of the sputtered atom
        Emax (float): maximum energy of the sputtered atom
        theta (float): azimuthal angle
        Eb (float): binding energy of the sputtered atom
    
    Returns:
        f(E,theat) (float): energy angular distribution function
    """

    f_total = ThompsonEnergyDistribution(E, Emax, Eb) * ThompsonAngularDistribution(theta)
    return np.where(E < Emax, f_total, 0)

def sample_Thompson_distribution(Emax, Eb, n_samples):
    """Sample random energies and angles from the Thompson energy and angular distribution functions

    Args:
        Emax (float): maximum energy of the sputtered atom
        Eb (float): binding energy of the sputtered atom
        n_samples (int): number of samples to generate
    
    Returns:
        E_samples (float): sampled energies
        theta_samples (float): sampled angles
    """
    # Sample random energies and angles from the Thompson energy and angular distribution functions
    E = np.linspace(0, Emax, 1000)
    theta = np.linspace(0, np.pi, 1000)
    f_E = ThompsonEnergyDistribution(E, Emax, Eb=Eb)
    f_theta = ThompsonAngularDistribution(theta)
    # Calculate the cumulative distribution functions for energy and angle
    CDF_E = np.cumsum(f_E) / np.sum(f_E)
    CDF_theta = np.cumsum(f_theta) / np.sum(f_theta)
    # Sample random values from the cumulative distribution functions
    E_samples = np.interp(np.random.rand(n_samples), CDF_E, E)
    theta_samples = np.interp(np.random.rand(n_samples), CDF_theta, theta)
    return E_samples, theta_samples


if __name__ == "__main__":
    from matplotlib import pyplot as plt
    # Example usage
    # Emax = 100  # eV
    # Eb = 10  # eV
    # b = 0.68
    # # sample just one energy and angle
    # n_samples = 1
    # E_samples, theta_samples = sample_Thompson_distribution(Emax, Eb, n_samples)
    # print(E_samples, theta_samples)


   # Plot reflection coefficient
    # target = "tungsten"
    # particle = "tungsten"
    # energy_eV = np.logspace(1.2, 3, 100)
    energy_eV = np.linspace(5, 1000, 100)
    # R_WB = np.zeros_like(energy_eV)

    # plt.subplot(1, 2, 1)
    # for i, energy in enumerate(energy_eV):
    #     R_WB[i] = reflectionWierzbickiBiersack(particle, target, energy)
    # # plt.figure()
    # plt.semilogx(energy_eV, R_WB)
    # plt.xlabel('Energy (eV)')
    # plt.ylabel('Reflection coefficient')


    # # Plot sputtering yield using Yamamura-Tawara model
    # plt.subplot(1, 2, 2)
    # energy_eV = np.linspace(20, 1000, 100)
    # S_Y = np.zeros_like(energy_eV)
    # S_Y = [sputteringYieldYamamura(particle, target,e) for e in energy_eV]
    # bohdansky_light_ion
    S_BO = np.zeros_like(energy_eV)
    # S_BO = [bohdansky_light_ion(particle, target,e) for e in energy_eV]

    # heavy_ion
    # SB_BO = np.zeros_like(energy_eV)
    # SB_BO = [bohdansky_heavy_ion(particle, target,e) for e in energy_eV]

    # plt.figure()
    # plt.semilogx(energy_eV, S_Y)
#    plt.semilogx(energy_eV, [bohdansky_light_ion("deuterium", "tungsten",e) for e in energy_eV], label="D on W")
#    plt.semilogx(energy_eV, [bohdansky_light_ion("helium", "tungsten",e) for e in energy_eV], label="He on W")
#    plt.semilogx(energy_eV, [bohdansky_light_ion("boron", "tungsten",e) for e in energy_eV], label="B on W")
    # plt.semilogx(energy_eV, [bohdansky_light_ion("oxygen", "tungsten",e) for e in energy_eV], label="O on W")
    # plt.semilogx(energy_eV, [bohdansky_light_ion("deuterium", "tungsten",e) for e in energy_eV], label="D on W")

# plt.semilogx(energy_eV, [bohdansky_light_ion("helium", "tungsten",e) for e in energy_eV], label="He on W")  
    # Get reflection coefficient reflectionWierzbickiBiersack(particle, target, energy)
    # plt.semilogx(energy_eV, [reflectionWierzbickiBiersack("deuterium", "tungsten",e) for e in energy_eV], label="D on W")
    # plt.semilogx(energy_eV, [reflectionWierzbickiBiersack("oxygen", "tungsten",e) for e in energy_eV], label="O on W")
#    plt.semilogx(energy_eV, [bohdansky_heavy_ion("tungsten", "tungsten",e) for e in energy_eV], label="W on W")
    # plt.semilogx(energy_eV, SB_BO, '--')
    
    for e in energy_eV:

            N = bohdansky_light_ion("deuterium", "tungsten",e) 
            D = bohdansky_light_ion("oxygen", "tungsten",e)
            if  D != 0:
                ratio = N/D
                print(e, ratio)
    print(ratio)
    # data = nc.Dataset("BCA/O_on_W.nc", "r")
    
    # rfyld = data['rfyld'][:]
    # spyld = data['spyld'][:]
    # E = data['nE'][:]
    # angle = data['nA'][:]

    # print('rflyd', rfyld[-1:,].shape)
    # print('E', E.shape)
    plt.legend()
    # plt.semilogx(E,rfyld[0])
    
        

    # plt.ylim(0, 1)
    plt.xlabel('Energy (eV)')
    plt.ylabel('Sputtering yield')
    plt.legend()
    plt.show()



    # 
    # # from materials import *
    # # materials = Materials()
    # target = "beryllium"
    # if hasattr(materials, target):
    #     target_material = getattr(materials, target)
    # print(target_material['Z'])
    # # target_material = getattr(str(deuterium), str(beryllium))
    # # # print(target_material)
