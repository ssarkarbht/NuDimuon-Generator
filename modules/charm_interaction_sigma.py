#!/bin/python

'''
Author : Sourav Sarkar
Date : February 12, 2023
Email: ssarkar1@ualberta.ca
Description : This script generates a charm hadron (meson or baryon)
    interaction cross section with an input medium.
'''

import numpy as np
import chromo as ch
import scipy.interpolate as ip


#Initialize a single chromo generator with the heaviest element
#Note: same chromo generator cannot be intialized multiple times
# multiple times in a script, in addtion once the generator is
# intialized, any heavier target element than the one used for
# initialization throws an error during cross-section computation
# Therefore, initialize the generator with the heaviest relevant
# element (i.e. Lead)

ev_kin = ch.kinematics.EventKinematics(321, (208,82),
                elab=1000.0*ch.constants.GeV)
generator = ch.models.DpmjetIII191(ev_kin)

#Define the CORSIKA parametrized charm meson/baryon cross-section
#with proton

def corsika_Mp(energy):
    ''' Takes the energy (in GeV) as inputs and outputs the
    cross section in mb. Note that the energy must be 1PeV
    or above for the correct parametrization.
    '''
    #Parametrized coefficients from CORSIKA
    CSK1 = 1.891
    CSK2 = 0.2095
    CSK3 = -2.157
    CSK4 = 1.263
    #check if the given energy is above the 1PeV threshold
    if isinstance(energy, np.ndarray):
        assert (energy>=1e6).all(), "Energy input has value(s) less than the threshold"
    else:
        assert energy>=1e6, "Energy input is less than the threshold"

    loge = np.log10(energy)

    xs = np.exp(CSK1 + CSK2*loge)+CSK3+CSK4*loge
    return xs


def corsika_Bp(energy):
    ''' Takes the energy (in GeV) as inputs and outputs the
    cross section in mb. Note that the energy must be 1PeV
    or above for the correct parametrization.
    '''
    #Parametrized coefficients from CORSIKA
    CSK1 = 2.269
    CSK2 = 0.207
    CSK3 = -0.9907
    CSK4 = 1.277
    #check if the given energy is above the 1PeV threshold
    if isinstance(energy, np.ndarray):
        assert (energy>=1e6).all(), "Energy input has value(s) less than the threshold"
    else:
        assert energy>=1e6, "Energy input is less than the threshold"

    loge = np.log10(energy)

    xs = np.exp(CSK1 + CSK2*loge)+CSK3+CSK4*loge
    return xs

def dpmjet_Ktarget(energy, nucleus=2212):
    global generator
    if isinstance(energy, np.ndarray):
        xsarr = np.zeros_like(energy)
        for i,e in enumerate(energy):
            newkin = ch.kinematics.EventKinematics(321, nucleus,
                    elab = e*ch.constants.GeV)
            xs = generator.cross_section(newkin)
            if isinstance(nucleus, int):
                xsarr[i]  = xs.total
            elif isinstance(nucleus, tuple):
                xsarr[i] = xs.inelastic
        return xsarr
    else:
        ev_kin = ch.kinematics.EventKinematics(321, nucleus,
                elab=energy*ch.constants.GeV)
        xs = generator.cross_section(ev_kin)
        if isinstance(nucleus, int):
            return xs.total
        elif isinstance(nucleus, tuple):
            return xs.inelastic

def sigma_Mp(energy):
    ''' This function computes the charm meson - proton
    cross section by scaling the DPMJETIII cross-section
    below 1PeV with the CORSIKA cross section above 1PeV
    '''
    assert isinstance(energy, np.ndarray), "The function takes numpy array input format"
    #compute the DPMJET and CORSIKA at 1PeV for scaling
    dpm = dpmjet_Ktarget(1e6)
    cor = corsika_Mp(1e6)
    alpha = cor/dpm

    #create an empty array for cross section
    xsarr = np.zeros_like(energy)

    #get the indices above and below 1 PeV
    above = np.where(energy>=1e6)[0]
    below = np.where(energy<1e6)[0]

    #compute the corsika above PeV cross-section
    xsarr[above] = corsika_Mp(energy[above])

    #compute the dpmjet below PeV cross-section
    xsarr[below] = alpha*dpmjet_Ktarget(energy[below])

    return xsarr


def sigma_Bp(energy):
    ''' This function computes the charm baryon - proton
    cross section by scaling the DPMJETIII cross-section
    below 1PeV with the CORSIKA cross section above 1PeV
    '''
    assert isinstance(energy, np.ndarray), "The function takes numpy array input format"
    #compute the DPMJET and CORSIKA at 1PeV for scaling
    dpm = dpmjet_Ktarget(1e6)
    cor = corsika_Bp(1e6)
    alpha = cor/dpm

    #create an empty array for cross section
    xsarr = np.zeros_like(energy)

    #get the indices above and below 1 PeV
    above = np.where(energy>=1e6)[0]
    below = np.where(energy<1e6)[0]

    #compute the corsika above PeV cross-section
    xsarr[above] = corsika_Bp(energy[above])

    #compute the dpmjet below PeV cross-section
    xsarr[below] = alpha*dpmjet_Ktarget(energy[below])

    return xsarr

def glauber_factor(energy, medium):
    assert isinstance(energy, np.ndarray), "The function takes numpy array input format"
    if medium==2212 or medium==(1,1):
        print ("No Glauber interpolation needed for proton/hydrogen.")
        return energy, np.ones_like(energy)
    beta = dpmjet_Ktarget(energy,nucleus=medium)/dpmjet_Ktarget(energy)
    return energy, beta

def charm_hadron_sigma(energy, charm_type, nucleus=2212):
    ''' Final function that calculates the charm hadron (meson or baryon)
    interaction cross section with the target nucleus.
    Parameters : 
        energy (np.ndarray) : input energy for cross section calculation
        charm_type (str) : Charm hadron type : 'Meson' or 'Baryon'
        nucleus (int, tuple) [optional] : Target nucleus type. Give interger
        for particle pdg code or tuple containing mass and atomic number of
        the nucleus.
    '''
    assert charm_type in ['Meson', 'Baryon'], f"Incompatible charm hadron type : {charm_type}"

    if charm_type == "Meson":
        barexs = sigma_Mp(energy)
    elif charm_type=="Baryon":
        barexs = sigma_Bp(energy)

    #compute the Glauber factor
    _,beta = glauber_factor(energy, nucleus)

    return beta*barexs
