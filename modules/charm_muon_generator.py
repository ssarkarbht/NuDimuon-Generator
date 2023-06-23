#!/bin/python

'''
Author: Sourav Sarkar
Date: June 15, 2022
Email: ssarkar1@ualberta.ca
Description: This script takes the pythia simulated
    charm hadron events and injects secondary muon
    sample from decay or interaction in ice. The 
    script also stores the probability of each event
    from charm prodcution cross-section and muon
    branching ratio.
'''

import numpy as np
import scipy.interpolate as ip
from numpy.random import default_rng
import sys, os

repo_dir = os.environ["DIMUON_REPO"]


sys.path.insert(0, '/data2/icecube/ssarkar/impy')
from impy.definitions import *
from impy.constants import *
from impy.kinematics import EventKinematics, CompositeTarget
from impy import impy_config, pdata


#Useful functions to be used in the main class
def convert_linlog(x, to='log'):
    ''' Converts the input to log or linear scale
    '''
    if to=='log':
        return np.log10(x)
    elif to=='lin':
        return 10**x

def get_decay_cdf(fname):
    '''Returns discrete cdf values at
    discrete energy fraction points
    '''
    data = np.loadtxt(fname, comments='#')
    xvals = data[0]
    yvals = data[1:]#for 300 hadron energy points
    #compute the average of the distributions
    #as the decay energy distribution for all energies
    #are alomst the same
    y_avg = np.average(yvals, axis=0)

    #get the bin width
    XBINS = np.linspace(0,1,101)#bins used in the dataset
    binwidth = np.diff(XBINS)[0]#all elements are same
    #compute the discrete cdf values
    cdf = np.cumsum(y_avg)*binwidth

    #add the last cdf point (1,1)
    cdf = np.append(cdf,1)
    xvals = np.append(xvals,1)

    return (xvals, cdf)

def logscale_cdf(x,y,binwidth):
    cdf = np.cumsum(y*x)*binwidth
    #normalize in case the logbin approx. is off
    cdf /= max(cdf)
    #add the last datapoint
    cdf = np.append(cdf, 1)
    x   = np.append(x, 1)
    return (x,cdf)

def get_interaction_cdf(fname):
    ''' Returns 2d interpolation function
    to compute cdf values from (at different
    hadron energies)
    '''
    #hadron energy simulation points in the dataset
    EBINS = np.logspace(2,8,300)
    #muon energy fraction bin edges computed in the
    #dataset
    XBINS = np.logspace(-8,0,101)
    #binwidth in logscale
    LOG_BIN_WIDTH = np.diff(np.linspace(-8,0,101))[0]

    data = np.loadtxt(fname, comments='#')
    xvals = data[0]
    yvals = data[1:]

    #prepare the 2d arrays for interpolation
    #array for discrete cdf points across the 2d grid
    cdf_arr = np.zeros((len(EBINS), len(XBINS)))
    for i in range(len(EBINS)):
        x, cdf = logscale_cdf(xvals, yvals[i], LOG_BIN_WIDTH)
        cdf_arr[i] = cdf
    #interpolation is performed in the log-log scale
    logevals = np.log10(EBINS)
    logxvals = np.log10(x)

    func = ip.RectBivariateSpline(logxvals, logevals,
            cdf_arr.T, kx=1, ky=1, s=0)
    return func

def br_intepolation(fname):
    data = np.loadtxt(fname, comments='#')
    x = data[:,0]
    y = data[:,1]/data[:,2]
    logx = np.log10(x)

    func = ip.interp1d(logx, y, kind='slinear')
    return func

#hadron-hadron interaction cross-section for charm hadrons
def sigma_meson(energy):
    loge = np.log10(energy)
    #get the fit constants
    PIP0 = 1.35485692
    PIP1 = -0.02237097
    PIP2 = 0.01398568
    KPP = 0.94044349
    DPP = 0.6864460702906747

    CSK1 = 1.891
    CSK2 = 0.2095
    CSK3 = -2.157
    CSK4 = 1.263

    if loge<6:
        y = DPP*10**(KPP*(PIP0+PIP1*loge+PIP2*loge**2))
    else:
        y = np.exp(CSK1 + CSK2*loge)+CSK3+CSK4*loge
    return y

def sigma_baryon(energy):
    loge = np.log10(energy)
    #get the fit constants
    PIP0 = 1.35485692
    PIP1 = -0.02237097
    PIP2 = 0.01398568
    KPP = 0.94044349
    LPP = 0.9600303037501511

    CSK1 = 2.269
    CSK2 = 0.207
    CSK3 = -0.9907
    CSK4 = 1.277

    if loge<6:
        y = LPP*10**(KPP*(PIP0+PIP1*loge+PIP2*loge**2))
    else:
        y = np.exp(CSK1 + CSK2*loge)+CSK3+CSK4*loge
    return y


def set_event_kinematics(proj, energy, target):
    if type(target)==int:
        ev_kin = EventKinematics(elab = energy * GeV,
            p1pdg = proj,
#            nuc2_prop = target)
            p2pdg = target)
    else:
        ev_kin = EventKinematics(elab = energy * GeV,
            p1pdg = proj,
            nuc2_prop = target)

    return ev_kin

def sigma_Dp(energy, generator):
    DPP = 0.6979928863374517
    CSK1 = 1.891
    CSK2 = 0.2095
    CSK3 = -2.157
    CSK4 = 1.263
    loge = np.log10(energy)

    if energy<1e6:
        generator.set_event_kinematics(set_event_kinematics(321, energy, 2212))
        xs = DPP*generator.sigma_inel()
    else:
        xs = np.exp(CSK1 + CSK2*loge)+CSK3+CSK4*loge
    return xs


def sigma_Lp(energy, generator):
    LPP = 0.9761791227127833
    CSK1 = 2.269
    CSK2 = 0.207
    CSK3 = -0.9907
    CSK4 = 1.277
    loge = np.log10(energy)

    if energy<1e6:
        generator.set_event_kinematics(set_event_kinematics(321, energy, 2212))
        xs = LPP*generator.sigma_inel()
    else:
        xs = np.exp(CSK1 + CSK2*loge)+CSK3+CSK4*loge
    return xs

######## Interpolated charm hadron - ice interaction cross section
dloc = '/data2/icecube/ssarkar/charm_muon_data/xs_fraction/charm_ice_xs.txt'
charm_ice_data = np.loadtxt(dloc, comments='#')
loge = np.log10(charm_ice_data[:,0])
logx1 = np.log10(charm_ice_data[:,1])
logx2 = np.log10(charm_ice_data[:,2])
dicell = ip.interp1d(loge, logx1, kind='slinear')
licell = ip.interp1d(loge, logx2, kind='slinear')
def sigma_Dice(energy):
    return 10**(dicell(np.log10(energy)))
def sigma_Lice(energy):
    return 10**(licell(np.log10(energy)))
###################################################################


class CharmMuonGenerator:
    '''
    Generates secondary charm muon from parametrized
    datasets (generated from PYTHIA and Sibyll2.3d (impy))
    and computes the probability of the event from
    charm production cross-section and muon branching
    ratio
    '''
    def __init__(self, decay_fractions, interaction_fractions,
            decay_br, interaction_br, charm_xs,
            seed):
        '''Intializes all the interpolation functions to be used
        in the sampling processes.
        Inputs: decay_fractions (list) : list of files containing
                the binned distributions of muon energy fractions
                (should have 11 files for 11 charmed hadrons)

                interaction_fractions (dict) : each key has list of 
                files (two keys: 'OX' and 'H' for interactions)
                containing the binned distributions of muon energy
                fractions from interaction with Oxygen anf Hydrogen

                seed (int) : Random seed to be used for sampling
        '''
        #--------------Energy Sampling Interpolations--------------
        #get the energy fraction cdfs to sample from
        self.decay_cdfs     = []
        self.decay_inv_cdfs = []
        for fname in decay_fractions:
            x,y = get_decay_cdf(fname)
            cdf = ip.interp1d(x,y,kind='linear',
                    bounds_error=False, fill_value='extrapolate')
            inv_cdf = ip.interp1d(y,x, kind='linear',
                    bounds_error=False, fill_value='extrapolate')

            self.decay_cdfs.append(cdf)
            self.decay_inv_cdfs.append(inv_cdf)

        interaction_Ox_2dcdfs     = []
        interaction_H_2dcdfs     = []
        for fname in interaction_fractions['OX']:
            func = get_interaction_cdf(fname)
            interaction_Ox_2dcdfs.append(func)
        for fname in interaction_fractions['H']:
            func = get_interaction_cdf(fname)
            interaction_H_2dcdfs.append(func)
        #store the functions in a dictionary
        self.interaction_2Dcdfs = {}
        self.interaction_2Dcdfs['OX'] = interaction_Ox_2dcdfs
        self.interaction_2Dcdfs['H']  = interaction_H_2dcdfs
        
        #--------------Random Number generator--------------------
        #initialize the ssample generator
        self.rng = default_rng(seed=seed)

        #--------------Muon Energy Threshold Setup----------------
        #set the muon energy threshold for sampling
        self.EMU_MIN = 10.0#GeV
        #minimum fraction value in the binned data
        XBINS_DECAY = np.linspace(0,1,101)
        self.XMIN_DECAY = 0.5*(XBINS_DECAY[:-1]+XBINS_DECAY[1:])[0]

        XBINS_INT = np.logspace(-8, 0, 101)
        self.XMIN_INT = 0.5*(XBINS_INT[:-1]+XBINS_INT[1:])[0]

        #--------------Decay/Interaction Length Setup-------------
        #Define the hadron masses and decay lifetime
        #taken from pdg tables
        #charm hadron mass in GeV
        self.MASS = [1.86966, 1.86966, 1.86484, 1.86484,
                1.96835, 1.96835, 2.46771, 2.46771,
                2.28646, 2.47044, 2.6952]
        #charm hadron mean decay lifetimes in seconds
        self.TAU = [1040e-15, 1040e-15, 410.1e-15, 410.1e-15,
                504e-15, 504e-15, 456e-15, 456e-15,
                202.4e-15, 153e-15, 268e-15]
        #speed of light in cm/s
        self.C = 2.99792458e10
        #define the nuclear interaction length in ice
        #LambdaI = 83.3#gm/cm^2
        #Rho     = 0.917#gm/cm^3
        #self.ILENGTH = LambdaI/Rho

        #ice molecular mass
        self.MICE = 18.02 #gm/mol
        #ice density
        self.RHO = 0.917 #gm/cm^3
        #avogadro number times mb to cm^2
        self.NA = 6.02214076e-4 #cm^2/mb/mol

        #--------------Muon BR Setup------------------------------
        #get interpolation function of decay BR
        self.decay_BR = []
        for fname in decay_br:
            func = br_intepolation(fname)
            self.decay_BR.append(func)

        #get interpolations for interaction BR
        self.interaction_BR = {'OX':[], 'H':[]}
        for fname in interaction_br['OX']:
            func = br_intepolation(fname)
            self.interaction_BR['OX'].append(func)
        for fname in interaction_br['H']:
            func = br_intepolation(fname)
            self.interaction_BR['H'].append(func)

        #--------------Charm Fractional XS Setup-----------------
        charm_data = np.loadtxt(charm_xs, comments='#')
        loge       = np.log10(charm_data[:,0])
        nu_xs      = charm_data[:,1]
        nubar_xs   = charm_data[:,2]

        func1 = ip.interp1d(loge, nu_xs, kind='slinear',
                bounds_error=False, fill_value='extrapolate')
        func2 = ip.interp1d(loge, nubar_xs, kind='slinear',
                bounds_error=False, fill_value='extrapolate')
        self.charm_xs = [func1, func2]


    def sample_decay_fraction(self, idx, energy, nsample=1):
        ''' Sample muon energy from charm hadron decays.
        Input: idx (int) : charm hadron identifier
               energy (float) : energy of the charm hadron (in GeV)
               (for the purpose of setting the sampling threshold)
        Returns: Sampled muon energy (floats)
        '''
        x_min = max(self.EMU_MIN/energy, self.XMIN_DECAY)
        assert x_min<=1., "Minimum energy fraction above the hadron energy"
        cdf_min = self.decay_cdfs[idx](x_min)

        cdf_sample = self.rng.random(int(nsample))*(1-cdf_min)+cdf_min
        x_sample = self.decay_inv_cdfs[idx](cdf_sample)
        return x_sample*energy
        #return x_sample, cdf_min

    def sample_interaction_fraction(self, target, idx, energy, nsample=1):
        ''' Sample muon energy from charm hadron interaction in Ice.
        Input: target (str) : target type to sample from ('OX' or 'H')
               idx (int) : charm hadron identifier
               energy (float) : energy of the charm hadron (in GeV)
        Returns: Sampled muon energy (floats)
        '''
        x_min = max(self.EMU_MIN/energy, self.XMIN_INT)
        assert x_min<=1., "Minimum energy fraction above the hadron energy"

        #get sample cdf points from 2d interpolation
        log_xmin = convert_linlog(x_min)
        grid_frac   = np.linspace(log_xmin, 0, 100)
        grid_energy = np.repeat(np.log10(energy), len(grid_frac))
        grid_cdf    = self.interaction_2Dcdfs[target][idx].ev(grid_frac, grid_energy)
        cdf_min = min(grid_cdf)

        #generate inverse cdf on-the-fly
        inv_cdf = ip.interp1d(grid_cdf, grid_frac, kind='linear',
                bounds_error=False, fill_value = 'extrapolate')
        
        cdf_sample = self.rng.random(int(nsample))*(1-cdf_min)+cdf_min
        x_sample = inv_cdf(cdf_sample)#in log-scale
        return 10**x_sample*energy
        #return x_sample, cdf_min

    def get_decay_interaction_length(self, idx, energy):
        '''Checking the decay and interaction length for
        each charm hadron
        '''
        DL = energy*self.TAU[idx]*self.C/self.MASS[idx]
        #get the interaction length sample
        if idx<6:#mesons
            #xs = sigma_meson(energy)
            #xs = sigma_Dp(energy, self.generator)
            xs = sigma_Dice(energy)
        elif idx>=6:#hadrons
            #xs = sigma_baryon(energy)
            #xs = sigma_Lp(energy, self.generator)
            xs = sigma_Lice(energy)
        IL = self.MICE/(self.RHO*self.NA*xs)

        return (DL, IL)

    def sample_decay_interaction(self, idx, energy, nsample=1):
        ''' Samples the fate of the Charm hadrons between
        decay and interaction in Ice.
        Input: idx (int) : Charm Hadron identifier
               energy (float) : Energy of the Charm hadron
        Returns: Interaction Type (int): 1 for decay, 2 for interaction
                 Target Type (str): 'OX' for Oxygen, 'H' for Hydrogen
        '''
        #get the decay length sample
        DL = energy*self.TAU[idx]*self.C/self.MASS[idx]
        cdf_sample = self.rng.random(nsample)
        #inverse cdf if decay length distribution
        Dx = -1*DL*np.log(1-cdf_sample)
        
        if energy<100:
            return 1, None
        #get the interaction length sample
        if idx<6:#mesons
            #xs = sigma_meson(energy)
            xs = sigma_Dice(energy)
        elif idx>=6:#hadrons
            #xs = sigma_baryon(energy)
            xs = sigma_Lice(energy)

        IL = self.MICE/(self.RHO*self.NA*xs)
        cdf_sample = self.rng.random(nsample)
        #inverse cdf of interaction length distribution
        Ix = -1*IL*np.log(1-cdf_sample)

        ###### DEBUG/TEST
        #Dx=0
        #Ix=1
        ######

        if Dx<=Ix: return 1, None
        elif Dx>Ix:
            #in case of interaction, sample the target as well
            #1: for hydrogen, 2: for Oxygen
            h2o = np.array([1,1,2])
            target = self.rng.choice(h2o)
            if target==1:
                return 2, 'H'
            elif target==2:
                return 2, 'OX'

    def sample_muon(self, idx, energy):
        '''This is a wrapper function that calls for the above
        sample functions to get secondary muon energy and other
        interaction info. (decay/interaction, target_type,
        fractional XS, Muon BR)

        Input : idx (int) : charm hadron identifier
                energy (float) :  charm hadron energy
        Returns : Muon energy (float)
                  Origin (int): 1:Decay, 2:Interaction
                  Target Type (str): 'OX':2, 'H':1, decay: None
                  Muon BR (float)
        '''
        #if charm hadron energy is below 10GeV, do nothing
        if energy<=10.: 
            print (f'\rCharm Hadron energy:{energy} below 10GeV, skipping Muon Generation',
                    end='\r')
            return None
        #if energy is below 100GeV, always decay be default
        elif energy>10. and energy<100.:
            br = self.decay_BR[idx](np.log10(energy))
            #if branching ratio is 0 or negative, do nothing
            if br<=0.0: return None
            muen = self.sample_decay_fraction(idx,energy)
            orig = 1
            ttype = None
            return (float(muen), orig, ttype, br)
        else:
            #decision on decay vs. interaction
            orig, ttype = self.sample_decay_interaction(idx, energy)
            if orig==1:
                muen = self.sample_decay_fraction(idx, energy)
                br = self.decay_BR[idx](np.log10(energy))
            elif orig==2:
                muen = self.sample_interaction_fraction(ttype, idx, energy)
                br = self.interaction_BR[ttype][idx](np.log10(energy))
            return (float(muen), orig, ttype, br)

    def compute_fractional_xs(self, energy, nutype):
        '''Computes the charm production fractional cross-section
        for the given neutrino energy and neutrino type
        '''
        if nutype == 1:
            return self.charm_xs[0](np.log10(energy))
        elif nutype == -1:
            return self.charm_xs[1](np.log10(energy))

    def inject_muon(self, idx, energy):
        '''Decides if there should be muon as decay product for
        events with more than one charm hadron (applicable to
        2nd,3rd... charm hadrons)
        '''
        if energy<=10.: return False
        orig, ttype = self.sample_decay_interaction(idx, energy)
        rand = self.rng.random(1)
        if orig==1:
            br = self.decay_BR[idx](np.log10(energy))
            if rand<=br: return True
            else: return False
        elif orig==2:
            br = self.interaction_BR[ttype][idx](np.log10(energy))
            if rand<=br: return True
            else: return False

if __name__=='__main__':
    pass
