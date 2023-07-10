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
from glob import glob
import json

#define the directory paths
repo_dir = os.environ["DIMUON_REPO"]
data_dir = repo_dir + "/data/"

#load the constants
cfile = data_dir + "constants_particles_materials/constants.json"
with open(cfile, 'r') as f:
    CONSTANT = json.load(f)

#particle list
pfile = data_dir + "constants_particles_materials/particles.json"
with open(pfile, 'r') as f:
    PARTICLE_DICT = json.load(f)
#get list of charm hadron pdg codes
pdglist = []
for key, val in PARTICLE_DICT.items():
    pdgid = val['pdgid']
    pdglist.append(pdgid)
#get a reverse dictionary for charm hadron names
names = list(PARTICLE_DICT.keys())
PDGNAMES = dict(zip(pdglist, names))
NAMESPDG = dict(zip(names, pdglist))

#material list
mfile = data_dir + "constants_particles_materials/medium_properties.json"
with open(mfile, "r") as f:
    MEDIUM = json.load(f)

#Charm Hadron - Nucleus Cross Section list
INTXS_DICT = {}
xsdir = data_dir + "charm_data/cross_sections/"
xsfiles = glob(xsdir+"*.dat")
for med in MEDIUM.keys():
    medname = "charm_"+med+"_sigma.dat"
    fname = xsdir+medname
    exist = fname in xsfiles
    if exist:
        INTXS_DICT[med] = fname
    else:
        print (f"WARNING! {medname} does not exist in the given location.")

class CharmSigma:
    ''' Load the charm interaction cross section data
    '''
    def __init__(self, x,y):
        logx = np.log10(x)
        logy = np.log10(y)
        self.func = ip.interp1d(logx, logy, kind='slinear')
    def __call__(self, energy):
        return 10**(self.func(np.log10(energy)))

def get_interaction_sigmas(xfile):
    ''' Returns the charm meson and baryon cross-section functions
    '''
    data = np.loadtxt(xfile, comments="#")
    return CharmSigma(data[:,0], data[:,1]), CharmSigma(data[:,0], data[:,2])


#Useful functions to be used in the main class
def convert_linlog(x, to='log'):
    ''' Converts the input to log or linear scale
    '''
    if to=='log':
        return np.log10(x)
    elif to=='lin':
        return 10**x

# Create Muon fractional energy (decay) CDF
# fractional energy in linear scale
def get_decay_cdf(data):
    '''Returns discrete cdf values at
    discrete energy fraction points
    '''
    xvals = data[0]
    yvals = data[1:]#for 350 hadron energy points
    #compute the average of the distributions
    #as the decay energy distribution for all energies
    #are alomst the same
    y_avg = np.average(yvals, axis=0)

    #get the bin width
    nbins = len(xvals)
    XBINS = np.linspace(0,1, nbins+1)#bins used in the dataset
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

def get_interaction_cdf(data, en_points):
    ''' Returns 2d interpolation function
    to compute cdf values from (at different
    hadron energies)
    '''
    #hadron energy simulation points in the dataset
    EBINS = en_points
    #muon energy fraction bin edges computed in the
    #dataset
    xvals = data[0]
    yvals = data[1:]
    nbins = len(xvals)
    XBINS = np.logspace(np.log10(1/EBINS[-1]),0,nbins+1)
    XLOGBINS = np.linspace(np.log10(1/EBINS[-1]),0,nbins+1)
    #binwidth in logscale
    LOG_BIN_WIDTH = np.diff(XLOGBINS)[0]

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
    return xvals, func

def br_interpolation(data):
    '''Builds the muon branching ratio for the 
    corresponding charm hadron and simulation type
    '''
    x = data[:,0]
    y = data[:,1]/data[:,2]
    logx = np.log10(x)
    func = ip.interp1d(logx, y, kind='slinear')
    return func

class CharmMuon:
    ''' Loads the muon energy fraction cdfs and BR
    '''
    def __init__(self, fname, simtype):
        self.data = np.load(fname)
        self.energy = self.data['EnergyPoints']
        self.type = simtype
        self.cdfs = {}
        self.inv_cdfs = {}
        self.brs = {}

    def build_brs(self, particles):
        ''' Takes the list of particles and builds the
        branchign ratio interpolator
        '''
        for part in particles:
            try:
                data = self.data[part+"_branchingRatio"]
            except:
                dummy = PDGNAMES[abs(NAMESPDG[part])]
                data = self.data[dummy+"_branchingRatio"]
                print (f"{part} data missing, using {dummy} instead")
            fn = br_interpolation(data)
            self.brs[part] = fn
        return None

    def build_decay_cdfs(self, particles):
        ''' Takes the muon energy fraction data and builds
        the interpolation functions for cdfs and inverse
        cdfs
        '''
        for part in particles:
            try:
                data = self.data[part+"_energyFraction"]
            except:
                dummy = PDGNAMES[abs(NAMESPDG[part])]
                data = self.data[dummy+"_energyFraction"]
                print (f"{part} data missing, using {dummy} instead")
            x,y = get_decay_cdf(data)
            cdf_fn = ip.interp1d(x,y,kind='linear',
                    bounds_error=False, fill_value='extrapolate')
            inv_cdf_fn = ip.interp1d(y,x, kind='linear',
                    bounds_error=False, fill_value='extrapolate')
            self.cdfs[part] = cdf_fn
            self.inv_cdfs[part] = inv_cdf_fn
        self.fraction = x
        return None

    def build_interaction_cdfs(self, particles):
        ''' Takes the muon energy fraction data and builds
        the 2D interpolation function for interaction cdfs
        '''
        for part in particles:
            try:
                data = self.data[part+"_energyFraction"]
            except:
                dummy = PDGNAMES[abs(NAMESPDG[part])]
                data = self.data[dummy+"_energyFraction"]
                print (f"{part} data missing, using {dummy} instead")
            x,cdf_fn = get_interaction_cdf(data, self.energy)

            self.cdfs[part] = cdf_fn
        self.fraction = x
        return None

    def build(self, particles):
        ''' This is a wrapper function that calls the decay/interaction
        cdf builders accrodingly
        '''
        self.build_brs(particles)
        if self.type == "interaction":
            self.build_interaction_cdfs(particles)
        elif self.type == "decay":
            self.build_decay_cdfs(particles)
        else:
            print (f"Unknown simulation type {self.type}, building nothing.")
        return None

class CharmMuonGenerator:
    '''
    Generates secondary charm muon from parametrized
    datasets (generated from PYTHIA and Sibyll2.3d (chromo))
    and computes the probability of the event from
    charm production cross-section and muon branching
    ratio
    '''
    def __init__(self, param_dir, prod_xsfile,
            seed, target='ice', mu_emin=10.0):
        '''Intializes all the interpolation functions to be used
        in the sampling processes.
        Inputs: param_dir (str) : directory path containing the
                parametric tables

                xsfile (str) : path to the charm hadron interaction
                cross section file

                seed (int) : Random seed to be used for sampling

                target (str) : Name of the target medium in which to
                sample the charm muons
        '''
        #--------------Decay Interpolations--------------
        decay_file = param_dir + "/decay_table.npz"
        self.decayMuon = CharmMuon(decay_file, 'decay')
        #build the interpolation functions
        self.decayMuon.build(NAMESPDG.keys())

        #--------------Interaction Interpolations--------
        self.interactionDict = {}
        self.atom_arr = []
        if MEDIUM[target]['composite'] is not None:
            for i, atom in enumerate(MEDIUM[target]['composite']):
                int_file = param_dir + "/" + atom + "_table.npz"
                intMuon = CharmMuon(int_file, 'interaction')
                intMuon.build(NAMESPDG.keys())
                self.interactionDict[atom] = intMuon
                self.atom_arr += [atom for _ in range(MEDIUM[target]['ratio'][i])]

        elif MEDIUM[target]['composite'] is None:
            int_file = param_dir + "/" + target + "_table.npz"
            intMuon = CharmMuon(int_file, 'interaction')
            intMuon.build(NAMESPDG.keys())
            self.interactionDict[target] = intMuon
            self.atom_arr += [target]
        #--------------Random Number generator--------------------
        #initialize the ssample generator
        self.rng = default_rng(seed=seed)

        #--------------Muon Energy Threshold Setup----------------
        #set the muon energy threshold for sampling
        self.EMU_MIN = mu_emin #GeV

        #--------------Decay/Interaction Length Setup-------------
        self.MASS = {}
        self.TAU = {}
        for key, val in PARTICLE_DICT.items():
            self.MASS[val['pdgid']] = val['mass']
            self.TAU[val['pdgid']] = val['lifetime']
        # Target properties
        self.PROP = MEDIUM[target]

        #--------------Charm cross-section Setup-----------------
        #interaction cross-section
        func1, func2 = get_interaction_sigmas(INTXS_DICT[target])
        self.charm_intxs = [func1, func2]

        #charm quark production fraction
        prod_data = np.loadtxt(prod_xsfile, comments='#')
        loge = np.log10(prod_data[:,0])
        nu_xs = prod_data[:,1]
        nubar_xs = prod_data[:,2]
        frac1 = ip.interp1d(loge, nu_xs, kind='slinear',
                bounds_error=False, fill_value='extrapolate')

        frac2 = ip.interp1d(loge, nubar_xs, kind='slinear',
                bounds_error=False, fill_value='extrapolate')
        self.charm_prodfrac = [frac1, frac2]

    def sample_decay_fraction(self, idx, energy, nsample=1):
        ''' Sample muon energy from charm hadron decays.
        Input: idx (int) : charm hadron pdgcode
               energy (float) : energy of the parent charm hadron (GeV)
        Returns: Sampled muon energy (floats)
        '''
        x_min = max(self.EMU_MIN/energy, self.decayMuon.fraction[0])
        assert x_min<=1., "Minimum energy fraction above the hadron energy"
        
        cdf_min = self.decayMuon.cdfs[PDGNAMES[idx]](x_min)
        cdf_sample = self.rng.random(int(nsample))*(1-cdf_min)+cdf_min
        x_sample = self.decayMuon.inv_cdfs[PDGNAMES[idx]](cdf_sample)

        return x_sample*energy

    def sample_interaction_fraction(self, target, idx, energy, nsample=1):
        ''' Sample outgoing muon energy from charm hadron interaction in 
        the given target atom.
        Input: target (str) : target atom name
               idx (int) : charm hadron pdg code
               energy (float) : energy of the charm hadron (in GeV)
        Returns: Sampled muon energy (floats)
        '''
        target_cdf = self.interactionDict[target]
        x_min = max(self.EMU_MIN/energy, target_cdf.fraction[0])
        assert x_min<=1., "Minimum energy fraction above the hadron energy"

        #get sample cdf points from 2d interpolation
        log_xmin = convert_linlog(x_min)
        #generate a 1D grid fro building inverse CDF for the given energy
        nbins = len(target_cdf.fraction)
        grid_frac   = np.linspace(log_xmin, 0, nbins)
        grid_energy = np.repeat(np.log10(energy), len(grid_frac))
        grid_cdf    = target_cdf.cdfs[PDGNAMES[idx]].ev(grid_frac, grid_energy)
        cdf_min = min(grid_cdf)

        #generate 1D inverse cdf on-the-fly
        inv_cdf = ip.interp1d(grid_cdf, grid_frac, kind='linear',
                bounds_error=False, fill_value = 'extrapolate')
        
        cdf_sample = self.rng.random(int(nsample))*(1-cdf_min)+cdf_min
        x_sample = inv_cdf(cdf_sample)#in log-scale
        return 10**x_sample*energy

    def get_decay_interaction_length(self, idx, energy):
        '''Computing the decay and interaction length for
        the given charm hadron and energy
        '''
        #Decay length
        DL = energy*self.TAU[idx]*CONSTANT['light_speed']['value']/self.MASS[idx]
        #get the interaction length sample
        if abs(idx)<1000:#mesons
            xs = self.charm_intxs[0](energy)
        elif abs(idx)>=1000:#baryons
            #xs = sigma_baryon(energy)
            #xs = sigma_Lp(energy, self.generator)
            xs = self.charm_intxs[1](energy)
        #Interaction length
        IL = self.PROP['mole_mass']/(self.PROP['density']*CONSTANT['avogadro']['value']*xs)

        return (DL, IL)

    def sample_decay_interaction(self, idx, energy, nsample=1):
        ''' Samples the fate of the Charm hadrons between
        decay and interaction in Ice.
        Input: idx (int) : Charm Hadron PDG CODE
               energy (float) : Energy of the Charm hadron
        Returns: Interaction Type (int): 1 for decay, 2 for interaction
                 Target Type (str): 'OX' for Oxygen, 'H' for Hydrogen
        '''
        #By default, the charm hadron decays for energy below 100GeV
        if energy<100.:
            return 1, None
        
        #get the decay and interaction lengths
        DL, IL = self.get_decay_interaction_length(idx, energy)

        cdf_sample = self.rng.random(nsample)
        #sampled length from the decay length distribution
        Dx = -1*DL*np.log(1-cdf_sample)
        
        cdf_sample = self.rng.random(nsample)
        #sampled length from the interaction length distribution
        Ix = -1*IL*np.log(1-cdf_sample)

        #Case : Charm Hadron Decays
        if Dx<=Ix: return 1, None
        #Case : Charm Hadron interactions
        elif Dx>Ix:
            #sample a random atom from the medium composition
            target = self.rng.choice(self.atom_arr)
            return 2, target

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
        if energy<=self.EMU_MIN: 
            print (f'\rCharm Hadron energy:{energy} below {self.EMU_MIN} GeV, skipping Muon Generation',
                    end='\r')
            return None

        #decision on decay vs. interaction
        orig, ttype = self.sample_decay_interaction(idx, energy)
        if orig==1:
            muen = self.sample_decay_fraction(idx, energy)
            br = self.decayMuon.brs[PDGNAMES[idx]](np.log10(energy))
        elif orig==2:
            muen = self.sample_interaction_fraction(ttype, idx, energy)
            br = self.interactionDict[ttype].brs[PDGNAMES[idx]](np.log10(energy))
        return (float(muen), orig, ttype, br)

    def compute_fractional_xs(self, energy, nutype):
        '''Computes the charm production fractional cross-section
        for the given neutrino energy and neutrino type
        '''
        if nutype == 1:
            return self.charm_prodfrac[0](np.log10(energy))
        elif nutype == -1:
            return self.charm_prodfrac[1](np.log10(energy))

    def inject_muon(self, idx, energy):
        '''Decides if there should be muon as decay product for
        events with more than one charm hadron (applicable to
        2nd,3rd... charm hadrons)
        '''
        if energy<=10.: return False
        orig, ttype = self.sample_decay_interaction(idx, energy)
        rand = self.rng.random(1)
        if orig==1:
            br = self.decayMuon.brs[PDGNAMES[idx]](np.log10(energy))
            if rand<=br: return True
            else: return False
        elif orig==2:
            br = self.interactionDict[ttype].brs[PDGNAMES[idx]](np.log10(energy))
            if rand<=br: return True
            else: return False

if __name__=='__main__':
    pass
