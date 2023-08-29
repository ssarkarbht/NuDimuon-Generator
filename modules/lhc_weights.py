#!/bin/python

'''
Author: Sourav Sarkar
Description: This script deals with the MC event weights
    for LHC neutrino collider experiments like FASERnu and
    SND@LHC
'''

import numpy as np
import scipy.interpolate as ip
import h5py as h5
import json
import os, sys

#define the directory paths
repo_dir = os.environ["DIMUON_REPO"]
data_dir = repo_dir + "/data/"

#Load constants
cfile = data_dir + "constants_particles_materials/constants.json"
with open(cfile, 'r') as f:
    CONSTANT = json.load(f)

#Load mediums
mfile = data_dir + "constants_particles_materials/medium_properties.json"
with open(mfile, 'r') as f:
    MEDIUM = json.load(f)

#>>>>>>>>>>>>>>>>>>>>>>> Cross-section data >>>>>>>>>>>>>
#get the files containing the cross-section data
xsdir = data_dir + 'cross_sections/'
# DIS Cross-sections
dis_dir = xsdir + 'dis_cross_section/DAT/'
nu_cc_iso = np.loadtxt(dis_dir+'total_nu_CC_iso_NLO_HERAPDF1.5NLO_EIG.dat',
                        skiprows=1)
nu_nc_iso = np.loadtxt(dis_dir+'total_nu_NC_iso_NLO_HERAPDF1.5NLO_EIG.dat',
                        skiprows=1)

nubar_cc_iso = np.loadtxt(dis_dir+'total_nubar_CC_iso_NLO_HERAPDF1.5NLO_EIG.dat',
                        skiprows=1)
nubar_nc_iso = np.loadtxt(dis_dir+'total_nubar_NC_iso_NLO_HERAPDF1.5NLO_EIG.dat',
                        skiprows=1)
dis_energies = np.linspace(1,12,111)

# Trident Cross-sections
tri_dir = xsdir + 'trident_cross_section/'
#following three are the trident total cross-section data (same for Nu/NuBar)
chdf_xs=np.loadtxt(tri_dir+'zb_CH+DF.csv',
                    delimiter=',',comments='#')
disq_xs=np.loadtxt(tri_dir+'zb_DIS_quark.csv',
                    delimiter=',',comments='#')
disp_xs=np.loadtxt(tri_dir+'zb_DIS_photon.csv',
                    delimiter=',',comments='#')
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# Build a generic class for cross sections
class SigmaCalc:
    ''' Generic class for cross-section calculation.
    '''
    def __init__(self, x,y,
            kind="cubic", scale="null_null", unit='mbTOcm2'):
        ''' Takes the log of the energy and cross section values (mb)
        to return an interpolated function.
        Input: x (ndarray): energy array
               y (ndarray): cross section array
               kind (str, optional): interpolation order
               scale (str, optional): change the input x,y scales before
                                performing interpolations.
        '''
        if scale=="log_null":
            x = np.log10(x)
        elif scale=="log_log":
            x = np.log10(x)
            y = np.log10(y)
        elif scale=="null_log":
            y = np.log10(y)

        self.function = ip.interp1d(x, y, kind=kind)
        self.scale = scale
        if unit is None:
            self.conv_factor = 1.0
        else:
            self.conv_factor = CONSTANT[unit]

    def __call__(self, energy):
        ''' Input: energy (double or ndarray): energy for cross-section
                                            evaluation
            Returns: Cross-section value (in mb)
        '''
        outscale = self.scale.split("_")[-1]
        if outscale=="null":
            return self.function(np.log10(energy))*self.conv_factor
        elif outscale=="log":
            return 10**(self.function(np.log10(energy)))*self.conv_factor
        else: assert False, "Unknown cross-section scaling type"

#================== Prepare the cross section functions =================
Nu_CC    = SigmaCalc(dis_energies, nu_cc_iso)
Nu_NC    = SigmaCalc(dis_energies, nu_nc_iso)
NuBar_CC = SigmaCalc(dis_energies, nubar_cc_iso)
NuBar_NC = SigmaCalc(dis_energies, nubar_nc_iso)

def dis_sigma(nuen, nutype, inttype='cc'):
    ''' takes the neutrino energy, nutype= 1 for nu, -1 for nubar
    and inttype='cc' or 'nc'. returns total cross-section in cm^2
    '''
    xsarr = np.zeros(len(nuen))
    if inttype=='cc':
        partid = np.where(nutype==1)[0]
        antipartid = np.where(nutype==-1)[0]
        xsarr[partid] = Nu_CC(nuen[partid])
        xsarr[antipartid] = NuBar_CC(nuen[antipartid])
    elif inttype=='nc':
        partid = np.where(nutype==1)[0]
        antipartid = np.where(nutype==-1)[0]
        xsarr[partid] = Nu_NC(nuen[partid])
        xsarr[antipartid] = NuBar_NC(nuen[antipartid])
    else:
        assert False, "Unknown attributes!"
    return xsarr

NuMu_CHDF = SigmaCalc(chdf_xs[:,0], chdf_xs[:,1],
                    kind='slinear', scale='log_log', unit=None)
NuMu_DISq = SigmaCalc(disq_xs[:,0], disq_xs[:,1],
                    kind='slinear', scale='log_log', unit=None)
NuMu_DISp = SigmaCalc(disp_xs[:,0], disp_xs[:,1],
                    kind='slinear', scale='log_log', unit=None)

def trident_sigma(nuen, medium='ice'):
    ''' Computes the total numu CC+NC cross section for
    oxygen atom. The function returns an average per atom
    cross section assuming water/ice medium.
    '''
    mprop = MEDIUM[medium]
    if mprop['composite'] is not None:
        num_nuc = sum(mprop['ratio'])
    else:
        num_nuc = 1.
    num_mass = mprop['mass_num']
    per_atom = num_mass/(num_nuc * 16.) #deriving from Oxygen Cross-section
    total = NuMu_CHDF(nuen)+NuMu_DISq(nuen)+NuMu_DISp(nuen)
    return total*nuen*per_atom


class LHCWeight:
    '''
    Calculates the fluxless MC event weights for
    neutrino collider like configuration of the
    detectors.
    '''
    def __init__(self, medium,
            detector_length,
            event_config):
        self.medium = medium
        self.med_prop = MEDIUM[medium]
        self.det_length = detector_length
        self.config = event_config

    def pgen(self, nuen):
        '''Calculates the MC event generation probability
        from the injected spectrum.
        '''
        index = self.config['gamma']
        emin = self.config['MinEnergy']
        emax = self.config['MaxEnergy']

        #in case of generating monoenergetic events
        if emin == emax:
            return 1.
        #treating E^-1 spectrum differently
        if index==1.:
            return nuen**(-1.*index)/(np.log10(emax)-np.log10(emin))
        else:
            return (nuen**(-1.*index)*(index-1.))/(emin**(1-index)-emax**(1-index))

    def interaction_prob(self, nuen, nutype):
        '''Calculates the MC event interaction probability
        inside the detector volume.
        '''
        totalXS = dis_sigma(nuen, nutype) * self.med_prop['mass_num']
        wint = (totalXS * self.med_prop['density'] * self.det_length *
                CONSTANT['avogadro_num']['value']) / self.med_prop['mole_mass']

        return wint

    def get_weight(self, nuen, nutype):
        int_prob = self.interaction_prob(nuen, nutype)
        gen_prob = self.pgen(nuen)

        #calculate the normalization factor
        event_num = self.config['event_number']
        #if the dataset is a mix of nu+nubar, half the normalization factor
        if len(self.config['nu_pdg'])==2:
            event_num *= 0.5
        
        return int_prob
        #return int_prob/(gen_prob * event_num)


class FluxCalc:
    '''Generic class for calculating fluxes.
    '''
    def __init__(self, fluxfile, fluxtype):
        '''Input:
        fluxfile (str) : compressed file containing numpy arrays
        fluxtype (str) : type of the flux to calculate
                (e.g. numu_total for total neutrino flux)
        '''
        fluxes = np.load(fluxfile)

        energy = fluxes['energies']
        flux = fluxes[fluxtype]

        #build the interpolation
        self.function = ip.interp1d(np.log10(energy),
                                    np.log10(flux),
                                    kind='linear',
                                    bounds_error = False,
                                    fill_value = 0.0)

    def __call__(self, nuen):
        return 10**(self.function(np.log10(nuen)))
