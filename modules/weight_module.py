#!/usr/bin/env python

'''
Author: Sourav Sarkar
Date: 7 Dec, 2020
Description: This script contains the python functions to calculate weight information
	for trident events. Parameters needed for calculating the weights come from
	parent LeptonInjector I3 files, nuSQuIDs and trident total cross-section.
'''

import numpy as np
from icecube import icetray, dataio, dataclasses
#import LeptonWeighter as lw
#import nuSQUIDSpy as nsq
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

def trident_sigma(nuen):
    ''' Computes the total numu CC+NC cross section for
    oxygen atom. The function returns an average per atom
    cross section assuming water/ice medium.
    '''
    per_atom = 18./(3.*16)
    total = NuMu_CHDF(nuen)+NuMu_DISq(nuen)+NuMu_DISp(nuen)
    return total*nuen*per_atom


class GenerateWeight(icetray.I3Module):
    '''
    I3Module class to extract weight information from initial LI file (hdf5)
    and lic file and calculate final weights and store it into the frames.
    '''
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        #Parameter for HDF5 filename
        self.AddParameter('LIFile',
                'Name of the initial Lepton Injector File',
                None)
        self.AddParameter('ConfigFile',
                'Name of the intial LI config file',
                None)
        #parameter for indicating weight calculation for trident/standard cc interaction
        self.AddParameter('InteractionType',# 1 for Trident interaction, 0 for Standard interaction
                        'Interaction type for weight calculation',# 2 for charm muon production
                        None)                                     # 3 for charm production
        self.AddParameter('CharmFile', #hdf5 file containing charm weights
                        'File containing Charm Muon Generation Weight',
                        None)

        self.AddParameter('WeightName', #Weight frame object name
                        'Weight Dictionary object name in the frame',
                        'I3MCWeightDict')

    def Configure(self):
        #Load the event properties from LI File
        li_fname = self.GetParameter('LIFile')
        self.evprops = h5.File(li_fname, 'r')['RangedInjector0']['properties']
        #Load the Simulation properties for LIC file
        cfg_fname = self.GetParameter('ConfigFile')
        with open(cfg_fname,'r') as f:
            propdict = json.load(f)
        self.simprop = propdict['LeptonInjector']
        #Interaction type
        self.inttype = self.GetParameter('InteractionType')

        #Load the charm weights in case of charm/charm muon production
        cname = self.GetParameter('CharmFile')
        if self.inttype in [2,3]:
            if cname is None:
                assert False, "CharmFile cannot be None for charm production events."
            #load the charm muon generation info
            self.charm_factor = h5.File(cname, 'r')['CharmWeights']
            #initialize the charm weight dictionary keys
            self.charm_keys = ['CharmHadronPDG', 'CharmHadronEnergy',
                    'CharmQuarkFraction', 'MuonBR', 'TotalFraction']

        #Initialize the WeightDictionary keys
        self.dict_keys = ['NEvents', 'Emin', 'Emax', 'ZenithMin', 'ZenithMax',
                'AzimuthMin', 'AzimuthMax', 'PowerlawIndex', 'InjectionRadius',
                'EndcapLength', 'PrimaryType', 'PrimaryEnergy', 'PrimaryZenith',
                'PrimaryAzimuth', 'ColumnDepth', 'TotalCrossection',
                'GenerationWeight', 'VolumeWeight', 'InteractionWeight',
                'OneWeight']
        
        #Initialize the event independent properties
        self.const_prop = [float(self.simprop['event_number']),
                self.simprop['MinEnergy'], self.simprop['MaxEnergy'],
                self.simprop['MinZenith']*np.pi/180., self.simprop['MaxZenith']*np.pi/180.,
                self.simprop['MinAzimuth']*np.pi/180., self.simprop['MaxAzimuth']*np.pi/180.,
                self.simprop['gamma'], self.simprop['injRadius'], self.simprop['endLength']]
        #Compute the weights for all events and store
        self.compute_all_weights()

        #Get the frame weight object name
        self.weightname = self.GetParameter('WeightName')

        #self.iterator = 0
        #self.nevents = self.simprop.events

    #Calculate the generation probability
    def pgen(self,nuen):
        #get the parameters
        self.index=self.simprop['gamma']
        self.emin=self.simprop['MinEnergy']
        self.emax=self.simprop['MaxEnergy']
        #case: monoenergetic events
        if self.emin==self.emax:
            return 1.
        #case: handle the integration diferently for E^-1 spectrum
        if self.index==1.:
            return nuen**(-1.*self.index)/(np.log10(self.emax)-np.log10(self.emin))
        #For all other powerlaw spectrum follow general calculation
        else:
            return (nuen**(-1.*self.index)*(self.index-1.))/(self.emin**(1-self.index)-self.emax**(1-self.index))

    #Calculate the injection area times solid angle
    def sim_volume(self):
        #get the parameters
        self.injrad=self.simprop['injRadius']*1e2 #convert from meter to cm
        self.zmax=self.simprop['MaxZenith']*np.pi/180.0 #convert from degree to rad
        self.zmin=self.simprop['MinZenith']*np.pi/180.0 #convert from degree to rad
        self.amin=self.simprop['MinAzimuth']*np.pi/180.0 #convert from degree to rad
        self.amax=self.simprop['MaxAzimuth']*np.pi/180.0 #convert from degree to rad
        #injection area
        injarea=np.pi*self.injrad**2
        #Solid angle used for simulation
        solidangle=(self.amax-self.amin)*(np.cos(self.zmin)-np.cos(self.zmax))
        return injarea*solidangle

    #Calculate the interaction weight within generation volume
    def interaction_weight(self,nuen, ptype, columndepth):
        #tridents
        if self.inttype==1:
            self.xs = trident_sigma(nuen) #in cm^2 unit
        #nu cc dis
        elif self.inttype==0:
            self.xs = dis_sigma(nuen, ptype) #in cm^2 unit
        #nu cc dis charm muon
        elif self.inttype in [2,3]:
            self.xs = dis_sigma(nuen, ptype)*self.charm_factor['totalXsFraction'] #in cm^2 unit
        else:
            assert False, f'Invalide Interaction type set: {self.inttype}'
        
        cd=columndepth * CONSTANT['avogadro_num']['value'] #in #targets/cm^-2 unit
        exp_factor=-1.*self.xs*cd
        val=1-np.exp(exp_factor)
        #if the exponent is too small and rounding error returns 0
        zero_idx = np.where(val==0)[0]
        val[zero_idx] = -1.*exp_factor[zero_idx]

        return val

    def compute_all_weights(self):
        #primary energy array
        p_en = self.evprops['totalEnergy']
        #total column depth
        col = self.evprops['totalColumnDepth']
        #primary particle type
        p_type = np.sign(self.evprops['initialType'])

        #compute the weights
        self.wgen = self.pgen(p_en)
        self.wvol = self.sim_volume()
        self.wint = self.interaction_weight(p_en, p_type, col)

        self.onew = self.wint*self.wvol/self.wgen
        return None

    #Function to process Q frames (storing weights)
    def DAQ(self,frame):
        #get the event index to extract all info from arrays
        event_idx  = frame['I3EventHeader'].event_id
        #get the event properties
        event_prop = self.evprops[event_idx]
        #store the relevant event properties
        event_param = [float(event_prop['initialType']),
                event_prop['totalEnergy'], event_prop['zenith'],
                event_prop['azimuth'], event_prop['totalColumnDepth']]
        #store the event weights
        event_weight = [self.xs[event_idx], self.wgen[event_idx],
                self.wvol, self.wint[event_idx],
                self.onew[event_idx]]
        #concatenate all the weight dictionary info
        wparam = self.const_prop + event_param + event_weight
        #create the key-value dictionary
        wdict = dict(zip(self.dict_keys, wparam))
        wdict_obj = dataclasses.I3MapStringDouble(wdict)

        #store into frame weight object
        frame.Put(self.weightname,wdict_obj)

        #in case of charm/charm muon production,
        #additionally store charm weights
        if self.inttype in [2,3]:
            wcharm = self.charm_factor[event_idx]
            charm_dict = dict(zip(self.charm_keys,
                [float(wcharm['hadronPDG']), wcharm['hadronEnergy'],
                    wcharm['xsFraction'], wcharm['branchingRatio'],
                    wcharm['totalXsFraction']]))
            cdict_obj = dataclasses.I3MapStringDouble(charm_dict)
            frame.Put('CharmWeightDict', cdict_obj)
        self.PushFrame(frame)

