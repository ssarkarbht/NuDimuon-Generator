#!/usr/bin/env python

'''
Author: Sourav Sarkar
Date: 7 Dec, 2020
Description: This script contains the python functions to calculate weight information
	for trident events. Parameters needed for calculating the weights come from
	parent LeptonInjector I3 files, nuSQuIDs and trident total cross-section.
'''

import numpy as np
from icecube import icetray, dataio, dataclasse
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
dis_dir = xs_dir + 'dis_cross_section/DAT/'
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
tri_dir = xs_dir + 'trident_cross_section/'
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
    if inttype=='cc' and nutype==1:
        return Nu_CC(nuen)
    elif inttype=='cc' and nutype==-1:
        return NuBar_CC(nuen)
    elif inttype='nc' and nutype==1:
        return Nu_NC(nuen)
    elif inttype='nc' and nutype==-1:
        return NuBar_NC(nuen)
    else:
        assert False, "Unknown attributes!"

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
                'Name of the initial Lepton Injector File',)
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

    def Configure(self):
        #Load the event properties from LI File
        li_fname = self.GetParameter('LIFile')
        self.evprops = h5.File(li_fname, 'r')['RangedInjector0']['properties']

        #Load the Simulation properties for LIC file
        cfg_fname = self.GetParameter('ConfigFile')
        self.simprop = lw.RangeSimulationDetails(cfg_fname)

        #Interaction type
        self.inttype = self.GetParameter('InteractionType')

        #Load the charm weights in case of charm/charm muon production
        cname = self.GetParameter('CharmFile')
        if self.inttype in [2,3]:
            if cname is None:
                assert False, "CharmFile cannot be None for charm production events."
            self.charm_factor = h5.File(cname, 'r')['CharmWeights']

        self.iterator = 0
        self.nevents = self.simprop.events

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
    def interaction_weight(self,nuen,columndepth):
        #tridents
        if self.inttype==1:
            self.xs = trident_sigma(nuen) #in cm^2 unit
        #nu cc dis
        elif self.inttype==0:
            self.xs = dis_sigma(nuen, self.ptype) #in cm^2 unit
        #nu cc dis charm muon
        elif self.inttype in [2,3]:
            self.xs = dis_sigma(nuen, self.ptype)*self.charmfactor['totalXsFraction'] #in cm^2 unit
        else:
            assert False, f'Invalide Interaction type set: {self.inttype}'
        
        cd=columndepth * CONSTANT['avogadro_num']['value'] #in #targets/cm^-2 unit
        exp_factor=-1.*self.xs*cd
        val=1-np.exp(exp_factor)
        if val==0.: #if the exponent is too small and rounding error returns 0
            return -1.*exp_factor
        return (1-np.exp(exp_factor))

    #calculate the neutrino propagation weight from atmosphere to the surface of generation volume
    def propagation_weight(self,nuen,zenith):
        numu=1
        units=nsq.Const()
        nuen=nuen*units.GeV
        nutype=int(0.5*(1-self.ptype))
        nu_state = self.nsq_atm.EvalFlavor(numu,np.cos(zenith),nuen,nutype)
        initial_flux=1.0e18*nuen**(-1.*1)##this spectral index should be the one with nsq comptation
        #initial_flux=1.0e18*nuen**(-1.*self.index)
        return nu_state/initial_flux

    #Function to process Q frames
    def DAQ(self,frame):
        #get the event dependent weight parameters
        eventprop=frame['EventProperties']
        #get event id to keep track of the weights in HDF5 file
        evtid=frame['I3EventHeader'].event_id
        #get primary neutrino pid to keep track of the weights in HDF5 file
        ppid=frame['I3MCTree'].primaries[0].id.minorID
        #extra factor for charm muon production
        if self.inttype==2:
            #self.charmfactor = frame['CharmMuonFactor'].value
            self.charmfactor = self.charm_factor[evtid]
        #get the neutrino/anti-neutrino type
        nuclass = [dataclasses.I3Particle.NuE,
                dataclasses.I3Particle.NuMu,
                dataclasses.I3Particle.NuTau]
        nubarclass = [dataclasses.I3Particle.NuEBar,
                dataclasses.I3Particle.NuMuBar, 
                dataclasses.I3Particle.NuTauBar]
        if eventprop.initialType in nuclass:
            self.ptype=1
        elif eventprop.initialType in nubarclass:
            self.ptype=-1
        else:
            print ("Check for bugs! Or upgrade particle type for weight dictionary")
        # Calculate generation probability
        wgen=self.pgen(eventprop.totalEnergy)
        # Calculate solid angle over which events are simulated
        wvol=self.sim_volume()
        # Calculate interaction probablity
        wint=self.interaction_weight(eventprop.totalEnergy,eventprop.totalColumnDepth)
        # Calculate propagation probability
        wprop=self.propagation_weight(eventprop.totalEnergy,eventprop.zenith)
        # Calculate 'Total Weight' (NuGen jargon)
        pint=wint*wprop
        # Calculate OneWeight
        oneweight=pint*wvol/wgen
        # Store all the weight parameters in the array container
        wparam=(evtid, ppid, self.index, self.emin, self.emax, self.injrad, self.zmin, self.zmax,
                self.amin, self.amax, eventprop.totalEnergy, eventprop.zenith, eventprop.azimuth,
                self.xs, eventprop.totalColumnDepth, eventprop.impactParameter, wprop, wint, pint, oneweight, 
                self.nevents, self.ptype)
        #self.weightarray[evtid]=wparam
        self.weightarray[self.iterator]=wparam
        self.PushFrame(frame)
        self.iterator+=1





class WeightAdd(icetray.I3Module):
	'''
	This I3Module class takes the hdf5 file containg the weight info and puts it in the I3File's
	WeightDictionary
	'''
	def __init__(self, context):
		icetray.I3Module.__init__(self, context)
		self.AddParameter('Filename', #HDF5 filename
				'Weight Dictionary filename',
				'nsq_propagation_weight.h5')

		self.AddParameter('Weightname', #Weight Object name
				'Weight Dictionary object name in the frame',
				'I3MCWeightDict')
		self.AddParameter('SkipEvents', #Events to skip weights for
				'List of event ids to skip',
				[])
		self.AddParameter('AddID', #If we want to add event ID object
				'Whether to Add SimEventID',
				True)
		self.AddParameter('MapID',
				'EventID or PrimaryID for mapping the weight array',
				'PID')#default is primary id, alternate: EID

	def Configure(self):
		h5filename = self.GetParameter('Filename')
		self.weightname = self.GetParameter('Weightname')
		self.skipevents=self.GetParameter('SkipEvents')
		self.addid = self.GetParameter('AddID')
		self.mapid = self.GetParameter('MapID')

		h5file = h5.File(h5filename,'r')
		self.weight_data = h5file['WeightDictionary']
		self.keys = ['SpectralIndex', 'Emin', 'Emax', 'InjectionRadius',
			'ZenithMin', 'ZenithMax', 'AzimuthMin', 'AzimuthMax',
			'PrimaryEnergy', 'PrimaryZenith', 'PrimaryAzimuth',
			'TotalCrossection', 'ColumnDepth', 'ImpactParamter','PropagationWeight', 
			'InteractionWeight', 'TotalWeight', 'OneWeight',
			'NEvents', 'NuType']
		#set the idx for pulling weights from dictionary
		if self.mapid=='PID':
			pp_ids = self.weight_data['primary_id']
		elif self.mapid=='EID':
			pp_ids = self.weight_data['event_id']
		weight_ids= np.arange(len(pp_ids))
		self.idx_map = dict(zip(pp_ids,weight_ids))
#		self.iterator=0
		#print (self.idx_map.keys())
		if self.addid: self.id_key = 'SimEventID'

	def DAQ(self, frame):
#		for i in self.skipevents:
#			if self.weight_data[self.iterator]['event_id']==i:
#				self.iterator+=1
#				print ("Skipping event with Event ID: ", self.weight_data[self.iterator][0])	
#		evt_weight=self.weight_data[self.iterator]
		if self.mapid=='PID':
			evt_weight=self.weight_data[self.idx_map[frame['I3MCTree'].primaries[0].id.minorID]]
		elif self.mapid=='EID':
			evt_weight=self.weight_data[self.idx_map[frame['I3EventHeader'].event_id]]

		evt_weight['NEvents'] -= len(self.skipevents)

		wdict=dict(zip(self.keys, evt_weight.tolist()[2:]))
		wdict_obj=dataclasses.I3MapStringDouble(wdict)
		frame.Put(self.weightname,wdict_obj)

		evtid = icetray.I3Int(int(evt_weight['event_id']))
		frame.Put(self.id_key, evtid)
		self.PushFrame(frame)
#		self.iterator += 1

class InsertID(icetray.I3Module):
	'''
	This I3Module class takes the hdf5 file containg the weight info (initial eventID) and puts it in the I3File's
	SimEventID object
	'''
	def __init__(self, context):
		icetray.I3Module.__init__(self, context)
		self.AddParameter('Filename', #HDF5 filename
				'Weight Dictionary filename')

	def Configure(self):
		h5filename = self.GetParameter('Filename')
		h5file = h5.File(h5filename,'r')
		self.weight_data = h5file['WeightDictionary']
		self.key = 'SimEventID'
		pp_ids = self.weight_data['primary_id']
		ev_ids = self.weight_data['event_id']
		self.idx_map = dict(zip(pp_ids,ev_ids))

	def DAQ(self, frame):
		evt_id=int(self.idx_map[frame['I3MCTree'].primaries[0].id.minorID])
		evtid = icetray.I3Int(evt_id)
		frame.Put(self.key, evtid)
		self.PushFrame(frame)


class Update_ID(icetray.I3Module):
	''' This class takes the weight hdf5 and updates the primary particle minor id
	motivation: for tridents, particle ids change when we create new mctree after
	calchep/madgraph event generation. where as the weight info is mainly constructed
	from lepton injector file that contains old mctree with different pids. so we
	need to update the ids in the weight dictionary for later precessing/analysis.
	'''
	def __init__(self, context):
		icetray.I3Module.__init__(self, context)
		self.AddParameter("HDF5Filename", "Name of HDF5 file containing Weight Dictionary",
				"sample.h5")
	def Configure(self):
		h5filename = self.GetParameter('HDF5Filename')
		self.h5file = h5.File(h5filename, 'a')
		evtid = self.h5file['WeightDictionary']['event_id']
		weight_idx = np.arange(len(evtid))
		self.idxmap = dict(zip(evtid,weight_idx))
		self.iterator = 0
	def DAQ(self, frame):
		pid = frame['I3MCTree'].primaries[0].id.minorID
		eid = frame['I3EventHeader'].event_id
		arr = self.h5file['WeightDictionary'][self.idxmap[eid]]
		arr['primary_id'] = pid
		self.h5file['WeightDictionary'][self.idxmap[eid]] = arr
		self.PushFrame(frame)
		self.iterator += 1
		if self.iterator==len(self.idxmap.keys()): self.h5file.close()
