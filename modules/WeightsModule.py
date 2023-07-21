#!/usr/bin/env python

'''
Author: Sourav Sarkar
Date: 7 Dec, 2020
Description: This script contains the python functions to calculate weight information
	for trident events. Parameters needed for calculating the weights come from
	parent LeptonInjector I3 files, nuSQuIDs and trident total cross-section.
'''

import numpy as np
from icecube import icetray, dataio, dataclasses, LeptonInjector
import nuSQUIDSpy as nsq
from scipy.interpolate import interp1d
import h5py as h5

#define the path for trident cross-section data
xsfilepath='/data/user/ssarkar/TridentProduction/simulation/WeightCalculation/cross-section_data/'

#Load the cross-section data:
#following three are the trident total cross-section data (same for Nu/NuBar)
chdf_xs=np.loadtxt(xsfilepath+'new_beacom_chdf.csv',delimiter=',',comments='#')
disq_xs=np.loadtxt(xsfilepath+'new_beacom_DIS_quark.csv',delimiter=',',comments='#')
disp_xs=np.loadtxt(xsfilepath+'new_beacom_DIS_photon.csv',delimiter=',',comments='#')
#following two are the standard CC total cross section
en_points = np.linspace(1,12,111)
cc_filepath = '/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms/'
nubar_cc = np.loadtxt(cc_filepath+'total_nubar_CC_iso_NLO_HERAPDF1.5NLO_EIG.dat',skiprows=1)
nu_cc    = np.loadtxt(cc_filepath+'total_nu_CC_iso_NLO_HERAPDF1.5NLO_EIG.dat', skiprows=1)
nubar_nc = np.loadtxt(cc_filepath+'total_nubar_NC_iso_NLO_HERAPDF1.5NLO_EIG.dat',skiprows=1)
nu_nc    = np.loadtxt(cc_filepath+'total_nu_NC_iso_NLO_HERAPDF1.5NLO_EIG.dat', skiprows=1)

nubar_cc_func = interp1d(en_points, nubar_cc, kind='cubic')
nu_cc_func    = interp1d(en_points, nu_cc, kind='cubic')

nubar_nc_func = interp1d(en_points, nubar_nc, kind='cubic')
nu_nc_func    = interp1d(en_points, nu_nc, kind='cubic')

#Function that calculates the total Standard CC xs from data as function of neutrino energy
def totalccxs(nuen, nutype):
	''' takes the neutrino energy and nutype= 1 for nu, -1 for nubar
	returns total cross-section in pb
	'''
	logen = np.log10(nuen)
	if nutype==-1:
		xsval = nubar_cc_func(logen)*1e9 #cross-section value in pb (1 pb=10^9 mb)
	elif nutype==1:
		xsval = nu_cc_func(logen)*1e9 #cross-section value in pb (1 pb=10^9 mb)
	else:
		assert False, f'Wrong nutype set: {nutype}'
	return xsval

def totalncxs(nuen, nutype):
	''' takes the neutrino energy and nutype= 1 for nu, -1 for nubar
	returns total cross-section in pb
	'''
	logen = np.log10(nuen)
	if nutype==-1:
		xsval = nubar_nc_func(logen)*1e9 #cross-section value in pb (1 pb=10^9 mb)
	elif nutype==1:
		xsval = nu_nc_func(logen)*1e9 #cross-section value in pb (1 pb=10^9 mb)
	else:
		assert False, f'Wrong nutype set: {nutype}'
	return xsval
#convert the data in log space for smooth interpolation
chdf_logen=np.log10(chdf_xs[:,0])
chdf_logxs=np.log10(chdf_xs[:,1])

disq_logen=np.log10(disq_xs[:,0])
disq_logxs=np.log10(disq_xs[:,1])

disp_logen=np.log10(disp_xs[:,0])
disp_logxs=np.log10(disp_xs[:,1])

#build the interpolation function
chdf_func=interp1d(chdf_logen,chdf_logxs,kind='slinear')
disq_func=interp1d(disq_logen,disq_logxs,kind='slinear')
disp_func=interp1d(disp_logen,disp_logxs,kind='slinear')


#Function that calculates the total trident cross-section from data
def totalxs(nuen):
	global chdf_func,disq_func,disp_func
	#calculate the (xs/en) from the interpolation functions
	logen=np.log10(nuen)
	chdf_val=10**chdf_func(logen)
	disq_val=10**disq_func(logen)
	disp_val=10**disp_func(logen)

	total_val=chdf_val+disq_val+disp_val
	xs=total_val*nuen*1e36 #cross-section value in pb (1 pb=10^-36 cm^2)
	return xs


class WeightExtract(icetray.I3Module):
	'''
	I3Module class to extract weight information from initial LI files
	and calculate final weights and store it into an hdf5 file 
	'''

	def __init__(self, context):
		icetray.I3Module.__init__(self, context)
		#Parameter for HDF5 filename
		self.AddParameter('WeightArray', #parameter for storing weight info
				'Container for weight information', #doc
				'warray') #default filename

		#paramter for input Simulation frame
		self.AddParameter('SFrame',
				'Input Simulation frame',
				'sframe')

		#paramter for input nuSQuIDs HDF5 file
		self.AddParameter('NSQFile', #file name containing neutrino evolved state
				'Input nuSQUIDs file',
				'nsq_propagation_weight.h5')

		#parameter for indicating weight calculation for trident/standard cc interaction
		self.AddParameter('InteractionType',# 1 for Trident interaction, 0 for Standard interaction
				'Interaction type for weight calculation',# 2 for charm muon production
				None)

		self.AddParameter('CharmFactor',#charm fractional cross-section for weight calculation
				'Charm Muon Fractional Cross-section',
				None)
	def Configure(self):
		#Load the array object
		self.weightarray = self.GetParameter('WeightArray')
		#Load the s frame object and event properties in it
		sframe   = self.GetParameter('SFrame')
		nsq_file = self.GetParameter('NSQFile')
		self.inttype = self.GetParameter('InteractionType')
		self.charm_factor = self.GetParameter('CharmFactor')

		self.simprop = sframe['LeptonInjectorProperties']
		self.nsq_atm = nsq.nuSQUIDSAtm(nsq_file)
		print ("nsqatm file loaded")

		self.conversion=1e-36 #cross-section conversion factor (from pb to cm^2)
		#self.density=0.92 #g/cm^3 denisty of ice for column depth calculation
		#self.proton_mass=1.6726219e-24 #gm
		self.NA = 6.022140857e+23 #Avogadro's number
		self.iterator = 0
		if self.inttype ==2:
			val = self.charm_factor.values()
			valarr = np.array(list(val))
			num_miss = len(np.where(valarr==0)[0])
			self.nevents = self.simprop.events - num_miss
		elif self.inttype!=2:
			self.nevents = self.simprop.events

		
		print ("Read All parameters!")
	#Calculate the generation probability
	def pgen(self,nuen):
		#get the parameters
		self.index=self.simprop.powerlawIndex
		self.emin=self.simprop.energyMinimum
		self.emax=self.simprop.energyMaximum
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
		self.injrad=self.simprop.injectionRadius*1e2 #convert from meter to cm
		self.zmax=self.simprop.zenithMaximum
		self.zmin=self.simprop.zenithMinimum
		self.amin=self.simprop.azimuthMinimum
		self.amax=self.simprop.azimuthMaximum
		#injection area
		injarea=np.pi*self.injrad**2
		#Solid angle used for simulation
		solidangle=(self.amax-self.amin)*(np.cos(self.zmin)-np.cos(self.zmax))
		return injarea*solidangle

	#Calculate the interaction weight within generation volume
	def interaction_weight(self,nuen,columndepth):
		#tridents
		if self.inttype==1:#divide by 16 to get iso cross section
			self.xs=totalxs(nuen)*self.conversion/16. #in cm^2 unit
		#nu cc dis
		elif self.inttype==0:
			self.xs=totalccxs(nuen, self.ptype)*self.conversion #in cm^2 unit
		#nu cc dis charm muon
		elif self.inttype==2:
			self.xs=totalccxs(nuen, self.ptype)*self.charmfactor*self.conversion #in cm^2 unit
		else:
			assert False, f'Invalide Interaction type set: {self.inttype}'
		cd=columndepth*self.NA #in #targets/cm^-2 unit
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
