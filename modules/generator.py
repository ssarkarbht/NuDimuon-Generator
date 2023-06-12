#!usr/bin/env python

'''
Author: Sourav Sarkar
Date: September 23, 2020
Description: Neutrino dimuon configuration module.
        Can be used for the following tasks:
        1. Sample neutrino energy using any given
        power law spectrum
        2. Sample Target (proton or neutron) accroding to
        the medium composition (e.g. ice - 10 proton, 8 neutron)
        3. 

'''

#import standard python libraries
import numpy as np
import h5py
import scipy.stats as st
from glob import glob
from optparse import OptionParser
from numpy.random import default_rng
import scipy.interpolate as ip
import scipy.integrate as it

#custom EPA approximation cross-section calculation module
import EPA_XS as epa

#IceCube Module for random number generators (optional)
try:
	from icecube import phys_services
except:
	print ("No icecube modules available, will be using default python modules.")

#for testing purpose
#import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.rc('figure',dpi=250)
#pltloc='/home/ssarkar/public_html/trident_plots/diff_distributions/'


#class to generate neutrino energy sample, interaction, nucleus type

class generator:
	'''
	This class generates the even generator configuration
	parameters for different neutrino trident interactions, using
	total cross-section results reported in [paper: Zhou, Beacom]
	and EPA method reportd in [paper: plastid]
	User also has the option to sample primary neutrino energy
	from a powerlaw energy spectrum flux using this module.

	All events generated based on this configuration assumes the
	interaction along Z axis (lhc geometry). Thus for a non-trivial
	geometry, user should rely on external detector simulation
	libraries (like LeptonInjector) or manual coordinate transformation
	of interacting particles post event generation.
	'''

	#Initialize with input paramter: random generator, cross-section data
	def __init__(self,rng_gen,xsfiles=None):
		self.rng=rng_gen
		if xsfiles is not None:
			#prepare the cross-section functions
			# Coherent + Diffractive XS
			chdf_xs=np.loadtxt(xsfiles[1],delimiter=',')
			logen=np.log10(chdf_xs[:,0])
			logxs=np.log10(chdf_xs[:,1])
			self.chdf_func=ip.interp1d(logen,logxs,kind='slinear')
			# DIS (quark mediated)
			disq_xs=np.loadtxt(xsfiles[2],delimiter=',')
			logen=np.log10(disq_xs[:,0])
			logxs=np.log10(disq_xs[:,1])
			self.disq_func=ip.interp1d(logen,logxs,kind='slinear')
			# DIS (photon mediated)
			disp_xs=np.loadtxt(xsfiles[3],delimiter=',')
			logen=np.log10(disp_xs[:,0])
			logxs=np.log10(disp_xs[:,1])
			self.disp_func=ip.interp1d(logen,logxs,kind='slinear')

	#If generating standalone (lhc geometry) events, sample neutrino energy based on given spectrum		
	#Parameters: Powerlaw spectral index, Minimum and Maximum energy, number of events (default 1)
	def sample_energy(self, index, emin, emax, n=1):

		#If only one energy allowed, return only possible energy (monoenergetic events)
		if emin==emax:
			return np.repeat(emin,n)
		#If powerlaw index is 1, sample uniformly in logspace
		if index==1.0:
			rand=np.log10(emin)+self.rng.random(size=n)*(np.log10(emax)-np.log10(emin))
			return 10**rand
		# For other powerlaw index, generate sample using imortance sampling from cdf
		else:
			rand=self.rng.random(size=n)
			en_prob=(1-rand)*emin**(1-index)+rand*emax**(1-index)
			return en_prob**(1/(1-index))

	#sample the type of interaction regime for given neutrino energy: (Coherent+Diffractive) vs. DIS
	def sample_interaction(self,en):
		#calculate the fraction of interaction cross-sections (coherent+diffractive/DIS)
		logen=np.log10(en)
		disxs=10**(self.disp_func(logen))+10**(self.disq_func(logen))
		chdxs=10**(self.chdf_func(logen))
		totalxs=chdxs+disxs
		chdf_frac=chdxs/totalxs
		dis_frac=disxs/totalxs

		if type(en)==float:
			num=1
		else:
			num=len(en)
		#generate random numbers uniformely between [0,1)
		rand=self.rng.random(size=num)
		interaction_arr=np.zeros(num)
		#randomely pick the interaction based on fractional probability
		int1=np.where(rand<=chdf_frac)
		int2=np.where(rand>chdf_frac)
		#interaction type=1 for coherent+diffractive, 2 for DIS
		interaction_arr[int1]=1.
		interaction_arr[int2]=2.
		#return (totalxs*en,interaction_arr)
		return interaction_arr

#sample the type of DIS interaction for a given neutrino energy: Photon vs. Quark mediated DIS interaction
	def sample_DIS(self, en):
		logen=np.log10(en)
		photon_xs=10**(self.disp_func(logen))
		quark_xs =10**(self.disq_func(logen))
		disxs=photon_xs+quark_xs
		photon_frac=photon_xs/disxs
		quark_frac =quark_xs/disxs
		if type(en)==float:
			num=1
		else:
			num=len(en)
		#generate random numbers uniformely between [0,1)
		rand=self.rng.random(size=num)
		interaction_arr=np.zeros(num)
		#randomely pick the interaction based on fractional probability
		int1=np.where(rand<photon_frac)
		int2=np.where(rand>=photon_frac)

		#interaction type=1 for photon initiated, 2 for quark initiated
		interaction_arr[int1]=1
		interaction_arr[int2]=2

		return interaction_arr

#sample nucleus type between proton and neutron
#ratio of proton and neutron in water is 10:8

	def sample_nucleus(self,n=1):
		nucleus_arr=np.zeros(n)
		num=self.rng.integers(18,size=n)
		int1=np.where(num<10)
		int2=np.where(num>9)
		#nucleus type=1 for proton, type=2 for neutron
		nucleus_arr[int1]=1
		nucleus_arr[int2]=2
		return nucleus_arr
'''
def converter_func(q,nuen):
	q2,val1=epa.WF(nuen,'OX')
	q2,val2=epa.WF(nuen,'H')
	q_arr=np.sqrt(q2).reshape(,len(q2))
	total_wf=val1+2.*val2
	prob=epa.epa_xs(4.*nuen*q_arr)*total_wf/q_arr
	prob_func=ip.interp1d(np.log10(q_arr[0,:]),np.log10(prob[0,:]))
	return 10**(func(np.log10(q)))
'''

#Following function loads EPA form factor weight data and returns the interpolation fucntion
def weight_function(wffiles):
	ox_data=np.loadtxt(wffiles[0])
	hy_data=np.loadtxt(wffiles[1])

	val1=ox_data[:,1]
	val2=hy_data[:,1]
	total_wf=val1+2.*val2#weight factor of oxygen plus 2 times the weight factor of hydrogen
	log_wf=ip.interp1d(np.log10(ox_data[:,0]),np.log10(total_wf))
	return log_wf

#
def converter_func(q,func):
	return 10**(func(np.log10(q)))

#class for generating interacting EPA photon energy from underlying distribution
#for a given neutrino energy

class EPAPhoton:
	def __init__(self,nuen,rng_gen,wfunc):

		self.rng=rng_gen
		#set the photon energy range for the incoming neutrino energy
		self.qmin=(0.0449/4./nuen)#GeV
		#self.qmax=1.
		self.qmax=1.295#GeV based on CTEQ14QED pdf
		q2arr=np.logspace(np.log10(self.qmin**2),np.log10(self.qmax**2),1000)
		self.q=np.sqrt(q2arr)

		#prepare the final pdf values (not normalized)
		wf_val=10**wfunc(np.log10(q2arr))
		self.prob=epa.epa_xs(4.*nuen*self.q)*wf_val/self.q

		#distances between the photon energy values
		logq=np.log10(self.q)
#		binwidth_log_avg=np.average(np.diff(logq))
#		bin_edge=0.5*(logq[:-1]+logq[1:])
#		bin_edge=np.append(bin_edge,bin_edge[-1]+binwidth_log_avg)
#		bin_edge=np.append(bin_edge[0]-binwidth_log_avg,bin_edge)
#		bin_edge=10**bin_edge
#		bin_diff=np.diff(bin_edge)

		cumsum=np.cumsum(self.prob*self.q)#*bin_diff
		self.cdf_val=cumsum/max(cumsum)

		self.log_inv_cdf=ip.interp1d(np.log10(self.cdf_val),np.log10(self.q))

#For testing purpose
#		plt.figure()
#		plt.plot(self.q,self.cdf_val)
#		plt.title('CDF for neutrino energy: %.2e' %nuen)
#		plt.xscale('log')
#		plt.yscale('log')
#		plt.xlabel("Photon energy (GeV)")
#		plt.savefig(pltloc+'cdf_check.pdf')

#interpolation function for extracting pdf
		self.prob_func=ip.interp1d(np.log10(self.q),np.log10(self.prob))

		#print ("Calculating the normalization factor for photon pdf")

#Calculation of the normalization factor for pdf

#		self.norm, err=it.quadrature(converter_func,self.qmin,self.qmax,args=(self.prob_func,))
		#self.norm, err=it.quad(converter_func,np.log10(self.qmin),np.log10(self.qmax),args=(self.prob_func,))
#		self.norm, err=it.quad(converter_func,self.qmin,self.qmax,args=(self.prob_func,))
		#self.norm, err=it.fixed_quad(converter_func,self.qmin,self.qmax,args=(self.prob_func,))
		#print ("Normalization+err is: ", self.norm, err)

	def pdf(self,q):
		return 10**self.prob_func(np.log10(q))#/self.norm

	def sampling(self,n):
#		u=self.rng.random(size=n)
		u=np.array([])
		for i in range(n):
			u=np.append(u,self.rng.uniform(0.,1.))
		logu=np.log10(u)
		logq=self.log_inv_cdf(logu)
		return (10**logq)
