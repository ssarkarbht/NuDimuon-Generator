#!/bin/python

'''
Author: Sourav Sarkar
Email: ssarkar1@ualberta.ca
Date: Feb 21, 2020
Desciption: The purpose of this script is to provide required modules
that runs appropriate functions and integation method to evaluate the 
final incoming neutrino-nucleon cross-section using the EPA model
for coherent and diffractive regime
'''

#User should have the option to choose different nuclear form factor models


import numpy as np
import scipy.interpolate as ip
from scipy.integrate import quad
import scipy.integrate as it
from glob import glob
import os

repo_dir = os.environ["DIMUON_REPO"]

#Plotting is for debugging and testing purposes
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('figure',dpi=350)

#################################---Segment for building EPA cross-section function----#########
#Loading the CM EPA cross-section data
xsloc = repo_dir + "/data/cross_sections/trident_cross_section/"
epa_data=np.loadtxt(xsloc + "EPA_CM_SM_XS_DATA.txt")
cmen=epa_data[:,0]
cmxs=epa_data[:,1]

#Defining the interpolation function in log scale (both in energy and cross-section)
epa_func=ip.UnivariateSpline(np.log10(cmen),np.log10(cmxs),s=0)

#define a wrapper function to spit out cross-section in pb
def epa_xs(s):
        lcmen=np.log10(np.sqrt(s))
        val=epa_func(lcmen)
        return 10**val

def epa_ana_xs(s):
	G=1.166e-5
	alpha=1/137.
	smin=0.044944
	xsa=G**2*alpha*s/np.pi**2/9*np.log(s/smin)*0.389379e9
	return xsa

##########################-----Form Factors------########################################
#Making the nuclear type library for form factor calculation
ntype_lib={'OX':(16.,8.),'H':(1.,1.), 
        'LAr':(40.,18.), 'Pb':(206.,82.),
        'TW':(184,74)}


#Coherent Form Factor:
#// Phys.Rev. D3 (1971) 2686-2706

#The function takes the photon momentum transfer squared and nucleus type
def FI_I(Q2,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	a = (5.076*(0.58+0.82*A**(1/3.)))**2
	fac = (1+a*Q2/12.)**2
	return (1/fac)**2/Q2

#Wood-Saxon form factor in the analytical form from Matheus Hoestert article
#
def FI_II(Q2,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	a = 0.523*5.068 #1fm=5.068GeV^-1
	r0 = 1.126*A**(1/3.)*5.068
	q=np.sqrt(Q2)
	val=(3*np.pi*a/(r0**2+np.pi*np.pi*a**2))*(np.pi*a*(1/np.tanh(np.pi*q*a))*np.sin(q*r0)-r0*np.cos(q*r0))/(q*r0*np.sinh(np.pi*q*a))
	return val**2/Q2

def logFI_II(lQ2,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	a = 0.523*5.068 #1fm=5.068GeV^-1
	r0 = 1.126*A**(1/3.)*5.068
	Q2=10**lQ2
	q=np.sqrt(Q2)
	val=(3*np.pi*a/(r0**2+np.pi*np.pi*a**2))*(np.pi*a*(1/np.tanh(np.pi*q*a))*np.sin(q*r0)-r0*np.cos(q*r0))/(q*r0*np.sinh(np.pi*q*a))
	return val**2


'''
q=np.logspace(-9,1,1000)
plt.figure()
plt.plot(q,q**2*FI_II(q**2,'OX'))
plt.plot(q,q**2*FI_I(q**2,'OX'))
#plt.plot(q,FI_I(q**2,'OX')/FI_I(q[0]**2,'OX'))
plt.yscale('log')
plt.xscale('log')
#plt.xlim([1e-2,7e-1])
plt.ylim([1e-5,1e2])
plt.savefig('testing.pdf')
'''

#Diffractive dipole form factor
#// 1612.05642

def Fdip(Q2,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	gdip = (1+Q2/0.71)**(-2)
	tau  = Q2/4.
	xi   = 4.7
	fac  = gdip*(1+tau*xi)/(1+tau)
	return fac**2/Q2

def logFdip(lQ2,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	Q2=10**lQ2
	gdip = (1+Q2/0.71)**(-2)
	tau  = Q2/4.
	xi   = 4.7
	fac  = gdip*(1+tau*xi)/(1+tau)
	return fac**2

##############################################################################################

#Final Integration function
#This function calculates the integration of the form factor term
#input argument takes the photon momentum transfer (not squared)
#Also, the function multplies the FF integration with the cross-section
#integrand for the CM energy integration


def FF_Int(s,nuen,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	q=s/2./nuen
	Q2 = q**2
	#Constant factor in the integral
	#For FI
	c1 = Z**2*1/137./np.pi
	#For Fdip
	c2 = Z*1/137./np.pi
	# Q boundary b/w Coherent and Diffractive regime
	Qlim = 0.217/A**(1/3.)
	#initializing the integrand values
	cy=dy=0.
	if q<Qlim:
#		cy=FI_II(Q2,ntype)
		cy, cerr = quad(FI_II,Q2,Qlim**2,args=(ntype,))
		dy, derr = quad(Fdip,Qlim**2,1.,args=(ntype,))
	else:
#		dy=Fdip(Q2,ntype)
		dy, derr = quad(Fdip,Q2,1.,args=(ntype,))
	return (c1*cy+c2*dy)*epa_xs(s)/s/2.
	#return (c1*cy+c2*dy)*epa_ana_xs(s)/s
	#return epa_xs(s)/s

def FF_UERR_Int(s,nuen,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	q=s/2./nuen
	Q2 = q**2
	#Constant factor in the integral
	#For FI
	c1 = Z**2*1/137./np.pi
	#For Fdip
	c2 = Z*1/137./np.pi
	# Q boundary b/w Coherent and Diffractive regime
	Qlim = 0.217/A**(1/3.)
	#initializing the integrand values
	cy=dy=0.
	cerr=derr=0.
	tot_err=0.
	if q<Qlim:
#		cy=FI_II(Q2,ntype)
		cy, cerr = quad(FI_II,Q2,Qlim**2,args=(ntype,))
		dy, derr = quad(Fdip,Qlim**2,1.,args=(ntype,))
		tot_err=np.sqrt((c1*cerr/cy)**2+(c2*derr/dy)**2)
	else:
#		dy=Fdip(Q2,ntype)
		dy, derr = quad(Fdip,Q2,1.,args=(ntype,))
		tot_err=(c2*derr/dy)
	return ((c1*cy+c2*dy)*(1+tot_err))*epa_xs(s/2)/s/2.

def FF_LERR_Int(s,nuen,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	q=s/2./nuen
	Q2 = q**2
	#Constant factor in the integral
	#For FI
	c1 = Z**2*1/137./np.pi
	#For Fdip
	c2 = Z*1/137./np.pi
	# Q boundary b/w Coherent and Diffractive regime
	Qlim = 0.217/A**(1/3.)
	#initializing the integrand values
	cy=dy=0.
	cerr=derr=0
	if q<Qlim:
#		cy=FI_II(Q2,ntype)
		cy, cerr = quad(FI_II,Q2,Qlim**2,args=(ntype,))
		dy, derr = quad(Fdip,Qlim**2,1.,args=(ntype,))
		tot_err=np.sqrt((c1*cerr/cy)**2+(c2*derr/dy)**2)
	else:
#		dy=Fdip(Q2,ntype)
		dy, derr = quad(Fdip,Q2,1.,args=(ntype,))
		tot_err=(c2*derr/dy)
	return ((c1*cy+c2*dy)*(1-tot_err))*epa_xs(s/2)/s/2.

def FFC(s,nuen,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	q=s/2./nuen
	Q2=q**2
	Q2lim=(0.217/A**(1/3.))**2
	if Q2<Q2lim:
		ffi,err=quad(FI_II,Q2,Q2lim,args=(ntype,))
	else:
		ffi=0.0
	return epa_xs(s)*ffi/s

def FFD(s,nuen,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	q=s/2./nuen
	Q2=q**2
	Q2lim=(0.217/A**(1/3.))**2
	Q2min=max(Q2,Q2lim)
	Q2max=1.69
	ffi,err=quad(Fdip,Q2min,Q2max,args=(ntype,))
	return epa_xs(s)*ffi/s

def WF(enu,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	Q2lim=(0.217/A**(1/3.))**2
	Q2min=(0.0449/4./enu)**2
	#tiny value added to the maximum allowed value so the lower and upper limit of the integrtal isn't zero
	Q2max=1.0+1e-10#1.69

	#constant factors:
	c1=Z**2*(1/137.)/np.pi
	c2=Z*(1/137.)/np.pi

	#Array of Q2
#	Q2=np.linspace(Q2min,1.0,1000)
	Q2=np.logspace(np.log10(Q2min),0.,1000)
	weightf=np.array([])
	for i in Q2:
		s=2.*enu*np.sqrt(i)
		chf=dff=0.
		#chf,err=quad(FI_II,i,Q2lim,args=(ntype,))
#		print (i,chf)
		if i<Q2lim:
			chf,err=quad(FI_II,i,Q2lim,args=(ntype,),limit=100,epsabs=1.5e-10)
			dff,err=quad(Fdip,Q2lim,Q2max,args=(ntype,))
#			chf,err=it.fixed_quad(FI_II,i,Q2lim,args=(ntype,))
		elif i>=Q2lim:
			dff,err=quad(Fdip,i,Q2max,args=(ntype,))
		#val=(c1*chf+c2*dff)/s
		val=(c1*chf+c2*dff)
		weightf=np.append(weightf,val)
	#print (max(weightf),min(weightf))
	return (Q2,weightf)


def new_WF(enu,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	Q2lim=(0.217/A**(1/3.))**2
	Q2min=(0.0449/4./enu)**2
	#tiny value added to the maximum allowed value so the lower and upper limit of the integrtal isn't zero
	Q2max=1.0+1e-10#1.69

	#constant factors:
	c1=Z**2*(1/137.)/np.pi
	c2=Z*(1/137.)/np.pi

	#Array of Q2
#	Q2=np.linspace(Q2min,1.0,1000)
#	Q2=np.logspace(np.log10(Q2min),0.,1000)
	logQ2=np.linspace(np.log10(Q2min),0.,1000)
	weightf=np.array([])
	for i in logQ2:
		chf=dff=0.
		Q2=10**i
		#chf,err=quad(FI_II,i,Q2lim,args=(ntype,))
#		print (i,chf)
		if Q2<Q2lim:
			chf,err=quad(logFI_II,i,np.log10(Q2lim),args=(ntype,),limit=100)
			dff,err=quad(logFdip,np.log10(Q2lim),np.log10(Q2max),args=(ntype,))
		elif Q2>=Q2lim:
			dff,err=quad(logFdip,i,np.log10(Q2max),args=(ntype,))
		#val=(c1*chf+c2*dff)/s
		val=(c1*chf+c2*dff)
		weightf=np.append(weightf,val)
	#print (max(weightf),min(weightf))
	return (10**(logQ2),weightf)

###Form Factor calculation independent of neutrino energy (input is given photon energy range)
def ranged_WF(emin,emax,ntype,n=1000):
        ninfo=ntype_lib.get(ntype)
        if not(ninfo):
                print ('Error: No know nucleus type ',ntype,' found in the library')
                print ('Available Nucleus types are: ', ntype_lib.keys())
                return 0
        else:
                A,Z=ninfo[0], ninfo[1]
        Q2lim=(0.217/A**(1/3.))**2
        Q2min=emin**2
        #tiny value added to the maximum allowed value so the lower and upper limit of the integrtal isn't zero
        epsilon=1e-10#tiny value to add to upper limit Q2max for last point of the numerical integration
        Q2max=emax**2#1.69

        #constant factors:
        c1=Z**2*(1/137.)/np.pi
        c2=Z*(1/137.)/np.pi

        #Array of Q2
#       Q2=np.linspace(Q2min,1.0,1000)
#       Q2=np.logspace(np.log10(Q2min),0.,1000)
        logQ2=np.linspace(np.log10(Q2min),np.log10(Q2max),n)
        weightf=np.array([])
        counter=1
        for i in logQ2:
                if counter%1000==0:
                        print ('%d number of points processed...' %counter)
                chf=dff=0.
                Q2=10**i
                #chf,err=quad(FI_II,i,Q2lim,args=(ntype,))
#               print (i,chf)
                if Q2<Q2lim:
                        chf,err=quad(logFI_II,i,np.log10(Q2lim),args=(ntype,),limit=100)
                        dff,err=quad(logFdip,np.log10(Q2lim),np.log10(Q2max+epsilon),args=(ntype,))
                elif Q2>=Q2lim:
                        dff,err=quad(logFdip,i,np.log10(Q2max+epsilon),args=(ntype,))
                #val=(c1*chf+c2*dff)/s
                val=(c1*chf+c2*dff)
                weightf=np.append(weightf,val)
                counter+=1
        #print (max(weightf),min(weightf))
        return (10**(logQ2),weightf)

