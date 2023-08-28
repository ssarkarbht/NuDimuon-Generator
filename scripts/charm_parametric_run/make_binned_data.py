#!/bin/python

'''
Author : Sourav Sarkar
Email : ssarkar1@ualberta.ca
Description : This script takes the grid simulation
    from pythia (for decay) or chromo (for interaction)
    and creates the binned histogram as the parametric
    tables to be used during the charm muon sampling.
'''

import numpy as np
from glob import glob
from optparse import OptionParser
import os,sys
import json
import h5py as h5

#add command-line arguments
parser = OptionParser()

parser.add_option("-t", "--Target", dest="TGT", type=str,
                default=None)
parser.add_option("-d", "--DataFolder", dest="DIR", type=str)
parser.add_option("-s", "--SimulationType", dest="SIM", type=str)
parser.add_option("-e", "--MinMuEnergy", dest="EMUMIN", type=float)

parser.add_option("-m", "--MinEnergy", dest="EMIN", type=float)
parser.add_option("-u", "--MaxEnergy", dest="EMAX", type=float)
parser.add_option("-n", "--BinNumber", dest="NBIN", type=int)
parser.add_option("-o", "--OutputDir", dest="OUT", type=str)

(options,args) = parser.parse_args()

# get the repo directory path
repo_dir = os.environ["DIMUON_REPO"]
dataloc = repo_dir + "/data/constants_particles_materials/"

# get the target materials and projectile properties
pfile = dataloc+"particles.json"
with open(pfile, 'r') as f:
    particle_dict = json.load(f)
#get list of charm hadron pdg codes
pdglist = []
for key, val in particle_dict.items():
    pdgid = val['pdgid']
    pdglist.append(pdgid)
#get a reverse dictionary for charm hadron names
names = list(particle_dict.keys())
pdgnames = dict(zip(pdglist, names))
namespdg = dict(zip(names, pdglist))

mfile = dataloc+"medium_properties.json"
with open(mfile, "r") as f:
    medium = json.load(f)

def make_histo(arr, energy, bins):
    ''' Creates the binned histogram data from
    simulated events.
    Input: arr (ndarray) : array of muon props
        energy (float): Charm hadron energy point
        bins (dict): three bins for histogramming
    '''
    mue_thr = options.EMUMIN/energy
    #get the index for non-zero muon production
    idxs = np.where(arr['pdgid']!=0)[0]
    #get the index for above threshold muon production
    thrs = np.where(arr['e_fraction']>=mue_thr)[0]

    #get total number of events w/ above thr. muons
    num = len(thrs)
    #get total number of simulated events
    dum = float(len(arr))

    #make histograms
    # Energy fraction
    fhist, _ = np.histogram(arr['e_fraction'][idxs],
            bins = bins['fbins'], density=True)
    # Opening angle
    ahist, _ = np.histogram(arr['theta'][idxs],
            bins = bins['abins'], density=True)
    # Muon Multiplicity
    mhist, _ = np.histogram(arr['multiplicity'],
            bins = bins['mbins'], density=True)

    return (num, dum, fhist, ahist, mhist)


#check if the sim type has the proper name
assert options.SIM in ['decay', 'interaction'], "Simulation type not decay or interaction!"
#parent directory
#datadir = options.DIR + "/" + options.SIM + "/"
datadir = options.DIR + "/"

#initiate empty disctionary of arrays for output
hist_dict = {}

#make a metadata array for the compressed file
metadata = np.array(['EnergyPoints : Parent hadron energies used for simulation',
    '<hadron>_<property> : first row is the bin centers of the histogrammed data',
    '<hadron>_branchingRatio : energy, # muons, # total evnets'])

#make the histogram bins
# fractional energy bins 
# (log-scale for interaction, lin-scale for decay)
if options.SIM == 'interaction':
    minfrac = np.log10(1./options.EMAX) #assuming lowest possible bins @1GeV Muons
    fbins = np.logspace(minfrac, 0., options.NBIN+1)
elif options.SIM == 'decay':
    fbins = np.linspace(0., 1., options.NBIN+1)
fbincen = 0.5*(fbins[:-1]+fbins[1:])

# opening angle bins (assuming highest opening angle @90 degrees)
abins = np.linspace(0., np.pi/2., options.NBIN+1)#rad
abincen = 0.5*(abins[:-1]+abins[1:])

# multiplicity bins (0 to 10 muons)
mbins = np.linspace(-0.5, 10.5, 12)
mbincen = 0.5*(mbins[:-1]+mbins[1:])

bin_dict = {'fbins':fbins, 'abins':abins, 'mbins':mbins}

#loop over the charm hadrons
for particle in names:
    if options.SIM=='decay':
        fname = datadir+particle+'.npz'
        try:
            data = np.load(fname)
            print (f"Histogramming {particle} decay...")
        except:
            dummy = pdgnames[abs(namespdg[particle])]
            print (f"Data for {particle} not found, replacing with {dummy}")
            fname = datadir+dummy+'.npz'
            data = np.load(fname)
    elif options.SIM=='interaction':
        fname = datadir + options.TGT + '_' + particle + '.npz'
        try:
            data = np.load(fname)
            print (f"Histogramming {particle} interaction...")
        except:
            dummy = pdgnames[abs(namespdg[particle])]
            print (f"Data for {particle} not found, replacing with {dummy}")
            fname = datadir + options.TGT + '_' +dummy+'.npz'
            data = np.load(fname)

    eninfo = data['energy']
    #build the charm hadron energy array
    energies = np.logspace(np.log10(eninfo[0]), 
                        np.log10(eninfo[1]),
                        int(eninfo[2]))

    #initialize empty arrays
    farr = np.zeros((len(energies)+1, len(fbincen)))
    farr[0] = fbincen

    aarr = np.zeros((len(energies)+1, len(abincen)))
    aarr[0] = abincen

    marr = np.zeros((len(energies)+1, len(mbincen)))
    marr[0] = mbincen

    br_arr = np.zeros((len(energies), 3))

    #loop over each energy points
    for idx, en in enumerate(energies):
        arr = data[str(idx)]
        num, dum, fhist, ahist, mhist = make_histo(arr, en, bin_dict)
        #start filling the arrays
        br_arr[idx] = (en, num, dum)
        farr[idx+1] = fhist
        aarr[idx+1] = ahist
        marr[idx+1] = mhist

    #store the arrays into parent dictionary
    hist_dict[particle+'_energyFraction'] = farr
    hist_dict[particle+'_openingAngle'] = aarr
    hist_dict[particle+'_multiplicity'] = marr
    hist_dict[particle+'_branchingRatio'] = br_arr

#save the energy points
hist_dict["EnergyPoints"] = energies
#save the histogrammed data into numpy compressed file
np.savez_compressed(options.OUT, metadata=metadata, **hist_dict)
