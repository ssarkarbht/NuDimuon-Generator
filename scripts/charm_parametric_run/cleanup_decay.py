#!/bin/python

import numpy as np
import json
from glob import glob
from optparse import OptionParser

#add command-line arguments
parser = OptionParser()

parser.add_option("-d", "--DataFolder", dest="DIR", type=str)
parser.add_option("-m", "--MinEnergy", dest="EMIN", type=float)
parser.add_option("-u", "--MaxEnergy", dest="EMAX", type=float)
parser.add_option("-n", "--NPoints", dest="NGRID", type=int)
parser.add_option("-o", "--OutputDir", dest="OUT", type=str)

(options,args) = parser.parse_args()

# get the repo directory path
repo_dir = os.environ["DIMUON_REPO"]
dataloc = repo_dir + "/data/constants_particles_materials/"

#dataloc="../data/constants_particles_materials/"

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

#simloc="/data2/icecube/ssarkar/charm_muon_data/decay/"
simloc = options.DIR

#build the energy array
#enarr = np.logspace(np.log10(options.EMIN),
#        np.log10(options.EMAX), options.NGRID)

#energy bin info to store in the output file
metadata = np.array([options.EMIN, options.EMAX, options.NGRID])

for pdg in pdglist:
    print (f"Cleaning up {pdgnames[pdg]} data")
    pname = pdgnames[pdg]
    files = glob(simloc+ "/" +pname+"_*.txt")
    #check if data for the current particle exist
    if len(files)==0:
        dummy = pdgnames[abs(namespdg[pname])]
        print (f"{pname} data not found, will be replacing with {dummy}")
        files = glob(simloc+ "/" +dummy+"_*.txt")
    #[backward incompatibility] store muon pdg code
    if pdg>0: mupdg = 13
    elif pdg<0: mupdg = -13
    #initialize the data dictionary to store in the output
    data_dict = {}
    data_dict["energy"] = metadata

    for fname in files:
        data = np.loadtxt(fname, comments="#", delimiter=",")
        mdata = data[:,0]
        fdata = data[:,1]
        adata = data[:,2]
        
        #store everything into structured array
        muon_data = np.zeros(len(mdata), dtype = [('pdgid', int),
                                        ('e_fraction', float),
                                        ('theta', float),
                                        ('multiplicity', int)])
        muon_data['pdgid'] = np.repeat(mupdg, len(mdata))
        muon_data['e_fraction'] = fdata
        muon_data['theta'] = adata
        muon_data['multiplicity'] = mdata

        #get the energy idx
        idx = fname.split('_')[-1].split('.')[0]

        data_dict[idx] = muon_data

    #store the data into numpy compressed file
    outfname = options.OUT + "/" + pname + ".npz"
    np.savez_compressed(outfname, **data_dict)
