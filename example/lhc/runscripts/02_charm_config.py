#!/bin/python

'''
Description : Example script to sample the target nucleons
    and store the per event configuration in the event
    hdf5 file (for lepton injector) as well as in a separate
    text file to be used by the PYTHIA8+DIRE script
'''

import numpy as np
import h5py as h5
from numpy.random import default_rng
from argparse import ArgumentParser
import generator as gn
import json

#get the hdf5 file
parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="file",
        type=str, required=True,
        help="Which file do you want to parse?")
parser.add_argument("-c", "--config", dest="cfg",
        type=str, required=True,
        help="Which json config file did you use to produce the event file?")

args = parser.parse_args()
#load the input file
filename = args.file
hf = h5.File(filename, "r")

#load the config info
cfile = args.cfg
with open(cfile) as f:
    config = json.load(f)

#get the generation details
if config['Experiment']=='telescope':
    econfig = config['LeptonInjector']
    #get the event information object
    evlist = hf["RangedInjector0"]["initial"]

elif config['Experiment']=='lhc':
    econfig = config['Settings']
    evlist = hf["InitialType"]

seed = econfig['random_seed']
#get the total number of events generated
evtnum = len(evlist)
#create unique event IDs
evt_ids = np.arange(evtnum)
#get the primary neutrino energy
evt_ens = evlist['Energy']
#get the primary neutrino pdg id
evt_pdg = evlist['ParticleType']

#creat a blank array 
data_arr = np.zeros(evtnum, dtype = [('event_id', int), ('nu_pdg', int),
    ('nu_energy', float), ('target', int)])

#generate random target nucleons for pythia settings
np_rng = default_rng(seed)
sampler = gn.generator(np_rng)
evt_tgt = sampler.sample_nucleus(n=evtnum)
hf.close()

data_arr['event_id'] = evt_ids
data_arr['nu_pdg'] = evt_pdg
data_arr['nu_energy'] = evt_ens
data_arr['target'] = evt_tgt

hf = h5.File(filename,'a')
try:
    del hf['PythiaTarget']
    print ("Will overwrite existing 'PythiaTarget' object.")
except:
    print ("Will be writing sampled target configuration into 'PythiaTarget'")

d = hf.create_dataset('PythiaTarget', data=data_arr)
hf.close()

#Also create a separate text file for passing to Pythia+DIRE script
#outfile = str(seed)+"_charmConfig.txt"
outfile = config["Settings"]["out02_filename"]

np.savetxt(outfile, data_arr, fmt='%d, %d, %.8e, %d', delimiter=',')

