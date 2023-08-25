#!/bin/python

'''
This script is used for generating the neutrino energies according to
a powerlaw spectrum specified in the configuration file.
'''

#import stuff
import numpy as np
import h5py as h5
import scipy.stats as st
import os
from optparse import OptionParser
import json

from generator import generator

#add command line input options
parser=OptionParser()
parser.add_option('-c', '--ConfigFile', dest='CFG', type=str)
#parse the arguments
(options, args) = parser.parse_args()

#load the configuration file
cfile = options.CFG
with open(cfile) as f:
    config = json.load(f)

#get the generation details
econfig = config['Standalone']

#set up the random seed details
randomGen = np.random.default_rng(econfig['random_seed'])

specgen = generator(randomGen)

event_energy = specgen.sample_energy(econfig['gamma'],
                            econfig['MinEnergy'],
                            econfig['MaxEnergy'],
                            n = econfig['event_number'])

#Prepare the preliminary event dataset

#Event IDs
event_idx = np.arange(econfig['event_number'])
#Neutrino PDGs
if len(econfig['nu_pdg'])==1:#if only neutrino OR anti-neutrino
    event_pdg = np.repeat(econfig['nu_pdg'][0],
                        econfig['event_number'])
elif len(econfig['nu_pdg'])==2:#if half-half (anti)neutrino
    event_pdg1 = np.repeat(econfig['nu_pdg'][0],
                    int(econfig['event_number']/2))
    event_pdg2 = np.repeat(econfig['nu_pdg'][1],
                    int(econfig['event_number']/2))

    event_pdg = np.append((event_pdg1, event_pdg2))

#create the output h5 file
h5f = h5.File(econfig['out_filename'], "w")
data_arr = np.zeros(econfig['event_number'],
        dtype = [('event_id', int),
            ('nu_pdg', int), ('nu_energy', float)])
data_arr['event_id'] = event_idx
data_arr['nu_pdg'] = event_pdg
data_arr['nu_energy'] = event_energy

h5f.create_dataset('InitialType', data=data_arr)
h5f.close()
