#!/bin/python

'''
Author : Sourav Sarkar
Email : ssarkar1@ualberta.ca
Description : This script performs the CHROMO
hadronic interaction simulation of charm hadrons
with target nucleus of the detector medium.
'''

import chromo as ch
import particle
from chromo.util import classproperty, Nuclei
import numpy as np
import os,sys
import json

# get the repo directory path
repo_dir = os.environ["DIMUON_REPO"]
dataloc = repo_dir + "/data/constants_particles_materials/"

# get the materials and constants properties
pfile = dataloc+"particles.json"
with open(pfile, 'r') as f:
    particle_dict = json.load(f)

#get list of charm hadron pdg codes
pdglist = []
projlist = []
for key, val in particle_dict.items():
    pdgid = val['pdgid']
    pdglist.append(pdgid)
    projlist.append(particle.PDGID(pdgid))

charm_hadrons = set(projlist)

class UpdatedSibyll23d(ch.models.Sibyll23d):
    '''Update the base sibyll model with new projectiles
    '''
    _projectiles = charm_hadrons
    _targets = Nuclei()

    @classproperty
    def projectiles(cls):
        return cls._projectiles

    #WARNING! following defintiion forces the class to allow
    #heavier targets, however sibyll backend fortran
    #run may fail due to out of limit error
    @classproperty
    def targets(cls):
        return cls._targets

#get a reverse dictionary for charm hadron names
def get_name(pdgid):
    names = list(particle_dict.keys())
    pdglist = [val['pdgid'] for _,val in particle_dict.items()]
    return dict(zip(pdglist, names))

#set the event kinematics based on projectile and target
def set_event_kinematics(proj, target, energy):
    ''' Create a new event kinematics instance for each
    simulated MC event energy grid point.
    '''
    evkin = ch.kinematics.EventKinematics(proj, target,
                elab = energy * ch.constants.GeV)
    return evkin

def get_muon_frac_theta(event):
    '''get kinematics of the final state muons in the
    particle event list
    '''
    pvals = event.pid[np.absolute(event.pid)==13]
    xvals = event.xlab[np.absolute(event.pid)==13]
    avals = event.theta[np.absolute(event.pid)==13]
    # if there are more than one muon choose the highest energy
    if len(xvals)>0:
        idx = np.argmax(xvals)
        return (pvals[idx], xvals[idx], avals[idx], len(xvals))
    return None

def run_get_events(generator, nevents, event_k):
    '''function that takes the model generator and
    simulates nevents number of events with event_k
    kinematics and returns the output muon properties
    '''
    muon_data = np.zeros(nevents, dtype = [("pdgid", int),
                                ("e_fraction", float),
                                ("theta", float),
                                ("multiplicity", int)])
    generator.kinematics = event_k
    for idx, event in enumerate(generator(nevents)):
        event = event.final_state()
        getval = get_muon_frac_theta(event)
        if getval is not None:
            muon_data[idx] = getval
    return muon_data


if __name__ == '__main__':
    #parse command-line arguments
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-t", "--Target", dest="TGT", type=str)
    parser.add_option("-p", "--Projectile", dest="PRJ", type=int)
    parser.add_option("-n", "--NEvents", dest="NEVE", type=int)
    parser.add_option("-l", "--Emin", dest="EMIN", type=float)
    parser.add_option("-u", "--Emax", dest="EMAX", type=float)
    parser.add_option("-g", "--GridPoints", dest="NGRID", type=int)
    parser.add_option("-o", "--OutfileName", dest="OUT", type=str)

    (options, args) = parser.parse_args()

    mfile = dataloc+"medium_properties.json"
    with open(mfile, "r") as f:
        medium = json.load(f)
    assert options.TGT in medium.keys(), f"Target {options.TGT} not in medium"

    target_prop = medium[options.TGT]
    target = (target_prop["mass_num"], target_prop["atom_num"])
    proj = options.PRJ
    assert proj in pdglist, f"Particle PDG {proj} not in the projectile list"

    #get the energy array
    enarr = np.logspace(np.log10(options.EMIN),
            np.log10(options.EMAX), options.NGRID)
    #energy grid info to save as metadata
    metadata = np.array([options.EMIN, options.EMAX, options.NGRID])

    #intialize the model
    init_kin = ch.kinematics.EventKinematics(proj, target,
            elab = 1000 * ch.constants.GeV)
    generator = UpdatedSibyll23d(init_kin)

    outfile = options.OUT
    array_dict = {}
    #loop over the energy points
    for i, energy in enumerate(enarr):
        evkin = set_event_kinematics(proj, target, energy)
        dataset = run_get_events(generator,
                                options.NEVE , evkin)
        array_dict[str(i)] = dataset

    #save all the arrays into npz file
    np.savez_compressed(outfile, energy=metadata, **array_dict)
