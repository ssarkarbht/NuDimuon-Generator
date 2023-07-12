#!/bin/python

import numpy as np
from glob import glob
import h5py as h5
from optparse import OptionParser
from math import *
from charm_muon_generator import *
import sys, os

def get_particle(line):
    ''' For a given particle line string
    extract the particle properties
    '''
    info   = line.split()
    pdg    = int(info[2])
    energy = float(info[3])
    theta  = float(info[4])
    phi    = float(info[5])

    return (pdg, energy, theta, phi)

class GenerateEvents(CharmMuonGenerator):
    ''' This class takes the input and output files,
    and the inputs for the inherited class for the
    charm muon generator and creates the final event
    list by generating secondary muons from the charm
    hadrons.
    '''
    def __init__(self, infile, outfile, paramdir, xsfile,
            seed, medium, ethreshold):
        super().__init__(paramdir, xsfile,
                seed, target=medium, mu_emin=ethreshold)

        # initialize the output file data structure
        self.hfile = h5.File(outfile, 'a')
        #if the file contains stuff, overwrite
        try:
            self.hgroup = self.hfile.create_group('EventParticleList')
        except:
            del self.hfile['EventParticleList']
            self.hgroup = self.hfile.create_group('EventParticleList')

        # load the input event file data
        f = open(infile, 'r')
        self.lines = np.asarray(f.readlines())
        # Read the lines with 'E' tag
        zeroid = np.char.find(self.lines, 'E')
        #select only the ids where 'E' tag is at the beginning
        #of the line
        self.idxs = np.where(zeroid==0)[0]

        #filterout the missing events (where there is no charm hadrons)
        start_id = self.idxs
        end_id   = np.append(self.idxs[1:], len(self.lines))
        #number of particles in each of the events
        self.npart    = end_id - start_id
        #Events with 4 or more lines have produced charm hadron
        self.filter = self.npart>=4
        #Check how many events are missing charm hadron
        nmiss = len(np.where(self.filter==False)[0])
        print (f"Missed {nmiss} events in {infile}")

        self.weight_dict = np.zeros(len(self.idxs), dtype=[('eventID', int),
                            ('nuEnergy', float),
                            ('hadronPDG', int),
                            ('hadronEnergy', float),
                            ('muonOrigin', int),
                            ('xsFraction', float),
                            ('branchingRatio', float),
                            ('totalXsFraction', float)])

    def get_event_info(self, idx):
        '''Get the initial information on the event
        '''
        event_info = self.lines[idx].split()
        nupdg = int(event_info[1])
        eid = event_info[2]
        nueng = float(event_info[3])
        return eid, nupdg, nueng

    def one_charm(self, idx, offset=0):
        '''If the event has only one charm hadron, call
        this function to generate charm muon.
        offset is used in the case for more than one charm
        hadron
        '''
        #initialize a blank list for particles
        part_list = []
        #get initial event info
        eid, nupdg, nuen = self.get_event_info(idx)

        #get the charm hadron
        charm = get_particle(self.lines[idx+3+offset])
        #generate a muon
        mu_sample = super().sample_muon(charm[0], charm[1])
        if mu_sample is None: return None, None

        #build the final particle list
        #tuple structure :
        # (tree_id, pdgid, energy, theta, phi, r/distance)
        #insert the primary incoming neutrino
        prim_nu = (1, nupdg, nuen, 0, 0, 0)
        part_list.append(prim_nu)

        #insert the primary outgoing muon (or other lepton)
        prim_mu = get_particle(self.lines[idx+1])
        part_list.append((2, prim_mu[0], prim_mu[1],
                    prim_mu[2], prim_mu[3], 0))

        #generate the secondary muon
        seco_mu = (2, int(np.sign(nupdg)*-13), mu_sample[0],
                    charm[2], charm[3], mu_sample[4])
        part_list.append(seco_mu)

        #add the rest energy as generic hadrons
        charm_quark = get_particle(self.lines[idx+2])
        #first hadron at the primary vertex
        en1 = charm_quark[1] - charm[1]
        hadron1 = (2, 9999999999, en1, charm_quark[2], charm_quark[3], 0)
        #second hadron at the charm hadron length
        en2 = charm[1] - mu_sample[0]
        hadron2 = (2, 9999999999, en2, charm[2], charm[3], mu_sample[4])

        #add the hadrons to the list
        part_list.append(hadron1)
        part_list.append(hadron2)

        #now fill the weight info
        xs = super().compute_fractional_xs(nuen, np.sign(nupdg))
        if mu_sample[2] is None: og = -1
        else: og = MEDIUM[mu_sample[2]]['atom_num']
        weights = (int(eid), nuen, charm[0], charm[1], og, xs,
                mu_sample[3], xs*mu_sample[3])

        return part_list, weights

    def many_charm(self, i, idx):
        ''' In case of more than one charm hadrons in the
        event list, call this function.
        '''
        #get the idx offset for the leading charm hadron in energy
        nhadrons = self.npart[i] - 3 #no. of charm hadrons
        offset = 0
        ch_en = -1.
        for j in range(nhadrons):
            chad = get_particle(self.lines[idx+3+j])
            if chad[1]>ch_en:
                offset = j
                ch_en = chad[1]

        #generate the usual particle list considering the leading
        #charm hadron
        part_list, weights = self.one_charm(idx, offset=offset)

        if part_list is None: return None, None
        #check if there's more muon to add from other charm hadrons
        new_muen = 0
        new_mulist = []
        for k in range(nhadrons):
            if k==offset: continue
            new_charm = get_particle(self.lines[idx+3+k])
            #check if we need to generate a muon for this charm hadron
            if super().inject_muon(new_charm[0], new_charm[1]):
                new_mu_sample = super().sample_muon(new_charm[0], new_charm[1])
                if new_mu_sample is not None:
                    new_mu = (2, 13, new_mu_sample[0], new_charm[2],
                            new_charm[3], new_mu_sample[4])
                    new_mulist.append(new_mu)
                    new_muen += new_mu[2]
        #if additional muon is generated, update the particle list
        if new_muen>0:
            assert len(new_mulist)>0, "Additional muon must exist!"
            #get the first hadron to update
            old_had = part_list[-2]
            new_had = (old_had[0], old_had[1], old_had[2]-new_muen,
                    old_had[3], old_had[4], old_had[5])
            #delete the old hadron
            del part_list[-2]
            #add the updated hadron
            part_list.append(new_had)
            for nmu in new_mulist:
                part_list.append(nmu)

        return part_list, weights

    def run_generator(self):
        ''' Wrapper function to iterate over all the events
        and add the output particles in the output file
        '''
        for i, idx in enumerate(self.idxs):
            if not(self.filter[i]): continue
            print (f'\rProcessing Event: {i}/{len(self.idxs)}', end='\r')
            if self.npart[i]==4:
                partlist, weights = self.one_charm(idx)
            elif self.npart[i]>4:
                partlist, weights = self.many_charm(i, idx)
            else: 
                print ('WARNING: No Charm Hadron in Eventlist Found!')
                continue
            #store the outputs
            if partlist is not None:
                key = str(weights[0])
                event_arr = np.zeros(len(partlist), dtype=[('ptype', int),
                ('pdg_id', int), ('energy', float), ('theta', float),
                ('phi', float), ('dist', float)])
                for ip, p in enumerate(partlist):
                    event_arr[ip] = p
                self.hgroup.create_dataset(key, data=event_arr)
                self.weight_dict[i] = weights
        #write the full weight dictionary, overwrite if exists
        try:
            self.hfile.create_dataset('CharmWeights', data=self.weight_dict)
        except:
            del self.hfile['CharmWeights']
            self.hfile.create_dataset('CharmWeights', data=self.weight_dict)
        #close the h5file
        self.hfile.close()
        return None

if __name__=='__main__':
    from optparse import OptionParser
    #add command-line arguments
    parser = OptionParser()
    # I/O
    parser.add_option("-i", "--InputFilename", dest="IFILE", type=str)
    parser.add_option("-o", "--OutputFilename", dest="OFILE", type=str)
    # Pregenerated Data
    parser.add_option("-p", "--ParamDir", dest="PDIR", type=str,
            default= data_dir+"charm_data/muon_tables")
    parser.add_option("-x", "--CrossSection", dest="CROSS", type=str,
            default= data_dir+"charm_data/cross_sections/CC_Charm_Fraction.txt")

    #Generator details
    parser.add_option("-s", "--RandomSeed", dest="SEED", type=int)
    parser.add_option("-m", "--TargetMedium", dest="MED", type=str)
    parser.add_option("-t", "--MuThreshold", dest="THR", type=float,
            default=1e1)

    (options,args) = parser.parse_args()

    gen_ev = GenerateEvents(options.IFILE, options.OFILE,
                        options.PDIR, options.CROSS, options.SEED,
                        options.MED, options.THR)
    gen_ev.run_generator()
    #Done.
