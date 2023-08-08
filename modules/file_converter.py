#!/usr/bin/env python

'''
Author: Sourav Sarkar
Date: 12 Sept, 2020
Description: This module takes the hdf5 file containing 
	list of merged events (input and output), events
	geometry and makes coordinate transformation
	(rotation and translation) of particle position
	and direction and puts it in an I3Frame.
'''

import numpy as np
import h5py as h5
import sys, os, json
try:
    from icecube import icetray, dataio, dataclasses
    from I3Tray import *
except:
    print ("IceTray environment not loaded.")
    sys.exit(1)

#Load the constant dictionary
#define the directory paths
repo_dir = os.environ["DIMUON_REPO"]
data_dir = repo_dir + "/data/"

#load the constants
cfile = data_dir + "constants_particles_materials/constants.json"
with open(cfile, 'r') as f:
    CONSTANT = json.load(f)

#Construct the rotation metrices
# theta rotation around x-axis
def x_rot(theta):
    '''Rotate theta amount around X-axis
    '''
    comp_x = np.array([1., 0., 0.])
    comp_y = np.array([0., np.cos(theta), np.sin(theta)])
    comp_z = np.array([0., -1.*np.sin(theta), np.cos(theta)])
    return np.vstack((comp_x,comp_y,comp_z))

# (pi/2 + phi) rotation around z-axis
def z_rot(phi):
    '''Rotate phi amount around Z-axis
    '''
    comp_x = np.array([np.sin(phi), np.cos(phi), 0.])
    comp_y = np.array([-1.*np.cos(phi), np.sin(phi), 0.])
    comp_z = np.array([0., 0., 1.])
    return np.vstack((comp_x,comp_y,comp_z))

#functions to convert zenith and azimuth into
#theta and phi

def zenith_to_theta(zenith):
    return np.pi-zenith
def azimuth_to_phi(azimuth):
    phi = np.pi + azimuth
    if isinstance(azimuth, np.ndarray):
        modidx = np.where(phi>=2*np.pi)[0]
        phi[modidx] = phi[modidx] - 2*np.pi
    else:
        if phi >= 2*np.pi: phi -= 2*np.pi
    return phi

class H5I3Converter(icetray.I3Module):
    '''
    This i3 class module takes the hdf5 file (kinematics+geometry)
    and converts events into I3MCTree
    '''
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)

        self.AddParameter('Filename', #name of the hdf5 file
				'HDF5 filename',
				'sample_file.h5')
        self.AddParameter('MCTreeName', #name of the I3MCTree object
				'Name of the output MCTree object',
				'I3MCTree')

    def Configure(self):
        #get the data
        h5filename = self.GetParameter('Filename')
        h5file = h5.File(h5filename,'r')
        self.geometry = h5file['EventGeometry']
        self.events   = h5file['EventParticleList']
        #get event IDs and total simulated events
        self.keys = list(self.events.keys())
        self.nevents = len(self.geometry.keys())
        #initialize the mctree object name to be injected in the frame
        self.mctreename = self.GetParameter('MCTreeName')
        # track the iteration of events and failed charm muon events
        self.iterator = 0
        self.failed_events = []

    def DAQ(self, frame):
        #print progress
        if (self.iterator+1)%1000==0:
            print (f"\rProcessing frame number: {self.iterator}", end="\r")
        # the key is not in event list, update failed event
        # and skip the event injection
        if str(self.iterator) not in self.keys:
            self.failed_events.append(self.iterator)
            self.iterator += 1
        else:#for all other events with charm muon
            mctree = self.create_mctree()
            frame.Put(self.mctreename, mctree)
            self.PushFrame(frame)
            self.iterator += 1
        #if the iterator hits the end idx, request suspension
        if self.iterator>=self.nevents: self.RequestSuspension()
		
    def create_mctree(self):
        '''Takes the particle and geomtry data of the event and 
        produces the I3MCTree object.
        '''
        #create the event idx from the iterator
        idx = str(self.iterator)
        #initialize an empty MCTree
        mctree = dataclasses.I3MCTree()
        #Load the Lepton Injector Geometry Info
        pos_x = self.geometry[idx][0]
        pos_y = self.geometry[idx][1]
        pos_z = self.geometry[idx][2]
        theta = self.geometry[idx][3]
        phi   = self.geometry[idx][4]
        #compute the rotational matrices
        xrot  = x_rot(theta)
        zrot  = z_rot(phi)
        #iterate over the list of particle in the event
        for p in self.events[idx]:
            particle = dataclasses.I3Particle()
            #get the particle geometry (along Z-axis)
            p_theta = p['theta']
            p_phi   = p['phi']
            ppx = np.sin(p_theta)*np.cos(p_phi)
            ppy = np.sin(p_theta)*np.sin(p_phi)
            ppz = np.cos(p_theta)
            pvec = np.array([ppx, ppy, ppz])[:, np.newaxis]
            # perform coordinate transformation
            px,py,pz = np.matmul(zrot,np.matmul(xrot,pvec)).flat
            # Map the generic hadron code from h5 to i3
            if p['pdg_id'] == 9999999999:
                #generic hadron pdg conversion
                particle.pdg_encoding = int(-2000001006)
            else:
                particle.pdg_encoding = int(p['pdg_id'])

            #add rest of the particle info
            particle.location_type_string = 'InIce'
            particle.dir = dataclasses.I3Direction(px, py, pz)
            # get the updated particle direction
            new_theta = particle.dir.theta
            new_phi   = particle.dir.phi
            # Compute the positional (+time) offset from the primary vertex
            offx = 0.0
            offy = 0.0
            offz = 0.0
            offt = 0.0
            if p['dist']!=0:
                offx = p['dist']*np.sin(new_theta)*np.cos(new_phi)
                offy = p['dist']*np.sin(new_theta)*np.sin(new_phi)
                offz = p['dist']*np.cos(new_theta)
                offt = (p['dist'] * icetray.I3Units.s
                        /CONSTANT['light_speed']['value']
                        /icetray.I3Units.cm)

            #update the particle position
            particle.pos = dataclasses.I3Position(pos_x+offx,
                    pos_y+offy, pos_z+offz)

            particle.energy = p['energy'] * icetray.I3Units.GeV
            particle.time = offt * icetray.I3Units.ns
            particle.length = float('Nan')

            if p['ptype']==1:
                particle.shape_string = 'Primary'
                particle.location_type_string = 'InIce'
                mctree.add_primary(particle)
                primary = particle
            elif p['ptype']==2:
                particle.shape_string = 'Null'
                particle.location_type_string = 'InIce'
                mctree.append_child(primary, particle)
            else: assert False, "Unknown particle type in HDF5 file."
        return mctree


class H5H5I3Converter(H5I3Converter):
    ''' This class inherits most of the methods
    from H5I3Converter with a modification to the
    input data structure streams.
    '''
    def __init__(self, context):
        ''' Add an additional input hdf5 file option
        to get the Lepton Injector event geometry from
        an additional file.
        '''
        super().__init__(context)
        self.AddParameter('LIFilename', #name of the initial LI file
                'Lepton Injector HDF5 filename',
                'li_file.h5')

    def Configure(self):
        ''' Overwrite the superclass method.
        '''
        #Get the particle kinematics from the first hdf5 file
        h5filename = super().GetParameter('Filename')
        h5file = h5.File(h5filename,'r')
        self.events   = h5file['EventParticleList']
        self.keys = list(self.events.keys())
        #get the geometry from the second hdf5 file
        liname = self.GetParameter('LIFilename')
        lifile = h5.File(liname, 'r')
        evprop = lifile['RangedInjector0']['properties']
        self.nevents = len(evprop)
        #update the geometry object from the new file
        geo_keys = np.arange(self.nevents).astype(str)
        geo_vals = np.hstack((evprop['x'].reshape(-1,1),
                        evprop['y'].reshape(-1,1),
                        evprop['z'].reshape(-1,1),
                        zenith_to_theta(evprop['zenith']).reshape(-1,1),
                        azimuth_to_phi(evprop['azimuth']).reshape(-1,1)))
        self.geometry = dict(zip(geo_keys, geo_vals))
        #initialize the mctree object name to be injected in the frame
        self.mctreename = super().GetParameter('MCTreeName')
        # track the iteration of events and failed charm muon events
        self.iterator = 0
        self.failed_events = []


