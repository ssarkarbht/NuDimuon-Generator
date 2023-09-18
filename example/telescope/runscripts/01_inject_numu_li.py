#!/bin/python

'''
Example script to generate energy and geometry of the 
injected neutrino interactions around a cylindrical
volume for given energy spectrum. Simulation details
can be modified in the configuration (json) file 
'''

#import stuff
import LeptonInjector as LI
from math import pi
import os
from optparse import OptionParser
import json

#add command line input options
parser=OptionParser()
parser.add_option('-c', '--ConfigFile', dest='CFG', type=str)

#parse the arguments
(options, args) = parser.parse_args()

#load the configuration file
cfile = options.CFG
with open(cfile) as f:
    config = json.load(f)

#load Lepton Injector configuration
lconfig = config['LeptonInjector']
#load the general event settings
sconfig = config['Settings']

# Now, we'll make a new injector for muon tracks 
#get the outgoing particle type
if lconfig['finalType_1']=='MuMinus':
    final_1 = LI.Particle.ParticleType.MuMinus
elif lconfig['finalType_1']=='MuPlus':
    final_1 = LI.Particle.ParticleType.MuPlus
#Check if the 2nd particle is hadron
assert lconfig['finalType_2']=='Hadrons', "2nd outgoing particle is not hadrons!"
final_2     = LI.Particle.ParticleType.Hadrons

the_injector = LI.Injector(sconfig['event_number'] ,
        final_1, final_2, 
        lconfig['CrossSectionDir']+'/dsdxdy_nu_CC_iso.fits',
        lconfig['CrossSectionDir']+'/sigma_nu_CC_iso.fits',
        lconfig['RangeMode'])

#degree to radian conversion factor
rad = pi/180.
# define some defaults 
minE        = sconfig['MinEnergy']     # [GeV]
maxE        = sconfig['MaxEnergy']   # [GeV]
gamma       = sconfig['gamma']
minZenith   = lconfig['MinZenith']*rad
maxZenith   = lconfig['MaxZenith']*rad
minAzimuth  = lconfig['MinAzimuth']*rad
maxAzimuth  = lconfig['MaxAzimuth']*rad

# construct the controller
controller  = LI.Controller( the_injector, minE, maxE, gamma, 
        minAzimuth, maxAzimuth, minZenith, maxZenith,
        lconfig["InjectionRadius"], lconfig['EndcapLength'])

# specify the output, earth model
path_to = lconfig['EarthModelDir']
controller.SetEarthModel("Planet", path_to)

#check if the output folder is defined in config file
outdir = sconfig['out_folder']
if outdir is None:
    #update the directory based on currect working directory location
    outdir = ""
else:
    #add a slash at the end if not present
    if outdir[-1]!='/': outdir += '/'

outfile = outdir + sconfig['out01_filename']
licfile = outdir + str(sconfig['random_seed']) + '_' + lconfig['LICFilename']

controller.NameOutfile(outfile)
controller.NameLicFile(licfile)
controller.setSeed(sconfig['random_seed'])
# run the simulation
controller.Execute()

