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

# Now, we'll make a new injector for muon tracks 
#get the outgoing particle type
if lconfig['finalType_1']=='MuMinus':
    final_1 = LI.Particle.ParticleType.MuMinus
elif lconfig['finalType_1']=='MuPlus':
    final_1 = LI.Particle.ParticleType.MuPlus

assert lconfig['finalType_2']=='Hadrons', "2nd outgoing particle is not hadrons!"
final_2     = LI.Particle.ParticleType.Hadrons

the_injector = LI.Injector( lconfig['event_number'] ,
        final_1, final_2, 
        lconfig['xs_folder']+'/test_xs.fits',
        lconfig['xs_folder']+'/test_xs_total.fits',
        lconfig['ranged_mode'])

deg = pi/180.

# define some defaults 
minE        = lconfig['MinEnergy']     # [GeV]
maxE        = lconfig['MaxEnergy']   # [GeV]
gamma       = lconfig['gamma']
minZenith   = lconfig['MinZenith']*deg
maxZenith   = lconfig['MaxZenith']*deg
minAzimuth  = lconfig['MinAzimuth']*deg
maxAzimuth  = lconfig['MaxAzimuth']*deg

# construct the controller
controller  = LI.Controller( the_injector, minE, maxE, gamma, minAzimuth, maxAzimuth, minZenith, maxZenith)

# specify the output, earth model
path_to = lconfig['earth_model']
controller.SetEarthModel("Planet", path_to)

#check if the output folder is defined in config file
outdir = lconfig['out_folder']
if outdir is None:
    #update the directory based on currect working directory location
    try:
        outdir = os.environ["WORKDIR"]
    except:
        assert False, "There's no output directory defined"
#add a slash at the end if not present
if outdir[-1]!='/': outdir += '/'

outfile = outdir + str(lconfig['random_seed']) + '_' + lconfig['out_filename']
licfile = outdir + str(lconfig['random_seed']) + '_' + lconfig['lw_filename']
controller.NameOutfile(outfile)
controller.NameLicFile(licfile)
controller.setSeed(lconfig['random_seed'])
# run the simulation
controller.Execute()

