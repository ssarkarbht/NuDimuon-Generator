#!/bin/python

'''
This is an example script showing how to run the weight generation
for charm dimuon events and store it into outout i3 file.
'''

from icecube import icetray
import json
import numpy as np
from I3Tray import *
from weight_module import *

from optparse import OptionParser
#add command-line arguments
parser = OptionParser()
# I/O
parser.add_option("-i", "--InputI3File", dest="I3FILE", type=str)
parser.add_option("-l", "--InputLIFile", dest="LIFILE", type=str)
parser.add_option("-c", "--InputCharmFile", dest="CMFILE", type=str)
parser.add_option("-f", "--ConfigFilename", dest="CFILE", type=str)
parser.add_option("-o", "--OutputFilename", dest="OFILE", type=str)
# Other arguments
#interaction type: 0: standard CC DIS
# 1: Trident dimuon, 2: charm dimuon, 3: charm hadron
parser.add_option("-t", "--Interaction", dest="ITYPE", type=int,
        default = 2)
parser.add_option("-w", "--WeightDictName", dest="WDICT", type=str,
        default="I3MCWeightDict")

(options, args) = parser.parse_args()

#================= Compute and Store Event weights =================
tray = I3Tray()

tray.AddModule("I3Reader", 'reader', Filename = options.I3FILE)
tray.AddModule(GenerateWeight, 'weightgen',
        LIFile=options.LIFILE, ConfigFile=options.CFILE,
        InteractionType=options.ITYPE, CharmFile=options.CMFILE,
        WeightName=options.WDICT)
tray.AddModule("I3Writer", "writer", Filename = options.OFILE)
tray.Execute()
tray.Finish()
