#!/bin/python

'''
An example script of running the file converter 
with two input hdf5 files and merging event info
into the output i3 file.
'''

from file_converter import *
from icecube import icetray
from I3Tray import *

from optparse import OptionParser
#add command-line arguments
parser = OptionParser()
# I/O
parser.add_option("-i", "--InputLIFile", dest="LIFILE", type=str)
parser.add_option("-f", "--InputCharmFile", dest="CMFILE", type=str)
parser.add_option("-o", "--OutputFilename", dest="OFILE", type=str)

#Generator details
parser.add_option("-s", "--RandomSeed", dest="SEED", type=int)
parser.add_option("-m", "--MCTreeName", dest="MCTR", type=str,
        default='I3MCTree')

(options,args) = parser.parse_args()

#Assign a unique run ID (random seed) for the file
RUNID = options.SEED
#start event ID with 0
event_id=0
def get_header(frame):
        global event_id, RUNID
        header          = dataclasses.I3EventHeader()
        header.event_id = event_id
        header.run_id   = RUNID
        frame["I3EventHeader"]=header
        event_id+=1


tray = I3Tray()

tray.AddModule("I3InfiniteSource", "source", stream=icetray.I3Frame.DAQ)

tray.AddModule(get_header, streams = [icetray.I3Frame.DAQ])

tray.AddModule(H5H5I3Converter, 'convert', Filename=options.CMFILE,
                                        LIFilename=options.LIFILE)

tray.AddModule("I3Writer","write_mctree", FileName=options.OFILE)

tray.Execute()
tray.Finish()

