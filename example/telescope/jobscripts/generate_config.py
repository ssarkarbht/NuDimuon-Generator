#!/bin/python

#Generate copies of the base config file with different random seeds

import json

#get the base config file
#NuMu
infile = "base_config_numu.json"

#NuMuBar
#infile = "base_config_numubar.json"

#number of batches
nbatch = 50
#The location to put the newly generated config files
outdir = "/data/icecube/ssarkar/dimuon_generator/rawdata/icecube_run/config_batches"

#Load the base configuration
with open(infile, "r") as f:
    config = json.load(f)

#get the base seed value
base_seed = config["Settings"]["random_seed"]
base_output = config["Settings"]["out04_filename"].split(".")

#update the iterative batch configs
for i in range(nbatch):
    config["Settings"]["random_seed"] = base_seed + i + 1
    seed_str = str(config["Settings"]["random_seed"])

    #define a seed specific config file name
    outfname = seed_str+"_config.json"
    config["WorkDir"] = seed_str

    #update the final output file name with seed specific value
    config["Settings"]["out04_filename"] = base_output[0]+"_"+\
            seed_str+"."+base_output[-1]
    #write the new config file
    with open(outdir+"/"+outfname, "w") as f:
        json.dump(config, f, indent=4)
    f.close()
