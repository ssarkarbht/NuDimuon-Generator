#!/bin/python

import json

infile = "base_config_numubar.json"
nbatch = 50
outdir = "/data/icecube/ssarkar/dimuon_generator/rawdata/gen_standalone/config_batches"

with open(infile, "r") as f:
    config = json.load(f)

base_seed = config["Standalone"]["random_seed"]

for i in range(nbatch):
    config["Standalone"]["random_seed"] = base_seed + i + 1
    seed_str = str(config["Standalone"]["random_seed"])
    outfname = seed_str+"_config.json"
    config["WorkDir"] = seed_str
    config["Standalone"]["out03_filename"] = "03_charm_"+seed_str+".txt"

    with open(outdir+"/"+outfname, "w") as f:
        json.dump(config, f, indent=4)
    f.close()
#print (config)
