#!/bin/bash

#This is a wrapper script to run all steps of the simulation chain in sequence

#load some helpful function to run
source bash_func.sh

#get the json config file and define full path for unambiguity
gconfig=$1
curr_dir=$(echo `pwd`)
gconfig=$(echo $curr_dir/$gconfig)

#=========== [Step: 1] Prepare the working directory ========
varname=$(get_workdir $gconfig)
export WORKDIR=$varname
colored_echo "$GREEN" "Setting working directory: $WORKDIR"

mkdir -p $WORKDIR
cd $WORKDIR

#=========== [Step: 2] Start lepton injector run ============
colored_echo "$YELLOW" "Starting Lepton Injector Run"

output=$(python3 test_python.py 2>&1)
#output=$

fout=$(echo "$output" | awk -F '<<LIOUT>>' '{print $2}')
echo $output
echo "Did I get the filename out?"
echo $fout
#python3 inject_muons.py -c $gconfig
#python3 fix_primary.py -f 
#This will create a bunch of files in the working location.

#=========== [Step: 3] Configure Charm Production ===========


#=========== [Step: 4] Run Pythia Generator =================


#=========== [Step: 5] Run Muon Sampler =====================

#=========== [Step: 6] Event Merger and Muon Sampler ========

#=========== [Step: 7] Event Weight Calculator ==============

