#!/bin/bash

#setup the repo paths
cd /data/icecube/ssarkar/dimuon_generator/repo/NuDimuon-Generator
source setup.sh

#get the json config file and define full path for unambiguity
gconfig=$1
dir_io=$2
cd $dir_io
curr_dir=$(echo `pwd`)

#get the working directory name
config_dir=/data/icecube/ssarkar/dimuon_generator/rawdata/gen_standalone/config_batches
dirname=$(get_workdir $config_dir/$gconfig)
mkdir -p $curr_dir/$dirname

cd $curr_dir/$dirname

#copy the scripts
cp $DIMUON_REPO/example/standalone/* .
#copy the config file
cp $config_dir/$gconfig .

#Start the runs
#Neutrino energy sampling
python3 01_inject_numu_standalone.py -c $gconfig
wait

#Charm configuration sampling
python3 02_charm_config.py -f 01_std.h5 -c $gconfig
wait

#Pythia run for charm hadron production
#get the output filename from config file
fvalue=$(jq '.Standalone.out03_filename' "$gconfig")
fvalue="${fvalue#\"}"
fvalue="${fvalue%\"}"

make dire08
wait
./dire08 nu_ccdis.cmnd 02_cnf.txt $fvalue
wait

#Transfer the final output
outdir=$(jq '.OutputDir' "$gconfig")
outdir="${outdir#\"}"
outdir="${outdir%\"}"

mv $fvalue $outdir/$fvalue

echo "Done."

