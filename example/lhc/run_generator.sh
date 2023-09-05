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

#[[STEP - 1]]
#Neutrino energy sampling
python3 01_inject_numu_standalone.py -c $gconfig
wait

#[[STEP - 2]]
#Charm configuration sampling
python3 02_charm_config.py -f 01_std.h5 -c $gconfig
wait

#[[STEP - 3]]
#Pythia run for charm hadron production
#get the output filename from config file
fvalue=$(jq '.Settings.out03_filename' "$gconfig")
fvalue="${fvalue#\"}"
fvalue="${fvalue%\"}"

make dire08
wait
./dire08 nu_ccdis.cmnd 02_cnf.txt $fvalue
wait


#[[STEP - 4]]
#Charm muon sampling from parametric tables
pyscript=$DIMUON_REPO/modules/get_charm_muons.py
#get the output filename from config file
ovalue=$(jq '.Settings.out04_filename' "$gconfig")
ovalue="${ovalue#\"}"
ovalue="${ovalue%\"}"

#get the parametric table directory
pvalue=$(jq '.CharmMuonSettings.ParamDir' "$gconfig")
pvalue="${pvalue#\"}"
pvalue="${pvalue%\"}"

#get the muon energy threshold
ethr=$(jq '.CharmMuonSettings.MuEThreshold' "$gconfig")
#get the random seed
seed=$(jq '.Settings.random_seed' "$gconfig")
#get the target medium
medval=$(jq '.Medium' "$gconfig")
medval="${medval#\"}"
medval="${medval%\"}"

python3 $pyscript -i $fvalue -o $ovalue -s $seed -p $pavlue -m $medval -t $ethr -f


#Transfer the final output
outdir=$(jq '.OutputDir' "$gconfig")
outdir="${outdir#\"}"
outdir="${outdir%\"}"

mv $ovalue $outdir/$ovalue

echo "Done."

