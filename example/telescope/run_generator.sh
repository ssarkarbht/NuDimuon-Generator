#!/bin/bash

#setup the repo paths
cd /data/icecube/ssarkar/dimuon_generator/repo/NuDimuon-Generator
source setup.sh

#get the json config file
gconfig=$1

#get the parent directory to start working
dir_io=$(jq '.ParentDir' "$gconfig")
dir_io="${dir_io#\"}"
dir_io="${dir_io%\"}"

#get the working directory name (not abosolute path)
dirname=$(get_workdir $gconfig)
mkdir -p $dir_io/$dirname

cd $dir_io/$dirname

#copy the scripts (in the runscripts directory)
cp $DIMUON_REPO/example/telescope/runscripts/* .
#copy the config file
cp $gconfig .

#Start the Lepton Injector runs
#[[STEP - 1.1]]
#Neutrino energy and cylindrical geometry sampling
python3 01_inject_numu_li.py -c $gconfig
wait

#[[STEP - 1.2]]
#fix the primary pdg code in the LI file
fvalue=$(jq '.Settings.out01_filename' "$gconfig")
fvalue="${fvalue#\"}"
fvalue="${fvalue%\"}"

python3 01_fix_primary.py

#Charm configuration sampling
python3 02_charm_config.py -f $fvalue -c $gconfig
wait

#[[STEP - 3]]
#get the input filename from config file
ivalue=$(jq '.Settings.out02_filename' "$gconfig")
ivalue="${ivalue#\"}"
ivalue="${ivalue%\"}"

#get the output filename from config file
fvalue=$(jq '.Settings.out03_filename' "$gconfig")
fvalue="${fvalue#\"}"
fvalue="${fvalue%\"}"

#compile the script if not already have executable
make dire08
wait
#Pythia run for charm hadron production
./dire08 nu_ccdis.cmnd $ivalue $fvalue
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

python3 $pyscript -i $fvalue -o $ovalue -s $seed -p $DIMUON_REPO/$pvalue -m $medval -t $ethr -f
wait

#Transfer the final output
outdir=$(jq '.OutputDir' "$gconfig")
outdir="${outdir#\"}"
outdir="${outdir%\"}"

mv $ovalue $outdir/$ovalue

#Clean the directory scripts copies
shopt -s extglob
rm !(*.txt|*.h5)


echo "Done."

