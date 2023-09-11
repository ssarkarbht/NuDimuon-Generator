#!/bin/bash

#config files
curr_dir=$(echo `pwd`)
cloc="/data/icecube/ssarkar/dimuon_generator/rawdata/lhc_run/config_batches"
#cd $cloc
#cfiles="*.json"
dagname="faserv_numu+numubar_2M.dag"
num=0
for i in $cloc/*.json
do
	echo "JOB job$num submitter.sub" >> $dagname
	echo "CATEGORY job$num GENERIC" >> $dagname
	echo "VARS job$num name=\"$num\" args=\"$i\"" >> $dagname
	echo "" >> $dagname
	num=$(($num+1))
done

#mv $dagname $curr_dir/$dagname
