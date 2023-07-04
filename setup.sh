#!/bin/bash

#------ set up path variables
# module directory
# data directory
# output file directory

echo "Setting up paths..."
#------ set up repository path
export DIMUON_REPO=$(echo `pwd`)
#repo parent directory path is required for some python module and data loading to work

#------ set up lhapdf paths
export LHAPDF_DATA_PATH=$DIMUON_REPO/data/pdfsets

# check if the pdfset folder exists, if they do don't download again
# if the folders or the tar files don't exist, download and unzip them

#------ make the pythia scripts
#using the make files compile the dire08.cc
#also, if needed to use the particle gun, compile standolone pythia particle gun script
echo Done.
