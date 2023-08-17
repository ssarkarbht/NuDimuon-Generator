#!/bin/bash

# Text color codes
RED="\033[0;31m"
GREEN="\033[0;32m"
YELLOW="\033[0;33m"
BLUE="\033[0;34m"
RESET="\033[0m"  # Reset text formatting

# Usage: colored_echo "color_code" "message"
colored_echo() {
  local color="$1"
  local message="$2"
  echo -e "${color}${message}${RESET}"
}

#colored_echo "$RED" "This is a red message."
#colored_echo "$GREEN" "This is a green message."
#colored_echo "$YELLOW" "This is a yellow message."
#colored_echo "$BLUE" "This is a blue message."

get_workdir() {
  local file="$1"
  value=$(jq '.WorkDir' "$file")
  if [[ "$value" == "null" ]]; then
          colored_echo "$RED" "Workdir not defined in config file" >&2
          colored_echo "$YELLOW" "Provide full path of a working directory: " >&2
          read workdir
          echo $workdir
  else
          value="${value#\"}"
          echo "${value%\"}"
  fi
}

#------ set up path variables
# module directory
# data directory
# output file directory


echo "Setting up paths..."
#------ set up repository path
export DIMUON_REPO=$(echo `pwd`)
#repo parent directory path is required for some python module and data loading to work

#add module directory to python path
export PYTHONPATH=$PYTHONPATH:$DIMUON_REPO/modules
#------ set up lhapdf paths
export LHAPDF_DATA_PATH=$DIMUON_REPO/data/pdfsets

# check if the pdfset folder exists, if they do don't download again
# if the folders or the tar files don't exist, download and unzip them

#------ make the pythia scripts
#using the make files compile the dire08.cc
#also, if needed to use the particle gun, compile standolone pythia particle gun script
echo Done.
