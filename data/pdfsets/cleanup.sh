#!/bin/bash

#specify the files that will not be deleted

keep=("setup_pdfs.sh" "cleanup.sh")

#remove rest of the files in the directory

find "." ! -name "${keep[0]}" ! -name "${keep[1]}" -exec rm -r {} \;
