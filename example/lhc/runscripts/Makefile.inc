# DIRE configuration file.
# Generated on Tue May 30 16:17:49 UTC 2023 with the user supplied options:
# --prefix=/usr/local
# --with-pythia8=/usr/local
# Configure string:
# ./configure --prefix=/usr/local --with-pythia8=/usr/local

# Install directory prefixes.
PREFIX_BIN=/usr/local/bin
PREFIX_INCLUDE=/usr/local/include
PREFIX_LIB=/usr/local/lib
PREFIX_SRC=/opt/DIRE/DIRE-2.002/src
PREFIX_SHARE=/usr/local/share/Dire

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=false
CXX=g++
CXX_COMMON=-O2 -std=c++98 -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so

# PYTHIA8 configuration.
PYTHIA8_USE=true
PYTHIA8_BIN=/usr/local/bin/
PYTHIA8_INCLUDE=/usr/local/include
PYTHIA8_LIB=/usr/local/lib
PYTHIA8_PREFIX=/usr/local
PYTHIA8_SHARE=/usr/local/share
PYTHIA8_EXAMPLES=/usr/local/share/Pythia8/examples

# MG5MES configuration.
MG5MES_USE=false
MG5MES_BIN=
MG5MES_INCLUDE=./
MG5MES_LIB=./

# OPENMP configuration.
OPENMP_USE=false
OPENMP_BIN=
OPENMP_INCLUDE=./
OPENMP_LIB=./
