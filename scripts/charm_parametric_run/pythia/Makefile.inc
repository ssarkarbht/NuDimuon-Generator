# PYTHIA configuration file.
# Generated on Tue May 30 16:14:56 UTC 2023 with the user supplied options:
# --prefix=/usr/local
# --with-hepmc2=/usr/local
# --with-hepmc2-include=/usr/local/include
# --with-lhapdf6=/usr/local
# --with-lhapdf6-include=/usr/local/include
# --with-gzip=/usr/local

# Install directory prefixes.
PREFIX_BIN=/usr/local/bin
PREFIX_INCLUDE=/usr/local/include
PREFIX_LIB=/usr/local/lib
PREFIX_SHARE=/usr/local/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=false
CXX=g++
CXX_COMMON=-O2  -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so

# EVTGEN configuration.
EVTGEN_USE=false
EVTGEN_BIN=
EVTGEN_INCLUDE=./
EVTGEN_LIB=./

# FASTJET3 configuration.
FASTJET3_USE=false
FASTJET3_BIN=
FASTJET3_INCLUDE=./
FASTJET3_LIB=./

# HEPMC2 configuration.
HEPMC2_USE=true
HEPMC2_BIN=/usr/local/bin/
HEPMC2_INCLUDE=/usr/local/include
HEPMC2_LIB=/usr/local/lib

# HEPMC3 configuration.
HEPMC3_USE=false
HEPMC3_BIN=
HEPMC3_INCLUDE=./
HEPMC3_LIB=./

# LHAPDF5 configuration.
LHAPDF5_USE=false
LHAPDF5_BIN=
LHAPDF5_INCLUDE=./
LHAPDF5_LIB=./

# LHAPDF6 configuration.
LHAPDF6_USE=true
LHAPDF6_BIN=/usr/local/bin/
LHAPDF6_INCLUDE=/usr/local/include
LHAPDF6_LIB=/usr/local/lib

# POWHEG configuration.
POWHEG_USE=false
POWHEG_BIN=
POWHEG_INCLUDE=./
POWHEG_LIB=./

# PROMC configuration.
PROMC_USE=false
PROMC_BIN=
PROMC_INCLUDE=./
PROMC_LIB=./

# ROOT configuration.
ROOT_USE=false
ROOT_BIN=
ROOT_INCLUDE=./
ROOT_LIB=./

# GZIP configuration.
GZIP_USE=true
GZIP_BIN=/usr/local//
GZIP_INCLUDE=/usr/local/include
GZIP_LIB=/usr/local/lib

# BOOST configuration.
BOOST_USE=false
BOOST_BIN=
BOOST_INCLUDE=./
BOOST_LIB=./

# PYTHON configuration.
PYTHON_USE=false
PYTHON_BIN=
PYTHON_INCLUDE=./
PYTHON_LIB=./
