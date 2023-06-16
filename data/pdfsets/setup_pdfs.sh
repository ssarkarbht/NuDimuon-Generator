#!/bin/bash

#go to the pdfset store directory
#CT14NLO
wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/CT14nlo.tar.gz
tar -xvzf CT14nlo.tar.gz
echo 13100 CT14nlo 1 >> pdfsets.index
#CT14lo
wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/CT14lo.tar.gz
tar -xvzf CT14lo.tar.gz
echo 13200 CT14lo 1 >> pdfsets.index
#CT18ANLO
wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/CT18ANLO.tar.gz
tar -xvzf CT18ANLO.tar.gz
echo 14600 CT18ANLO 1 >> pdfsets.index
#CT18NLO
wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/CT18NLO.tar.gz
tar -xvzf CT18NLO.tar.gz
echo 14400 CT18NLO 1 >> pdfsets.index



