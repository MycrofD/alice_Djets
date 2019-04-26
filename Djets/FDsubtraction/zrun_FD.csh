#!/bin/bash
#---------------------
# A.Mohanty
# Utrecht University
# auro.mohanty@cern.ch
#---------------------

# B-FD subtraction
#---------------------

dataFile=$1
dataAnalysisFile=$2
OUT=$3
simDir=$OUT/SimFiles/BFeedDown
comMatrixFile=$OUT/ResponseMatrix/combineMatrixFD
outSpectraDir=$OUT/FDsubtraction

root -l -b -q subtractFD_z.C'("'$dataFile'","'$dataAnalysisFile'","'$simDir'","'$comMatrixFile'","'$outSpectraDir'")'
