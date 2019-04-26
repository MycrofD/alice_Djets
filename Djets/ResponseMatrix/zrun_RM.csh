#!/bin/bash
#---------------------
# A.Mohanty
# Utrecht University
# auro.mohanty@cern.ch
#---------------------

# Response Matrix
#---------------------


isPrompt=$1
dataFile=$2
outDirRM=$3

root -l -b -q DetRM_z.C'('$isPrompt',"'$dataFile'","'$outDirRM'")'
root -l -b -q combineRM_z.C'('$isPrompt',"'$outDirRM'")'
