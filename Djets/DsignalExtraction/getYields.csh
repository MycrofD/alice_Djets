#!/bin/bash

# B. Trzeciak (Utrecht University)
#-- script to run D-jets yield extraction

##### data output settings
dataFile=$1
isEff=$2
effFile=$3
isRefl=$4
refFile=$5
ispostfix=$6
postfix=$7
dirOut=$8
isMoreFile=$9
prod=${10}
################################################
############### D-jet signal
################################################

root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'")' #> outfile 2>&1
