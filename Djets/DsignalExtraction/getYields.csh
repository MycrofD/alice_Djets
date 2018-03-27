#!/bin/bash

# B. Trzeciak -- script to run D-jets analysis chain for the efficiency-corrected signal extraction

##### data output settings
dataFile=$1
lhcprod=$2
prod=$4
effFile=$6
refFile=$8
postfix=${10}
dirOut=${11}

################################################
############### D-jet signal
################################################

root -l -b -q signalExtraction_SB.C'("'$dataFile'", "'$lhcprod'",'$3', "'$prod'",'$5',"'$effFile'",'$7',"'$refFile'",1,'$9',"'$postfix'","'$dirOut'")' #> outfile 2>&1
