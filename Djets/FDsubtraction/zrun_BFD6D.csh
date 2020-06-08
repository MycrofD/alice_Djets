#!/bin/bash
#---------------------
# A.Mohanty
# Utrecht University
# auro.mohanty@cern.ch
#---------------------

# B-FD subtraction
#---------------------
dataInFileDir=$1
outDir=$2
simDir=$3
dataAnalysisFile=$4
MCAnalysisFile=$5
isprefix=$6
ispostfix=$7
DjetEff=$8
efffile=$9

root -l  subtractFD_zjet.C'( "'$dataInFileDir'","'$outDir'","'$simDir'","'$dataAnalysisFile'","'$MCAnalysisFile'",'$isprefix','$ispostfix','$DjetEff',"'$efffile'")'
