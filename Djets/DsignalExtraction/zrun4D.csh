#! /bin/bash

#--------------------------------------------------
#  Author A.Mohanty
#  Utrecht University
#  auro.mohanty@cern.ch
#--------------------------------------------------

doBFD=0
doUnfold=1

R=$1
EOS_local=$2
#/eos/user/a/amohanty/
OUT=${EOS_local}/media/jackbauer/data/z_out/R_0$R
#OUT=${OUT}_25
OUT=${OUT}_finaltry

##################################### THINGS TO CHANGE
flagCUT=0
cutfileNo=0
#####################################
source zrun_settings.csh
######### for FD
signalExtractionDir=$OUT/signalExtraction/plots
FDsubtractionDir=$OUT/FDsubtraction
simDir=$OUT/SimFiles/BFeedDown
isprefix=0
ispostfix=0
DjetEff=1
efffile=$OUT/efficiency/DjetEff_prompt_jetpt
######### other unfolding parameters
detRMFilePrompt=$OUT/ResponseMatrix/DetMatrix_prompt
bkgRMFile=$OUT/ResponseMatrix/DetMatrix_
unfoldingDir=$OUT/unfolding
regPar=5
isPrior=0
priorType=1 #1 to 8
isBkgRM=0
isFDUpSys=0
isFDDownSys=0

isPrompt=1
listName="FD"

## B FEED DOWN=================================
if [ $doBFD -eq 1 ]; then
 cd ../FDsubtraction
  root -l  subtractFD_zjet.C'( "'$signalExtractionDir'","'$FDsubtractionDir'","'$simDir'","'$data'","'$effFile'",'$isprefix','$ispostfix','$DjetEff',"'$efffile'")'
 cd ../DsignalExtraction
fi
## B FEED DOWN ends============================

## UNFOLDING=================================
if [ $doUnfold -eq 1 ]; then
 cd ../unfolding
  root -l  unfold_Bayeszjet.C'("'$listName'",'$isPrompt','$ispostfix','$isprefix',"'$effFile'","'$FDsubtractionDir'","'$detRMFilePrompt'","'$bkgRMFile'","'$unfoldingDir'",'$regPar',0,'$priorType','$isBkgRM','$isFDUpSys','$isFDDownSys')'

 if [ $isPrior -eq 1 ]; then
  for thing in 1 2 3 4 5 6 7 8; do #these are priors
   root -l  unfold_Bayeszjet.C'("'$listName'",'$isPrompt','$ispostfix','$isprefix',"'$effFile'","'$FDsubtractionDir'","'$detRMFilePrompt'","'$bkgRMFile'","'$unfoldingDir'",'$regPar','$isPrior','$priorType','$isBkgRM','$isFDUpSys','$isFDDownSys')'
  done
 fi
 cd ../DsignalExtraction
fi
## UNFOLDING ends=============================
