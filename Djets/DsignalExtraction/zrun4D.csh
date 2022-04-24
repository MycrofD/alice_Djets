#! /bin/bash

#--------------------------------------------------
#  Author A.Mohanty
#  Utrecht University
#  auro.mohanty@cern.ch
#--------------------------------------------------

R=$1
OUT=${EOS_local}/media/jackbauer/data/z_out/R_0$R
#####################################
source zrun_settings.csh
######### for FD
signalExtractionDir=$OUT/signalExtraction/plots
FDsubtractionDir=$OUT/FDsubtraction
simDir=$OUT/SimFiles/BFeedDown
DjetEff=1
prompteff=$OUT/efficiency/DjetEff_prompt_jetpt
######### other unfolding parameters
unfoldingDir=$OUT/unfolding
regPar=5
isPrior=0
priorType=1 #1 to 8
isFDUpSys=0
isFDDownSys=0
NTrials=10

## B FEED DOWN=================================
if [ $doBFD -eq 1 ]; then
 cd ../FDsubtraction
  root -l  subtractFD_zjet.C'( "'$listName'","'$signalExtractionDir'","'$FDsubtractionDir'","'$simDir'","'$data'","'$effFile'",'$isprefix','$ispostfix','$DjetEff',"'$prompteff'")'
 cd ../DsignalExtraction
fi
## B FEED DOWN ends============================

## UNFOLDING=================================
if [ $doUnfold -eq 1 ]; then
 cd ../unfolding
  root -l  unfold_Bayeszjet.C'("'$listName'",'$ispostfix','$isprefix','$DjetEff',"'$prompteff'","'$effFile'","'$FDsubtractionDir'","'$unfoldingDir'",'$regPar',0,'$priorType','$isFDUpSys','$isFDDownSys','$NTrials')'

 if [ $isPrior -eq 1 ]; then
  for thing in 1 2 3 4 5 6 7 8; do #these are priors
   root -l  unfold_Bayeszjet.C'("'$listName'",'$ispostfix','$isprefix','$DjetEff',"'$prompteff'","'$effFile'","'$FDsubtractionDir'","'$unfoldingDir'",'$regPar','$isPrior','$priorType','$isFDUpSys','$isFDDownSys','$NTrials')'
  done
 fi
 cd ../DsignalExtraction
fi
## UNFOLDING ends=============================
