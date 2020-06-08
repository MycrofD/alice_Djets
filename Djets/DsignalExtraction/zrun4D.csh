#! /bin/bash

#--------------------------------------------------
#  Author A.Mohanty
#  Utrecht University
#  auro.mohanty@cern.ch
#--------------------------------------------------

doBFD=1
doUnfold=0

R=$1
EOS_local=$2
#/eos/user/a/amohanty/
OUT=${EOS_local}/media/jackbauer/data/z_out/R_0$R
#OUT=${OUT}_25
OUT=${OUT}_finaltry

##################################### THINGS TO CHANGE
flagCUT=0
cutfileNo=0
#if [ $flagCUT -eq 1 ]; then
#OUT=${OUT}_cutsysDeDe
#fi 

if [ $flagCUT -eq 1 ]; then
 data=${dataCUT2}
 effFile=${effCUT2}
 if [ $cutfileNo -eq 1 ]; then
  OUT=${OUT}_cutsysLL1
 elif [ $cutfileNo -eq 2 ]; then 
  OUT=${OUT}_cutsysLL2
 elif [ $cutfileNo -eq 3 ]; then 
  OUT=${OUT}_cutsysLL3
 elif [ $cutfileNo -eq 4 ]; then 
  OUT=${OUT}_cutsysTT2
 elif [ $cutfileNo -eq 0 ]; then 
  OUT=${OUT}_cutsysDeDe
# elif [ $cutfileNo -eq 5 ]; then
#  OUT=${OUT}_cutsysT3
 fi
fi 

#####################################
dataInFileDir=$OUT/signalExtraction/plots
outDirFD=$OUT/FDsubtraction
simDir=$OUT/SimFiles/BFeedDown
isprefix=0
ispostfix=0
DjetEff=1
efffile=$OUT/efficiency/DjetEff_prompt_jetpt
#####################################
dataUnfoldInDir=$OUT/FDsubtraction
detRMFilePrompt=$OUT/ResponseMatrix/DetMatrix_prompt
bkgRMFile=$OUT/ResponseMatrix/DetMatrix_
unfoldingDirOut=$OUT/unfolding
regPar=5
isPrior=0
priorType=1 #1 to 8
isBkgRM=0
isFDUpSys=0
isFDDownSys=0

source zrun_settings.csh

if [ $flagCUT -eq 1 ]; then
 data=${dataCUT2}
 effFile=$effCUT2
fi

## B FEED DOWN=================================
if [ $doBFD -eq 1 ]; then
 cd ../FDsubtraction
  #bash zrun_BFD6D.csh $dataInFileDir $outDirFD $simDir $data $effFile $isprefix $ispostfix $DjetEff $efffile
  root -l  subtractFD_zjet.C'( "'$dataInFileDir'","'$outDirFD'","'$simDir'","'$data'","'$effFile'",'$isprefix','$ispostfix','$DjetEff',"'$efffile'")'
 cd ../DsignalExtraction
fi

## B FEED DOWN ends=================================
## UNFOLDING=================================
if [ $doUnfold -eq 1 ]; then
 cd ../unfolding
  bash zrun_unfold4D.csh $dataUnfoldInDir $detRMFilePrompt $bkgRMFile $unfoldingDirOut $regPar $isPrior $priorType $isBkgRM $isFDUpSys $isFDDownSys $effFile
 cd ../DsignalExtraction
 
 if [ $isPrior -eq 1 ]; then
 for thing in 1 2 3 4 5 6 7 8; do #these are priors
   cd ../unfolding
     bash zrun_unfold4D.csh $dataUnfoldInDir $detRMFilePrompt $bkgRMFile $unfoldingDirOut $regPar 1 $thing $isBkgRM $isFDUpSys $isFDDownSys $effFile
   cd ../DsignalExtraction
 done
 fi
fi
## UNFOLDING ends=================================
