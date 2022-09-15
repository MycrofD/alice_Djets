#!/bin/bash

# B. Trzeciak (Utrecht University)
# A. Mohanty (Utrecht University)
#-- script to run D-jets yield extraction and raw yield systematics

source ../run_settings.csh

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
isprefix=${11}
################################################
############### D-jet signal
################################################

root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix')' #> outfile 2>&1
#################################################
################ Cut systematics fixing sigma to data
#################################################
#if[ $flagCUT -eq 1 ]; then
#  root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix',5,1)' #> outfile 2>&1
#fi

################################################
############### Raw Yield systematics
################################################
dorawsysSB=$flagSBSig #for signal SB range variation

dorawsysMT=$flagMulti #for multi-trial
################################################
############### :->Signal SB ranges
################################################
boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=0;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=2 #default

# Side-Band Signal ranges
if [ $dorawsysSB -eq 1 ]; then
boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=0;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=2                                                                                                                   
ctry=0
  for sigmaWindow in 2 3; do
    for SBout in 9 8; do
      for SBint in 4 3.5 4.5; do
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]
      done
    done
  done
fi

################################################
############### :->Multi-trial
################################################
#forlist1:(4) boundSigma=0 (free sigma); boundSigma=1 (sigma fixed to MC),fsigmafactor=1,1.1,0.9
#forlist2:(2) mass=0,1
#forlist3:(3) bkg=0,1,2
#forlist4:(3) lowlim=1.71,1.70,1.72
#forlist5:(3) upplim=2.10,2.09,2.11
#forlist6:(3) massrebin=2,4,0.5
#########################################################
if [ $dorawsysMT -eq 1 ]; then
ctry=1
fsigmafactor=1

  for boundSigma in 0; do
    for bkgType in 0 1 2; do
      for fixedMass in 0 1; do
        for minfSys in 1.71 1.72 1.70; do
          for maxfSys in 2.1 2.09 2.11; do
            for fMassBinWidthFactor in 2 4 1; do
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysMT')'
    ctry=$[ctry+1]
            done #massbinwidth
          done #maxfsys
        done #minfsys
      done #fixedMass
    done #bkgtype
  done #boundsigma

  for boundSigma in 1; do
  for fsigmafactor in 1 1.1 0.9; do
    for bkgType in 0 1 2; do
      for fixedMass in 0 1; do
        for minfSys in 1.71 1.72 1.70; do
          for maxfSys in 2.1 2.09 2.11; do
            for fMassBinWidthFactor in 2 4 1; do
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysMT')'
    ctry=$[ctry+1]
            done #massbinwidth
          done #maxfsys
        done #minfsys
      done #fixedMass
    done #bkgtype
  done #fsigmafactor
  done #boundsigma

fi
########################################################
