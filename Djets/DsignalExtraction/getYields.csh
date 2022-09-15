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

ctry=1
# 2, (3,8)
sigmaWindow=2;SBout=8;SBint=3
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]
# 2, (4,10)
sigmaWindow=2;SBout=10;SBint=4
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]
# 2, (4,12)
sigmaWindow=2;SBout=12;SBint=4
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]
# 2, (4,15)
sigmaWindow=2;SBout=15;SBint=4
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]
# 3, (4,9)
sigmaWindow=3;SBout=9;SBint=4
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]
# 3, (4,12)
sigmaWindow=3;SBout=12;SBint=4
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]

fi


#ctry=0
#  for sigmaWindow in 2 3; do
#    for SBout in 9 8; do
#      for SBint in 4 3.5 4.5; do
#    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
#ctry=$[ctry+1]
#      done
#    done
#  done
#fi

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

# (1)*2*2*2*2*2=32
  for boundSigma in 0; do
    for bkgType in 0 1; do
      for fixedMass in 0 1; do
        for minfSys in 1.72 1.74; do
          for maxfSys in 2.00 2.03; do
            for fMassBinWidthFactor in 2 4; do
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysMT')'
    ctry=$[ctry+1]
            done #massbinwidth
          done #maxfsys
        done #minfsys
      done #fixedMass
    done #bkgtype
  done #boundsigma


# (1*1*)2*2*2*2*2=32
  for boundSigma in 1; do
  for fsigmafactor in 1; do
    for bkgType in 0 1; do
      for fixedMass in 0 1; do
        for minfSys in 1.72 1.74; do
          for maxfSys in 2.00 2.03; do
            for fMassBinWidthFactor in 2 4; do
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysMT')'
    ctry=$[ctry+1]
            done #massbinwidth
          done #maxfsys
        done #minfsys
      done #fixedMass
    done #bkgtype
  done #fsigmafactor
  done #boundsigma

# (1*)*2*2*1*2*2*2=32
  for boundSigma in 1; do
  for fsigmafactor in 1.1 0.9; do
    for bkgType in 0 1; do
      for fixedMass in 0; do
        for minfSys in 1.72 1.74; do
          for maxfSys in 2.00 2.03; do
            for fMassBinWidthFactor in 2 4; do
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
## for compatibility with 7 TeV and 13 TeV:
#• free width σ, free mean mD0 ,
#• fixed width σ = σMC, free mean mD0 ,
#• fixed width σ = 0.85 · σMC, free mean mD0 ,
#• fixed width σ = 1.15 · σMC, free mean mD0 ,
#• free width σ, fixed mean mD0 = mPDG,
#• fixed width σ = σMC, fixed mean mD0 = mPDG.

#and following variations were considered on the fitting procedure
#• background functions: exponential, linear,
#• lower limit of fit range: 1.72, 1.74 GeV/c2,
#• upper limit of fit range: 2.00, 2.03 GeV/c2,
#• mass bin width: 5 MeV/c2, 10 MeV/c2

# So, which numbers do they belong to in our trials 
# free width σ, free mean mD0 - 
#     boundSigma=0, fixedMass=0
# fixed width σ = σMC, free mean mD0
#     boundSigma=1, fsigmafactor=1, fixedMass=0
# fixed width σ = 0.85 (0.9) · σMC, free mean mD0
#     boundSigma=1, fsigmafactor=0.9, fixedMass=0
# fixed width σ = 1.15 (1.1) · σMC, free mean mD0 ,
#     boundSigma=1, fsigmafactor=1.1, fixedMass=0
# free width σ, fixed mean mD0 = mPDG,
#     boundSigma=0, fixedMass=1
# fixed width σ = σMC, fixed mean mD0 = mPDG
#     boundSigma=1, fsigmafactor=1, fixedMass=1
#
# what is missing from sigma mass combination?
# 4 sigma, 2 mass
# free s, free m- OK
# free s, fixed m-OK
# MC s, free m-OK
# MC s, fixed m-OK
# 1.1 s, free m-OK
# 1.1 s, fixed m
# 0.9 s, free m-OK
# 0.9 s, fixed m
# So, 162 free sigma OK.
# Next 162: MC sigma OK.
# From third 162, 162/3, /2
# From fourth 162, 162/3, /2
#

# background functions: exponential, linear,
#     bkgType [0, 1]
# lower limit of fit range: 1.72, 1.74 GeV/c2,
#     minfSys in 
# upper limit of fit range: 2.00, 2.03 GeV/c2,
# mass bin width: 5 MeV/c2, 10 MeV/c2
#     fMassBinWidthFactor in [2, 1]





# background: polynomials
# (648/4) = 162 for each sigma
# bkgtype: polynomial. divide by 3, remove third component.
# index number is 1 less than plotted trial number 
# 162/3 = 54
# index number in [108, 161], [270, 323], [432, 485], [594, 647]
#
# mass bin width: remove 2.5 MeV, use only 5, 10 MeV. 
# i.e. remove 1 from [2,4,1] from the for loop
# e.g. index number, remove 2 from 0,1,2. Remove 5 from 3, 4, 5 etc.
# index number in [3x-1] for x in range(1,648/3=216)
#
