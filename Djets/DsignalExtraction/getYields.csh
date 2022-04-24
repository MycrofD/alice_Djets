#!/bin/bash

# B. Trzeciak (Utrecht University)
# A. Mohanty (Utrecht University)
#-- script to run D-jets yield extraction and raw yield systematics

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

################################################
############### Raw Yield systematics
################################################
dorawsysSB=0 #for signal SB range variation

dorawsysMT=0 #for multi-trial
################################################
############### :->Signal SB ranges
################################################
boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=0;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=2 #default
if [ $dorawsysSB -eq 1 ]; then
    ctry=1;fsigmaSignal=3;fbkgll=-9;fbkglh=-4;fbkgrl=4;fbkgrh=9
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'
  
    ctry=2;fsigmaSignal=2;fbkgll=-8;fbkglh=-4;fbkgrl=4;fbkgrh=8
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'
  
    ctry=3;fsigmaSignal=3;fbkgll=-8;fbkglh=-4;fbkgrl=4;fbkgrh=8
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'
  
    ctry=4;fsigmaSignal=2;fbkgll=-7;fbkglh=-4;fbkgrl=4;fbkgrh=7
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'
  
    ctry=5;fsigmaSignal=3;fbkgll=-7;fbkglh=-4;fbkgrl=4;fbkgrh=7
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'

    ctry=6;fsigmaSignal=2;fbkgll=-9;fbkglh=-4.5;fbkgrl=4.5;fbkgrh=9
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'

    ctry=7;fsigmaSignal=3;fbkgll=-9;fbkglh=-4.5;fbkgrl=4.5;fbkgrh=9
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'
  
    ctry=8;fsigmaSignal=2;fbkgll=-8;fbkglh=-4.5;fbkgrl=4.5;fbkgrh=8
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'

    ctry=9;fsigmaSignal=3;fbkgll=-8;fbkglh=-4.5;fbkgrl=4.5;fbkgrh=8
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'
  
    ctry=10;fsigmaSignal=2;fbkgll=-7;fbkglh=-4.5;fbkgrl=4.5;fbkgrh=7
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'
  
    ctry=11;fsigmaSignal=3;fbkgll=-7;fbkglh=-4.5;fbkgrl=4.5;fbkgrh=7
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'
  
    ctry=12;fsigmaSignal=2;fbkgll=-9;fbkglh=-3.5;fbkgrl=3.5;fbkgrh=9
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'

    ctry=13;fsigmaSignal=3;fbkgll=-9;fbkglh=-3.5;fbkgrl=3.5;fbkgrh=9
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'
  
    ctry=14;fsigmaSignal=2;fbkgll=-8;fbkglh=-3.5;fbkgrl=3.5;fbkgrh=8
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'

    ctry=15;fsigmaSignal=3;fbkgll=-8;fbkglh=-3.5;fbkgrl=3.5;fbkgrh=8
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'
  
    ctry=16;fsigmaSignal=2;fbkgll=-7;fbkglh=-3.5;fbkgrl=3.5;fbkgrh=7
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'
  
    ctry=17;fsigmaSignal=3;fbkgll=-7;fbkglh=-3.5;fbkgrl=3.5;fbkgrh=7
    root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysSB','$fsigmaSignal','$fbkgll','$fbkglh','$fbkgrl','$fbkgrh')'

fi
  
################################################
############### :->Multi-trial
################################################
#trials=20
###parameters, default
#boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=0;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=2
#boundSigma=1;fsigmafactor=1;fixedMass=0;bkgType=0;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=2
#boundSigma=1;fsigmafactor=1.15;fixedMass=0;bkgType=0;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=2
#boundSigma=1;fsigmafactor=0.85;fixedMass=0;bkgType=0;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=2
#boundSigma=0;fsigmafactor=1;fixedMass=1;bkgType=0;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=2
#boundSigma=1;fsigmafactor=1;fixedMass=1;bkgType=0;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=2
#boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=1;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=2
#boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=2;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=2
#boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=0;minfSys=1.72;maxfSys=2.1;fMassBinWidthFactor=2
#boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=0;minfSys=1.70;maxfSys=2.1;fMassBinWidthFactor=2
#boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=0;minfSys=1.70;maxfSys=2.09;fMassBinWidthFactor=2
#boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=0;minfSys=1.70;maxfSys=2.11;fMassBinWidthFactor=2
#boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=0;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=4
#boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=0;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=1

## variation 0: default
## variation 1: fixed sigma=sigmaMC
## variation 2: fixed sigma=sigmaMC*1.15
## variation 3: fixed sigma=sigmaMC*0.85
## variation 4: free sigma, fixed mass
## variation 5: fixed sigma, fixed mass
## variation 6: background type 1
## variation 7: background type 2
## variation 8: lower limit fit range 1.72
## variation 9: lower limit fit range 1.70
## variation10: upper limit fit range 2.09
## variation11: upper limit fit range 2.11
## variation12: mass rebin times 2*2
## variation13: mass rebin times 2*0.5

#forlist1:(4) boundSigma=0; boundSigma=1,fsigmafactor=1,1.1,0.9 
#forlist2:(2) mass=0,1
#forlist3:(3) bkg=0,1,2
#forlist4:(3) lowlim=1.71,1.70,1.72
#forlist5:(3) upplim=2.10,2.09,2.11
#forlist6:(3) massrebin=2,4,0.5

#if [ $dorawsysMT -eq 1 ]; then
#ctry=1
#fsigmafactor=1
#
#for fixedMass in 0 1; do
#  for bkgType in 0 1 2; do
#    for minfSys in 1.71 1.72 1.70; do
#      for maxfSys in 2.1 2.09 2.11; do
#        for fMassBinWidthFactor in 2 4 1; do
#          for boundSigma in 0; do
#fsigmafactor=1
#root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysMT')'
#ctry=$[ctry+1]
#          done
#          for boundSigma in 1; do
#            #for fsigmafactor in 1 1.1 1.15 0.9 0.85; do
#            for fsigmafactor in 1 1.1 0.9; do
#root -l -b -q signalExtraction_SB.C'("'$dataFile'", '$isEff', "'$effFile'",'$isRefl', "'$refFile'", '$ispostfix', "'$postfix'", "'$dirOut'", 1, '$isMoreFile',"'$prod'",'$isprefix','$boundSigma','$fsigmafactor','$fixedMass','$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$dorawsysMT')'
#ctry=$[ctry+1]
#            done
#          done
#        done
#      done    
#    done
#  done
#done
#fi
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
