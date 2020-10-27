#!/bin/bash

#--------------------------------------------------
#  Author A.Mohanty
#  Utrecht University
#  auro.mohanty@cern.ch
#--------------------------------------------------
#-- script to run D-jets z_||,ch analysis

defaultrun=1
multitrial=0
sbsigsys=0

zrun_file=zrun4.csh
unoriginal=unoriginal_bash
cat $unoriginal > zrun.csh
cat $zrun_file >> zrun.csh

#EOS_local=/eos/user/a/amohanty/
EOS_local=

## FLAGS for all analysis steps
oneD=1    #SUPREME FLAG FOR 1D
#--------------------------
flagEff=1 #1 Getting efficiencies and MC sigmas in different jetpt intervals for all Dpt bins
flagRef=1 #2 Reflections for different jetpt intervals for all Dpt bins
flagSBs=1 #3 Side Band subtraction method
flagSim=1 #4 Simulation for non-prompt and prompt D-jets

# 2D settings
doBFD=0
doUnfold=0

OUT=${EOS_local}/media/jackbauer/data/z_out/R_0
OUTsuffix=_finaltry

# Running for all required bins of jet pT
#----------------------------------------- #{5.0, 7.0, 10.0, 15.0, 36.0, 5.0, 15.0, 30.0, 15.0, 50.0, 10.0, 16.0, 36.0, 3.0, 5.0};
if [ $defaultrun -eq 1 ]; then
#for R in 2; do
for R in 4 6; do
    if [ $oneD -eq 1 ]; then
	    #for thing in 16 1 2 3 9 6 18 19; do
	    #for thing in 1 2 3 9 18 19; do
	    #for thing in 26 1 2 3 9 24; do
	    #for thing in 1 2 3 9 16; do
	    #for thing in 1 2 3 9; do
	    for thing in 1 2 3 9 24; do
	    #for thing in 1 ; do
	    	bash zrun.csh $thing $R ${OUT}$R${OUTsuffix} $flagEff $flagRef $flagSBs $flagSim
	    done
    fi
	bash zrun4D.csh $R ${OUT}$R${OUTsuffix} $doBFD $doUnfold
done
#fi
## Multi trial systematics
##------------------------
elif [ $multitrial -eq 1 ] || [ $sbsigsys -eq 1 ]; then
boundSigma=$1
fsigmafactor=$2
fixedMass=$3
bkgType=$4
minfSys=$5
maxfSys=$6
fMassBinWidthFactor=$7
sysNum=$8
    
    for R in 2; do
        for thing in 1 2 3 9; do
        
            if [ $sbsigsys -eq 0 ]; then
            	bash zrun.csh $thing $R ${OUT}$R${OUTsuffix} 0 0 0 0 $boundSigma $fsigmafactor $fixedMass $bkgType $minfSys $maxfSys $fMassBinWidthFactor $sysNum 
            elif [ $sbsigsys -eq 1 ]; then
            sigmaSignal=$9
            sigmaBkgll=${10}
            sigmaBkglh=${11}
            sigmaBkgrl=${12}
            sigmaBkgrh=${13}
            	bash zrun.csh $thing $R ${OUT}$R${OUTsuffix} 0 0 0 0 $boundSigma $fsigmafactor $fixedMass $bkgType $minfSys $maxfSys $fMassBinWidthFactor $sysNum $sigmaSignal $sigmaBkgll $sigmaBkglh $sigmaBkgrl $sigmaBkgrh
            fi
        done
    done
#fi #elif [ $multitrial -eq 1 ] || [ $sbsigsys -eq 1 ]; then
#fi #defaultrun if else end
#----------------------------
fi
