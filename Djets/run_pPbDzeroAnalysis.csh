#!/bin/bash

currDir=`pwd`

ptbinning=$1
jetpttruemin=$2
jetpttruemax=$3
jetptmeasmin=$4
jetptmeasmax=$5

outputdirectorybase=$HOME/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_
outputdirectory=$HOME/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts
outputdirectorySignal=Default
if [ $ptbinning -eq 0 ]; then
outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_ppbinning
fi
if [ $ptbinning -eq 1 ]; then
outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_finebinning
fi
if [ $ptbinning -eq 2 ]; then
outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_PbPbbinning
fi
if [ $ptbinning -eq 3 ]; then
outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_binning2
fi


lhcprod=_LHC16R03
efficiencyfile=AnalysisResults_fast_R03_D0MCPythia_default.root
detRMpromptfile=AnalysisResults_fast_R03_D0MCPythia_default.root
detRMnonpromptfile=AnalysisResults_fast_R03_D0MCPythia_default.root

ispostfix=0
postfix=Cut
ispostfixFD=1
postfixFD=FD

bkgRMtype=$6
unfType=$7 #bayes
regPar=$8
isPrior=$9
priorType=${10}

# systematics
doBkg=${11}
doPrior=${12}
doBkgPrior=${13}
doCutVar=${14}
doJESSys=${15}
doFDSys=${16}
doSystematics=${17}
doRawCutVar=${18}
doRawSpectra=${19}
doSignal=${20}

unfCutSys=0

################################################
############### Default
################################################

if [ $doRawSpectra -eq 1 ]; then
 ./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 0 0 0 1
 exit 1
fi


if [ $doSignal -eq 1 ]; then
./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 0 0 0 0
fi

if [ $doFDSys -eq 1 ]; then
	./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 1 1 0 0

	./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 1 0 1 0
fi


################################################
############### Background variations
################################################

if [ $doBkg -gt 0 ]; then

for bkg in `seq 1 1 11`; do


		bkgRMtype=$bkg
	./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 0 0 0 0

done

fi
bkgRMtype=0

################################################
############### Prior variations
################################################

if [ $doPrior -gt 0 ]; then

for prior in `seq 0 1 8`; do

	isPrior=1
	priorType=$prior
	./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 0 0 0 0

done

fi
isPrior=0
priorType=0

################################################
############### Background and prior variations
################################################

if [ $doBkgPrior -gt 0 ]; then

for bkg in `seq 0 1 11`; do

	bkgRMtype=0
	unfType=0 #bayes
	regPar=3
	isPrior=0
	priorType=0

		bkgRMtype=$bkg
	./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 0 0 0 0

	for prior in `seq 0 1 8`; do
		isPrior=1
		priorType=$prior
		./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 0 0 0 0
	done

done

fi

################################################
############### JES systematics
################################################

if [ $doJESSys -eq 1 ]; then

	./run_JES.csh $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectory

fi

################################################
############### Cut Variation
################################################

if [ $doCutVar -eq 1 ]; then

	./run_cut.csh $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectorybase $doRawSpectra

fi

################################################
############### Get systematic uncertanties
################################################


if [ $doSystematics -eq 1 ]; then
	cd systematics
		root -l -b -q JetSpectrumSys.C'('$regPar',"'$outputdirectory'","'$outputdirectorySignal'",1)'
		root -l -b -q cutsSystematics.C'('$regPar',"'$outputdirectorybase'","'$outputdirectorySignal'",'$doRawCutVar',1,1)'
		root -l -b -q cutsSystematics.C'('$regPar',"'$outputdirectorybase'","'$outputdirectorySignal'",'$doRawCutVar',0,1)'
	cd $currDir
fi

exit 1
