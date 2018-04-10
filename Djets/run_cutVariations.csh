#!/bin/bash

ptbinning=$3
jetpttruemin=$4
jetpttruemax=$5
jetptmeasmin=$6
jetptmeasmax=$7

doCuts=$2
outputdirectorySignal=${13}
out=${14}

doRawSpectra=${15}

base=CutVarBase
outputdirectory=$out$base
if [ $doCuts -gt 0 ]; then
	cut=Cut$1
	outputdirectory=$out$cut
fi

lhcprod=_LHC16R03_CutVariation
efficiencyfile=AnalysisResults_fast_R03_D0MCPythia_CutVariation_1.root
detRMpromptfile=AnalysisResults_fast_R03_D0MCPythia_CutVariation_1.root
detRMnonpromptfile=AnalysisResults_fast_R03_D0MCPythia_CutVariation_1.root


ispostfix=0
postfix=Cut
if [ $doCuts -gt 0 ]; then
	ispostfix=1
	postfix=$1
fi


ispostfixFD=1
postfixFD=FD
if [ $doCuts -gt 0 ]; then
	ispostfixFD=1
	postfixFD=$1FD
fi


bkgRMtype=$8
unfType=$9 #bayes
regPar=${10}
isPrior=${11}
priorType=${12}



################################################
############### Cut variation
################################################
if [ $doRawSpectra -eq 0 ]; then
  ./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 0 1 0 0 0 0
fi

if [ $doRawSpectra -eq 1 ]; then
  ./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 0 1 0 0 0 1
fi

exit 1
