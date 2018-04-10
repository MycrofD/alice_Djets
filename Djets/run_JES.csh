#!/bin/bash

ptbinning=$1
jetpttruemin=$2
jetpttruemax=$3
jetptmeasmin=$4
jetptmeasmax=$5

bkgRMtype=$6
unfType=$7 #bayes
regPar=${8}
isPrior=${9}
priorType=${10}

outputdirectorySignal=${11}\_JES

#outputdirectory=DzeroR03_RefDPt3PythiaEff_BaseCuts
outputdirectory=${12}

lhcprod=_LHC16R03
efficiencyfile=AnalysisResults_fast_R03_D0MCPythia_default.root
detRMpromptfile=AnalysisResults_fast_R03_D0MCPythia_default.root
detRMnonpromptfile=AnalysisResults_fast_R03_D0MCPythia_default.root

ispostfix=0
postfix=Cut
ispostfixFD=1
postfixFD=FD


################################################
############### Default
################################################

./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 0 0 0 0 0 0

################################################
############### JES systematics
################################################


	detRMpromptfile=AnalysisResults_fast_R03_D0MCPythia_JES96_1.root
	detRMnonpromptfile=AnalysisResults_fast_R03_D0MCPythia_JES96_1.root
	eff=96
	outputdirectorySignal1=$outputdirectorySignal$eff
	./run.csh $outputdirectory $outputdirectorySignal1 $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 0 0 0 0 0 0

	detRMpromptfile=AnalysisResults_fast_R03_D0MCPythia_JES95_1.root
	detRMnonpromptfile=AnalysisResults_fast_R03_D0MCPythia_JES95_1.root
	eff=95
	outputdirectorySignal1=$outputdirectorySignal$eff
	./run.csh $outputdirectory $outputdirectorySignal1 $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 0 0 0 0 0 0

	detRMpromptfile=AnalysisResults_fast_R03_D0MCPythia_JES90_1.root
	detRMnonpromptfile=AnalysisResults_fast_R03_D0MCPythia_JES90_1.root
	eff=90
	outputdirectorySignal1=$outputdirectorySignal$eff
	./run.csh $outputdirectory $outputdirectorySignal1 $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 0 0 0 0 0 0

	exit 1
