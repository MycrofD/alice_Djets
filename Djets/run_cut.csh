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

outputdirectorySignal=${11}
outputdirectorybase=${12}

doRawSpectra=${13}

./run_cutVariations.csh cut 0 $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectorybase $doRawSpectra
./run_cutVariations.csh L0 1 $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectorybase $doRawSpectra
./run_cutVariations.csh L0 1 $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectorybase $doRawSpectra
./run_cutVariations.csh L1 1 $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectorybase $doRawSpectra
./run_cutVariations.csh L2 1 $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectorybase $doRawSpectra
./run_cutVariations.csh L3 1 $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectorybase $doRawSpectra
./run_cutVariations.csh T0 1 $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectorybase $doRawSpectra
./run_cutVariations.csh T1 1 $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectorybase $doRawSpectra
./run_cutVariations.csh T2 1 $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectorybase $doRawSpectra
