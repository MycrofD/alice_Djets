#!/bin/bash

# B. Trzeciak -- script to run D-jets analysis chain for the efficiency-corrected signal extraction

currDir=`pwd`

dataFileDir=/home/basia/Work/alice/analysis/pPb_run2/outData/
jetSignalDir=/home/basia/Work/alice/analysis/pPb_run2/DstarSignal

effFileDir=/home/basia/Work/alice/analysis/pPb_run2/outMC
effDir=/home/basia/Work/alice/analysis/pPb_run2/Efficiencies

prod=$1

##### data output settings 
dataFile=$2
outDataResult=$3
jetSignalOut=$jetSignalDir/$outDataResult

mkdir -p $jetSignalOut/plots

#### efficiency output seetings
effOutPrompt=$effFileDir/$4
effFileDirOut=$effDir/$5

mkdir -p $effFileDirOut

##### postfix for tree (for a cut variation mostly)
if [ $7 -gt 0 ]
then
	postfix=$8
fi

################################################
############### efficiency
################################################

if [ $9 -eq 0 ]

cd $effDir
jetptmin=2
jetptmax=50

#### get prompt efficiency
if [ $7 -eq 0 ]
then
	./getEfficiency.csh 1 $effOutPrompt $effFileDirOut 0
else
	./getEfficiency.csh 1 $effOutPrompt $effFileDirOut 1 $postfix
fi

effPrompt=$effFileDirOut/DjetEff_prompt_jetpt${jetptmin}_${jetptmax}.root

fi

cd $currDir

################################################
############### D-jet signal
################################################


# already with D bins defined
cd $jetSignalDir
## get efficiency corrected yields, the 6th argument is if eff corrected
if [ $7 -eq 0 ]
then
	root -l -b -q plotsJetPtSpectrum_sideBand.C'("'$effPrompt'", "'$prod'", "'$dataFileDir'", "'$dataFile'", "'$outDataResult'",'$6','$9','$7')' > outfile 2>&1
else
	root -l -b -q plotsJetPtSpectrum_sideBand.C'("'$effPrompt'", "'$prod'", "'$dataFileDir'", "'$dataFile'", "'$outDataResult'",'$6','$9','$7',"'$postfix'")'  > outfile 2>&1
fi

dataEffCorrOutFile=$jetSignalOut/JetPtSpectra_SB_${prod}_eff_ptD3.root

cd $currDir




