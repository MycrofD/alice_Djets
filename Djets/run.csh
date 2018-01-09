#!/bin/bash

# B. Trzeciak -- script to run D-jets analysis chain

outdirectory=$HOME/Work/alice/analysis/testOutDzeroRef/
mkdir -p $outdirectory

Dmeson=0 #0: D0, 1: D*
isRefl=1 #reflections for D0

if [ $Dmeson -eq 0 ]; then
  cat configDzero.h > config.h
else
  cat configDstar.h > config.h
fi

echo >> config.h
echo "TString OUTDIRECTORY=\"$outdirectory\";"  >> config.h

isSignal=1
isEffPrompt=0
isEffNonPrompt=0

currDir=`pwd`
signalDir=DsignalExtraction
effDir=efficiency

sigalDirOut=$outdirectory/signalExtraction  #don't change this, should be hardcoded in the macro
effDirOut=$outdirectory/efficiency          #don't change this, should be hardcoded in the macro

if [ $isEffPrompt -eq 0 ]; then
  signalFile=$sigalDirOut/JetPtSpectra_SB_noEff.root
else
  signalFile=$sigalDirOut/JetPtSpectra_SB_eff.root
fi

effFilePrompt=$effDirOut/DjetEff_prompt.root
effFileNonPrompt=$effDirOut/DjetEff_nonPrompt.root


################################################
############### D-jet efficiency
################################################
if [ $isEffPrompt -gt 0 ]; then
  cd $effDir
  ./getEfficiency.csh 1 $effOutPrompt $effFileDirOut 1 $postfix
  cd $currDir
fi

if [ $isEffNonPrompt -gt 0 ]; then
  cd $effDir
  cd $currDir
fi

################################################
############### D-jet signal
################################################
if [ $isSignal -gt 0 ]; then
cd $signalDir

#pattern: $analysisFile_$lhcprod$prod.root
analysisfile=/home/basia/Work/alice/analysis/pPb_run2/outData/cuts1605/AnalysisResults
analysisfile=/home/basia/Work/alice/analysis/pPb_run2/D0jet/outData/AnalysisResults
lhcprod=woSDD
lhcprod=LHC16woSDDFAST
isMoreFiles=0
prod=kl
reflfile=/home/basia/Work/alice/analysis/pPb_run2/D0jet/outMC/reflections/reflections_fitted_DoubleGaus.root
ispostfix=0
postfix=Cut

# 1: analysis output file;
# 2: lhcprod;
# 3: if more than one file
# 4: prod (if more than one file);
# 5: if eff corrected;
# 6: efficiency file
# 7: if to fit reflection for D0
# 8: name of file with reflection templates
# 9: if output list with postfix;
# 10: postfix
./getYields.csh $analysisfile $lhcprod $isMoreFiles $prod $isEffPrompt $effFilePrompt $isRefl $reflfile $ispostfix $postfix

cd $currDir

fi

################################################
###############
################################################
