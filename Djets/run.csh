#!/bin/bash

# B. Trzeciak -- script to run D-jets analysis chain

############### set up the output directory name (the whole analysis chain output will be saved there)
outdirectory=$HOME/Work/alice/analysis/testOutDzero/
outdirectory=$HOME/Work/alice/analysis/testOutDzeroRef/
mkdir -p $outdirectory

Dmeson=0 #0: D0, 1: D*
isRefl=1 #reflections for D0

############### configure pT bins, etc., in the corresponding Dstar/Dzero config files
if [ $Dmeson -eq 0 ]; then
  cat configDzero.h > config.h
else
  cat configDstar.h > config.h
fi

echo >> config.h
echo "TString OUTDIRECTORY=\"$outdirectory\";"  >> config.h

############### analysis flags
isSignal=1
isEffPrompt=0
isEffNonPrompt=0

################################################
############### D-jet signal config
################################################
outfiledir=$HOME/Work/alice/analysis/pPb_run2/outData/cuts1605
outfiledir=$HOME/Work/alice/analysis/pPb_run2/D0jet/outData
lprod=woSDD
lprod=LHC16woSDDFAST
isMoreFiles=0
prod=kl #if more than one file to analyse

################################################
############### in/out directories config
################################################
currDir=`pwd`
signalDir=DsignalExtraction
effDir=efficiency

signalDirOut=$outdirectory/signalExtraction  #don't change this, should be hardcoded in the macro
effDirOut=$outdirectory/efficiency          #don't change this, should be hardcoded in the macro

if [ $isEffPrompt -eq 0 ]; then
  signalFile=$signalDirOut/JetPtSpectra_SB_noEff.root
else
  signalFile=$signalDirOut/JetPtSpectra_SB_eff.root
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

#pattern: $analysisFile_$lhcprod$production.root
analysisfile=$outfiledir/AnalysisResults
lhcprod=$lprod
lhcprod=$lprod
production=$prod #if more than one file to analyse
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
./getYields.csh $analysisfile $lhcprod $isMoreFiles $production $isEffPrompt $effFilePrompt $isRefl $reflfile $ispostfix $postfix

cd $currDir

fi

################################################
###############
################################################
