#!/bin/bash
#
#//-----------------------------------------------------------------------
#//  Author B.Trzeciak
#//  Utrecht University
#//  barbara.antonina.trzeciak@cern.ch
#//-----------------------------------------------------------------------
#-- script to run D-jets analysis
# examples from D0-jet p-Pb analysis

# ============= set up config =================

Dmeson=0 # 0: D0, 1: D*
if [ $Dmeson -eq 0 ]; then
  conffile=configDzero_pp.h
else
  conffile=configDstar.h
fi

############## RooUnfold path
roounfoldpwd=$HOME/ALICE_HeavyFlavour/RooUnfold-1.1.1/libRooUnfold

currDir=`pwd`
#outputdirectorybase=$HOME/Work/alice/analysis/pp5TeV/D0jet/results_APW/DzeroR03_prelim
outputdirectorybase=$HOME/Work/alice/analysis/pp5TeV/D0jet/results_APW/DzeroR03_
outputdirectory=${outputdirectorybase}pPbCuts
#outputdirectory=${outputdirectorybase}cuts2
outputdirectorySignal=Default

# ========== input directories
anaoutfiledir=$HOME/Work/alice/analysis/pp5TeV/D0jet/outData
MCoutfiledir=$HOME/Work/alice/analysis/pp5TeV/D0jet/outMC
refoutfiledir=$HOME/Work/alice/analysis/pp5TeV/D0jet/outMC/reflections
#simFilesDir=$HOME/Work/alice/analysis/pp5TeV/D0jet/outMC/POWHEGSimulations/fastSim_pp5TeV # POWHEG sim
simFilesDir=/home/jackbauer/ALICE_HeavyFlavour/work/Djets/out/outMC/all
bkgRMDir=$currDir/ResponseMatrix/BkgRM03

# ========== file names
analysisdatafile=AnalysisResults_437.root
#analysisdatafile=AnalysisResults_LHC17pq_FASTwoSDD.root
#analysisdatafile=AnalysisResults_503_R04.root
#analysisdatafile=AnalysisResults_504_R06.root
isMoreFiles=0                                     # 1 if there are more input files to be read
prod=kl                                           # if there are more input files to be read
reflfile=reflectionTemplates_pp.root
effMCfile=AnalysisResults_454.root
#effMCfile=AnalysisResults_fast_R03_D0MC_def.root
#effMCfile=AnalysisResults_634_pp5TeV_z.root
#effMCfile=AnalysisResults_642_pp5TeV_z.root
#effMCfile=AnalysisResults_683_ppMC_R06.root
MCfile=$effMCfile #AnalysisResults_fast_R03_D0MC_def.root
#MCfile=AnalysisResults_fast_R03_D0MCPythia_JES96_1.root
efficiencyfile=$effMCfile
detRMpromptfile=$MCfile
detRMnonpromptfile=$MCfile
isprefix=0                                       #  if prefix is one, postfix must be one.
ispostfix=0                                       # if container in the analysis output file has different name than default you set here if and what is the postfix, this is set up in the signal, efficiency and RM extraction macros
postfix=HP4 #cut2
ispostfixFD=0                                     # if container in the analysis output file has different name than default you set here if and what is the postfix, for the FD part wagons are usually configured with additional "FD" string in the container name, you should adjust this to yours configuration
#postfixFD=FD
postfixFD=HP4 #FDcut2

isRefl=1
isBkgRM=0

######## !!! POWHEG simulations config
nSimFilesB=10                                    # have to correspond to number of files defined in the config file
nSimFilesC=9                                    # have to correspond to number of files defined in the config file

unfType=$1
regPar=$2
isPrior=$3
priorType=$4
bkgRMtype=$5
doRawSpectra=$6

############### bkg fluctuations matrix
bkgRMfileName=tmp.root

############### RUN ###############
./run.csh $outputdirectory $outputdirectorySignal $conffile $Dmeson $anaoutfiledir $analysisdatafile $refoutfiledir $reflfile $isRefl $isMoreFiles $prod $ispostfix $postfix $ispostfixFD $postfixFD $MCoutfiledir $efficiencyfile $detRMpromptfile $detRMnonpromptfile $isBkgRM $bkgRMtype $bkgRMDir $bkgRMfileName $simFilesDir $unfType $regPar $isPrior $priorType 0 0 $doRawSpectra $nSimFilesB $nSimFilesC $roounfoldpwd $isprefix


################################################
############### Additional
################################################

#### systematics  -- these are additional settings - in order to use them you need to configure this and (if needed) further scripts
doCutVar=${7}
doRawCutVar=${8}
doJESSys=${9}
doFDSys=${10}
doSystematics=${11}

################################################
############### B FD systematics -- this will work w/o additinal modifications
################################################
if [ $doFDSys -eq 1 ]; then

  ./run.csh $outputdirectory $outputdirectorySignal $conffile $Dmeson $anaoutfiledir $analysisdatafile $refoutfiledir $reflfile $isRefl $isMoreFiles $prod $ispostfix $postfix $ispostfixFD $postfixFD $MCoutfiledir $efficiencyfile $detRMpromptfile $detRMnonpromptfile $isBkgRM $bkgRMtype $bkgRMDir $bkgRMfileName $simFilesDir $unfType $regPar $isPrior $priorType 1 0 0 $nSimFilesB $nSimFilesC $roounfoldpwd $isprefix

  ./run.csh $outputdirectory $outputdirectorySignal $conffile $Dmeson $anaoutfiledir $analysisdatafile $refoutfiledir $reflfile $isRefl $isMoreFiles $prod $ispostfix $postfix $ispostfixFD $postfixFD $MCoutfiledir $efficiencyfile $detRMpromptfile $detRMnonpromptfile $isBkgRM $bkgRMtype $bkgRMDir $bkgRMfileName $simFilesDir $unfType $regPar $isPrior $priorType 0 1 0 $nSimFilesB $nSimFilesC $roounfoldpwd $isprefix

fi

################################################
############### JES systematics -- configure below names for the responce matrices
################################################

if [ $doJESSys -eq 1 ]; then

    detRMpromptfile1=AnalysisResults_fast_R03_D0MCPythia_JES96_1.root
 # 	detRMnonpromptfile1=AnalysisResults_fast_R03_D0MCPythia_JES96_1.root
detRMnonpromptfile1=AnalysisResults_fast_R03_D0MC_def.root
  	eff=95
  	outputdirectorySignal1=$outputdirectorySignal$eff

    ./run.csh $outputdirectory $outputdirectorySignal1 $conffile $Dmeson $anaoutfiledir $analysisdatafile $refoutfiledir $reflfile $isRefl $isMoreFiles $prod $ispostfix $postfix $ispostfixFD $postfixFD $MCoutfiledir $efficiencyfile $detRMpromptfile1 $detRMnonpromptfile1 $isBkgRM $bkgRMtype $bkgRMDir $bkgRMfileName $simFilesDir $unfType $regPar $isPrior $priorType 0 0 0 $nSimFilesB $nSimFilesC $roounfoldpwd $isprefix

#    detRMpromptfile1=AnalysisResults_fast_R03_D0MCPythia_JES95_1.root
#    detRMnonpromptfile1=AnalysisResults_fast_R03_D0MCPythia_JES95_1.root
#    eff=95
#    outputdirectorySignal1=$outputdirectorySignal$eff
#
#    ./run.csh $outputdirectory $outputdirectorySignal1 $conffile $Dmeson $anaoutfiledir $analysisdatafile $refoutfiledir $reflfile $isRefl $isMoreFiles $prod $ispostfix $postfix $ispostfixFD $postfixFD $MCoutfiledir $efficiencyfile $detRMpromptfile1 $detRMnonpromptfile1 $isBkgRM $bkgRMtype $bkgRMDir $bkgRMfileName $simFilesDir $unfType $regPar $isPrior $priorType 0 0 0 $nSimFilesB $nSimFilesC $roounfoldpwd $isprefix
#
#    detRMpromptfile1=AnalysisResults_fast_R03_D0MCPythia_JES90_1.root
#    detRMnonpromptfile1=AnalysisResults_fast_R03_D0MCPythia_JES90_1.root
#    eff=90
#    outputdirectorySignal1=$outputdirectorySignal$eff
#
#    ./run.csh $outputdirectory $outputdirectorySignal1 $conffile $Dmeson $anaoutfiledir $analysisdatafile $refoutfiledir $reflfile $isRefl $isMoreFiles $prod $ispostfix $postfix $ispostfixFD $postfixFD $MCoutfiledir $efficiencyfile $detRMpromptfile1 $detRMnonpromptfile1 $isBkgRM $bkgRMtype $bkgRMDir $bkgRMfileName $simFilesDir $unfType $regPar $isPrior $priorType 0 0 0 $nSimFilesB $nSimFilesC $roounfoldpwd $isprefix

fi

################################################
############### Cut Variation -- configure this !!!
################################################

if [ $doCutVar -eq 1 ]; then

  # ========== file names
  analysisdatafileCut=AnalysisResults_LHC16R03_CutVariation.root
  efficiencyfileCut=AnalysisResults_fast_R03_D0MCPythia_CutVariation_1.root
  detRMpromptfileCut=AnalysisResults_fast_R03_D0MCPythia_CutVariation_1.root
  detRMnonpromptfileCut=AnalysisResults_fast_R03_D0MCPythia_CutVariation_1.root
  isReflCut=0
  outputdirectorySignalCut=Default

  ispostfixCut=0
  ispostfixFDCut=1
  postfixCut=Cut
  postfixFDCut=FD
  outputdirectoryCut=${outputdirectorybase}Cuts_VarDef
	./run.csh $outputdirectoryCut $outputdirectorySignalCut $conffile $Dmeson $anaoutfiledir $analysisdatafileCut $refoutfiledir $reflfile $isReflCut $isMoreFiles $prod $ispostfixCut $postfixCut $ispostfixFDCut $postfixFDCut $MCoutfiledir $efficiencyfileCut $detRMpromptfileCut $detRMnonpromptfileCut $isBkgRM $bkgRMtype $bkgRMDir $bkgRMfileName $simFilesDir $unfType $regPar $isPrior $priorType 0 0 $doRawSpectra $nSimFilesB $nSimFilesC $roounfoldpwd $isprefix

arraySize
declare -a cutarray=( L0 L1 L2 L3 T0 T1 T2 ) #declare postfixes

for i in "${cutarray[@]}"
do
  ispostfixCut=1
  ispostfixFDCut=1
  postfixCut=${cutarray[$i]}
  postfixFDCut=${cutarray[$i]}FD
  outputdirectoryCut=${outputdirectorybase}Cuts_${cutarray[$i]}
  ./run.csh $outputdirectoryCut $outputdirectorySignalCut $conffile $Dmeson $anaoutfiledir $analysisdatafileCut $refoutfiledir $reflfile $isReflCut $isMoreFiles $prod $ispostfixCut $postfixCut $ispostfixFDCut $postfixFDCut $MCoutfiledir $efficiencyfileCut $detRMpromptfileCut $detRMnonpromptfileCut $isBkgRM $bkgRMtype $bkgRMDir $bkgRMfileName $simFilesDir $unfType $regPar $isPrior $priorType 0 0 $doRawSpectra $nSimFilesB $nSimFilesC $roounfoldpwd $isprefix
done

  # cut variation systematics
#  cd systematics
  #after unfolding
#  root -l -b -q cutsSystematics.C'('$regPar',"'$outputdirectorybase'","'$outputdirectorySignal'",'$doRawCutVar',1,1)'
  #before unfolding
#  root -l -b -q cutsSystematics.C'('$regPar',"'$outputdirectorybase'","'$outputdirectorySignal'",'$doRawCutVar',0,1)'
  cd $currDir

fi

exit 1

################################################
############### Get systematic uncertanties
################################################

#if [ $doSystematics -eq 1 ]; then
#	cd systematics
#		root -l -b -q JetSpectrumSys.C'('$regPar',"'$outputdirectory'","'$outputdirectorySignal'",1)'
#	cd $currDir
#fi

exit 1
