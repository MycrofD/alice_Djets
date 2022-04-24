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
outputdirectorybase=$HOME/Work/alice/analysis/pp5TeV/D0jet/results_APW/DzeroR06_
outputdirectory=${outputdirectorybase}QM19Cuts_Dpt3
#outputdirectory=${outputdirectorybase}pPbCuts
#outputdirectory=${outputdirectorybase}cuts2
outputdirectorySignal=Default

# ========== input directories
anaoutfiledir=$HOME/Work/alice/analysis/pp5TeV/D0jet/outData
MCoutfiledir=$HOME/Work/alice/analysis/pp5TeV/D0jet/outMC
refoutfiledir=$HOME/Work/alice/analysis/pp5TeV/D0jet/outMC/reflections
#simFilesDir=$HOME/Work/alice/analysis/pp5TeV/D0jet/outMC/POWHEGSimulations/fastSim_pp5TeV # POWHEG sim
#simFilesDir=/home/jackbauer/ALICE_HeavyFlavour/work/Djets/out/outMC/allR
simFilesDir=/home/jackbauer/ALICE_HeavyFlavour/work/Djets/out/outMC/all
bkgRMDir=$currDir/ResponseMatrix/BkgRM03

# ========== file names
###### R02
#analysisdatafile=AnalysisResults_632_R02.root #R02
#effMCfile=AnalysisResults_868_R02ppMC.root #R02 
###### R04
#analysisdatafileCut=CutSys/AnalysisResults_R04_693_cutsys.rtoot #R04
#efficiencyfileCut=CutSys/AnalysisResults_803_R04ppMCcuts.root #R04
#analysisdatafile=AnalysisResults_503_R04.root #R04
#effMCfile=CutSys/AnalysisResults_R04_693_cutsys.root #R04 
##### R06
analysisdatafileCut=CutSys/AnalysisResults_587_R06ppcuts.root #R06
efficiencyfileCut=CutSys/AnalysisResults_811_R06_cutsys2.root #R06
analysisdatafile=AnalysisResults_504_R06.root
effMCfile=AnalysisResults_871_R06ppMC.root #R06

detRMpromptfile1=newMC_JES/AnalysisResults_R06_721.root  # R06 JES
detRMnonpromptfile1=newMC_JES/AnalysisResults_R06_721.root  # R06 JES

#analysisdatafile=AnalysisResults_437.root #R03 ?
#analysisdatafile=AnalysisResults_LHC17pq_FASTwoSDD.root
#analysisdatafile=AnalysisResults_503_R04.root
isMoreFiles=0                                     # 1 if there are more input files to be read
prod=kl                                           # if there are more input files to be read
reflfile=reflectionTemplates_pp.root

#effMCfile=AnalysisResults_454.root #R03 ?
#effMCfile=AnalysisResults_fast_R03_D0MC_def.root
#effMCfile=AnalysisResults_634_pp5TeV_z.root
#effMCfile=AnalysisResults_642_pp5TeV_z.root
#effMCfile=AnalysisResults_683_ppMC_R06.root

#effMCfile=AnalysisResults_870_R04ppMC.root #R04
MCfile=$effMCfile #AnalysisResults_fast_R03_D0MC_def.root
#MCfile=AnalysisResults_fast_R03_D0MCPythia_JES96_1.root
efficiencyfile=$effMCfile
detRMpromptfile=$MCfile
detRMnonpromptfile=$MCfile
isprefix=0                                       #  if prefix is one, postfix must be one.
ispostfix=0                                       # if container in the analysis output file has different name than default you set here if and what is the postfix, this is set up in the signal, efficiency and RM extraction macros
postfix=HP2 #cut2
ispostfixFD=$ispostfix #0                                     # if container in the analysis output file has different name than default you set here if and what is the postfix, for the FD part wagons are usually configured with additional "FD" string in the container name, you should adjust this to yours configuration
#postfixFD=FD
postfixFD=$postfix #FDcut2
####

isRefl=1
isBkgRM=0
######## !!! POWHEG simulations config
nSimFilesB=10                                    # have to correspond to number of files defined in the config file. Def 10
nSimFilesC=9                                    # have to correspond to number of files defined in the config file. Def 9

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


#    detRMpromptfile1=AnalysisResults_fast_R03_D0MCPythia_JES96_1.root
 # 	detRMnonpromptfile1=AnalysisResults_fast_R03_D0MCPythia_JES96_1.root
#detRMnonpromptfile1=AnalysisResults_fast_R03_D0MC_def.root
  	eff=96
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
  ## R03
  #analysisdatafileCut=AnalysisResults_LHC16R03_CutVariation.root
  #efficiencyfileCut=AnalysisResults_fast_R03_D0MCPythia_CutVariation_1.root
  # R06
  #analysisdatafileCut=CutSys/AnalysisResults_587_R06ppcuts.root 
  #efficiencyfileCut=CutSys/AnalysisResults_811_R06_cutsys2.root
  detRMpromptfileCut=$efficiencyfileCut #AnalysisResults_fast_R03_D0MCPythia_CutVariation_1.root
  detRMnonpromptfileCut=$efficiencyfileCut #AnalysisResults_fast_R03_D0MCPythia_CutVariation_1.root
  isReflCut=0
  outputdirectorySignalCut=Default

  ispostfixCut=$ispostfix 
  ispostfixFDCut=$ispostfixFD 
  postfixCut=$postfix 
  postfixFDCut=$postfixFD 
  outputdirectoryCut=${outputdirectorybase}Cuts_VarDef
	./run.csh $outputdirectoryCut $outputdirectorySignalCut $conffile $Dmeson $anaoutfiledir $analysisdatafileCut $refoutfiledir $reflfile $isReflCut $isMoreFiles $prod $ispostfixCut $postfixCut $ispostfixFDCut $postfixFDCut $MCoutfiledir $efficiencyfileCut $detRMpromptfileCut $detRMnonpromptfileCut $isBkgRM $bkgRMtype $bkgRMDir $bkgRMfileName $simFilesDir $unfType $regPar $isPrior $priorType 0 0 $doRawSpectra $nSimFilesB $nSimFilesC $roounfoldpwd $isprefix

arraySize
#declare -a cutarray=( L0 L1 L2 L3 T0 T1 T2 ) #declare postfixes
declare -a cutarray=( HP2 ) #HP3 HP5 HP6 ) #declare postfixes

for i in "${cutarray[@]}"
do
  ispostfixCut=1
  isprefixCut=1
  ispostfixFDCut=1
  postfixCut=${cutarray[$i]}
  postfixFDCut=${cutarray[$i]}
  outputdirectoryCut=${outputdirectorybase}Cuts_${cutarray[$i]}
  ./run.csh $outputdirectoryCut $outputdirectorySignalCut $conffile $Dmeson $anaoutfiledir $analysisdatafileCut $refoutfiledir $reflfile $isReflCut $isMoreFiles $prod $ispostfixCut $postfixCut $ispostfixFDCut $postfixFDCut $MCoutfiledir $efficiencyfileCut $detRMpromptfileCut $detRMnonpromptfileCut $isBkgRM $bkgRMtype $bkgRMDir $bkgRMfileName $simFilesDir $unfType $regPar $isPrior $priorType 0 0 $doRawSpectra $nSimFilesB $nSimFilesC $roounfoldpwd $isprefixCut
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
