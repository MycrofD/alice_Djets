#! /bin/bash

#--------------------------------------------------
#  Author A.Mohanty
#  Utrecht University
#  auro.mohanty@cern.ch
#--------------------------------------------------
#source zrun.csh

doBFD=1
doUnfold=0

R=$1
EOS_local=$2
#/eos/user/a/amohanty/
OUT=${EOS_local}/media/jackbauer/data/z_out/R_0$R
#OUT=${OUT}_36
#OUT=${OUT}_35
#OUT=${OUT}_34
#OUT=${OUT}_25
#OUT=${OUT}_25UnfCheck
OUT=${OUT}_finaltry

##################################### THINGS TO CHANGE
flagCUT=0
cutfileNo=0
#if [ $flagCUT -eq 1 ]; then
#OUT=${OUT}_cutsysDeDe
#fi 

if [ $flagCUT -eq 1 ]; then
 data=${dataCUT2}
 effFile=${effCUT2}
 if [ $cutfileNo -eq 1 ]; then
  OUT=${OUT}_cutsysLL1
 elif [ $cutfileNo -eq 2 ]; then 
  OUT=${OUT}_cutsysLL2
 elif [ $cutfileNo -eq 3 ]; then 
  OUT=${OUT}_cutsysLL3
 elif [ $cutfileNo -eq 4 ]; then 
  OUT=${OUT}_cutsysTT2
 elif [ $cutfileNo -eq 0 ]; then 
  OUT=${OUT}_cutsysDeDe
# elif [ $cutfileNo -eq 5 ]; then
#  OUT=${OUT}_cutsysT3
 fi
fi 

#####################################
dataUnfoldInDir=$OUT/FDsubtraction
detRMFilePrompt=$OUT/ResponseMatrix/DetMatrix_prompt
bkgRMFile=$OUT/ResponseMatrix/DetMatrix_
unfoldingDirOut=$OUT/unfolding
regPar=5
isPrior=0
priorType=1 #1 to 8
isBkgRM=0
isFDUpSys=0
isFDDownSys=0


if [ $R -eq 2 ]; then
 data=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_632_R02.root 
 effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_868_R02ppMC.root                                                                                                       
elif [ $R -eq 3 ]; then
 data=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/trial_437.root
 #effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_634_pp5TeV_z.root
 effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_454.root                                                                                                               
elif [ $R -eq 4 ]; then
 data=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_503_R04.root
 #effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_642_pp5TeV_z.root # default
 effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/CutSys/AnalysisResults_R04_693_cutsys.root 
 effFileJES=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/newMC_JES/AnalysisResults_R04_JES4_722.root #JES 4%                                                                                 
  
 dataCUT=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/CutSys/AnalysisResults_R04_515_cutsys.root 
 effCUT=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/CutSys/AnalysisResults_R04_693_cutsys.root 
 dataCUT2=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/CutSys/AnalysisResults_584_R04ppcuts.root  
 effCUT2=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/CutSys/AnalysisResults_803_R04ppMCcuts.root                                                                                            
elif [ $R -eq 6 ]; then
 data=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_504_R06.root
 effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_683_ppMC_R06.root
 effFileJES=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/newMC_JES/AnalysisResults_R06_721.root #JES 4%                                                                                      
  
 dataCUT=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/CutSys/AnalysisResults_R06_522_cutsys.root 
 effCUT=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet_pp5TeV_z/cutsystematics_z/AnalysisResultsR06MC_cuts719.root 
 dataCUT2=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/CutSys/AnalysisResults_587_R06ppcuts.root  
 effCUT2=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/CutSys/AnalysisResults_811_R06_cutsys2.root                                                                                            
fi

if [ $flagCUT -eq 1 ]; then
 data=${dataCUT2}
 effFile=$effCUT2
fi

## B FEED DOWN=================================
if [ $doBFD -eq 1 ]; then
 cd ../FDsubtraction
  bash zrun_BFD5D.csh 
 cd ../DsignalExtraction
fi

## B FEED DOWN ends=================================
## UNFOLDING=================================
if [ $doUnfold -eq 1 ]; then
 cd ../unfolding
  bash zrun_unfold4D.csh $dataUnfoldInDir $detRMFilePrompt $bkgRMFile $unfoldingDirOut $regPar $isPrior $priorType $isBkgRM $isFDUpSys $isFDDownSys $effFile
 cd ../DsignalExtraction
 
 if [ $isPrior -eq 1 ]; then
 for thing in 1 2 3 4 5 6 7 8; do #these are priors
   cd ../unfolding
     bash zrun_unfold4D.csh $dataUnfoldInDir $detRMFilePrompt $bkgRMFile $unfoldingDirOut $regPar 1 $thing $isBkgRM $isFDUpSys $isFDDownSys $effFile
   cd ../DsignalExtraction
 done
 fi
fi
## UNFOLDING ends=================================
