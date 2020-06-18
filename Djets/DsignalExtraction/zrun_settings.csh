#! /bin/bash

##################################### MAIN ANALYSISRESULTS.ROOT files
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
##################################### THINGS TO CHANGE
isprefix=0
ispostfix=0
##################################### THINGS TO CHANGE
# remember to change ND (number of Dmesons in config file)
flagCUT=0 # 
cutfileNo=0 #2,3,5,6
listName="" #SQ2 3 5 6
if [ $flagCUT -eq 1 ]; then
 data=${dataCUT}
 effFile=${effCUT}
 if [ $cutfileNo -eq 2 ]; then
  OUT=${OUT}_cutsL2
  listName="SQ2"
 elif [ $cutfileNo -eq 3 ]; then
  OUT=${OUT}_cutsL3
  listName="SQ3"
 elif [ $cutfileNo -eq 5 ]; then
  OUT=${OUT}_cutsT2
  listName="SQ5"
 elif [ $cutfileNo -eq 6 ]; then
  OUT=${OUT}_cutsT3
  listName="SQ6"
 fi
fi
 
##################################### THINGS TO CHANGE
