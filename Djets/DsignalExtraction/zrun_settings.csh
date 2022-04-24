#! /bin/bash

##################################### MAIN ANALYSISRESULTS.ROOT files
if [ $R -eq 2 ]; then
 data=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_744_R02_pp_5cuts.root 
 effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_1059_R02_ppMC_5cuts.root
 effFileJES=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/newMC_JES/AnalysisResults_1060_R02_JESMC.root #JES 4%
elif [ $R -eq 3 ]; then
 #data=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_757_R03_pp_5cuts.root 
 #effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_1071_R03_ppMC_5cuts.root
 #data_pPb=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/test/AnalysisResults_465_pPb.root
 #data_PbPb=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/test/AnalysisResults_40_PbPb.root
 data=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_769_R03_pp_TPC.root
 effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_1090_R03_ppMC_TPC.root
elif [ $R -eq 4 ]; then
 data=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_745_R04_pp_5cuts.root 
 effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_1061_R04_ppMC_5cuts.root 
 effFileJES=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/newMC_JES/AnalysisResults_R04_JES4_722.root #JES 4%
 data_pp13=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/test/AnalysisResults_pp13TeV.root #13 TeV
elif [ $R -eq 6 ]; then
 data=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_751_R06_pp_5cuts.root 
 effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_1062_R06_ppMC_5cuts.root 
 effFileJES=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/newMC_JES/AnalysisResults_R06_721.root #JES 4%
fi
##################################### GENERAL STUFF
#doBFD=0
#doUnfold=0
#OUT=${OUT}_finaltry #OUT is in zrun4 and zrun4D
##################################### general, changed later for CUTSYS
listName="" #SQ2 3 5 6 
isprefix=0 #
ispostfix=0 #
##################################### RAW SYS multi trial
flagMulti=0
flagSBSig=1
##################################### JES SYS
flagJES=0
##################################### CUT SYS
# remember to change ND (number of Dmesons in config file) if needed
flagCUT=0
cutfileNo=0 #2,3,5,6
##################################### THINGS TO CHANGE
if [ $flagMulti -eq 1 ]; then
    OUT=${OUT}/RawSys_Multi
elif [ $flagSBSig -eq 1 ]; then
    OUT=${OUT}/RawSys_SBSig
elif [ $flagJES -eq 1 ]; then
    OUT=${OUT}/JES
    effFile=$effFileJES
elif [ $flagCUT -eq 1 ]; then
    isprefix=1
    ispostfix=1
    if [ $cutfileNo -eq 2 ]; then
        #OUT=${OUT}_L2
        OUT=${OUT}/SQ2
        listName="SQ2"
    elif [ $cutfileNo -eq 3 ]; then
        #OUT=${OUT}_L3
        OUT=${OUT}/SQ3
        listName="SQ3"
    elif [ $cutfileNo -eq 5 ]; then
        #OUT=${OUT}_T2
        OUT=${OUT}/SQ5
        listName="SQ5"
    elif [ $cutfileNo -eq 6 ]; then
        #OUT=${OUT}_T3
        OUT=${OUT}/SQ6
        listName="SQ6"
    fi
    #OUT=${OUT}/cutfix_data
    #OUT=${OUT}_cutfix_MC
fi
##################################### THINGS TO CHANGE
