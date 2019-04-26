#!/bin/bash
# ---------------------
# A.Mohanty
# Utrecht University
# auro.mohanty@cern.ch
# ---------------------

# Getting simulations for prompt and non-prompt
# ---------------------------------------------
OUT=$1

isBsim=$7            # switch this flag on if you haven't prepared output of the simulations yet with current efficiencies, it takes time -- if the simulation output directory is empty the simulations will be run anyway
isCsim=1            # switch this flag on if you haven't prepared output of the simulations yet with current efficiencies, it takes time -- if the simulation output directory is empty the simulations will be run anyway
currDir=`pwd`
MainDir=/home/jackbauer/ALICE_HeavyFlavour/work/Djets/alice_Djets/Djets
SimDir=$MainDir/POWHEGSim
#outdirBase=/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet_pp5TeV_z/results
nSimFilesB=$2     # have to correspond to number of files defined in the config file
nSimFilesC=$3     # have to correspond to number of files defined in the config file
BFDSimDirOut=$OUT/SimFiles/BFeedDown #$outdirBase/Simulations/BFeedDown
CSimDirOut=$OUT/SimFiles/Prompt
simFilesDir=/home/jackbauer/ALICE_HeavyFlavour/work/Djets/out/outMC/beauty/
simFilesDirC=/home/jackbauer/ALICE_HeavyFlavour/work/Djets/out/outMC/charm/
#effFilePrompt=/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results/DzeroR03_pPbCuts/Default/efficiency/DjetEff_prompt_jetpt5_50.root
#effFileNonPrompt=/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results/DzeroR03_pPbCuts/Default/efficiency/DjetEff_nonPrompt_jetpt5_50.root
effFilePrompt=$4
effFileNonPrompt=$5
isEff=$6
############### B feed-down simulations
################################################

if [ ! -d "$BFDSimDirOut" ] || [ $isBsim -eq 1 ]; then

cd $SimDir
# 1: number of simulation files (as defined in the config file)
# 2: simulation file directory;
# 3: quark type: 0: charm, 1: beauty
# 4: zfrac:: D-jet z spectrum: 1, D-meson pT spectrum: 0
# 5: if D meson pT cut applied (for D-jet pT spectrum case), if yes the lower and upper values from the D-meson pT bins from the config file are taken
# 6: if efficiency applied (for B simulations, ratio of non-prompt/prompt efficiency)
# 7: prompt efficiency file
# 8: non-prompt efficiency file
# 9: directory for the output files

# ./doGetSimOut_z.csh $nSimFilesB $simFilesDir 1 1 1 1 $effFilePrompt $effFileNonPrompt $BFDSimDirOut
 ./doGetSimOut_z.csh $nSimFilesB $simFilesDir 1 1 1 $isEff $effFilePrompt $effFileNonPrompt $BFDSimDirOut
#plot all the variations

elif [ ! -d "$CSimDirOut" ] || [ $isCsim -eq 1 ]; then
cd $SimDir
 ./doGetSimOut_z.csh $nSimFilesC $simFilesDirC 0 1 1 $isEff $effFilePrompt $effFileNonPrompt $CSimDirOut

cd $currDir
fi
