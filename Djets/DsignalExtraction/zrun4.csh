#! /bin/bash

#--------------------------------------------------
#  Author A.Mohanty
#  Utrecht University
#  auro.mohanty@cern.ch
#--------------------------------------------------
#-- script to run D-jets z_||,ch analysis
#-- (development stage)

#- Instructions.
#---------------
# Executing 'bash zrun_main.csh' on your terminal in AliPhysics environment,
# in turn, runs this file 'zrun.csh'

EOS_local=$3
R=$2		# Jet radius, fed as an argument from zrun_main.csh 
#/eos/user/a/amohanty/

## Writing the configDzero_pp.h file: Setting up D pT bins
##--------------------------------------------------------
bin_zjet=$1 # change the bin number you want to use here.
        # 0 means no cut on jet pt
        # 1 through x mean all bins of jet pt
        # The arguments here are passed from zrun_main.csh
#----- Flags for all analysis steps
flagEff=1 #1 Getting efficiencies and MC sigmas in different jetpt intervals for all Dpt bins
flagRef=1 #2 Reflections for different jetpt intervals for all Dpt bins
flagSBs=1 #3 Side Band subtraction method
flagSim=1 #4 Simulation for non-prompt and prompt D-jets
#-----
# change flagCUT, flagJES  in zrun_settings.csh
#-----
finerunfold=0
boundSigma=6 #also see bin_zjet=24	# if needed to fit certain Dpt bins with a bounded sigma: sigma +/- some fraction of this sigma
if [ $bin_zjet -eq 24 ]; then #2-5
 boundSigma=0
fi
fBsimN=10 #11	# number of non-prompt sim files
fCsimN=1 #9#11	# number of prompt sim files
unoriginal_tag=unoriginal	# adding a tag that says the next file zrun.csh is not the original file to be edited, rather this is
conffile_z1=configDzero_ppz1.h	# first half of config file
conffile_z2=configDzero_ppz2.h	# second half of config file
cat $unoriginal_tag > configDzero_ppz.h
cat $conffile_z1 >> configDzero_ppz.h
echo "const double  fRpar = 0.${R};" >> configDzero_ppz.h	# writing after first half of config file
echo "const int  Rpar = ${R};" >> configDzero_ppz.h		# 
echo "const int fBsimN = ${fBsimN};" >> configDzero_ppz.h	#
echo "const int fCsimN = ${fCsimN};" >> configDzero_ppz.h	#
echo "int zjetbin = ${bin_zjet};" >> configDzero_ppz.h		#

# defining Dpt intervals and, measured and unfolding true bins for different jetpt intervals
if [ $bin_zjet -eq 26 ]; then  # full jet bin		# before it was 14
 echo "double        fptbinsDA[] = {1,2,3,4,5,6,7,8,10,12,16,24,36};" >> configDzero_ppz.h
 echo " double fdplotmin = 2000, dataplotmax = 200000;" >> configDzero_ppz.h
elif [ $bin_zjet -eq 24 ]; then #2-5
 echo "double        fptbinsDA[] = {2,3,4,5};" >> configDzero_ppz.h
 echo " double fdplotmin = 2000, dataplotmax = 200000;" >> configDzero_ppz.h
elif [ $bin_zjet -eq 16 ]; then #3-5
 echo "double        fptbinsDA[] = {2,3,4,5};" >> configDzero_ppz.h
 echo " double fdplotmin = 2000, dataplotmax = 200000;" >> configDzero_ppz.h
elif [ $bin_zjet -eq 1 ]; then #5-7
 echo "double        fptbinsDA[] = {2,3,4,5,6,7};" >> configDzero_ppz.h
 echo " double fdplotmin = 2000, dataplotmax = 200000;" >> configDzero_ppz.h
elif [ $bin_zjet -eq 2 ]; then #7-10
 if [ $R -eq 2 ]; then
  echo "double        fptbinsDA[] = {4,5,6,7,8,10};" >> configDzero_ppz.h # report from z = 0.4
 else
 echo "double        fptbinsDA[] = {3,4,5,6,7,8,10};" >> configDzero_ppz.h # report from z = 0.4
 fi
 echo " double fdplotmin = 800, dataplotmax = 100000;" >> configDzero_ppz.h
elif [ $bin_zjet -eq 3 ]; then #10-15
 echo "double        fptbinsDA[] = {5,6,7,8,10,12,15};" >> configDzero_ppz.h
 echo " double fdplotmin = 100, dataplotmax = 40000;" >> configDzero_ppz.h
elif [ $bin_zjet -eq 6 ]; then #5-15
 echo "double        fptbinsDA[] = {2,3,4,5,6,7,8,10,12,15};" >> configDzero_ppz.h
 echo " double fdplotmin = 1000, dataplotmax = 200000;" >> configDzero_ppz.h
elif [ $bin_zjet -eq 9 ]; then #15-50 ...
 if [ $R -eq 2 ]; then
  #echo "double        fptbinsDA[] = {8,10,12,16,24,36};" >> configDzero_ppz.h
  echo "double        fptbinsDA[] = {10,12,16,24,36};" >> configDzero_ppz.h
 else 
  echo "double        fptbinsDA[] = {5,8,10,12,16,24,36};" >> configDzero_ppz.h
 fi 
 echo " double fdplotmin = 20, dataplotmax = 20000;" >> configDzero_ppz.h
fi
cat $conffile_z2 >> configDzero_ppz.h

# Output directory
# Remember to change this also in zrun4D.csh
OUT=${EOS_local}/media/jackbauer/data/z_out/R_0$R
OUT=${OUT}_finaltry
if [ $boundSigma -eq 1 ]; then
 OUT=${OUT}_boundSigma1
elif [ $boundSigma -eq 2 ]; then
 OUT=${OUT}_boundSigmaall10
elif [ $boundSigma -eq 3 ]; then
 OUT=${OUT}_boundSigmaall20
fi
## Getting data and MC files
##-----------------------------------------------------
source zrun_settings.csh

## Getting Efficiencies for each bin of jet pt interval
##-----------------------------------------------------------
outDir=$OUT/efficiency
if [ $flagEff -eq 1 ]; then
 cd ../efficiency
 root -l -b -q DjetEfficiency_z.C'(1, "'$effFile'","'$outDir'", 0, '$ispostfix', "'$listName'", '$isprefix')'
 root -l -b -q DjetEfficiency_z.C'(0, "'$effFile'","'$outDir'", 0, '$ispostfix', "'$listName'", '$isprefix')'
 cd ../DsignalExtraction
fi
## Reflections
##------------
outRefl=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/reflections
listNameRef=""
kl=""
if [ $flagRef -eq 1 ]; then
 root -l -b -q signalExtraction_zrefl.C'("'$effFile'","'$outRefl'",0,"'$listNameRef'",0,"'$kl'")'
fi
## Signal Extraction Side Bands, with/without efficiency correction
##---------------------------------------------------------
#prompteff=/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results/DzeroR03_pPbCuts/Default/efficiency/DjetEff_prompt_jetpt5_50.root
prompteff=$outDir/DjetEff_prompt_jetpt 
refFile=$outRefl/reflectionTemplates_pp_
if [ $flagRef -eq 1 ]; then
 out=$OUT/signalExtraction
elif [ $flagRef -eq 0 ]; then
 out=$OUT/signalExtraction_noRef
fi
save=1
isMoreFiles=0
prod=kl
saveDir=Z0to102

if [ $flagSBs -eq 1 ]; then
 root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix', "'$saveDir'", '$boundSigma')'
fi
## B-feed down simulation
##-----------------------
##/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/POWHEGSimulations/fastSim_pp5TeV
effFilePrompt=$prompteff
effFileNonPrompt=$outDir/DjetEff_nonPrompt_jetpt
isBsim=1
if [ $flagSim -eq 1 ]; then
 cd ../POWHEGSim
 bash zrun_sim.csh $OUT $fBsimN $fCsimN $effFilePrompt $effFileNonPrompt $flagEff $isBsim
 cd ../DsignalExtraction
fi

###################################################################
#
#               THE END 
#
###################################################################
###################################################################

##for running custom cxx macros from AliPhysics
##gInterpreter->ProcessLine(".x AliHFInvMassFitterDJET.cxx++g");
#if [ $flagSBs -eq 1 ]; then
#root -l -b << EOF
#gInterpreter->ProcessLine(".x AliHFInvMassFitterDJET.cxx++g");
#.L signalExtraction_SBz.C
#signalExtraction_SBz("$data", $flagEff, "$prompteff", $flagRef, "$refFile", $ispostfix, "$listName", "$out", $save, $isMoreFiles, "$prod", $isprefix, "$saveDir", $boundSigma)
#EOF
#fi

