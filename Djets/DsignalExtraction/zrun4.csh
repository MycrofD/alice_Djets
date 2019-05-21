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
#/eos/user/a/amohanty/

## Writing the configDzero_pp.h file: Setting up D pT bins
##--------------------------------------------------------
bin_zjet=$1 # change the bin number you want to use here.
		# 0 means no cut on jet pt
		# 1 through x mean all bins of jet pt
		# The arguments here are passed from zrun_main.csh
#----- Flags for all analysis steps
flagEff=0 #1 Getting efficiencies and MC sigmas in different jetpt intervals for all Dpt bins
flagRef=0 #2 Reflections for different jetpt intervals for all Dpt bins
flagSBs=0 #3 Side Band subtraction method
flagSim=0 #4 Simulation for non-prompt and prompt D-jets
flagRes=1 #5 Response matrix
flagBFD=0 #6 B-feed down subtraction
flagUnf=1 #7 Unfolding
#-----
flagRawSys=0 #8 Raw Yield Systematics
#-----
finerunfold=0
boundSigma=0	# if needed to fit certain Dpt bins with a bounded sigma: sigma +/- some fraction of this sigma
R=$2		# Jet radius, fed as an argument from zrun_main.csh 
fBsimN=11	# number of non-prompt sim files
fCsimN=9	# number of prompt sim files
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
if [ $bin_zjet -eq 14 ]; then  # full jet bin			 
 echo "double        fptbinsDA[] = {1,2,3,4,5,6,7,8,10,12,16,24,36};" >> configDzero_ppz.h
 echo " const int     fptbinsZFinalN = 6;" >> configDzero_ppz.h
 echo " double        fptbinsZFinalA[fptbinsZFinalN+1] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.02};" >> configDzero_ppz.h
 echo " double fdplotmin = 2000, dataplotmax = 200000;" >> configDzero_ppz.h
elif [ $bin_zjet -eq 1 ]; then #5-7
 echo "double        fptbinsDA[] = {2,3,4,5,6,7};" >> configDzero_ppz.h
 echo " const int     fptbinsZFinalN = 6;" >> configDzero_ppz.h
 echo " double        fptbinsZFinalA[fptbinsZFinalN+1] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.02};" >> configDzero_ppz.h
 echo " double fdplotmin = 2000, dataplotmax = 200000;" >> configDzero_ppz.h
elif [ $bin_zjet -eq 2 ]; then #7-10
 echo "double        fptbinsDA[] = {3,4,5,6,7,8,10};" >> configDzero_ppz.h # report from z = 0.4
 echo " const int     fptbinsZFinalN = 6;" >> configDzero_ppz.h
 echo " double        fptbinsZFinalA[fptbinsZFinalN+1] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.02};" >> configDzero_ppz.h
 echo " double fdplotmin = 800, dataplotmax = 100000;" >> configDzero_ppz.h
elif [ $bin_zjet -eq 3 ]; then #10-15
 echo "double        fptbinsDA[] = {5,6,7,8,10,12,15};" >> configDzero_ppz.h
 if [ $finerunfold -eq 0 ]; then #10-15
  echo " const int     fptbinsZFinalN = 5;" >> configDzero_ppz.h
  echo " double        fptbinsZFinalA[fptbinsZFinalN+1] = {0.4, 0.6, 0.7, 0.8, 0.9, 1.02};" >> configDzero_ppz.h
 elif [ $finerunfold -eq 1 ]; then #10-15
  echo " const int     fptbinsZFinalN = 6;" >> configDzero_ppz.h
  echo " double        fptbinsZFinalA[fptbinsZFinalN+1] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.02};" >> configDzero_ppz.h
 fi
 echo " double fdplotmin = 100, dataplotmax = 40000;" >> configDzero_ppz.h
elif [ $bin_zjet -eq 6 ]; then #5-15
 echo "double        fptbinsDA[] = {2,3,4,5,6,7,8,10,12,15};" >> configDzero_ppz.h
 echo " const int     fptbinsZFinalN = 6;" >> configDzero_ppz.h
 echo " double        fptbinsZFinalA[fptbinsZFinalN+1] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.02};" >> configDzero_ppz.h
 echo " double fdplotmin = 1000, dataplotmax = 200000;" >> configDzero_ppz.h
elif [ $bin_zjet -eq 9 ]; then #15-50 ...
 echo "double        fptbinsDA[] = {5,8,10,12,16,24,36};" >> configDzero_ppz.h
 
 echo " const int     fptbinsZFinalN = 6;" >> configDzero_ppz.h
 echo " double        fptbinsZFinalA[fptbinsZFinalN+1] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.02};" >> configDzero_ppz.h
 echo " double fdplotmin = 20, dataplotmax = 20000;" >> configDzero_ppz.h
fi
cat $conffile_z2 >> configDzero_ppz.h

# Output directory
OUT=${EOS_local}/media/jackbauer/data/z_out/R_0$R
if [ $boundSigma -eq 1 ]; then
 OUT=${OUT}_boundSigma1
elif [ $boundSigma -eq 2 ]; then
 OUT=${OUT}_boundSigmaall10
elif [ $boundSigma -eq 3 ]; then
 OUT=${OUT}_boundSigmaall20
fi
## Getting Efficiencies for each bin of jet pt interval
##-----------------------------------------------------
if [ $R -eq 3 ]; then
 data=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/trial_437.root
 effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_634_pp5TeV_z.root
elif [ $R -eq 4 ]; then
 data=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_503_R04.root
 effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_642_pp5TeV_z.root
elif [ $R -eq 6 ]; then
 data=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_504_R06.root
 effFile=${EOS_local}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_683_ppMC_R06.root
fi
isPrompt=1
outDir=$OUT/efficiency
postfix=0
listName=""
isprefix=0
if [ $flagEff -eq 1 ]; then
 cd ../efficiency
 root -l -b -q DjetEfficiency_z.C'(1, "'$effFile'","'$outDir'", 0, '$postfix', "'$listName'", '$isprefix')'
 root -l -b -q DjetEfficiency_z.C'(0, "'$effFile'","'$outDir'", 0, '$postfix', "'$listName'", '$isprefix')'
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
isEff=1
#efffile=/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results/DzeroR03_pPbCuts/Default/efficiency/DjetEff_prompt_jetpt5_50.root
efffile=$outDir/DjetEff_prompt_jetpt
isRef=$flagRef
refFile=$outRefl/reflectionTemplates_pp_
postfix=0
listName=Cut
if [ $flagRef -eq 1 ]; then
out=$OUT/signalExtraction
elif [ $flagRef -eq 0 ]; then
out=$OUT/signalExtraction_noRef
fi
save=1
isMoreFiles=0
prod=kl
isprefix=0
saveDir=Z0to102

##for running custom cxx macros from AliPhysics
##gInterpreter->ProcessLine(".x AliHFInvMassFitterDJET.cxx++g");
#if [ $flagSBs -eq 1 ]; then
#root -l -b << EOF
#gInterpreter->ProcessLine(".x AliHFInvMassFitterDJET.cxx++g");
#.L signalExtraction_SBz.C
#signalExtraction_SBz("$data", $isEff, "$efffile", $isRef, "$refFile", $postfix, "$listName", "$out", $save, $isMoreFiles, "$prod", $isprefix, "$saveDir", $boundSigma)
#EOF
#fi

if [ $flagSBs -eq 1 ]; then
 root -l -b -q signalExtraction_SBz.C'("'$data'", '$isEff', "'$efffile'", '$isRef', "'$refFile'", '$postfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix', "'$saveDir'", '$boundSigma')'
fi
## B-feed down simulation
##-----------------------
##/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/POWHEGSimulations/fastSim_pp5TeV
effFilePrompt=$efffile
effFileNonPrompt=$outDir/DjetEff_nonPrompt_jetpt
isBsim=1
if [ $flagSim -eq 1 ]; then
 cd ../POWHEGSim
 bash zrun_sim.csh $OUT $fBsimN $fCsimN $effFilePrompt $effFileNonPrompt $isEff $isBsim
 cd ../DsignalExtraction
fi
## Response
##------------------------
#dataFile=/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_414_z.root
dataFile=$effFile
outDirRM=$OUT/ResponseMatrix
if [ $flagRes -eq 1 ]; then
 cd ../ResponseMatrix
 bash zrun_RM.csh 0 $dataFile $outDirRM
 bash zrun_RM.csh 1 $dataFile $outDirRM
 cd ../DsignalExtraction
fi

## B-feed down subtraction
##------------------------
dataFileFDdir=$out/plots/$saveDir
if [ $flagBFD -eq 1 ]; then
 cd ../FDsubtraction
 bash zrun_FD.csh $dataFileFDdir $data $OUT
 cd ../DsignalExtraction
fi

## Unfolding
##----------
#dataFileFDdir=$out/plots/$saveDir
dataUnfoldInDir=$OUT/FDsubtraction
detRMFilePrompt=$OUT/ResponseMatrix/DetMatrix_prompt
bkgRMFile=$OUT/ResponseMatrix/DetMatrix_
 if [ $finerunfold -eq 0 ]; then #10-15
unfoldingDirOut=$OUT/unfolding
 elif [ $finerunfold -eq 1 ]; then #10-15
unfoldingDirOut=$OUT/unfolding_finer
 fi
regPar=5
isPrior=0
priorType=0 
isBkgRM=0
isFDUpSys=0 
isFDDownSys=0

if [ $flagUnf -eq 1 ]; then
 cd ../unfolding
 bash zrun_unfold.csh $dataUnfoldInDir $detRMFilePrompt $bkgRMFile $unfoldingDirOut $regPar $isPrior $priorType $isBkgRM $isFDUpSys $isFDDownSys
 cd ../DsignalExtraction
fi


#####--- Systematics ---#####
## Raw Yield Systematics
##---------------------------
#if [ $flagRawSys -eq 1 ]; then
# cd ../systematics/YieldExtraction
# bash zrun_rawsys.csh 
# cd ../../DsignalExtraction
#fi
