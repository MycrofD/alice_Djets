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
## Writing the configDzero_ppz.h file: Setting up D pT bins
##--------------------------------------------------------
bin_zjet=$1 # The arguments here are passed from zrun_main.csh
R=$2        # Jet radius, fed as an argument from zrun_main.csh 
OUT=$3      # Output directory
##----- Flags for all analysis steps
flagEff=$4 #1 Getting efficiencies and MC sigmas in different jetpt intervals for all Dpt bins
flagRef=$5 #2 Reflections for different jetpt intervals for all Dpt bins
flagSBs=$6 #3 Side Band subtraction method
flagSim=$7 #4 Simulation for non-prompt and prompt D-jets
#-----
#-----
# change flagCUT, flagJES  in zrun_settings.csh
finerunfold=0
boundSigma=0 #also see bin_zjet=24	# if needed to fit certain Dpt bins with a bounded sigma: sigma +/- some fraction of this sigma
#if [ $R -eq 2 ]; then
#  if [ $bin_zjet -eq 2 ]; then
#    boundSigma=5
#  fi
#fi
## boundSigma is 5 below for cut systematics
if [ $bin_zjet -eq 24 ]; then #2-5
 boundSigma=0
fi
fBsimN=10 #11	# number of non-prompt sim files
fCsimN=12 #9#11	# number of prompt sim files: 9 pow+pyt, 1 pyth6, 1 pyth8, 1 pyth 8-2
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

##### different Dpt ranges for different R
if [ $R -eq 2 ]; then
 echo "double fDptRangesA[] = {2,2,4,5,10};//Dpt mins in jetpt bins, 36 is max in last bin">> configDzero_ppz.h #
else
 echo "double fDptRangesA[] = {2,2,3,5,5};//Dpt mins in jetpt bins, 36 is max in last bin">> configDzero_ppz.h #
fi
##### defining Dpt intervals and, measured and unfolding true bins for different jetpt intervals
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

## Getting data and MC files SOURCING THE SETTINGS
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
outRefl=/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/reflections
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

if [ $flagCUT -eq 1 ]; then
    boundSigma=5
fi

if [ $flagSBs -eq 1 ]; then
 root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix', "'$saveDir'", '$boundSigma')'
fi
# Multi trial
if [ $flagMulti -eq 1 ]; then
out=$OUT/signalExtraction
ctry=1
fsigmafactor=1

# (1)*2*2*2*2*2=32
  for boundSigma in 0; do
    for bkgType in 0 1; do
      for fixedMass in 0 1; do
        for minfSys in 1.72 1.74; do
          for maxfSys in 2.00 2.03; do
            for fMassBinWidthFactor in 2 4; do

root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix', "'$saveDir'", '$boundSigma', '$fsigmafactor', '$fixedMass', '$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$flagMulti')'

ctry=$[ctry+1]

            done
          done
        done
      done
    done
  done

# (1*1*)2*2*2*2*2=32
  for boundSigma in 1; do
  for fsigmafactor in 1; do
    for bkgType in 0 1; do
      for fixedMass in 0 1; do
        for minfSys in 1.72 1.74; do
          for maxfSys in 2.00 2.03; do
            for fMassBinWidthFactor in 2 4; do

root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix', "'$saveDir'", '$boundSigma', '$fsigmafactor', '$fixedMass', '$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$flagMulti')'
ctry=$[ctry+1]

            done
          done
        done
      done
    done
  done
  done

# (1*)*2*2*1*2*2*2=32
  for boundSigma in 1; do
  for fsigmafactor in 1.1 0.9; do
    for bkgType in 0 1; do
      for fixedMass in 0; do
        for minfSys in 1.72 1.74; do
          for maxfSys in 2.00 2.03; do
            for fMassBinWidthFactor in 2 4; do

root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix', "'$saveDir'", '$boundSigma', '$fsigmafactor', '$fixedMass', '$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$flagMulti')'
ctry=$[ctry+1]

            done
          done
        done
      done
    done
  done
  done


# Side-Band Signal ranges
elif [ $flagSBSig -eq 1 ]; then
boundSigma=0;fsigmafactor=1;fixedMass=0;bkgType=0;minfSys=1.71;maxfSys=2.1;fMassBinWidthFactor=2
out=$OUT/signalExtraction
    root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix', "'$saveDir'", '$boundSigma')'
ctry=1
#  for sigmaWindow in 2 3; do
#    for SBout in 9 8; do
#      for SBint in 4 3.5 4.5; do
#root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix',"'$saveDir'", '$boundSigma', '$fsigmafactor', '$fixedMass', '$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$flagSBSig','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
#ctry=$[ctry+1]
#      done
#    done
#  done
#fi

# 2, (3,8)
sigmaWindow=2;SBout=8;SBint=3
    root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix',"'$saveDir'", '$boundSigma', '$fsigmafactor', '$fixedMass', '$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$flagSBSig','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]
# 2, (4,10)
sigmaWindow=2;SBout=10;SBint=4
    root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix',"'$saveDir'", '$boundSigma', '$fsigmafactor', '$fixedMass', '$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$flagSBSig','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]
# 2, (4,12)
sigmaWindow=2;SBout=12;SBint=4
    root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix',"'$saveDir'", '$boundSigma', '$fsigmafactor', '$fixedMass', '$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$flagSBSig','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]
# 2, (4,15)
sigmaWindow=2;SBout=15;SBint=4
    root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix',"'$saveDir'", '$boundSigma', '$fsigmafactor', '$fixedMass', '$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$flagSBSig','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]
# 3, (4,9)
sigmaWindow=3;SBout=9;SBint=4
    root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix',"'$saveDir'", '$boundSigma', '$fsigmafactor', '$fixedMass', '$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$flagSBSig','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]
# 3, (4,12)
sigmaWindow=3;SBout=12;SBint=4
    root -l -b -q signalExtraction_SBz.C'("'$data'", '$flagEff', "'$prompteff'", '$flagRef', "'$refFile'", '$ispostfix', "'$listName'", "'$out'", '$save', '$isMoreFiles', "'$prod'", '$isprefix',"'$saveDir'", '$boundSigma', '$fsigmafactor', '$fixedMass', '$bkgType','$minfSys','$maxfSys','$fMassBinWidthFactor','$ctry','$flagSBSig','$sigmaWindow','-$SBout','-$SBint','$SBint','$SBout')'
ctry=$[ctry+1]

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

