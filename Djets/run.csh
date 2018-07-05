#!/bin/bash

# B. Trzeciak -- script to run D-jets analysis chain

###############
### Prepare:
# config Dzero.h or config Dstar.h
# output root file from the data processing
# output root file for the efficiencies from MC processing
# output root file for the detector responce matrix from MC processing
# B feed-down simulation outputs, with the correct D meson pt
# bkg. response matrix if needed (p-Pb, Pb-Pb)
# D pT bins are hard coded also in efficiency and B feed-down extraction macros
###############

outputdirectory=$1
outputdirectorySignal=$2
lhcprod=$3
efficiencyfile=$4
detRMpromptfile=$5
detRMnonpromptfile=$6

############### set up the output directory name (the whole analysis chain output will be saved there)
outdirBase=$outputdirectory
outdir=$outputdirectorySignal
currDir=`pwd`

############### analysis flags (you can switch some of the flags to run on a part of the analysis chain)
Dmeson=0            #0: D0, 1: D*
isSignal=1
isRefl=1            #reflections for D0
isEffPrompt=1
isEffNonPrompt=1
isDetRMPrompt=1
isDetRMNonPrompt=1
isUnfolding=1
isFinalSpectra=1

isFDSub=1
isFDSys=${21}
isFDUpSys=${22}
isFDDownSys=${23}
isBkgRM=1
bkgRMtype=$7         # which bkg RM, 0 is the default one: RandCones_BkgM_Djet5Excl.root
isCsim=0              # switch one (provide below path (simFilesDir) to your raw simulation files, and file names in the config header file), swith this flag on if you haven't got output of the simulation yet, it takes time
isBsim=0            # switch one (provide below path (simFilesDir) to your raw simulation files, and file names in the config header file), swith this flag on if you haven't got output of the simulation yet with a current efficiencies, it takes time
isBsimNoEff=0       # switch one if you want non-prompt simulations without scaling for non-prompt/prompt efficiency

ispostfix=${12}
postfix=${13}
ispostfixFD=${14}
postfixFD=${15}

ptbinning=${16}
jetpttruemin=${17}
jetptmeasmin=${18}

isDefaultAn=${19}
doCutVar=${20}
isRawSpectra=${24}

isOnlyFinal=0
if [ $isOnlyFinal -eq 1 ]; then
  isSignal=0
  isRefl=0            #reflections for D0
  isEffPrompt=0
  isEffNonPrompt=0
  isDetRMPrompt=0
  isDetRMNonPrompt=0
  isUnfolding=0
fi

############### unfolding config
unfType=$8 #0:bayes, 1:SVD
#regPar=4 #regularization parameter for unfolding
regPar=$9 #regularization parameter for unfolding
isPrior=${10}
priorType=${11}

############### bkg fluctuations matrix
case $bkgRMtype in
0)
  bkgRMfileName=RandCones_BkgM_Djet5Excl.root
  ;;
1)
  bkgRMfileName=RandCones_BkgM_Djet5Excl_Dlead.root
  ;;
2)
  bkgRMfileName=RandCones_BkgM_Djet5Perp.root
  ;;
3)
  bkgRMfileName=RandCones_BkgM_Djet5Perp_Dlead.root
  ;;
4)
  bkgRMfileName=RandCones_BkgM_Djet10Excl.root
  ;;
5)
  bkgRMfileName=RandCones_BkgM_Djet10Excl_Dlead.root
  ;;
6)
  bkgRMfileName=RandCones_BkgM_Djet10Perp.root
  ;;
7)
  bkgRMfileName=RandCones_BkgM_Djet10Perp_Dlead.root
  ;;
8)
  bkgRMfileName=RandCones_BkgM_Alljet5Excl.root
  ;;
9)
  bkgRMfileName=RandCones_BkgM_Alljet5Perp.root
  ;;
10)
  bkgRMfileName=RandCones_BkgM_Alljet10Excl.root
  ;;
11)
  bkgRMfileName=RandCones_BkgM_Alljet10Perp.root
  ;;
*)
  echo "!!! Wrong bkg RM typ !!!"
  exit 1
  ;;
esac

if [ $bkgRMtype -gt 0 ]; then
    outdir=$outdir\_BkgM$bkgRMtype
fi

################################################
############### in/out directories config
################################################

############### code dirs
signalDir=DsignalExtraction
effDir=efficiency
FDDir=FDsubtraction
RMDir=ResponseMatrix
SimDir=POWHEGSim
unfoldingDir=unfolding
finalDir=finalSpectra

############### out dirs
signalDirOut=$outdirBase/$outdir/signalExtraction
effDirOut=$outdirBase/$outdir/efficiency
FDSubDirOut=$outdirBase/$outdir/FDsubtraction
RMDirOut=$outdirBase/$outdir/ResponseMatrix
unfoldingDirOut=$outdirBase/$outdir/unfolding
if [ $unfType -eq  0 ]; then
  unfoldingDirOut=$unfoldingDirOut\_Bayes\_
fi
if [ $unfType -eq  1 ]; then
  unfoldingDirOut=$unfoldingDirOut\_SVD\_
fi
unfoldingDirOut=$unfoldingDirOut$regPar
if [ $isPrior -eq 1 ]; then
  unfoldingDirOut=$unfoldingDirOut\_priorType$priorType
fi

if [ $isFDUpSys -eq 1 ]; then
  unfoldingDirOut=$unfoldingDirOut\_FDsysUp
fi
if [ $isFDDownSys -eq 1 ]; then
  unfoldingDirOut=$unfoldingDirOut\_FDsysDown
fi

finalDirOut=$unfoldingDirOut/finalSpectra
BFDSimDirOut=$outdirBase/Simulations/BFeedDown
PromptSimDirOut=$outdirBase/Simulations/Prompt


############### D-jet signal config
anaoutfiledir=$HOME/Work/alice/analysis/pPb_run2/D0jet/outData
analysisfile=$anaoutfiledir/AnalysisResults
isMoreFiles=0
production=kl #if more than one file to analyse
analysisDataFile=$anaoutfiledir/AnalysisResults$lhcprod.root
reflfile=$HOME/Work/alice/analysis/pPb_run2/D0jet/outMC/reflections/reflections_fitted_DoubleGaus.root

############### efficiency config
MCoutfiledir=$HOME/Work/alice/analysis/pPb_run2/D0jet/outMC
#efffile=$MCoutfiledir/AnalysisResults_fast_D0MCHijing_SMQcorr2.root
efffile=$MCoutfiledir/$efficiencyfile
jetpteffmin=5
jetpteffmax=50
recoPt=0   #if numinator is cut on gen or reco jet pT

############### response matrices config
bkgRMDir=$currDir/ResponseMatrix/BkgRM03
bkgRMFile=$bkgRMDir/$bkgRMfileName
#detRMPrompt=$MCoutfiledir/AnalysisResults_fast_D0MCPythia_SMQcorr2.root
#detRMNonPrompt=$MCoutfiledir/AnalysisResults_fast_D0MCPythia_SMQcorr2.root
detRMPrompt=$MCoutfiledir/$detRMpromptfile
detRMNonPrompt=$MCoutfiledir/$detRMnonpromptfile

############### POWHEG simulations
simFilesDir=/media/basia/Disk2/Work/Djets/POWHEGSimulations/fastSim_pPb5TeV
#simFilesDir=/home/basia/Work/alice/analysis/fastSim_pPb5TeV/files #old simulations
nSimFilesB=9           # have to correspond to what defined in the config file -1
nSimFilesC=8            # have to correspond to what defined in the config file -1


if [ $isRawSpectra -eq 1 ]; then
  isEffPrompt=0
fi

############### configure pT bins, etc., in the corresponding Dstar/Dzero config files
if [ $Dmeson -eq 0 ]; then
  cat configDzero.h > config.h
else
  cat configDstar.h > config.h
fi


#outdirectory=$outdir
mkdir -p $outdirBase/$outdir
echo >> config.h
echo "TString OUTDIRECTORY=\"$outdirBase\";"  >> config.h

if [ $isDefaultAn -eq 1 ]; then
  echo "const int ND = 4;"  >> config.h
  echo "const int NDMC = 3;"  >> config.h
else
  echo "const int ND = 2;"  >> config.h
  echo "const int NDMC = 2;"  >> config.h
fi

# pp binning
if [ $ptbinning -eq 0 ]; then

  case $jetpttruemin in
  3)
    echo "const int fptbinsJetTrueN = 9;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 3,4,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
    #echo "const int fptbinsJetTrueN = 8;" >> config.h
    #echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 3,4,5,6,8,10,14,20,30 };" >> config.h
    #;;
  4)
    echo "const int fptbinsJetTrueN = 8;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 4,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
  5)
    echo "const int fptbinsJetTrueN = 7;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 5,6,8,10,14,20,30,50 };" >> config.h
    ;;
    #echo "const int fptbinsJetTrueN = 6;" >> config.h
    #echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 5,6,8,10,14,20,30 };" >> config.h
    #;;
  *)
    echo "!!! Wrong value of min true pT  !!!"
    exit 1
    ;;
  esac

  case $jetptmeasmin in
  3)
    echo "const int fptbinsJetMeasN = 9;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 3,4,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
    #echo "const int fptbinsJetMeasN = 8;" >> config.h
    #echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 3,4,5,6,8,10,14,20,30 };" >> config.h
    #;;
  4)
    echo "const int fptbinsJetMeasN = 8;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 4,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
  5)
    echo "const int fptbinsJetMeasN = 7;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 5,6,8,10,14,20,30,50 };" >> config.h
    ;;
    #echo "const int fptbinsJetMeasN = 6;" >> config.h
    #echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 5,6,8,10,14,20,30 };" >> config.h
    #;;
  *)
    echo "!!! Wrong value of min true pT  !!!"
    exit 1
    ;;
  esac

fi

# pp binning 2
if [ $ptbinning -eq 8 ]; then

  case $jetpttruemin in
  3)
    echo "const int fptbinsJetTrueN = 8;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 3,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
    #echo "const int fptbinsJetTrueN = 8;" >> config.h
    #echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 3,4,5,6,8,10,14,20,30 };" >> config.h
    #;;
  4)
    echo "const int fptbinsJetTrueN = 8;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 4,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
  5)
    echo "const int fptbinsJetTrueN = 7;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 5,6,8,10,14,20,30,50 };" >> config.h
    ;;
    #echo "const int fptbinsJetTrueN = 6;" >> config.h
    #echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 5,6,8,10,14,20,30 };" >> config.h
    #;;
  *)
    echo "!!! Wrong value of min true pT  !!!"
    exit 1
    ;;
  esac

  case $jetptmeasmin in
  3)
    echo "const int fptbinsJetMeasN = 8;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 3,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
    #echo "const int fptbinsJetMeasN = 8;" >> config.h
    #echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 3,4,5,6,8,10,14,20,30 };" >> config.h
    #;;
  4)
    echo "const int fptbinsJetMeasN = 8;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 4,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
  5)
    echo "const int fptbinsJetMeasN = 7;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 5,6,8,10,14,20,30,50 };" >> config.h
    ;;
    #echo "const int fptbinsJetMeasN = 6;" >> config.h
    #echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 5,6,8,10,14,20,30 };" >> config.h
    #;;
  *)
    echo "!!! Wrong value of min true pT  !!!"
    exit 1
    ;;
  esac

fi


if [ $ptbinning -eq 5 ]; then

  case $jetpttruemin in
  3)
    echo "const int fptbinsJetTrueN = 16;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 3,4,5,6,8,10,12,14,16,18,20,22,24,25,30,40,50 };" >> config.h
    ;;
  4)
    echo "const int fptbinsJetTrueN = 8;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 4,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
  5)
    echo "const int fptbinsJetTrueN = 7;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 5,6,8,10,14,20,30,50 };" >> config.h
    ;;

  *)
    echo "!!! Wrong value of min true pT  !!!"
    exit 1
    ;;
  esac

  case $jetptmeasmin in
  3)
    echo "const int fptbinsJetMeasN = 16;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 3,4,5,6,8,10,12,14,16,18,20,22,24,25,30,40,50 };" >> config.h
    ;;
  4)
    echo "const int fptbinsJetMeasN = 8;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 4,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
  5)
    echo "const int fptbinsJetMeasN = 7;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 5,6,8,10,14,20,30,50 };" >> config.h
    ;;
  *)
    echo "!!! Wrong value of min true pT  !!!"
    exit 1
    ;;
  esac

fi
# Pb-Pb binning
if [ $ptbinning -eq 2 ]; then

  case $jetpttruemin in
  3)
    echo "const int fptbinsJetTrueN = 7;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 3,5,10,15,20,25,35,50 };" >> config.h
    ;;
  5)
    echo "const int fptbinsJetTrueN = 6;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 5,10,15,20,25,35,50 };" >> config.h
    ;;

  *)
    echo "!!! Wrong value of min true pT  !!!"
    exit 1
    ;;
  esac

  case $jetptmeasmin in
  3)
    echo "const int fptbinsJetMeasN = 7;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 3,5,10,15,20,25,35,50 };" >> config.h
    ;;
  5)
    echo "const int fptbinsJetMeasN = 6;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] ={ 5,10,15,20,25,35,50 };" >> config.h
    ;;
  *)
    echo "!!! Wrong value of min true pT  !!!"
    exit 1
    ;;
  esac

fi

if [ $ptbinning -eq 1 ]; then

  case $jetpttruemin in
  3)
    echo "const int fptbinsJetTrueN = 11;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 3,4,5,6,8,10,12,16,20,24,30,50 };" >> config.h
    ;;
  4)
    echo "const int fptbinsJetTrueN = 10;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 4,5,6,8,10,12,16,20,24,30,50 };" >> config.h
    ;;
  5)
    echo "const int fptbinsJetTrueN = 9;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 5,6,8,10,12,16,20,24,30,50 };" >> config.h
    ;;

  *)
    echo "!!! Wrong value of min true pT  !!!"
    exit 1
    ;;
  esac

  case $jetptmeasmin in
  3)
    echo "const int fptbinsJetMeasN = 11;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 3,4,5,6,8,10,12,16,20,24,30,50 };" >> config.h
    ;;
  4)
    echo "const int fptbinsJetMeasN = 10;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 4,5,6,8,10,12,16,20,24,30,50 };" >> config.h
    ;;
  5)
    echo "const int fptbinsJetMeasN = 9;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] ={ 5,6,8,10,12,16,20,24,30,50 };" >> config.h
    ;;
  *)
    echo "!!! Wrong value of min true pT  !!!"
    exit 1
    ;;
  esac

fi
#binning 2
if [ $ptbinning -eq 3 ]; then

 case $jetpttruemin in
 3)
   echo "const int fptbinsJetTrueN = 7;" >> config.h
   echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 3,5,8,12,16,22,28,40 };" >> config.h
   ;;
 5)
   echo "const int fptbinsJetTrueN = 6;" >> config.h
   echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 5,8,12,16,22,28,40 };" >> config.h
   ;;

 *)
   echo "!!! Wrong value of min true pT  !!!"
   exit 1
   ;;
 esac

 case $jetptmeasmin in
 3)
   echo "const int fptbinsJetMeasN = 7;" >> config.h
   echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 3,5,8,12,16,22,28,40 };" >> config.h
   ;;
 5)
   echo "const int fptbinsJetMeasN = 6;" >> config.h
   echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] ={ 5,8,12,16,22,28,40 };" >> config.h
   ;;
 *)
   echo "!!! Wrong value of min true pT  !!!"
   exit 1
   ;;
 esac

fi
#binning 3
if [ $ptbinning -eq 4 ]; then

 case $jetpttruemin in
 3)
   echo "const int fptbinsJetTrueN = 10;" >> config.h
   echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 3,4,6,9,11,14,17,23,28,35,50 };" >> config.h
   ;;
 5)
   echo "const int fptbinsJetTrueN = 6;" >> config.h
   echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 5,8,12,16,22,28,40 };" >> config.h
   ;;

 *)
   echo "!!! Wrong value of min true pT  !!!"
   exit 1
   ;;
 esac

 case $jetptmeasmin in
 3)
   echo "const int fptbinsJetMeasN = 10;" >> config.h
   echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 3,4,6,9,11,14,17,23,28,35,50 };" >> config.h
   ;;
 5)
   echo "const int fptbinsJetMeasN = 6;" >> config.h
   echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] ={ 5,8,12,16,22,28,40 };" >> config.h
   ;;
 *)
   echo "!!! Wrong value of min true pT  !!!"
   exit 1
   ;;
 esac

fi

#binning 4
if [ $ptbinning -eq 6 ]; then

 case $jetpttruemin in
 3)
   echo "const int fptbinsJetTrueN = 8;" >> config.h
   echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 3,5,7,10,13,18,25,32,50 };" >> config.h
   ;;
 5)
   echo "const int fptbinsJetTrueN = 6;" >> config.h
   echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 5,8,12,16,22,28,40 };" >> config.h
   ;;

 *)
   echo "!!! Wrong value of min true pT  !!!"
   exit 1
   ;;
 esac

 case $jetptmeasmin in
 3)
   echo "const int fptbinsJetMeasN = 8;" >> config.h
   echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 3,5,7,10,13,18,25,32,50 };" >> config.h
   ;;
 5)
   echo "const int fptbinsJetMeasN = 6;" >> config.h
   echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] ={ 5,8,12,16,22,28,40 };" >> config.h
   ;;
 *)
   echo "!!! Wrong value of min true pT  !!!"
   exit 1
   ;;
 esac

fi

if [ $ptbinning -eq 7 ]; then

  case $jetpttruemin in
  3)
    echo "const int fptbinsJetTrueN = 8;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 3,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
  4)
    echo "const int fptbinsJetTrueN = 8;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 4,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
  5)
    echo "const int fptbinsJetTrueN = 7;" >> config.h
    echo "double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 5,6,8,10,14,20,30,50 };" >> config.h
    ;;

  *)
    echo "!!! Wrong value of min true pT  !!!"
    exit 1
    ;;
  esac

  case $jetptmeasmin in
  3)
    echo "const int fptbinsJetMeasN = 8;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 3,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
  4)
    echo "const int fptbinsJetMeasN = 8;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 4,5,6,8,10,14,20,30,50 };" >> config.h
    ;;
  5)
    echo "const int fptbinsJetMeasN = 7;" >> config.h
    echo "double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 5,6,8,10,14,20,30,50 };" >> config.h
    ;;
  *)
    echo "!!! Wrong value of min true pT  !!!"
    exit 1
    ;;
  esac

fi



################################################
############### D-jet efficiency
################################################
cd $effDir
if [ $isEffPrompt -eq 1 ]; then
    aliroot -l -b -q DjetEfficiency.C'(1,"'$efffile'","'$effDirOut'",'$jetpteffmin','$jetpteffmax','$recoPt','$ispostfix',"'$postfix'")'
fi

if [ $isEffNonPrompt -eq 1 ]; then
    aliroot -l -b -q DjetEfficiency.C'(0,"'$efffile'","'$effDirOut'",'$jetpteffmin','$jetpteffmax','$recoPt','$ispostfixFD',"'$postfixFD'")'
fi

effFilePrompt=$effDirOut/DjetEff_prompt_jetpt$jetpteffmin\_$jetpteffmax.root
effFileNonPrompt=$effDirOut/DjetEff_nonPrompt_jetpt$jetpteffmin\_$jetpteffmax.root

if [[ $isEffPrompt -eq 1 && $isEffNonPrompt -eq 1 ]]; then
#plot prompt and non-prompt D-jet efficiencies, last two arguments: plot min and plot max
  aliroot -l -b -q drawEff.C'("'$effFilePrompt'","'$effFileNonPrompt'","'$effDirOut'",'$jetpteffmin','$jetpteffmax',2,36)'
fi

cd $currDir


################################################
############### D-jet signal
################################################


############### run
if [ $isSignal -eq 1 ]; then
cd $signalDir
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
if [ $isEffPrompt -eq 1 ]; then
  if [ ! -f $effFilePrompt ]; then
    echo "!!! Prompt efficiency file: \"$effFilePrompt\"  does not exist !!!! "
    exit 1
  fi
fi

./getYields.csh  $analysisfile $lhcprod $isMoreFiles $production $isEffPrompt $effFilePrompt $isRefl $reflfile $ispostfix $postfix $signalDirOut

if [ $isEffPrompt -eq 0 ]; then
  signalFile=$signalDirOut/JetPtSpectra_SB_noEff.root
else
  signalFile=$signalDirOut/JetPtSpectra_SB_eff.root
fi

cd $currDir
fi

if [ ! -f $signalFile ]; then
  echo "!!! D-jet signal extraction file: \"$signalFile\"  does not exist !!!! "
  exit 1
fi


if [ $isRawSpectra -eq 1 ]; then
  exit 1
fi

################################################
############### Detector Response Matrix
################################################
cd $RMDir

if [ $isDetRMPrompt -eq 1 ]; then
    aliroot -l -b -q DetRM.C'(1,"'$detRMPrompt'","'$RMDirOut'",'$ispostfix',"'$postfix'")'
fi

if [ $isDetRMNonPrompt -eq 1 ]; then
  aliroot -l -b -q DetRM.C'(0,"'$detRMNonPrompt'","'$RMDirOut'",'$ispostfixFD',"'$postfixFD'")'
fi

detRMFilePrompt=$RMDirOut/DetMatrix_prompt.root
detRMFileNonPrompt=$RMDirOut/DetMatrix_nonPrompt.root

cd $currDir

################################################
############### B feed-down simulations
################################################

if [ ! -d "$BFDSimDirOut" ] || [ $isBsim -eq 1 ]; then

cd $SimDir
# 1: number of simulation files (as defined in the config file)
# 2: simulation file directory;
# 3: quark type: 0: charm, 1: beauty
# 4: D-jet pT spectrum: 1, D-meson pT spectrum: 0
# 5: if D meson pT cut applied (for D-jet pT spectrum case), if yes the lower and upper values from the D-meson pT bins from the config file are taken
# 6: if efficiency applied (for B simulations, ratio of non-prompt/prompt efficiency)
# 7: prompt efficiency file
# 8: non-prompt efficiency file
# 9: direcotry for the output files

 ./doGetSimOut.csh $nSimFilesB $simFilesDir 1 1 1 1 $effFilePrompt $effFileNonPrompt $BFDSimDirOut
#plot all the variations
  aliroot -l -b -q plotSimSpectra.C'(1,1,1,1,"'$BFDSimDirOut'","'$BFDSimDirOut'")'

cd $currDir
fi

if [ $isBsimNoEff -eq 1 ]; then
cd $SimDir

  #no eff scaling
  if [ $isBsimNoEff -eq 1 ]; then
    ./doGetSimOut.csh $nSimFilesB $simFilesDir 1 1 1 0 $effFilePrompt $effFileNonPrompt $BFDSimDirOut
    aliroot -l -b -q plotSimSpectra.C'(1,1,1,0,"'$BFDSimDirOut'","'$BFDSimDirOut'")'
  fi
cd $currDir
fi


################################################
############### B feed-down subtraction
################################################

if [ $isFDSub -eq 1 ]; then
  if [ ! -f $detRMFileNonPrompt ]; then
    echo "!!! Non-prompt det RM file: \"$detRMFileNonPrompt\"  does not exist !!!! "
    exit 1
  fi

  cd $RMDir
  aliroot -l -b -q combineRM.C'(0,"'$RMDirOut'","'$detRMFileNonPrompt'",'$isBkgRM',"'$bkgRMFile'")'
  combRMFDFile=$RMDirOut/combineMatrixFD.root
  if [ ! -f $combRMFDFile ]; then
    echo "!!! Non-prompt combined RM file: \"$combRMFDFile\"  does not exist !!!! "
    exit 1
  fi
  cd $currDir

  cd $FDDir
  aliroot -l -b -q subtractFD.C'("'$signalFile'","'$analysisDataFile'","'$BFDSimDirOut'","'$combRMFDFile'","'$FDSubDirOut'")'

  cd $currDir
fi

if [ $isFDSub -eq 1 ]; then
  signalBCorrFile=$FDSubDirOut/JetPtSpectrum_FDsub.root
fi
if [ $isFDSub -eq 0 ]; then
  signalBCorrFile=$signalFile
fi


################################################
############### Unfolding
################################################
if [ $isUnfolding -eq 1 ]; then

  if [ ! -f $detRMFilePrompt ]; then
    echo "!!! Non-prompt det RM file: \"$detRMFilePrompt\"  does not exist !!!! "
    exit 1
  fi

  cd $RMDir
    aliroot -l -b -q combineRM.C'(1,"'$RMDirOut'","'$detRMFilePrompt'",'$isBkgRM',"'$bkgRMFile'")'
    combRMFile=$RMDirOut/combineMatrix.root
  cd $currDir

  cd $unfoldingDir
      if [ $unfType -eq  0 ]; then
          aliroot -l -b -q unfold_Bayes.C'("'$signalBCorrFile'","'$detRMFilePrompt'","'$bkgRMFile'","'$unfoldingDirOut'",'$regPar','$isPrior','$priorType','$isBkgRM','$isFDUpSys','$isFDDownSys')'
      fi
      if [ $unfType -eq  1 ]; then
          aliroot -l -b -q unfold_SVD.C'("'$signalBCorrFile'","'$detRMFilePrompt'","'$bkgRMFile'","'$unfoldingDirOut'",'$regPar','$isPrior','$priorType','$isBkgRM','$isFDUpSys','$isFDDownSys')'
      fi


  cd $currDir

fi

cd $currDir

if [ $isUnfolding -eq 1 ]; then
  signalUnfoldedFile=$unfoldingDirOut/unfoldedSpectrum_unfoldedJetSpectrum.root
fi
if [ $isUnfolding -eq 0 ]; then
  signalUnfoldedFile=$signalBCorrFile
fi


################################################
############### prompt D-jet simulations
################################################

if [ ! -d "$PromptSimDirOut" ] || [ $isCsim -eq 1 ]; then

cd $SimDir
# 1: number of simulation files (as defined in the config file)
# 2: simulation file directory;
# 3: quark type: 0: charm, 1: beauty
# 4: D-jet pT spectrum: 1, D-meson pT spectrum: 0
# 5: if D meson pT cut applied (for D-jet pT spectrum case), if yes the lower and upper values from the D-meson pT bins from the config file are taken
# 6: if efficiency applied (for B simulations, ratio of non-prompt/prompt efficiency)
# 7: prompt efficiency file
# 8: non-prompt efficiency file
# 9: direcotry for the output files
  ./doGetSimOut.csh $nSimFilesC $simFilesDir 0 1 1 0 $effFilePrompt $effFileNonPrompt $PromptSimDirOut
  #plot all the variations
  aliroot -l -b -q plotSimSpectra.C'(0,1,1,0,"'$PromptSimDirOut'","'$PromptSimDirOut'")'

cd $currDir
fi

################################################
############### Final spectra
################################################

if [ $isFinalSpectra -eq 1 ]; then
sysDir=$outdirBase/$outdir/systematics
  cd $finalDir
  aliroot -l -b -q finalJetSpectra.C'("'$signalUnfoldedFile'","'$analysisDataFile'","'$PromptSimDirOut'","'$finalDirOut'","'$sysDir'")'

  cd $currDir
fi

################################################
############### END
################################################
