#!/bin/bash

currDir=`pwd`

ptbinning=$1
jetpttruemin=$2
jetpttruemax=$3
jetptmeasmin=$4
jetptmeasmax=$5

lhcprod=_LHC16R03
efficiencyfile=AnalysisResults_fast_R03_D0MCPythia_default.root
detRMpromptfile=AnalysisResults_fast_R03_D0MCPythia_default.root
detRMnonpromptfile=AnalysisResults_fast_R03_D0MCPythia_default.root


outputdirectorybase=$HOME/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_
#outputdirectorybase=$HOME/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt324PythiaEff_
outputdirectory=${outputdirectorybase}BaseCuts
outputdirectorySignal=Default

if [ $ptbinning -eq 0 ]; then
  outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_ppbinning
#outputdirectorySignal=Refl15_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_ppbinning
#outputdirectorySignal=S2sigma_SB59_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_ppbinning
#outputdirectorySignal=SB58_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_ppbinning
fi
if [ $ptbinning -eq 2 ]; then
outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_PbPbbinning
#outputdirectorySignal=Refl15_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_PbPbbinning
#outputdirectorySignal=S2sigma_SB59_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_PbPbbinning
#outputdirectorySignal=SB58_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_PbPbbinning
fi
if [ $ptbinning -eq 8 ]; then
  outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_ppbinning2
fi
if [ $ptbinning -eq 1 ]; then
outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_finebinning
fi
if [ $ptbinning -eq 3 ]; then
outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_binning2
fi
if [ $ptbinning -eq 4 ]; then
outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_binning3
fi
if [ $ptbinning -eq 6 ]; then
outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_binning4
fi
if [ $ptbinning -eq 7 ]; then
outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_binning5
fi
if [ $ptbinning -eq 5 ]; then
outputdirectorySignal=Default_jetMeas$jetptmeasmin\_$jetptmeasmax\_jetTrue$jetpttruemin\_$jetpttruemax\_finebinning5
fi

ispostfix=0
postfix=Cut
ispostfixFD=1
postfixFD=FD

bkgRMtype=$6
unfType=$7 #bayes
regPar=$8
isPrior=$9
priorType=${10}

# systematics
doBkg=${11}
doPrior=${12}
doBkgPrior=${13}
doCutVar=${14}
doJESSys=${15}
doFDSys=${16}
doSystematics=${17}
doRawCutVar=${18}
doRawSpectra=${19}
doSignal=${20}

unfCutSys=0

################################################
############### Default
################################################

  ################################################
  ############### without efficiency correction
  ################################################
  if [ $doRawSpectra -eq 1 ]; then
 ./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 0 0 0 1

     ################################################
     ############### Cut Variation
     ################################################

 if [ $doCutVar -eq 1 ]; then
 	./run_cut.csh $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectorybase $doRawSpectra
fi

exit 1
fi

  ################################################
  ############### full signal, with efficiency correction
  ################################################

if [ $doSignal -eq 1 ]; then
./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 0 0 0 0
fi

    ################################################
    ############### B FD systematics
    ################################################
if [ $doFDSys -eq 1 ]; then
	./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 1 1 0 0

	./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 1 0 1 0

fi

################################################
############### Cut Variation
################################################

if [ $doCutVar -eq 1 ]; then

	./run_cut.csh $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectorybase $doRawSpectra

  # cut variation systematics
  cd systematics
  #after unfolding
  root -l -b -q cutsSystematics.C'('$regPar',"'$outputdirectorybase'","'$outputdirectorySignal'",'$doRawCutVar',1,1)'
  #before unfolding
  root -l -b -q cutsSystematics.C'('$regPar',"'$outputdirectorybase'","'$outputdirectorySignal'",'$doRawCutVar',0,1)'
  cd $currDir
fi

################################################
############### JES systematics
################################################

if [ $doJESSys -eq 1 ]; then

	./run_JES.csh $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $outputdirectorySignal $outputdirectory

fi

################################################
############### Background variations
################################################

if [ $doBkg -gt 0 ]; then

for bkg in `seq 1 1 11`; do

		bkgRMtype=$bkg
	./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 0 0 0 0

done

fi
bkgRMtype=0

################################################
############### Prior variations
################################################

if [ $doPrior -gt 0 ]; then

for prior in `seq 0 1 8`; do

	isPrior=1
	priorType=$prior
	./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 0 0 0 0

done

fi
isPrior=0
priorType=0

################################################
############### Background and prior variations
################################################

if [ $doBkgPrior -gt 0 ]; then

for bkg in `seq 0 1 11`; do

	bkgRMtype=0
	unfType=0 #bayes
	regPar=3
	isPrior=0
	priorType=0

		bkgRMtype=$bkg
	./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 0 0 0 0

	for prior in `seq 0 1 8`; do
		isPrior=1
		priorType=$prior
		./run.csh $outputdirectory $outputdirectorySignal $lhcprod $efficiencyfile $detRMpromptfile $detRMnonpromptfile $bkgRMtype $unfType $regPar $isPrior $priorType $ispostfix $postfix $ispostfixFD $postfixFD $ptbinning $jetpttruemin $jetptmeasmin 1 0 0 0 0 0
	done

done


fi

################################################
############### Get systematic uncertanties
################################################

if [ $doSystematics -eq 1 ]; then
	cd systematics
		root -l -b -q JetSpectrumSys.C'('$regPar',"'$outputdirectory'","'$outputdirectorySignal'",1)'
	cd $currDir
fi

exit 1
