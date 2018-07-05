#!/bin/bash

################################################
############### FULL with sys
################################################
doreg4=0

utype=2 #0:Bayes, 1: SVD
prior1=3
prior2=4

ptbinning=2 			# 0: pp bining 1: more fine binning 2: Pb-Pb binning 5: more fine pp binning

dbkg=0
dprior=0
dcutvar=1
dJESsys=0
dFDsys=0
drawcutvar=1
dsys=0

################################################
############### ptmeas: 3-50, pttrue: 3-50, reg=3
################################################
jetptmeasmin=3
jetptmeasmax=50
jetpttruemin=3
jetpttruemax=50

bkgRMtype=0
unfType=$utype #0: bayes. 1: SVD
regPar=$prior1
isPrior=0
priorType=0

doBkg=$dbkg
doPrior=$dprior
doBkgPrior=0
doCutVar=$dcutvar
doJESSys=$dJESsys
doFDSys=$dFDsys
doSys=$dsys
doRawCutVarSys=$drawcutvar

############### Raw Signal only
doRawSignal=1
doSignal=0
./run_pPbDzeroAnalysis.csh  $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $doBkg $doPrior $doBkgPrior $doCutVar $doJESSys $doFDSys $doSys $doRawCutVarSys $doRawSignal $doSignal

############### Default
doRawSignal=0
doSignal=1
./run_pPbDzeroAnalysis.csh  $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $doBkg $doPrior $doBkgPrior $doCutVar $doJESSys $doFDSys $doSys $doRawCutVarSys $doRawSignal $doSignal

################################################
############### ptmeas: 3-50, pttrue: 5-50, reg=3
################################################
jetptmeasmin=3
jetptmeasmax=50
jetpttruemin=5
jetpttruemax=50

bkgRMtype=0
unfType=$utype
regPar=$prior1
isPrior=0
priorType=0

doBkg=$dbkg
doPrior=$dprior
doBkgPrior=0
doCutVar=$dcutvar
doJESSys=$dJESsys
doFDSys=$dFDsys
doSys=$dsys
doRawCutVarSys=$drawcutvar


############### Raw Signal only
doRawSignal=1
doSignal=0
./run_pPbDzeroAnalysis.csh  $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $doBkg $doPrior $doBkgPrior $doCutVar $doJESSys $doFDSys $doSys $doRawCutVarSys $doRawSignal $doSignal

############### Default
doRawSignal=0
doSignal=1
./run_pPbDzeroAnalysis.csh  $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $doBkg $doPrior $doBkgPrior $doCutVar $doJESSys $doFDSys $doSys $doRawCutVarSys $doRawSignal $doSignal


exit 1
################################################
############### FULL with sys
################################################
doreg4=0

utype=0 #0:Bayes, 1: SVD
prior1=3
prior2=4

ptbinning=2 			# 0: pp bining 1: more fine binning 2: Pb-Pb binning 5: more fine pp binning

dbkg=0
dprior=0
dcutvar=1
dJESsys=0
dFDsys=0
drawcutvar=1
dsys=0

################################################
############### ptmeas: 3-50, pttrue: 3-50, reg=3
################################################
jetptmeasmin=3
jetptmeasmax=50
jetpttruemin=3
jetpttruemax=50

bkgRMtype=0
unfType=$utype #0: bayes. 1: SVD
regPar=$prior1
isPrior=0
priorType=0

doBkg=$dbkg
doPrior=$dprior
doBkgPrior=0
doCutVar=$dcutvar
doJESSys=$dJESsys
doFDSys=$dFDsys
doSys=$dsys
doRawCutVarSys=$drawcutvar

############### Raw Signal only
doRawSignal=1
doSignal=0
./run_pPbDzeroAnalysis.csh  $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $doBkg $doPrior $doBkgPrior $doCutVar $doJESSys $doFDSys $doSys $doRawCutVarSys $doRawSignal $doSignal

############### Default
doRawSignal=0
doSignal=1
./run_pPbDzeroAnalysis.csh  $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $doBkg $doPrior $doBkgPrior $doCutVar $doJESSys $doFDSys $doSys $doRawCutVarSys $doRawSignal $doSignal

################################################
############### ptmeas: 3-50, pttrue: 5-50, reg=3
################################################
jetptmeasmin=3
jetptmeasmax=50
jetpttruemin=5
jetpttruemax=50

bkgRMtype=0
unfType=$utype
regPar=$prior1
isPrior=0
priorType=0

doBkg=$dbkg
doPrior=$dprior
doBkgPrior=0
doCutVar=$dcutvar
doJESSys=$dJESsys
doFDSys=$dFDsys
doSys=$dsys
doRawCutVarSys=$drawcutvar


############### Raw Signal only
doRawSignal=1
doSignal=0
./run_pPbDzeroAnalysis.csh  $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $doBkg $doPrior $doBkgPrior $doCutVar $doJESSys $doFDSys $doSys $doRawCutVarSys $doRawSignal $doSignal

############### Default
doRawSignal=0
doSignal=1
./run_pPbDzeroAnalysis.csh  $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $doBkg $doPrior $doBkgPrior $doCutVar $doJESSys $doFDSys $doSys $doRawCutVarSys $doRawSignal $doSignal


if [ $doreg4 -eq 1 ]; then

################################################
############### ptmeas: 3-50, pttrue: 3-50, reg=4
################################################

jetptmeasmin=3
jetptmeasmax=50
jetpttruemin=3
jetpttruemax=50

bkgRMtype=0
unfType=$utype
regPar=$prior2
isPrior=0
priorType=0

doBkg=$dbkg
doPrior=$dprior
doBkgPrior=0
doCutVar=$dcutvar
doJESSys=$dJESsys
doFDSys=$dFDsys
doSys=$dsys
doRawCutVarSys=$drawcutvar


############### Raw Signal only
doRawSignal=1
doSignal=0
./run_pPbDzeroAnalysis.csh  $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $doBkg $doPrior $doBkgPrior $doCutVar $doJESSys $doFDSys $doSys $doRawCutVarSys $doRawSignal $doSignal

############### Default
doRawSignal=0
doSignal=1
./run_pPbDzeroAnalysis.csh  $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $doBkg $doPrior $doBkgPrior $doCutVar $doJESSys $doFDSys $doSys $doRawCutVarSys $doRawSignal $doSignal


################################################
############### ptmeas: 3-50, pttrue: 5-50, reg=4
################################################
jetptmeasmin=3
jetptmeasmax=50
jetpttruemin=5
jetpttruemax=50

bkgRMtype=0
unfType=$utype
regPar=$prior2
isPrior=0
priorType=0

doBkg=$dbkg
doPrior=$dprior
doBkgPrior=0
doCutVar=$dcutvar
doJESSys=$dJESsys
doFDSys=$dFDsys
doSys=$dsys
doRawCutVarSys=$drawcutvar


############### Raw Signal only
doRawSignal=1
doSignal=0
./run_pPbDzeroAnalysis.csh  $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $doBkg $doPrior $doBkgPrior $doCutVar $doJESSys $doFDSys $doSys $doRawCutVarSys $doRawSignal $doSignal

############### Default
doRawSignal=0
doSignal=1
./run_pPbDzeroAnalysis.csh  $ptbinning $jetpttruemin $jetpttruemax $jetptmeasmin $jetptmeasmax $bkgRMtype $unfType $regPar $isPrior $priorType $doBkg $doPrior $doBkgPrior $doCutVar $doJESSys $doFDSys $doSys $doRawCutVarSys $doRawSignal $doSignal

fi

exit 1
