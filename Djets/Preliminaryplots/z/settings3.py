# this is settings 3

import ROOT
import array
import numpy as np
from settings2 import *

# data histogram
hU = []
hU.append(RFILE.Get('UnfProjectX_'+str(jetbin)))

histData = hU[0].Clone("")
sethistoType(histData,ROOT.kRed+2,1.7,20,3)

if xsec == 0:
    histData.Scale(1/hU[0].Integral())
    histData.Scale(1,"width")
    histData.GetYaxis().SetRangeUser(-0.4,5)
    histData.GetYaxis().SetTitle("Probability")
elif xsec == 1:
    histData.Scale(datascaling)
    histData.Scale(1/dy)
    histData.Scale(1,"width")
    histData.GetYaxis().SetRangeUser(-0.1*hU[0].GetMaximum(),hU[0].GetMaximum()*3)
    histData.GetYaxis().SetTitle("#frac{d^{3}#sigma}{dz d p_{T} d#eta}")

# Thoery predictions
# powheg+pythia6
hT1 = TheoryPred(THFILES, powpyt6names, xsec, jetbin)
# pythia6
hT2 = TheoryPred(THFILES2, pyt6names, xsec, jetbin)
# pythia8
#hT3 = TheoryPred(THFILES3, pyt8names, xsec, jetbin)
##----------
hTR0 = histRatio(hU, histData)
##----------
[yerr1do, yvalues1, yerr1up] = SysRange(fptbinsZN,powpyt6names,hT1)
hTR1 = histRatio(hT1, histData)
##----------
[yerr1doR, yvalues1R, yerr1upR] = SysRange(fptbinsZN,powpyt6names,hTR1)
hTR2 = histRatio(hT2, histData)
##----------
#sethistoType(hTR0, colors[1], 1.7, 20, 0)

sethistoType(hT1[0], colors[0], 1.7, 4, 0)
sethistoType(hTR1[0], colors[0], 1.7, 4, 0)
sethistoType(hT2[0], ROOT.kGreen+2, 2, 33, 0)
sethistoType(hTR2[0], ROOT.kGreen+2, 2, 33, 0)
val=[]
for i in range(6): val.append(histData.GetBinContent(i+1))
yvalues2[5]=2.2*(max(val[:]))



#savedFile = ROOT.TFile("finalXsection_R"+str()+str(),"recreate")
savedFile = ROOT.TFile("finalXsection_R"+str(R)+str(jetbin)+".root","recreate")
hT1[0].Write()
hT1[0].SetName("central_powPyth6")
histData.Write()
histData.SetName("histData")
savedFile.Close()
