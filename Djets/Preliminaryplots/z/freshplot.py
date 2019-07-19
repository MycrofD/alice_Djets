## Z X-sec and Prob
import ROOT
import array
import rootpy
import scipy as sp
#import matplotlib.pyplot as plt
import numpy as np
from settings2 import *
from settings3 import *

from rootpy.io import root_open

##settings to be transferred back
########################## CANVAS
c1 = ROOT.TCanvas("c1","",200,10,900,800)
########################## PAD 1
pad1 = ROOT.TPad("","",0,0.3,1,1.0)
pad1.SetMargin(0.13,1,0,1) #left right bottom top
pad1.Draw()
pad1.cd()
mgT = ROOT.TMultiGraph()
mgT.SetTitle("")
#mgT.GetYaxis().SetTitle(ytitletop)
histData.SetTitle("Data")
hT1[0].SetTitle("POWHEG+PYTHIA6")
hT2[0].SetTitle("PYTHIA6")

#histData.Draw()

gr0 = TGRAPH(numbincenters,xvalues1,yvalues2,xerr1,yerr1do,yerr1up,xedges,ROOT.kWhite)
gr0.SetLineWidth(0)
gr0.SetTitle(jetbintitle)
gr1 = TGRAPH(numbincenters,xvalues1,yvalues1,xerr1,yerr1do,yerr1up,xedges,colors[0])
gr1.SetTitle("Syst. Unc. (Theory)")
mgT.Add(gr0)
mgT.Add(gr1)
#histData.GetYaxis().SetRangeUser(0,3)
mgT.SetTitle(";;"+ytitletop)
mgT.GetYaxis().SetTitleSize(ytitletopsize)
mgT.GetYaxis().SetLabelSize(0.04)
mgT.GetXaxis().SetRangeUser(0.4,1.02)
mgT.GetYaxis().SetTitleOffset(1.1)
mgT.Draw("AZ5")
histData.Draw("same")
hT1[0].Draw("same")
if xsec == 0: hT2[0].Draw("same")
pad1.BuildLegend(0.15,0.65,0.5,0.85)
########################## PAD 2
c1.cd()
pad2 = ROOT.TPad("","",0,0.05,1,0.3)
pad2.SetMargin(0.13,1,0.33,0) #left right bottom top
pad2.Draw()
pad2.cd()
mgB = ROOT.TMultiGraph()
mgB.SetTitle("")

gr1b = TGRAPH(numbincenters,xvalues1,yvalues1R,xerr1,yerr1doR,yerr1upR,xedges,colors[0])
mgB.Add(gr1b)
mgB.SetTitle(";"+xtitle+";"+ ytitlebot)
mgB.GetYaxis().SetRangeUser(0,3)
mgB.GetYaxis().SetTitleSize(0.15)
mgB.GetYaxis().SetTitleOffset(0.22)
mgB.GetXaxis().SetTitleSize(0.15)
mgB.GetXaxis().SetTitleOffset(0.82)
mgB.GetYaxis().SetLabelSize(0.1)
mgB.GetXaxis().SetLabelSize(0.1)
mgB.Draw("AZ5")
hTR1[0].Draw("same")
hTR2[0].Draw("same")
########################## BACK to CANVAS
c1.Draw()
c1.SaveAs(str(xsec)+"_R"+str(R)+"_jet"+str(jetbin)+".pdf")
c1.SaveAs(str(xsec)+"_R"+str(R)+"_jet"+str(jetbin)+".png")
#############
#
