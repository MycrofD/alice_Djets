import ROOT
import array
import rootpy
import scipy as sp
#import matplotlib.pyplot as plt
import numpy as np
from settings import *

from rootpy.io import root_open

def sethistoType(hist, color, msize, mstyle, lwidth):
    hist.SetMarkerColor(color)
    hist.SetLineColor(color)
    hist.SetMarkerSize(msize)
    hist.SetMarkerStyle(mstyle)
    hist.SetLineWidth(lwidth)
    return hist

# Theory predictions
#-------------------
def TheoryPred(THFILE, powpyt6names, xsec):
    #return a list of theory histograms depending on xsec is needed or probability
    hT = []
    for files in range(len(powpyt6names)):
        h = THFILE[files+i*len(powpyt6names)].Get('hPt')
        if xsec == 0:
            h.Scale(1/h.Integral())
        elif xsec == 1:
            h.Scale(simScaling)
            h.Scale(1/dy)
        h.Scale(1,"width")
        hT.append(h)
    return hT

# Theory systematics
#------------------
fptbinsZN = 6
def SysRange(fptbinsZN, powpyt6names, hlist):
    # powpyt6names only provides the length of the list of theory simulations
    # hlist is the list of histograms whose systematics will be provided
    # returns list of up an low deviations
    sysup, syscent, sysdo, cSim = [], [], [], len(powpyt6names)
    for ibinsys in range(fptbinsZN):
        minval = maxval = centval = hlist[0].GetBinContent(ibinsys+1)
        for kCsim in range(1,cSim):
            maxval = max(maxval,hlist[kCsim].GetBinContent(ibinsys+1))
            minval = min(minval,hlist[kCsim].GetBinContent(ibinsys+1))
        sysup.append(maxval - centval)
        syscent.append(centval)
        sysdo.append(centval - minval)
    return [sysup, syscent, sysdo]

# Graph settings
#------------------
def TGRAPH(bins, histcenters, syscent, exlh, sysdo, sysup, graphcolor, fillstyle, histedges, histsetting, ytitle, xtitle):
    #      bins, histXcenters, histYcenters, xerr, xerr, yerrdown, yerrup, color, fillstyle, edges of hist, base histogram
    graph = ROOT.TGraphAsymmErrors(bins, array.array('d',histcenters), array.array('d',syscent), array.array('d',exlh), array.array('d',exlh), array.array('d',sysdo), array.array('d',sysup))
    graph.SetFillColor(graphcolor)
    graph.SetFillStyle(fillstyle)
    graphhist = ROOT.TH1F("","",fptbinsZN,array.array('d',histedges))
    graphhist.SetMinimum(histsetting.GetMaximum()*(-0.04))
    graphhist.SetMaximum(histsetting.GetMaximum()*1.2)

    graphhist.GetYaxis().SetTitle(ytitle)
    graphhist.GetXaxis().SetTitle(xtitle)
#    graphhist.SetTitleSize(0.0,"XYZ")
#    graphhist.SetTitleOffset(1.6,"XYZ")
    graphhist.SetTitle("")
    graphhist.SetMinimum(0)

    graph.SetHistogram(graphhist)
    return [graph, graphhist]

# Colors
#----------------
colors2use = [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kYellow+2, ROOT.kOrange+2, ROOT.kBlack, ROOT.kBlue+10, ROOT.kBlue-2, ROOT.kRed-2, ROOT.kGreen-2]

ROOT.gStyle.SetOptStat(000)
# unfolded data
r3file2D = root_open('/media/jackbauer/data/z_out/R_03_35/unfolding/Bayes/alljetz2D/unfold2DoutFile.root','read')
r4file2D = root_open('/media/jackbauer/data/z_out/R_04_35/unfolding/Bayes/alljetz2D/unfold2DoutFile.root','read')
r6file2D = root_open('/media/jackbauer/data/z_out/R_06_35/unfolding/Bayes/alljetz2D/unfold2DoutFile.root','read')

# Data Train
data3=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_437.root","read")
data4=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_503_R04.root","read")
data6=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_504_R06.root","read")

# theory simulations
jetRs = ['3','4','6']
jetnames = ['5_7','7_10','10_15','15_50']
dptnames = ['2_7','3_10','5_15','5_36']
powpyt6names = ['central', 'F1R2_1536594271','F1R05_1536598175','F2R1_1535916012','F2R2_1535895146','F05R1_1536604800','F05R05_1535894261','m13_1536595965','m17_1536655729']

# Scaling settings
simScaling = 0.5
BRDzero = 0.0389
sigma_in = 50.77
dataLum3 = data3.Get("PWG3_D2H_DmesonsForJetCorrelationsMBN0").Get("NormalizationCounter").GetNEventsForNorm()/sigma_in
dataLum4 = data4.Get("PWG3_D2H_DmesonsForJetCorrelationsMBN0").Get("NormalizationCounter").GetNEventsForNorm()/sigma_in
dataLum6 = data6.Get("PWG3_D2H_DmesonsForJetCorrelationsMBN0").Get("NormalizationCounter").GetNEventsForNorm()/sigma_in

    
# Settings for plotting: Which R and (xsec or prob)
#------------------------------------------------
xsec = 1 # if xsec = 0 it is probability, else xsec
R = 3

# List of theory sim files for each radius
Tr3file2D = []
Tr4file2D = []
Tr6file2D = []

for jetR in range(len(jetRs)):
    for jdname in range(len(jetnames)):
        for powpyt6name in powpyt6names:
            filename = ROOT.TFile('/media/jackbauer/data/z_out/R_0'+jetRs[jetR]+'_35/SimFiles/Prompt/Jetpt'+jetnames[jdname]+'/Z_AnalysisResults_FastSim_powheg+pythia6_charm_'+powpyt6name+'_Jetpt'+jetnames[jdname]+'_Dpt'+dptnames[jdname]+'_effScaled_Dzero.root','read')
            if jetR == 0: Tr3file2D.append(filename)
            elif jetR == 1: Tr4file2D.append(filename)
            elif jetR == 2: Tr6file2D.append(filename)

# -------------------
if R == 3:
    dataLum = dataLum3
    RFILE = r3file2D
    THFILE = Tr3file2D
elif R == 4:
    dataLum = dataLum4
    RFILE = r4file2D
    THFILE = Tr4file2D
elif R == 6:
    dataLum = dataLum6
    RFILE = r6file2D
    THFILE = Tr6file2D

datascaling = 1/(2*BRDzero*dataLum)
dy = 2*((9-R)/10.0)




# Unfolded spectra
hU=[]
#jetbins=4
title = ["5-7", "7-10", "10-15", "15-50"]
#c1 = ROOT.TCanvas('c1', '', 200, 10, 2800, 500 )
#c1.Divide(jetbins,1)
c1 = ROOT.TCanvas('c1', '', 200, 10, 800, 800 )
#c1.SetLeftMargin(4.5)

leg3=[]
start = 1
for i in range(jetbins):
    hU.append(RFILE.Get('UnfProjectX_'+str(i+start)))
    sethistoType(hU[i],ROOT.kBlue+2,1.7,21,4)
jetrange = 1
# now for loop is inactive, "iterates" over one value only
for i in range(jetrange-1,jetrange):
    pad1 = ROOT.TPad("","",0,0.3,1,1.0)
    pad1.SetBottomMargin(0.01)
    pad1.SetRightMargin(1.1)
    pad1.SetLeftMargin(0.13)
    pad1.Draw()
    pad1.cd()

    hU[i].SetStats(0)
    
#    hU[i].SetTitle("Jet pt: "+title[i]+" GeV")
#    hU[i].GetXaxis().SetTitle("z")
    if xsec == 0:
        hU[i].Scale(1/hU[i].Integral())
        hU[i].Scale(1,"width")
        hU[i].GetYaxis().SetRangeUser(-0.4,5)
        hU[i].GetYaxis().SetTitle("Probability")
    elif xsec == 1:
        hU[i].Scale(datascaling)
        hU[i].Scale(1/dy)
        hU[i].Scale(1,"width")
        hU[i].GetYaxis().SetRangeUser(-0.1*hU[i].GetMaximum(),hU[i].GetMaximum()*3)
        hU[i].GetYaxis().SetTitle("#frac{d^{3}#sigma}{dz d p_{T} d#eta}")
    hU[i].SetTitleOffset(1.2,"XYZ")
    hU[i].Draw()

    # pad1 settings contd...
    hU[i].GetYaxis().SetLabelSize(0.03)
    axis = ROOT.TGaxis(-5,20,-5,220,20,220,510,"")
    axis.SetLabelFont(43)
    axis.SetLabelSize(15)
    axis.Draw()

    # Thoery predictions
    hT = TheoryPred(THFILE, powpyt6names, xsec)

    ##----------
    [sysup, syscent, sysdo] = SysRange(fptbinsZN,powpyt6names,hT)
    zval = [0.45,0.55,0.65,0.75,0.85,0.96]
    zedges = [0.4,0.5,0.6,0.7,0.8,0.9,1.02]
    exlh = [0.05,0.05,0.05,0.05,0.05,0.06]

    ytitle = "#frac{d^{3}#sigma}{dz_{||}^{ch} dp_{T} d#it{#eta}}"
 #   ytitle = "d^{3}#sigma/dz_{||}^{ch} dp_{T} d#it{#eta}"

    [graph, graphhist] = TGRAPH(fptbinsZN, zval, syscent, exlh, sysdo, sysup, colors2use[0], 3004, zedges, hU[i], ytitle, "")
    graph.Draw("a2same")
    hU[i].Draw("same")
    hT[0].Draw("same")

    ##----------
    # Lower plot settings
    c1.cd()
    pad2 = ROOT.TPad("","",0,0.05,1,0.3)
    pad2.SetTopMargin(0.015)
    pad2.SetBottomMargin(0.4)
    pad2.SetRightMargin(1.1)
    pad2.SetLeftMargin(0.13)
    pad2.Draw()
    pad2.cd()

    mg = ROOT.TMultiGraph()

    hTR = []
    for item in range(len(hT)):
        hratio = hT[item].Clone("hratio")
        hratio.Divide(hU[i])
        hTR.append(hratio)

    ##----------
    hURatio = hU[i].Clone("hURatio")
    hURatio.Divide(hU[i])
    hURatio.Draw()
    [sysupR, syscentR, sysdoR] = SysRange(fptbinsZN,powpyt6names,hTR)
    #xtitle = "#it{z}_{||}^{ch} = #frac{#vec{p}_{ch jet}.#vec{p}_{D}}{#vec{p}_{ch jet}.#vec{p}_{ch jet}}"
    xtitle = "#it{z}_{||}^{ch}"
    [graphR, graphhistR] = TGRAPH(fptbinsZN, zval, syscentR, exlh, sysdoR, sysupR, colors2use[0], 3004, zedges, hURatio, "Theory/data", xtitle)
    graphhistR.SetTitleSize(0.1,"X")
    graphhistR.SetTitleOffset(1.55,"X")
    graphhistR.SetMinimum(-0.5)
    graphhistR.SetMaximum(5)

    graphhistR.SetTitleSize(0.1,"Y")
    graphhistR.GetYaxis().SetTitle("Theory/Data")
    graphhistR.SetTitleOffset(0.3,"Y")
    graphhistR.GetYaxis().SetLabelFont(43)
    graphhistR.GetYaxis().SetLabelSize(20)

    graphhistR.GetXaxis().SetLabelSize(20)
    graphhistR.GetXaxis().SetLabelFont(43)

    #[graphR2, graphhistR2] = TGRAPH(fptbinsZN, zval, syscentR, exlh, sysdoR, sysupR, colors2use[1], 3005, zedges, hURatio, "Theory/data", xtitle)
    #graphR.Draw("a2same")
    mg.Add(graphR)
    #mg.Add(graphR2)
    mg.Draw("AC")
    #hURatio.Draw("same")
    #hTR[0].Draw("same")

ROOT.gStyle.SetOptStat(000)
c1.Draw()

c1.SaveAs("dataVtheory.png")

print(dy)
print(datascaling)
print(simScaling)
