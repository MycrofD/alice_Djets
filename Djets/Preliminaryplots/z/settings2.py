# this is settings 2

import ROOT
import array
import numpy as np

# Settings for plotting: Which R and (xsec or prob)
#------------------------------------------------
xsec = 0 # if xsec = 0 it is probability, else xsec
R = 3
#jetbins = 4 # actually 5
jetbin = 4 # 0,1,2,3,4 don't use 0
if xsec == 0:
    ytitletop = "Probability"
elif xsec == 1:
    ytitletop="#frac{d^{3}#sigma}{dz_{||}^{ch} dp_{T} d#it{#eta}}"

jetbintitle = "R = 0."+str(R)
if jetbin == 1:
    jetbintitle += ", 5<Jet p_{T}< 7 GeV/#it{c}"
elif jetbin == 2:
    jetbintitle += ", 7<Jet p_{T}<10 GeV/#it{c}"
elif jetbin == 3:
    jetbintitle += ", 10<Jet p_{T}<15 GeV/#it{c}"
elif jetbin == 4:
    jetbintitle += ", 15<Jet p_{T}<50 GeV/#it{c}"

yvalues2 = [0,0,0,0,0,5]
#if xsec == 0:
#    yvalues2 = [0,0,0,0,0,5]
#elif xsec == 1:
#    if jetbin == 1: a = 0.1
#    elif jetbin == 2: a = 0.04
#    elif jetbin == 3: a = 0.02
#    elif jetbin == 4: a = 0.01
#    yvalues2 = [0,0,0,0,0,a]
#############
def sethistoType(hist, color, msize, mstyle, lwidth):
    hist.SetMarkerColor(color)
    hist.SetLineColor(color)
    hist.SetMarkerSize(msize)
    hist.SetMarkerStyle(mstyle)
    hist.SetLineWidth(lwidth)
    return hist

def TGRAPH(n,xvalues1,yvalues1,xerr1,yerr1do,yerr1up,xedges,color):
    graph = ROOT.TGraphAsymmErrors(n,array.array('d',xvalues1),array.array('d',yvalues1),array.array('d',xerr1),array.array('d',xerr1),array.array('d',yerr1do),array.array('d',yerr1up))
    graph.SetLineColor(color)
    graph.SetLineWidth(2) #1504
    graph.SetFillStyle(3005)
    return graph

def hist2list(hist):
    ybins = []
    for i in range(hist.GetNbinsX()):
        ybins.append(hist.GetBinContent(i+1))
    return ybins

def histRatio(hT,hUhist):
    hTR = []
    for item in range(len(hT)):
        hratio = hT[item].Clone("hratio")
        hratio.Divide(hUhist)
        hTR.append(hratio)

    return hTR
#-------------------
def TheoryPred(THFILES, powpyt6names, xsec,jetbin):
    #return a list of theory histograms depending on xsec is needed or probability
    hT = []
    for files in range(len(powpyt6names)):
        h = THFILES[files+(jetbin-1)*len(powpyt6names)].Get('hPt')
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
    return [sysdo, syscent, sysup]

#plot settings
# Colors
#----------------
colors = [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kYellow+2, ROOT.kOrange+2, ROOT.kBlack, ROOT.kBlue+10, ROOT.kBlue-2, ROOT.kRed-2, ROOT.kGreen-2]


ytitletopsize=0.05
ytitlebot="Theory/Data"
xtitle = "#it{z}_{||}^{ch}"
numbincenters=6

xvalues1 = [0.45,0.55,0.65,0.75,0.85,0.96]
yvalues1 = []
yerr1, yerr2, xerr1= [],[],[]
xedges=[0.4,0.5,0.6,0.7,0.8,0.9,1.02]

xerr1=[0.05,0.05,0.05,0.05,0.05,0.06]
for i in range(numbincenters):
    yvalues1.append(10*np.sin(xvalues1[i]))
    yerr1.append(1.5)
    yerr2.append(1.5)



# ROOT FILES, datafiles
# unfolded data
r3file2D = ROOT.TFile('/media/jackbauer/data/z_out/R_03_35/unfolding/Bayes/alljetz2D/unfold2DoutFile.root','read')
r4file2D = ROOT.TFile('/media/jackbauer/data/z_out/R_04_35/unfolding/Bayes/alljetz2D/unfold2DoutFile.root','read')
r6file2D = ROOT.TFile('/media/jackbauer/data/z_out/R_06_35/unfolding/Bayes/alljetz2D/unfold2DoutFile.root','read')

# Data Train
data3=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_437.root","read")
data4=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_503_R04.root","read")
data6=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_504_R06.root","read")

# theory simulations
jetRs = ['3','4','6']
jetnames = ['5_7','7_10','10_15','15_50']
dptnames = ['2_7','3_10','5_15','5_36']
powpyt6names = ['central', 'F1R2_1536594271','F1R05_1536598175','F2R1_1535916012','F2R2_1535895146','F05R1_1536604800','F05R05_1535894261','m13_1536595965','m17_1536655729']
pyt6names = ['1552224706']
pyt8names = ['1552144258']
powpytDIJETnames = ['central', 'F1R2_1536594271','F1R05_1536598175','F2R1_1535916012','F2R2_1535895146','F05R1_1536604800','F05R05_1535894261','m13_1536595965','m17_1536655729']

# Scaling settings
simScaling = 0.5
BRDzero = 0.0389
sigma_in = 50.77
dataLum3 = data3.Get("PWG3_D2H_DmesonsForJetCorrelationsMBN0").Get("NormalizationCounter").GetNEventsForNorm()/sigma_in
dataLum4 = data4.Get("PWG3_D2H_DmesonsForJetCorrelationsMBN0").Get("NormalizationCounter").GetNEventsForNorm()/sigma_in
dataLum6 = data6.Get("PWG3_D2H_DmesonsForJetCorrelationsMBN0").Get("NormalizationCounter").GetNEventsForNorm()/sigma_in

# List of theory sim files for each radius
Tr3file2D, Tr3filePYT6, Tr3filePYT8 = [], [], []
Tr4file2D, Tr4filePYT6, Tr4filePYT8 = [], [], []
Tr6file2D, Tr6filePYT6, Tr6filePYT8 = [], [], []

#Pow+Pyt6
for jetR in range(len(jetRs)):
    for jdname in range(len(jetnames)):
        for powpyt6name in powpyt6names:
            filename = ROOT.TFile('/media/jackbauer/data/z_out/R_0'+jetRs[jetR]+'_35/SimFiles/Prompt/Jetpt'+jetnames[jdname]+'/Z_AnalysisResults_FastSim_powheg+pythia6_charm_'+powpyt6name+'_Jetpt'+jetnames[jdname]+'_Dpt'+dptnames[jdname]+'_effScaled_Dzero.root','read')
            if jetR == 0: Tr3file2D.append(filename)
            elif jetR == 1: Tr4file2D.append(filename)
            elif jetR == 2: Tr6file2D.append(filename)

#Pyt6
for jetR in range(len(jetRs)):
    for jdname in range(len(jetnames)):
        for pyt6name in pyt6names:
            filename = ROOT.TFile('/media/jackbauer/data/z_out/R_0'+jetRs[jetR]+'_35/SimFiles/Prompt/Jetpt'+jetnames[jdname]+'/Z_AnalysisResults_FastSim_pythia6_charm_'+pyt6name+'_Jetpt'+jetnames[jdname]+'_Dpt'+dptnames[jdname]+'_effScaled_Dzero.root','read')
            if jetR == 0: Tr3filePYT6.append(filename)
            elif jetR == 1: Tr4filePYT6.append(filename)
            elif jetR == 2: Tr6filePYT6.append(filename)

#Pyt8
for jetR in range(len(jetRs)):
    for jdname in range(len(jetnames)):
        for pyt8name in pyt8names:
            filename = ROOT.TFile('/media/jackbauer/data/z_out/R_0'+jetRs[jetR]+'_35/SimFiles/Prompt/Jetpt'+jetnames[jdname]+'/Z_AnalysisResults_FastSim_pythia8_charm_'+pyt8name+'_Jetpt'+jetnames[jdname]+'_Dpt'+dptnames[jdname]+'_effScaled_Dzero.root','read')
            if jetR == 0: Tr3filePYT8.append(filename)
            elif jetR == 1: Tr4filePYT8.append(filename)
            elif jetR == 2: Tr6filePYT8.append(filename)

# -------------------
if R == 3:
    dataLum = dataLum3
    RFILE = r3file2D
    THFILES = Tr3file2D
    THFILES2 = Tr3filePYT6
    THFILES3 = Tr3filePYT8
elif R == 4:
    dataLum = dataLum4
    RFILE = r4file2D
    THFILES = Tr4file2D
    THFILES2 = Tr4filePYT6
    THFILES3 = Tr4filePYT8
elif R == 6:
    dataLum = dataLum6
    RFILE = r6file2D
    THFILES = Tr6file2D
    THFILES2 = Tr6filePYT6
    THFILES3 = Tr6filePYT8

datascaling = 1/(2*BRDzero*dataLum)
dy = 2*((9-R)/10.0)

