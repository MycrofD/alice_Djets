import ROOT as RT 
import ROOT 
import array
##import rootpy
import scipy as sp
import numpy as np

RT.gStyle.SetErrorX(0) # to remove the horizontal error bars

#######################################################################1
#########          FUNCTIONS/SETTINGS USED  ###########################
#######################################################################
######## USED IN DATA HISTOS
def sethistoType(hist, color, msize, mstyle, lwidth):
    hist.SetMarkerColor(color);hist.SetLineColor(color);hist.SetMarkerSize(msize);hist.SetMarkerStyle(mstyle);hist.SetLineWidth(lwidth)
    return hist

######### USED IN THEORY HISTOS
def TheoryPred(THFILES, powpyt6names, xsec,jetbin):
    #return a list of theory histograms depending on xsec is needed or probability
    hT = []
    for files in range(len(powpyt6names)):
        h = THFILES[files+(jetbin-1)*len(powpyt6names)].Get('hPt')
        h = h.Rebin(fptbinsZN,'h',array.array('d',fptbinsZlh))
        if xsec == 0:
            h.Scale(1/h.Integral())
        elif xsec == 1:
            h.Scale(simScaling);h.Scale(1/dy)
        h.Scale(1,'width')
        hT.append(h)
    return hT

######### USED IN THEORY HISTOS
def histRatio(hT,hUhist):
    hTR = []
    for item in range(len(hT)):
        hratio = hT[item].Clone('hratio')
        hratio.Divide(hUhist)
        hTR.append(hratio)
    return hTR

######### USED IN THEORY HISTOS SYSTEMATICS
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

######### USED FOR ASYMMETRICAL ERRORS
def TGRAPH(n,xvalues1,yvalues1,xerr1,yerr1do,yerr1up,xedges,color):
    graph = ROOT.TGraphAsymmErrors(n,array.array('d',xvalues1),array.array('d',yvalues1),array.array('d',xerr1),array.array('d',xerr1),array.array('d',yerr1do),array.array('d',yerr1up))
    graph.SetLineColor(color)
    graph.SetLineWidth(2) #1504
    #graph.SetFillStyle(3005)
    return graph

#######################################################################2
#########          CENTRAL SETTINGS         ###########################
#######################################################################

xsec = 1    ## if xsec == 0, it is probability, else xsec
R = 4       ## JET RADIUS
jetbin = 4  ## JET PT interval: 0,1,2,3,4 don't use 0
if R == 4:
    yups = [2.2,2.3,2.2,2.2] ## x-sec upper limit, this times the max
    ylimitRs = [[0.01,1.5],[0.05,1.5],[0.21,1.7],[0.1,2.7]]
    ylimits = [[-0.15,5.6], [-0.15,5.6], [-0.15,5.6], [-0.15,5.6]] ## Probability upper limit
    yupsR = [0.55,0.55,0.55,0.55] ## prob upperlimit ratio, 1 +/- this
elif R == 6: 
    yups = [2.2,2.2,2.2,2.8]
    ylimitRs = [[0.01,1.5],[0.01,1.5],[0.01,1.5],[0.1,2.7]]
    ylimits = [[-0.15,4.9], [-0.15,4.9], [-0.15,4.9], [-0.15,4.9]] ## Probability upper limit
    yupsR = [0.71,0.71,0.51,1.35] ## prob upperlimit ratio, 1 +/- this

ylimitR = ylimitRs[jetbin-1]
ylimit = ylimits[jetbin-1]

YUP = yups[jetbin-1]
YUPR = yupsR[jetbin-1]

############## Y axis title
#ytitletops = ['Probability','#frac{d^{3}#sigma}{dz_{||}^{ch} dp_{T} d#it{#eta}}']
ytitletops = ['Probability Density','#frac{d^{2}#it{#sigma}}{d#it{z}_{||}^{ch} d#it{#eta}} (mb)']
ytitletop = ytitletops[xsec]
ytitletopsize=ylabeltopsize=0.05 #ylabeltopsize=0.05
ytopoffset=1.50
if xsec == 0:
    ytopoffset = 1.2

ytitlebot='theory / data'
ytitlebotsize=ylabelbotsize=0.09 #=0.09
ybotoffset = 0.69

xlabelsize=0.085
xtitlesize=0.13
############## X axis title
xtitle = '#it{z}_{||}^{ch}'
############## jetpt bin legend
jetbintitles = ['5 < #it{p}_{T,ch jet} < 7 GeV/#it{c}', '7 < #it{p}_{T,ch jet} < 10 GeV/#it{c}','10 < #it{p}_{T,ch jet} < 15 GeV/#it{c}','15 < #it{p}_{T,ch jet} < 50 GeV/#it{c}']
jetbintitle = jetbintitles[jetbin-1]
############## 
fptbinsZN = numbincenters = 5
yvalues1 = []
xvalues1 = [0.5,0.65,0.75,0.85,0.95]
xedges = fptbinsZlh = [0.4,0.6,0.7,0.8,0.9,1.0] 
xerr1 = [0.1,0.05,0.05,0.05,0.05]
yerr1 = [1.5,1.5,1.5,1.5,1.5] ## dummy values
Dup = [7,10,15,36]
Ddo = [2,3,5,5]
############## Colors
colors = [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kYellow+2, ROOT.kOrange+2, ROOT.kBlack, ROOT.kBlue+10, ROOT.kBlue-2, ROOT.kRed-2, ROOT.kGreen-2]
#######################################################################3
#########          DATA    SETTINGS         ###########################
#######################################################################
######### FILES contatining unfolded spectra
dirname = '_25'
EOS=''#'/eos/user/a/amohanty/'
r4file2D = RT.TFile(EOS+'/media/jackbauer/data/z_out/R_04'+dirname+'/unfolding/Bayes/alljetz2D/unfold2DoutFileAPW.root','read')
r6file2D = RT.TFile(EOS+'/media/jackbauer/data/z_out/R_06'+dirname+'/unfolding/Bayes/alljetz2D/unfold2DoutFileAPW.root','read')
######### DATA Train files
data4=RT.TFile(EOS+'/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_503_R04.root','read')
data6=RT.TFile(EOS+'/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_504_R06.root','read')
# Scaling settings
simScaling = 0.5
BRDzero = 0.0389
sigma_in = 50.77
dataLum4 = data4.Get('PWG3_D2H_DmesonsForJetCorrelationsMBN0').Get('NormalizationCounter').GetNEventsForNorm()/sigma_in
dataLum6 = data6.Get('PWG3_D2H_DmesonsForJetCorrelationsMBN0').Get('NormalizationCounter').GetNEventsForNorm()/sigma_in
#######################################################################4
#########          THEORY  SETTINGS         ###########################
#######################################################################
jetRs = ['4','6']
jetnames = ['5_7','7_10','10_15','15_50']
dptnames = ['2_7','3_10','5_15','5_36']
powpyt6names = ['central', 'F1R2_1536594271','F1R05_1536598175','F2R1_1535916012','F2R2_1535895146','F05R1_1536604800','F05R05_1535894261','m13_1536595965','m17_1536655729']
Tr3file2D, Tr4file2D, Tr6file2D = [], [], [] ###### List of theory sum files for each radius
#Pow+Pyt6
for jetR in range(len(jetRs)):
    for jdname in range(len(jetnames)):
        for powpyt6name in powpyt6names:
            filename = RT.TFile('/media/jackbauer/data/z_out/R_0'+jetRs[jetR]+'_25Sim/SimFiles/Prompt/Jetpt'+jetnames[jdname]+'/Z_AnalysisResults_FastSim_powheg+pythia6_charm_'+powpyt6name+'_Jetpt'+jetnames[jdname]+'_Dpt'+dptnames[jdname]+'_Dzero.root','read')
            #if jetR == 0: Tr3file2D.append(filename)
            if jetR == 0: Tr4file2D.append(filename)
            elif jetR == 1: Tr6file2D.append(filename)

#######################################################################5
#########          FINAL DATA and THEORY  SETTINGS         ############
#######################################################################
if R == 4:
    dataLum = dataLum4
    RFILE = r4file2D
    THFILES = Tr4file2D
elif R == 6:
    dataLum = dataLum6
    RFILE = r6file2D
    THFILES = Tr6file2D 
print('data lum',dataLum)
datascaling = 1/(2*BRDzero*dataLum)
dy = 2*((9-R)/10.0)
print('dy: ',dy)
#######################################################################6
#########          FINAL DATA and THEORY HISTOGRAMS        ############
#######################################################################
######### DATA HISTO
#ylimit = [0.0001,1]

hU = [RFILE.Get('UnfProjectX_'+str(jetbin))]
histData = hU[0].Clone('hU')
#histData.GetXaxis().SetRangeUser(xedges[0],xedges[-1])
sethistoType(histData,RT.kRed+2,1.3,20,2)

if xsec == 0:
    #ylimit = [-0.1,4.4]
    #ylimitR = [0.3,1.95]
    histData = histData.Rebin(fptbinsZN,'histData',array.array('d',fptbinsZlh))
    #histData.Scale(1/histData.Integral());
    histData.Scale(1/hU[0].Integral());
    histData.Scale(1,"width");
    histData.GetYaxis().SetRangeUser(ylimit[0],ylimit[-1])
elif xsec == 1:
    #ylimitR = [0.1,2.7]
    histData = histData.Rebin(fptbinsZN,'histData',array.array('d',fptbinsZlh))
    histData.Scale(datascaling);
    histData.Scale(1/dy);histData.Scale(1,'width');
    #ylimit = [0.01*histData.GetMinimum(),100*histData.GetMaximum()]#histData.GetYaxis().SetRangeUser(-0.1*hU[0].GetMaximum(),hU[0].GetMaximum()*3)
    ylimit = [-histData.GetMaximum()/5,histData.GetMaximum()*YUP]#histData.GetYaxis().SetRangeUser(-0.1*hU[0].GetMaximum(),hU[0].GetMaximum()*3)
######### DATA SYSTEMATICS
yvalues2 = []
yerr2do = []
yerr2up = []
syspercentUp = np.array([85,95,95,95,95])*0.01
syspercentDo = np.array([85,95,95,95,95])*0.01

if R == 4:
    if jetbin == 1:
        syspercentUp = np.array([15,11, 8, 5, 9])*0.01
        syspercentDo = np.array([17,12, 9, 6, 9])*0.01
    elif jetbin == 2:
        syspercentUp = np.array([17,13, 7, 7, 5])*0.01
        syspercentDo = np.array([23,13, 8, 7, 6])*0.01
    elif jetbin == 3:
        syspercentUp = np.array([14,11, 7, 6, 7])*0.01
        syspercentDo = np.array([19,11, 8, 7, 7])*0.01
    elif jetbin == 4:
        syspercentUp = np.array([18, 7,10,10,11])*0.01
        syspercentDo = np.array([27, 7,10,11,11])*0.01
elif R == 6:
    if jetbin == 1:
        syspercentUp = np.array([12, 8, 7, 6,12])*0.01
        syspercentDo = np.array([12, 9, 7, 7,12])*0.01
    elif jetbin == 2:
        syspercentUp = np.array([ 9, 7, 7, 7, 8])*0.01
        syspercentDo = np.array([11, 7, 7, 7, 9])*0.01
    elif jetbin == 3:
        syspercentUp = np.array([13,11, 6, 8, 9])*0.01
        syspercentDo = np.array([14,11, 6, 8,10])*0.01
    elif jetbin == 4:
        syspercentUp = np.array([12, 8, 7,14,16])*0.01
        syspercentDo = np.array([14, 8, 8,14,16])*0.01
normScale = (2.1**2 + BRDzero**2)*0.0001
if xsec == 1:
    syspercentUp = np.sqrt(syspercentUp**2 + normScale)
    syspercentDo = np.sqrt(syspercentDo**2 + normScale)
for i in range(histData.GetNbinsX()):
    yerr2up.append(syspercentUp[i]*histData.GetBinContent(i+1))
    yerr2do.append(syspercentDo[i]*histData.GetBinContent(i+1))
    yvalues2.append(histData.GetBinContent(i+1))

######### THEORY HISTO
# powheg+pythia6
hT1 = TheoryPred(THFILES, powpyt6names, xsec, jetbin)

print(len(Tr6file2D))
#### data data ratio
hTR0 = histRatio([histData],histData)

#### theory data ratio
[yerr1do, yvalues1, yerr1up] = SysRange(fptbinsZN,powpyt6names,hT1)
hTR1 = histRatio(hT1, histData)
[yerr1doR, yvalues1R, yerr1upR] = SysRange(fptbinsZN,powpyt6names,hTR1)
yvalues2R = np.array([1,1,1,1,1])

print(type(yerr1do))

sethistoType(hT1[0], colors[0], 1.3, 4, 0)
sethistoType(hTR1[0], colors[0], 1.3, 4, 0)
#######################################################################7
#########          CANVAS  SETTINGS         ###########################
#######################################################################

#DEBUGGING
print(hT1[0].GetBinContent(1)*2+
      hT1[0].GetBinContent(2)+
      hT1[0].GetBinContent(3)+
      hT1[0].GetBinContent(4)+
      hT1[0].GetBinContent(5)
    )
##############
######### READYING THE CANVAS
#c1 = RT.TCanvas('c1','c1',800,20,700,700)
c1 = RT.TCanvas('c1','c1',800,20,700,700)
#l1 = RT.TLegend(0.5,0.07,0.8,0.27,'','NB NDC')
l1 = RT.TLegend(0.61,0.63,0.87,0.87,'','NB NDC')
l1.SetBorderSize(0)
p1 = RT.TPaveText(0.18,0.57,0.52,0.9,'NB NDC')

horizonMargin=0.022; linezero = RT.TLine(0.4+horizonMargin, 0, 1-horizonMargin, 0);
linezero.SetLineColor(RT.kWhite)

p1.SetBorderSize(0);
p1.SetFillStyle(0);
p1.SetTextAlign(13);
p1.SetTextFont(43);
p1.SetTextSize(18);

p1.AddText('ALICE Preliminary')
p1.AddText('pp, #sqrt{#it{s}} = 5.02 TeV')
p1.AddText('charged jets, anti-#it{k}_{T}, #it{R}=0.'+str(R)+', |#it{#eta}_{lab}^{jet}|<0.'+str(9-R))
#p1.AddText('#it{R} = 0.'+str(R)+', |#it{#eta}_{lab}^{jet}| < 0.'+str(9-R))
#p1.AddText('with D^{0} ,  '+str(Ddo[jetbin-1])+' < #it{p}_{T,D} < '+str(Dup[jetbin-1])+' GeV/#it{c}')
p1.AddText('with D^{0}, '+str(Ddo[jetbin-1])+' < #it{p}_{T,D^{0}} < '+str(Dup[jetbin-1])+' GeV/#it{c}')
p1.AddText(jetbintitle)
#p1.AddText(' < #it{p}_{T,ch jet} <  GeV/#it{c}')
#p1.AddText(jetbintitle)
#p1.AddEntry()

######### READYING THE FIRST PAD
pad1 = RT.TPad('p1','p1',0.01,0.38,1,1)   ####      Defining the first pad
pad1.SetMargin(0.16,1,0,1)              ####      Setting margin
pad1.Draw()                             ####      Drawing the first pad
pad1.cd()                               ####      Entering the first pad
#if xsec == 1:
#    pad1.cd().SetLogy()                     ####      
#pad1.cd().SetLogy()                     ####      
pad1.SetTickx(1)
pad1.SetTicky(1)

######### MULTIGRAPH
mgT = RT.TMultiGraph()
mgT.SetTitle(";;"+ytitletop)
#mgT.GetYaxis().SetTitleSize(0.95)
#mgT.GetYaxis().SetTitleOffset(0.1)
histData.SetTitle('')
hT1[0].SetTitle('')

######### TGraphAsymmErrors
gr1 = TGRAPH(numbincenters,xvalues1,yvalues1,xerr1,yerr1do,yerr1up,xedges,colors[0])
gr1.SetTitle('')
gr1.SetFillStyle(0)

#cid = RT.TColor.GetColor('#990000')
cidsys = RT.TColor.GetColor('#cccccc')
gr2 = TGRAPH(numbincenters,xvalues1,yvalues2,xerr1,yerr2do,yerr2up,xedges,cidsys)
gr2.SetFillStyle(1001)
#transcol = RT.GetColorTransparent(cidsys,0.35)
transcol=cidsys
gr2.SetFillColorAlpha(transcol,0.7)
gr2.SetLineWidth(0) #1504
gr2.SetLineColorAlpha(cidsys,0) #1504

mgT.Add(gr2)
mgT.Add(gr1)
mgT.GetYaxis().SetRangeUser(ylimit[0],ylimit[-1])
#mgT.GetYaxis().SetRangeUser(histData.GetMaximum()*0.07,histData.GetMaximum()*2.3)
#mgT.GetYaxis().SetRangeUser(-0.001,histData.GetMaximum()*2.8)
#mgT.GetYaxis().SetRangeUser(-histData.GetMaximum()/5,histData.GetMaximum()*2.8)
#print(histData.GetMaximum())
mgT.GetYaxis().SetTitleSize(ytitletopsize)
mgT.GetYaxis().SetLabelSize(ylabeltopsize)
mgT.GetYaxis().SetTitleOffset(ytopoffset)

mgT.Draw('AZ5')
histData.Draw('same')
hT1[0].Draw('same')

entry1 = l1.AddEntry(histData,'Data','p')
entry1.SetFillColor(cidsys)
entry2 = l1.AddEntry(gr2,'Sys. Unc. (data)','f')
entry2.SetFillColor(cidsys)
#entry2.SetFillStyle(1001)
entry3 = l1.AddEntry(hT1[0],'POWHEG+PYTHIA6','p')
entry4 = l1.AddEntry(gr1,'Sys. Unc. (theory)','f')
entry4.SetFillStyle(0)
l1.Draw('same')
p1.Draw('same')
linezero.Draw('same')

######### READYING THE SECOND PAD
c1.cd()                                 ####      Going back to the canvas
pad2 = RT.TPad('p2','p2',0.01,0.05,1,0.38)  ####      Defining the second pad
pad2.SetMargin(0.16,1,0.27, 0)          ####      Setting margin          ##  , bottom, top
pad2.Draw()                             ####      Drawing the second pad
pad2.cd()                               ####      Entering the second pad
pad2.SetGridy()
pad2.SetTicky(1)
pad2.SetTickx(1)
mgB = RT.TMultiGraph()
mgB.SetTitle('')

gr1b = TGRAPH(numbincenters,xvalues1,yvalues1R,xerr1,yerr1doR,yerr1upR,xedges,colors[0])
gr1b.SetFillStyle(0)
gr2b = TGRAPH(numbincenters,xvalues1,yvalues2R,xerr1,syspercentDo,syspercentUp,xedges,colors[0])
gr2b.SetFillStyle(1001)
gr2b.SetFillColorAlpha(cidsys,0.70)
gr2b.SetLineWidth(0) #1504
gr2b.SetLineColorAlpha(cidsys,0) #1504

mgB.Add(gr2b) #data systematics
mgB.Add(gr1b) #theory systematics
mgB.SetTitle(';'+xtitle+';'+ytitlebot)
if xsec == 1:
    mgB.GetYaxis().SetRangeUser(ylimitR[0],ylimitR[-1])
if xsec == 0:
    mgB.GetYaxis().SetRangeUser(1-YUPR,1+YUPR)
mgB.GetYaxis().SetTitleSize(ytitlebotsize)
mgB.GetYaxis().SetTitleOffset(ybotoffset)
mgB.GetXaxis().SetTitleSize(xtitlesize)
mgB.GetXaxis().SetTitleOffset(0.82)
mgB.GetYaxis().SetLabelSize(ylabelbotsize)
mgB.GetXaxis().SetLabelSize(xlabelsize)
mgB.Draw('AZ5')
hTR1[0].Draw('same') #Theory central
hTR0[0].Draw('same')  #Data central

######### DRAWING THE CANVAS
c1.Update()
c1.Draw()
c1.SaveAs('QM_JIRA_linear_29Oct/pdf/'+str(R)+'_'+str(xsec)+'_'+str(jetbin)+'.pdf')
c1.SaveAs('QM_JIRA_linear_29Oct/eps/'+str(R)+'_'+str(xsec)+'_'+str(jetbin)+'.eps')
c1.SaveAs('QM_JIRA_linear_29Oct/png/'+str(R)+'_'+str(xsec)+'_'+str(jetbin)+'.png')

######### SAVING IN A ROOT FILE
yerr2doROOT = histData.Clone('yerr2doROOT')
yerr2upROOT = histData.Clone('yerr2upROOT')
for i in range(len(yerr2do)):
    yerr2doROOT.SetBinContent(i+1,yerr2do[i]);yerr2upROOT.SetBinContent(i+1,yerr2up[i])

savefile = RT.TFile(str(R)+'_'+str(xsec)+'_'+str(jetbin)+'.root', 'RECREATE')
histData.Write('histData')
hT1[0].Write('hT0')
hT1[1].Write('hT1')
hT1[2].Write('hT2')
hT1[3].Write('hT3')
hT1[4].Write('hT4')
hT1[5].Write('hT5')
hT1[6].Write('hT6')
hT1[7].Write('hT7')
hT1[8].Write('hT8')
yerr2doROOT.Write('yerr2doROOT')
yerr2upROOT.Write('yerr2upROOT')
savefile.Close()

input()
