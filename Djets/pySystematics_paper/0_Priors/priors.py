import os, os.path, sys
import matplotlib.pyplot as plt
import ROOT as RT
import rootpy as rp
import numpy as np
import scipy as sp
from rootpy.io import root_open
import array

from matplotlib import colors as mcolors
##----------------------------------------------------------------
from ROOT import TCanvas, TLegend, TLine
##----------------------------------------------------------------
#--------------------------------------------------
if len(sys.argv)==1:
    print("   === Usage example: python priors.py R jetbin ")
    print("   === e.g.: python priors.py 02 1")
    exit()

# reading bin contents
def BinValues(myhist):
    c = []
    bins = myhist.GetXaxis().GetNbins()
    for i in range(1,bins+1):
        bincenter = myhist.GetBinCenter(i)
        binnum = myhist.GetXaxis().FindBin(bincenter)
        c.append(myhist.GetBinContent(binnum))
    return np.array(c)

# variation ratios
def SimpleVarRatio(hist, defaulthist): # needs array of arrays, i.e. array of TH1 histograms
    m = np.array(BinValues(hist) / BinValues(defaulthist))
    return m

# FinalRMS()
def FinalRMS(histarrays, defaulthist): # needs array of arrays, i.e. array of TH1 histograms
    squaresum, trials = 0, len(histarrays)
    for i in range(trials):
        squaresum += ( BinValues(histarrays[i]) - BinValues(defaulthist) )**2
    meansquare = squaresum/(trials)
    rms = np.sqrt(meansquare)/BinValues(defaulthist)
    print(trials,rms)
    return rms*100

#--------------------------------------------------
fptbinsZN = 5
fptbinsZA = [0.5,0.65,0.75,0.85,0.95]
fptbinsZlo = [0.4,0.6,0.7,0.8,0.9]
fptbinsZhi = [0.6,0.7,0.8,0.9,1.0]
fptbinsZlh = [0.4,0.6,0.7,0.8,0.9,1.02]

Rtitle = str(sys.argv[1]) #'02'
jetbin = int(sys.argv[2]) #1
jetbintitle=['5-7 GeV', '7-10 GeV', '10-15 GeV', '15-50 GeV']
jetbinname=['5_7', '7_10', '10_15', '15_50']

#R='04_34'
R=Rtitle+'_finaltry'
priorType=['','1','2','3','4','5','6','7','8']

RTColors = [RT.kRed+2, RT.kGreen+2, RT.kBlue+2, RT.kOrange+2, RT.kViolet+2, RT.kYellow+2, RT.kCyan+2, RT.kAzure+2, RT.kMagenta-6]
#PriorNames = ["Def MC", "z*j", "1/z * 1/j", "1/z", "z", "1/j", "j", "z/j", "j/z"] #old
PriorNames = ["Def MC", "1+(z-j/70)/2", "1-(z-j/70)/2", "1+(2*z-j/70)/3", "1-(2*z-j/70)/3", "1+(3*z-j/70)/4", "1-(3*z-j/70)/4", "1+(1+z-2*j/70)/3", "1-(1+z-2*j/70)/3"]#new

datafile, hist, priorhists, priorfile = [], [], [], []
for i in range(len(priorType)):
    datafile.append(
        RT.TFile("/media/jackbauer/data/z_out/R_%s/unfolding/Bayes/alljetz2D/unfold2DoutFileAPW%s.root"%(R,priorType[i]) ))
    priorfile.append(
        RT.TFile("/media/jackbauer/data/z_out/R_%s/unfolding/Bayes/alljetz2D/responseobject%s.root"%(R,priorType[i])) )
    print(priorfile[i])
    exec("hist.append(datafile[i].UnfProjectX_%d.Clone())"%(jetbin+1))

#######
RT.gStyle.SetOptStat(0000)
############## PRIORS Plot
############## -----------------------------
cp = TCanvas("cp","cp",900,600)
cp.Divide(3,3)
cp.SetLogz()
for i in range(9):
    exec("hp%s = priorfile[i].hGenRebin.Clone('hp%s')"%(i,i))
    cp.cd(i+1)
    cp.cd(i+1).SetLogz()
    exec("hp%s.SetTitle('%s: '+PriorNames[i])"%(i,i))
    exec("hp%s.GetXaxis().SetTitle('z')"%(i))
    exec("hp%s.GetYaxis().SetTitle('jet pt')"%(i))
    exec("hp%s.Draw('surf1')"%(i))
cp.Update()
wait = input()
############## PRIORS V DATA Plot
############## -----------------------------
exec("dataprior = datafile[0].hDataProjectX_%s"%(jetbin+1))
dataprior = dataprior.Rebin(fptbinsZN,'dataprior',array.array('d',fptbinsZlh))
dataprior.Scale(1,'width')
dataprior.Scale(1/dataprior.Integral())
dataprior.SetMarkerStyle(20)
dataprior.SetLineWidth(3)
cpd = TCanvas("cpd","cpd",900,600)
lpd = TLegend(0.6,0.1,0.9,0.4)
cpd.SetLogy()
for i in range(len(priorType)):
    exec("hpD%s = priorfile[i].hGenRebin.Clone('hpD%s')"%(i,i))
    exec("hpD%s.GetXaxis().SetRangeUser(0.4,1.02)"%(i))
    exec("hpd%s = hpD%s.ProjectionX('hpd%s',jetbin+1,jetbin+1,'E')"%(i,i,i)) 
    exec("hpd%s.Scale(1,'width')"%(i))
    exec("hpd%s.Scale(1/hpd%s.Integral())"%(i,i))
    exec("hpd%s.SetLineColor(RTColors[i])"%(i))
    exec("hpd%s.SetMarkerColor(RTColors[i])"%(i))
    exec("hpd%s.SetLineWidth(2)"%(i))
    exec("hpd%s.SetMarkerSize(1)"%(i))
    if (i==0):
        exec("hpd%s.Draw()"%(i))
        exec("hpd%s.SetTitle('Data v Priors')"%(i))
        exec("hpd%s.SetMarkerStyle(23)"%(i))
    else:
        exec("hpd%s.Draw('same')"%(i))
        exec("hpd%s.SetMarkerStyle(21)"%(i))
    exec("lpd.AddEntry(hpd%s,PriorNames[i],'l')"%(i))
dataprior.Draw('same')
lpd.AddEntry(dataprior,'Data','l')
lpd.Draw('same')
cpd.Update()

############## PRIORS V DATA Ratio Plot
############## -----------------------------
cpdR = TCanvas("cpdR","cpdR",500,300)
lpdR = TLegend(0.1,0.4,0.3,0.9)
lR = TLine(0.4, 1, 1.02, 1)
for i in range(len(priorType)):
    exec("hpdR%s = hpd%s.Clone('hpdR%s')"%(i,i,i))
    exec("hpdR%s.Divide(dataprior)"%(i))
    if i == 0:
        exec("hpdR%s.Draw()"%(i))
        exec("hpdR%s.SetTitle('%s, %s')"%(i,Rtitle, jetbintitle[jetbin]))
        exec("hpdR%s.GetYaxis().SetTitle('Priors V. Data')"%(i))
        exec("hpdR%s.SetTitleSize(0.06,'Y')"%(i))
        exec("hpdR%s.GetYaxis().SetTitleOffset(0.5)"%(i))
    else:
        exec("hpdR%s.Draw('same')"%(i))
    exec("lpdR.AddEntry(hpdR%s,PriorNames[i],'l')"%(i))
lR.Draw('same')
#lpdR.Draw('same')
cpdR.Update()

############## PRIORS V DefPrior Ratio Plot
############## -----------------------------
cpdR2 = TCanvas("cpdR2","cpdR2",500,300)
lR2 = TLine(0.4, 1, 1.02, 1)
exec("hpdRQ%s = hpd%s.Clone('hpdRQ%s')"%(0,0,0))
for i in range(1,len(priorType)):
    exec("hpdRQ%s = hpd%s.Clone()"%(i,i))
    exec("hpdRQ%s.Divide(hpdRQ0)"%(i))
    if i == 1:
        exec("hpdRQ%s.Draw()"%(i))
        exec("hpdRQ%s.SetTitle('%s, %s')"%(i,Rtitle, jetbintitle[jetbin]))
        exec("hpdRQ%s.GetYaxis().SetTitle('Priors V. Def MC')"%(i))
        exec("hpdRQ%s.SetTitleSize(0.06,'Y')"%(i))
        exec("hpdRQ%s.GetYaxis().SetTitleOffset(0.5)"%(i))
    else:
        exec("hpdRQ%s.Draw('same')"%(i))
lR2.Draw('same')
cpdR2.Update()
wait = input()
print("eeyes")
#exit()
#print("nnnes")


############## TCANVAS showing different prior-unfolded spectra ratio
############## ------------------------------------------------------
c1 = TCanvas("c1","ex",800,900)
h0 = hist[0].Clone('h0');
for i in range(len(priorType)-1):
    exec("h%s = hist[%d].Clone('h%s');"%(i+1,i+1,i+1))
    exec("h%s.SetLineColor(RTColors[%d]);"%(i+1,i+1))
    exec("h%s.SetMarkerColor(RTColors[%d]);"%(i+1,i+1))
    exec("h%s.SetLineWidth(2);"%(i+1))
    exec("h%s.SetMarkerStyle(21);"%(i+1))
    exec("h%s.SetMarkerSize(1);"%(i+1))
    exec("h%s.Divide(h0);"%(i+1))
    if (i==0): 
        exec("h%d.Draw();"%(i+1))
    else:
        exec("h%d.Draw('same');"%(i+1))
    exec("h%d.GetXaxis().SetRangeUser(0.4,1.0)"%(i+1))
c1.Update()
#wait = input("PRESS ENTER TO CONTINUE.")

############# RATIO PLOT
############# --------------------------------------------------------
plt.figure()
plt.hlines(SimpleVarRatio(hist[1],hist[0]),fptbinsZlo,fptbinsZhi,label=PriorNames[1],colors='red')
plt.hlines(SimpleVarRatio(hist[2],hist[0]),fptbinsZlo,fptbinsZhi,label=PriorNames[2],colors='blue')
plt.hlines(SimpleVarRatio(hist[3],hist[0]),fptbinsZlo,fptbinsZhi,label=PriorNames[3],colors='green')
plt.hlines(SimpleVarRatio(hist[4],hist[0]),fptbinsZlo,fptbinsZhi,label=PriorNames[4],colors='orange')
plt.hlines(SimpleVarRatio(hist[5],hist[0]),fptbinsZlo,fptbinsZhi,label=PriorNames[5],colors='yellowgreen')
plt.hlines(SimpleVarRatio(hist[6],hist[0]),fptbinsZlo,fptbinsZhi,label=PriorNames[6],colors='purple')
plt.hlines(SimpleVarRatio(hist[7],hist[0]),fptbinsZlo,fptbinsZhi,label=PriorNames[7],colors='peru')
plt.hlines(SimpleVarRatio(hist[8],hist[0]),fptbinsZlo,fptbinsZhi,label=PriorNames[8],colors='indigo')
plt.legend(loc='upper right')
plt.plot((0.4,1.0),(1,1),'k.-.')
plt.xlabel('$z_{||}^{ch}$')
plt.ylabel('Deviation')
plt.title('R=%s, jet $p_{T}$: %s'%(Rtitle, jetbintitle[jetbin]))
axes = plt.gca()
axes.set_xlim(0.4,1.0)
axes.set_ylim(0.97,1.10)
#plt.show()
#plt.yticks(np.arange(0.9, 1.1, 0.01))
plt.grid()
plt.savefig('unfoldingPriors%s%s.pdf'%(Rtitle,jetbinname[jetbin]))
plt.savefig('unfoldingPriors%s%s.png'%(Rtitle,jetbinname[jetbin]))
wait = input("PRESS ENTER TO CONTINUE.")

############## RMS PLOT
############## ----------------------------------------------------------------------
histarray0 = np.array(hist)
plt.figure()
plt.hlines(FinalRMS(hist[1:],hist[0]),fptbinsZlo,fptbinsZhi,label='',colors='red')
for a,b in zip(fptbinsZA,FinalRMS(hist[1:],hist[0])):
    plt.text(a,b, '%.2f'%(abs(b)))
#plt.legend(loc='upper right')
#plt.plot((0.4,1.0),(1,1),'k.-.')
plt.xlabel('$z_{||}^{ch}$')
plt.ylabel('RMS in %')
plt.title('R=%s, jet $p_{T}$: %s'%(Rtitle, jetbintitle[jetbin]))
plt.yticks(np.arange(0, 6, 1))
plt.grid()
axes = plt.gca()
axes.set_xlim(0.4,1.0)
axes.set_ylim(0,6)
#plt.show()
plt.savefig('unfoldingPriorsRMS%s%s.pdf'%(Rtitle,jetbinname[jetbin]))
plt.savefig('unfoldingPriorsRMS%s%s.png'%(Rtitle,jetbinname[jetbin]))

