import os, os.path, sys
import matplotlib.pyplot as plt
import ROOT as RT
import rootpy as rp
import numpy as np
import scipy as sp
from rootpy.io import root_open
from root_numpy import hist2array

import array
from matplotlib import colors as mcolors
##----------------------------------------------------------------
from ROOT import TCanvas
##----------------------------------------------------------------
#--------------------------------------------------

if len(sys.argv)==1:
    print("   === Usage example: python bfdsys.py R jetbin ")
    print("   === e.g.: python bfdsys.py 02 1")
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
    meansquare = squaresum/(trials-1)
    rms = np.sqrt(meansquare)/BinValues(defaulthist)
    print(trials,rms)
    return rms*100

#--------------------------------------------------
binlowedges = np.array([0.0,0.4,0.6,0.7,0.8,0.9])
binupedges = np.array([0.4,0.6,0.7,0.8,0.9,1.0])

fptbinsZN = 6
fptbinsZA = [0.5,0.65,0.75,0.85,0.95]
fptbinsZlo = [0.4,0.6,0.7,0.8,0.9]
fptbinsZhi = [0.6,0.7,0.8,0.9,1.0]
fptbinsZlh = [0.2,0.4,0.6,0.7,0.8,0.9,1.02]


jetbin = int(sys.argv[2]) # 0 1 2 3
Rtitle = str(sys.argv[1]) #'04'
jetbintitle=['5-7 GeV', '7-10 GeV', '10-15 GeV', '15-50 GeV']
jetbinname=['5_7', '7_10', '10_15', '15_50']

#R='04_34'
R=Rtitle+'_25'
if Rtitle == '02':
    R=Rtitle+'_finaltry'
else:
    R=Rtitle+'_25'

bfdsys=['','Up','Do']

datafile, hist = [], []
#for i in range(len(bfdsys)):
for i in range(1):
    datafile.append(
    #    RT.TFile('/media/jackbauer/data/z_out/R_%s/FDsubtraction/Jetbin_%s/plots/JetPtSpectrum_FDsub.root'%(R,jetbinname[jetbin])))
    #exec("hist.append(  datafile[i].hData_binned_sub.Clone())")
    #exec("hist.append(datafile[i].hData_binned_sub_up.Clone())")
    #exec("hist.append(datafile[i].hData_binned_sub_down.Clone())")
        #RT.TFile("/media/jackbauer/data/z_out/R_%s/unfolding/Bayes/alljetz2D/unfold2DoutFile2350APW.root"%(R) ))
        RT.TFile("/media/jackbauer/data/z_out/R_%s/unfolding/Bayes/alljetz2D/unfold2DoutFileAPW.root"%(R) ))
    exec("hist.append(  datafile[i].UnfProjectX_%s.Clone())"%(jetbin+1))
    exec("hist.append(datafile[i].UnfProjectXUp_%s.Clone())"%(jetbin+1))
    exec("hist.append(datafile[i].UnfProjectXDo_%s.Clone())"%(jetbin+1))

#--------------------------------------------------
c0 = RT.TCanvas('c0','c0',800,600)
hspec0 = hist[0].Rebin(fptbinsZN,'hspec0',array.array('d',fptbinsZlh))  #Clone('hspec0')
hspec1 = hist[1].Rebin(fptbinsZN,'hspec1',array.array('d',fptbinsZlh))  #Clone('hspec1')
hspec2 = hist[2].Rebin(fptbinsZN,'hspec2',array.array('d',fptbinsZlh))  #Clone('hspec1')
hspec0.Scale(1,'width');hspec1.Scale(1,'width');hspec2.Scale(1,'width')
hspec1.SetLineColor(RT.kRed+2);hspec2.SetLineColor(RT.kGreen+2);
hspec1.SetMarkerColor(RT.kRed+2);hspec2.SetMarkerColor(RT.kGreen+2);

hspec0.GetXaxis().SetRangeUser(0.4,1.02)
hspec0.SetTitle('R='+Rtitle+': '+jetbintitle[jetbin])
hspec0.Draw();hspec1.Draw('same');hspec2.Draw('same')

c0.Update()
c0.SaveAs('bfdYield%s%s.pdf'%(Rtitle,jetbinname[jetbin]))
c0.SaveAs('bfdYield%s%s.png'%(Rtitle,jetbinname[jetbin]))
##--------------------------------------------------
#c1 = TCanvas("c1","ex",800,900)
#h0 = hist[0].Clone('h0')
#h1 = hist[1].Clone('h1')
#h2 = hist[2].Clone('h2')
#h1.SetLineColor(RT.kRed+2);h2.SetLineColor(RT.kGreen+2);
#h1.SetMarkerColor(RT.kRed+2);h2.SetMarkerColor(RT.kGreen+2);
#h1.SetLineWidth(2);h2.SetLineWidth(2);
#h1.SetMarkerStyle(21);h2.SetMarkerStyle(21);
#h1.SetMarkerSize(1);h2.SetMarkerSize(1);
#h1.Divide(h0)
#h2.Divide(h0)
#h1.Draw()
#h2.Draw("same")
#h1.GetXaxis().SetRangeUser(0.4,1.0)
#c1.Update()
##help()
#wait = input("PRESS ENTER TO CONTINUE.")

#RATIO PLOT
#----------
plt.figure()
#plt.hlines(SimpleVarRatio(hist[1],hist[0]),binlowedges,binupedges,label='FD Up',colors='red')
#plt.hlines(SimpleVarRatio(hist[2],hist[0]),binlowedges,binupedges,label='FD Down',colors='blue')
plt.hlines(hist2array(hspec1)/hist2array(hspec0),binlowedges,binupedges,label='FD Up',colors='red')
plt.hlines(hist2array(hspec2)/hist2array(hspec0),binlowedges,binupedges,label='FD Down',colors='blue')
plt.legend(loc='upper right')
plt.plot((0.4,1.0),(1,1),'k.-.')
plt.xlabel('$z_{||}^{ch}$')
plt.ylabel('Deviation from central')
plt.title('R=%s, jet $p_{T}$: %s'%(Rtitle, jetbintitle[jetbin]))
axes = plt.gca()
axes.set_xlim(0.4,1.0)
axes.set_ylim(0.4,1.6)
#plt.show()
plt.savefig('bfdsys%s%s.pdf'%(Rtitle,jetbinname[jetbin]))
plt.savefig('bfdsys%s%s.png'%(Rtitle,jetbinname[jetbin]))
#
print(((hist2array(hspec1)[1:])/(hist2array(hspec0)[1:]))-1)
print(1-((hist2array(hspec2)[1:])/(hist2array(hspec0)[1:])))
input()

## Stat Unc
##---------
plt.figure()
#statUnc = np.array([hspec0.GetBinError(1),hspec0.GetBinError(2),hspec0.GetBinError(3),hspec0.GetBinError(4),hspec0.GetBinError(5),])/np.array([hspec0.GetBinContent(1),hspec0.GetBinContent(2),hspec0.GetBinContent(3),hspec0.GetBinContent(4),hspec0.GetBinContent(5)])
#plt.hlines(statUnc,binlowedges,binupedges,label='stat unc.',colors='red')
##plt.legend(loc='upper right')
#plt.xlabel('$z_{||}^{ch}$')
#plt.ylabel('Relative stat unc. in %')
#plt.show()

statUnc = np.array([hspec0.GetBinError(1),hspec0.GetBinError(2),hspec0.GetBinError(3),hspec0.GetBinError(4),hspec0.GetBinError(5),]) 
cont = np.array([hspec0.GetBinContent(1),hspec0.GetBinContent(2),hspec0.GetBinContent(3),hspec0.GetBinContent(4),hspec0.GetBinContent(5)])
print(statUnc)
print(cont)
##RMS PLOT
##--------
#histarray0 = np.array(hist)
#plt.figure()
#plt.hlines(FinalRMS(hist,hist[0]),binlowedges,binupedges,label='',colors='red')
##plt.hlines(RMS(hist[2],hist[0]),binlowedges,binupedges,label='Gentler prior',colors='blue')
##plt.legend(loc='upper right')
##plt.plot((0.4,1.0),(1,1),'k.-.')
#plt.xlabel('$z_{||}^{ch}$')
#plt.ylabel('RMS in %')
#plt.title('R=%s, jet $p_{T}$: %s'%(Rtitle, jetbintitle[jetbin]))
##plt.yticks(np.arange(0, 6, 1))
#plt.grid()
#axes = plt.gca()
#axes.set_xlim(0.4,1.0)
#axes.set_ylim(0,4)
##plt.show()
#plt.savefig('RegParamsRMS%s%s.pdf'%(Rtitle,jetbinname[jetbin]))
#plt.savefig('RegParamsRMS%s%s.png'%(Rtitle,jetbinname[jetbin]))
