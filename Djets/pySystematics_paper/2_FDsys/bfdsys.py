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
def setHistoDetails(hh, color, Mstyle, Msize, Lwidth, Lstyle):
       # setHistoDetails(hhdata_sub  ,2,20,0.8,1,1);
       # setHistoDetails(hhdata_subup,2,20,0  ,1,2);
       # setHistoDetails(hhdata_subdo,2,20,0  ,1,2);
    hh.SetMarkerColor(color);
    hh.SetMarkerStyle(Mstyle);
    hh.SetLineColor(color);
    hh.SetLineWidth(Lwidth);
    hh.SetMarkerSize(Msize);
    hh.SetLineStyle(Lstyle);
    hh.SetTitle("");
    hh.GetXaxis().SetTitle("#it{z}_{||}");
    return hh
def setErrorsZero1D(hh):
    for i in range(1,hh.GetNbinsX()+1):
        hh.SetBinError(i,0)
    return hh


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

fptbinsZN = 5
fptbinsZA = [0.5,0.65,0.75,0.85,0.95]
fptbinsZlo = [0.4,0.6,0.7,0.8,0.9]
fptbinsZhi = [0.6,0.7,0.8,0.9,1.0]
fptbinsZlh = [0.4,0.6,0.7,0.8,0.9,1.02]


jetbin = int(sys.argv[2]) # 0 1 2 3
Rtitle = str(sys.argv[1]) #'04'
jetbintitle=['5-7 GeV', '7-10 GeV', '10-15 GeV', '15-50 GeV']
jetbinname=['5_7', '7_10', '10_15', '15_50']

#R='04_34'
R=Rtitle+'_finaltry'

bfdsys=['','Up','Do']

hist = []
datafile=RT.TFile("/media/jackbauer/data/z_out/R_%s/FDsubtraction/outFD_1.root"%(R) )
exec("hist.append(datafile.hsub_c%s.Clone('h_ce'))"%(jetbin+2))
exec("hist.append(datafile.hsub_u%s.Clone('h_up'))"%(jetbin+2))
exec("hist.append(datafile.hsub_d%s.Clone('h_do'))"%(jetbin+2))
#datafile=RT.TFile("/media/jackbauer/data/z_out/R_%s/unfolding/Bayes/alljetz2D/unfold2DoutFileAPW.root"%(R) )
#exec("hist.append(datafile.UnfProjectX_%s.Clone('h_ce'))"%(jetbin+1))
#exec("hist.append(datafile.UnfProjectXUp_%s.Clone('h_up'))"%(jetbin+1))
#exec("hist.append(datafile.UnfProjectXDo_%s.Clone('h_up'))"%(jetbin+1))

#--------------------------------------------------
c0 = RT.TCanvas('c0','c0',800,600)
c0.SetLogy()
hspec0 = hist[0].Rebin(fptbinsZN,'hspec0',array.array('d',fptbinsZlh))  #Clone('hspec0')
hspec1 = hist[1].Rebin(fptbinsZN,'hspec1',array.array('d',fptbinsZlh))  #Clone('hspec1')
hspec2 = hist[2].Rebin(fptbinsZN,'hspec2',array.array('d',fptbinsZlh))  #Clone('hspec1')
hspec0.Scale(1,'width');hspec1.Scale(1,'width');hspec2.Scale(1,'width')
hspec0=setHistoDetails(hspec0,RT.kRed +2,20,0.8,1,1);
hspec1=setHistoDetails(hspec1,RT.kBlue+2,20,0.8,1,2);hspec1=setErrorsZero1D(hspec1)
hspec2=setHistoDetails(hspec2,RT.kBlue+2,20,0.8,1,2);hspec2=setErrorsZero1D(hspec2)

hspec0.GetXaxis().SetRangeUser(0.4,1.02)
hspec0.GetYaxis().SetRangeUser(0.5*hspec2.GetMinimum(), 1.5*hspec1.GetMaximum())
hspec0.SetTitle('R='+Rtitle+': '+jetbintitle[jetbin])
hspec0.Draw();hspec1.Draw('same');hspec2.Draw('same')

c0.Update()
c0.SaveAs('bfdYield%s%s.pdf'%(Rtitle,jetbinname[jetbin]))
c0.SaveAs('bfdYield%s%s.png'%(Rtitle,jetbinname[jetbin]))
##--------------------------------------------------
#RATIO PLOT
#----------
plt.figure()
plt.hlines(hist2array(hspec1)/hist2array(hspec0),fptbinsZlo,fptbinsZhi,label='FD Up',colors='red')
for a,b in zip(fptbinsZA,hist2array(hspec1)/hist2array(hspec0)):
    plt.text(a,b, '%.2f'%(abs(b-1)))
plt.hlines(hist2array(hspec2)/hist2array(hspec0),fptbinsZlo,fptbinsZhi,label='FD Down',colors='blue')
for a,b in zip(fptbinsZA,hist2array(hspec2)/hist2array(hspec0)):
    plt.text(a,b, '%.2f'%(abs(b-1)))
plt.legend(loc='upper right')
plt.plot((0.4,1.0),(1,1),'k.-.')
plt.xlabel('$z_{||}^{ch}$')
plt.ylabel('Deviation from central')
plt.title('R=%s, jet $p_{T}$: %s'%(Rtitle, jetbintitle[jetbin]))
axes = plt.gca()
axes.set_xlim(0.4,1.0)
#axes.set_ylim(0.4,1.4)
axes.set_ylim(0,2)
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
statUnc = np.array([hspec0.GetBinError(1),hspec0.GetBinError(2),hspec0.GetBinError(3),hspec0.GetBinError(4),hspec0.GetBinError(5),]) 
cont = np.array([hspec0.GetBinContent(1),hspec0.GetBinContent(2),hspec0.GetBinContent(3),hspec0.GetBinContent(4),hspec0.GetBinContent(5)])
print(statUnc)
print(cont)
