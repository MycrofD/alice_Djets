import os, os.path, sys
import matplotlib.pyplot as plt
import ROOT as RT
import rootpy as rp
import numpy as np
import scipy as sp
from rootpy.io import root_open

from matplotlib import colors as mcolors
##----------------------------------------------------------------
from ROOT import TCanvas
##----------------------------------------------------------------
#--------------------------------------------------
if len(sys.argv)==1:
    print("   === Usage example: python iterations2.py R jetbin ")
    print("   === e.g.: python iterations2.py 02 1")
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

Rtitle = str(sys.argv[1]) #'02'
jetbin = int(sys.argv[2]) #1
jetbintitle=['5-7 GeV', '7-10 GeV', '10-15 GeV', '15-50 GeV']
jetbinname=['5_7', '7_10', '10_15', '15_50']

R=Rtitle+'_34'
#R='04_34'
if Rtitle == '02':
    R=Rtitle+'_finaltry'
else:
    R=Rtitle+'_34'

regbayes=['5','4','6']

hist = []
datafile= RT.TFile("/media/jackbauer/data/z_out/R_%s/unfolding/Bayes/alljetz2D/unfold2DoutFileAPW.root"%(R) )
exec("hist.append(datafile.UnfProjectX_%s.Clone())"%(jetbin+1))
exec("hist.append(datafile.UnfProjectXIterPre_%s.Clone())"%(jetbin+1))
exec("hist.append(datafile.UnfProjectXIterPost_%s.Clone())"%(jetbin+1))


c1 = TCanvas("c1","ex",800,900)
h0 = hist[0].Clone()
h1 = hist[1].Clone()
h2 = hist[2].Clone()
h1.SetLineColor(RT.kRed+2);h2.SetLineColor(RT.kGreen+2);
h1.SetMarkerColor(RT.kRed+2);h2.SetMarkerColor(RT.kGreen+2);
h1.SetLineWidth(2);h2.SetLineWidth(2);
h1.SetMarkerStyle(21);h2.SetMarkerStyle(21);
h1.SetMarkerSize(1);h2.SetMarkerSize(1);
h1.Divide(h0)
h2.Divide(h0)
h1.Draw()
h2.Draw("same")
h1.GetXaxis().SetRangeUser(0.4,1.0)
c1.Update()
#help()
wait = input("PRESS ENTER TO CONTINUE.")

##RATIO PLOT
#----------
plt.figure()
plt.hlines(SimpleVarRatio(hist[1],hist[0]),binlowedges,binupedges,label='Bayes Iter = 4',colors='red')
plt.hlines(SimpleVarRatio(hist[2],hist[0]),binlowedges,binupedges,label='Bayes Iter = 6',colors='blue')
plt.legend(loc='upper right')
plt.plot((0.4,1.0),(1,1),'k.-.')
plt.xlabel('$z_{||}^{ch}$')
plt.ylabel('Deviation from default Bayes Iter = 5')
plt.title('R=%s, jet $p_{T}$: %s'%(Rtitle, jetbintitle[jetbin]))
axes = plt.gca()
axes.set_xlim(0.4,1.0)
axes.set_ylim(0.97,1.25)
#plt.show()
plt.savefig('RegParams%s%s.pdf'%(Rtitle,jetbinname[jetbin]))
plt.savefig('RegParams%s%s.png'%(Rtitle,jetbinname[jetbin]))

##RMS PLOT
##--------
histarray0 = np.array(hist)
plt.figure()
plt.hlines(FinalRMS(hist,hist[0]),binlowedges,binupedges,label='',colors='red')
#plt.hlines(RMS(hist[2],hist[0]),binlowedges,binupedges,label='Gentler prior',colors='blue')
#plt.legend(loc='upper right')
#plt.plot((0.4,1.0),(1,1),'k.-.')
plt.xlabel('$z_{||}^{ch}$')
plt.ylabel('RMS in %')
plt.title('R=%s, jet $p_{T}$: %s'%(Rtitle, jetbintitle[jetbin]))
#plt.yticks(np.arange(0, 6, 1))
plt.grid()
axes = plt.gca()
axes.set_xlim(0.4,1.0)
axes.set_ylim(0,6)
#plt.show()
plt.savefig('RegParamsRMS%s%s.pdf'%(Rtitle,jetbinname[jetbin]))
plt.savefig('RegParamsRMS%s%s.png'%(Rtitle,jetbinname[jetbin]))
