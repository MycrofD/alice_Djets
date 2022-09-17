import os, os.path
import array
import matplotlib.pyplot as plt
import ROOT as RT
import rootpy as rp
import numpy as np
import scipy as sp
#from rootpy.io import root_open
from root_numpy import hist2array

from matplotlib import colors as mcolors
import random
##----------------------------------
##----------------------------------
#--------------------------------------------------
# variable/constant values used
fptbinsZN = 6
fptbinsZA = [0.5,0.65,0.75,0.85,0.95]
fptbinsZlo = [0.4,0.6,0.7,0.8,0.9]
fptbinsZhi = [0.6,0.7,0.8,0.9,1.0]
fptbinsZlh = [0.2,0.4,0.6,0.7,0.8,0.9,1.02]
pyColor = ['red','chocolate','olive','darkorange','goldenrod','yellowgreen','indianred','cadetblue','dodgerblue','purple','crimson']
RTcolors = [RT.kGreen+1, RT.kRed+2, RT.kGreen+2, RT.kBlue+2, RT.kOrange+2, RT.kViolet+1, RT.kYellow+1,RT.kRed+1, RT.kOrange+1]
#--------------------------------------------------
#--------------------------------
# CONFIG SETTINGS
#--------------------------------
jetbinSlNo=0
R='06'
Rtitle=R+'_25'
jetbinname = ["5_7","7_10","10_15","15_50"]
djetbin = [5,6,6,6]
dptbins = djetbin[jetbinSlNo]
jetbin = jetbinname[jetbinSlNo]
jetbintitle=['5-7 GeV', '7-10 GeV', '10-15 GeV', '15-50 GeV']

#--------------------------------------------------
file1 = RT.TFile('/media/jackbauer/data/z_out/R_%s/unfolding/Bayes/alljetz2D/unfold2DoutFile2350APW.root'%(Rtitle),'read')
file2 = RT.TFile('/media/jackbauer/data/z_out/R_%s/unfolding/Bayes/alljetz2D/unfold2DoutFileJES4APW.root'%(Rtitle),'read')
files = [file1, file2]

exec("hist1 = file1.UnfProjectX_%s.Clone('h1')" %(jetbinSlNo+1))
exec("hist2 = file2.UnfProjectX_%s.Clone('h2')" %(jetbinSlNo+1))

h1Reb = hist1.Rebin(fptbinsZN,'',array.array('d',fptbinsZlh))
h2Reb = hist2.Rebin(fptbinsZN,'',array.array('d',fptbinsZlh))
h1Reb.Scale(1,'width')
h2Reb.Scale(1,'width')

hists = [h1Reb, h2Reb]

##---------FIGURE 1 Ratio
c1 = RT.TCanvas('c1','c1',800,900)
l1 = RT.TLegend(0.1,0.7,0.3,0.9)

hists[0].SetMarkerColor(RT.kRed)
hists[0].SetLineColor(RT.kRed)
hists[1].SetMarkerColor(RT.kGreen+2)
hists[1].SetLineColor(RT.kGreen+2)
hists[0].SetTitle('R'+R+': '+jetbintitle[jetbinSlNo])
hists[0].GetXaxis().SetTitle('z')
hists[0].GetYaxis().SetTitle('N')
hists[0].Draw()
hists[1].Draw('same')
l1.AddEntry(hists[0],'Default')
l1.AddEntry(hists[1],'4% tracks lost')
l1.Draw('same')
c1.Update()
c1.SaveAs('yieldJES_R'+Rtitle+str(jetbinSlNo)+'.pdf')
c1.SaveAs('yieldJES_R'+Rtitle+str(jetbinSlNo)+'.png')
wait=input()

##---------FIGURE 2 Ratio
plt.figure()
plt.hlines(  hist2array(hists[1])[1:] / hist2array(hists[0])[1:] , fptbinsZlo, fptbinsZhi, colors='red')
plt.plot(  (0.4,1.0) , (1.0,1.0), 'k.-.')
plt.ylabel('Ratio')
plt.xlabel('$z_{||}^{ch}$')
plt.title('R=%s, jet $p_{T}$: %s'%(R, jetbintitle[jetbinSlNo]))

#plt.show()
plt.savefig('JES_'+Rtitle+str(jetbinSlNo)+'.pdf')
plt.savefig('JES_'+Rtitle+str(jetbinSlNo)+'.png')
print(hist2array(hists[1])[1:])
print(hist2array(hists[0])[1:])
