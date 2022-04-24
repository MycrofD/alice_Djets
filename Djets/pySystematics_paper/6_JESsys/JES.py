import os, os.path, sys
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

if len(sys.argv)==1:
    print("   === Usage example: python file.py R jetbin ")
    print("   === e.g.: python file.py 02 0")
    exit()

#RT.gStyle.SetOptStat(0000)
##----------------------------------
##----------------------------------
#--------------------------------------------------
# variable/constant values used
fptbinsZN = 5
fptbinsZA = [0.5,0.65,0.75,0.85,0.95]
fptbinsZlo = [0.4,0.6,0.7,0.8,0.9]
fptbinsZhi = [0.6,0.7,0.8,0.9,1.0]
fptbinsZlh = [0.4,0.6,0.7,0.8,0.9,1.02]
pyColor = ['red','chocolate','olive','darkorange','goldenrod','yellowgreen','indianred','cadetblue','dodgerblue','purple','crimson']
RTcolors = [RT.kGreen+1, RT.kRed+2, RT.kGreen+2, RT.kBlue+2, RT.kOrange+2, RT.kViolet+1, RT.kYellow+1,RT.kRed+1, RT.kOrange+1]
#--------------------------------------------------
#--------------------------------
# CONFIG SETTINGS
#--------------------------------
jetbinSlNo=int(sys.argv[2]) # 0 1 2 3
R=str(sys.argv[1]) #'04'
Rtitle=R+'_finaltry'
jetbinname = ["5_7","7_10","10_15","15_50"]
jetbin = jetbinname[jetbinSlNo]
jetbintitle=['5-7 GeV', '7-10 GeV', '10-15 GeV', '15-50 GeV']

#--------------------------------------------------
file1 = RT.TFile('/media/jackbauer/data/z_out/R_%s/unfolding/Bayes/alljetz2D/unfold2DoutFileAPW.root'%(Rtitle),'read')
file2 = RT.TFile('/media/jackbauer/data/z_out/R_%s_JES/unfolding/Bayes/alljetz2D/unfold2DoutFileAPW.root'%(Rtitle),'read')
files = [file1, file2]

exec("hist1 = file1.UnfProjectX_%s.Clone('h1')" %(jetbinSlNo+1))
exec("hist2 = file2.UnfProjectX_%s.Clone('h2')" %(jetbinSlNo+1))

h1Reb = hist1.Rebin(fptbinsZN,'',array.array('d',fptbinsZlh))
h2Reb = hist2.Rebin(fptbinsZN,'',array.array('d',fptbinsZlh))
h1Reb.Scale(1,'width')
h2Reb.Scale(1,'width')

hists = [h1Reb, h2Reb]

##---------FIGURE 1 Overlap
c1 = RT.TCanvas('c1','c1',700,600)
l1 = RT.TLegend(0.1,0.8,0.3,0.9)
c1.SetLogy()

hists[0].SetMarkerColor(RT.kRed)
hists[0].SetLineColor(RT.kRed)
hists[1].SetMarkerColor(RT.kGreen+2)
hists[1].SetLineColor(RT.kGreen+2)
hists[0].SetTitle('R'+R+': '+jetbintitle[jetbinSlNo])
hists[0].GetXaxis().SetTitle('z')
hists[0].GetYaxis().SetTitle('yield')
hists[0].Draw()
hists[1].Draw('same')
l1.AddEntry(hists[0],'100% tracks, default')
l1.AddEntry(hists[1],'96% tracks')
l1.Draw('same')
c1.Update()
c1.SaveAs('yieldJES_R'+Rtitle+str(jetbinSlNo)+'.pdf')
c1.SaveAs('yieldJES_R'+Rtitle+str(jetbinSlNo)+'.png')
wait=input()

##---------FIGURE 2 Ratio
plt.figure()
plt.hlines(  hist2array(hists[1])[0:] / hist2array(hists[0])[0:] , fptbinsZlo, fptbinsZhi, colors='red')
for a,b in zip(fptbinsZA,hist2array(hists[1])/hist2array(hists[0])):
    plt.text(a,b, '%.2f'%(abs(b-1)))
plt.plot(  (0.4,1.0) , (1.0,1.0), 'k.-.')
plt.ylabel('Ratio: 96% tracks / 100% tracks')
plt.xlabel('$z_{||}^{ch}$')
plt.title('R=%s, jet $p_{T}$: %s'%(R, jetbintitle[jetbinSlNo]))

#plt.show()
plt.savefig('JES_'+Rtitle+str(jetbinSlNo)+'.pdf')
plt.savefig('JES_'+Rtitle+str(jetbinSlNo)+'.png')
print(hist2array(hists[1]))
print(hist2array(hists[0]))
print(type(hist2array(hists[1])/hist2array(hists[0])))
