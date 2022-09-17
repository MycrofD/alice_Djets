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
##----------------------------------
if len(sys.argv)==1:
    print("   === Usage example: python SigSBsys.py R jetbin ")
    print("   === e.g.: python SigSBsys.py 02 1")
    exit()

## FUNCTIONS

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
jetbinSlNo=int(sys.argv[2])#3
Rtitle=str(sys.argv[1])#'06'
jetbinname = ["5_7","7_10","10_15","15_50"]
djetbin = [5,6,6,6]
dptbins = djetbin[jetbinSlNo]
jetbin = jetbinname[jetbinSlNo]
jetbintitle=['5-7 GeV', '7-10 GeV', '10-15 GeV', '15-50 GeV']

#--------------------------------------------------

fileIndices = [0,1,2,6,7]
fileLabels = ['2 sig, 4-9 SB', '2 sig, 4-7 SB', '2 sig, 4-8 SB', '2 sig, 5-7 SB', '2 sig, 5-8 SB', '2 sig, 5-9 SB', '3 sig, 4-8 SB', '3 sig, 4-9 SB', '3 sig, 5-9 SB']
#---accessing files and histos----
#fileIndices = [0,1,2,3,4,5,6,7,8]
#fileLabels = ['2 sig, 4-7 SB', '2 sig, 4-8 SB', '2 sig, 5-7 SB', '2 sig, 5-8 SB', '2 sig, 5-9 SB', '3 sig, 4-8 SB', '3 sig, 4-9 SB', '3 sig, 5-9 SB']
files = []
hists = []
for i in range(len(fileIndices)):
    f = RT.TFile('/media/jackbauer/data/z_out/R_%s_QM_sys/signalExtractionSBranges/plots/Z0to102_jetbin_%s/JetPtSpectra_SB_eff%s.root' %(Rtitle,jetbin,fileIndices[i]), 'read')
    files.append(f)
    h = files[i].hjetptspectrumReb.Clone('h'+str(i))
    hists.append(h)


#--------------------------------------------------FIGURE 1: The histos
RMShists = []
c1 = RT.TCanvas('c1','c1',800,600)
for i in range(0,len(hists)):
    hi = hists[i].Rebin(fptbinsZN,"",array.array("d",fptbinsZlh))
    hi.Scale(1,'width')
    RMShists.append(hi)
    RMShists[i].SetMarkerColor(RTcolors[i])
    RMShists[i].SetLineColor(RTcolors[i])
    RMShists[i].GetXaxis().SetRangeUser(0.4,1.02)
    if i == 0:
        RMShists[i].Draw()
        RMShists[i].SetTitle('R=' + Rtitle + ', ' +jetbintitle[jetbinSlNo])
    else:
        RMShists[i].Draw('same')

c1.Update()
c1.SaveAs('histos'+Rtitle+jetbinname[jetbinSlNo]+'.pdf')
c1.SaveAs('histos'+Rtitle+jetbinname[jetbinSlNo]+'.png')

#--------------------------------------------------
#--------------------------------------------------FIGURE 2: The ratios
def varRatio(RMShists):
    varhists = []
    for i in range(1,len(hists)):
        k = hist2array(RMShists[i])/hist2array(RMShists[0])
        varhists.append(k[1:])

    return np.array(varhists)

print(varRatio(RMShists))

#--------------------------------------------------
plt.figure()
for i in range(len(varRatio(RMShists))):
    plt.hlines(varRatio(RMShists)[i], fptbinsZlo, fptbinsZhi,colors=pyColor[i],label=fileLabels[fileIndices[i+1]])

plt.plot((0.4,1.0),(1.0,1.0),'k.-.')
plt.legend(loc='upper right')
axes = plt.gca()
axes.set_xlim(0.4,1.0)
axes.set_ylim(0.8,1.2)
plt.ylabel('Ratio')
plt.xlabel('$z_{||}^{ch}$')
plt.title('R=%s, jet $p_{T}$: %s'%(Rtitle, jetbintitle[jetbinSlNo]))

#plt.show()
plt.savefig('SigSBratio%s%s.pdf'%(Rtitle,jetbinname[jetbinSlNo]))
plt.savefig('SigSBratio%s%s.png'%(Rtitle,jetbinname[jetbinSlNo]))
#--------------------------------------------------
#--------------------------------------------------FIGURE 2: The RMS
def FinalRMS2(RMShists):
    RMSval = []
    sumval = 0
    for i in range(1,len(RMShists)):
        sumval += (  ( hist2array(RMShists[i])-hist2array(RMShists[0]) )/(hist2array(RMShists[0])) )**2  /  ( len(RMShists)-1  )


    return np.array(np.sqrt(sumval))[1:]*100

print(FinalRMS2(RMShists))

plt.figure()
plt.hlines(FinalRMS2(RMShists), fptbinsZlo, fptbinsZhi)
axes = plt.gca()
axes.set_xlim(0.4,1.0)
axes.set_ylim(0,9)
plt.grid()
plt.ylabel('RMS in %')
plt.xlabel('$z_{||}^{ch}$')
plt.title('R=%s, jet $p_{T}$: %s'%(Rtitle, jetbintitle[jetbinSlNo]))

plt.savefig('SigSBrms%s%s.pdf'%(Rtitle,jetbinname[jetbinSlNo]))
plt.savefig('SigSBrms%s%s.png'%(Rtitle,jetbinname[jetbinSlNo]))
#plt.show()
#--------------------------------------------------
wait = input()


