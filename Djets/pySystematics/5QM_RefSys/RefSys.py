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
# Run just the SB_z.C again but now with RefSys ON, and changing both percentage and index accordingly.
##----------------------------------
if len(sys.argv)==1:
    print("   === Usage example: python RefSys.py R jetbin ")
    print("   === e.g.: python RefSys.py 02 1")
    exit()

##----------------------------------
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
jetbinSlNo=int(sys.argv[2])
R=str(sys.argv[1])#'06'

if R == '02':
    Rtitle=R+'_finaltry'
else:
    Rtitle=R+'_25'

jetbinname = ["5_7","7_10","10_15","15_50"]
djetbin = [5,6,6,6]
dptbins = djetbin[jetbinSlNo]
jetbin = jetbinname[jetbinSlNo]
jetbintitle=['5-7 GeV', '7-10 GeV', '10-15 GeV', '15-50 GeV']
#--------------------------------------------------
file0=RT.TFile('/media/jackbauer/data/z_out/R_%s/signalExtraction/plots/Z0to102_jetbin_%s/JetPtSpectra_SB_eff.root'%(Rtitle,jetbin),'read')
file1=RT.TFile('/media/jackbauer/data/z_out/R_%s/signalExtraction/plots/Z0to102_jetbin_%s/JetPtSpectra_SB_effrefSys0.root'%(Rtitle,jetbin),'read')
file2=RT.TFile('/media/jackbauer/data/z_out/R_%s/signalExtraction/plots/Z0to102_jetbin_%s/JetPtSpectra_SB_effrefSys1.root'%(Rtitle,jetbin),'read')

hist0 = file0.hjetptspectrumReb.Clone('h0')
hist1 = file1.hjetptspectrumReb.Clone('h1')
hist2 = file2.hjetptspectrumReb.Clone('h2')

h0 = hist0.Rebin(fptbinsZN, '', array.array('d',fptbinsZlh))
h1 = hist1.Rebin(fptbinsZN, '', array.array('d',fptbinsZlh))
h2 = hist2.Rebin(fptbinsZN, '', array.array('d',fptbinsZlh))

h0.Scale(1,'width')
h1.Scale(1,'width')
h2.Scale(1,'width')

#---------------------------------------------------
plt.figure()
plt.hlines( (hist2array(h1))[1:] / (hist2array(h0))[1:], fptbinsZlo, fptbinsZhi, colors='red', label='+50% ref')
plt.hlines( (hist2array(h2))[1:] / (hist2array(h0))[1:], fptbinsZlo, fptbinsZhi, colors='blue', label='-50% ref')
plt.legend(loc='lower left')

axes = plt.gca()
axes.set_xlim(0.4,1.0)
#axes.set_ylim(0,9)
plt.grid()
plt.ylabel('Ratio in %')
plt.xlabel('$z_{||}^{ch}$')
plt.title('R=%s, jet $p_{T}$: %s'%(R, jetbintitle[jetbinSlNo]))
plt.savefig('RefSysRatio%s%s.pdf'%(Rtitle,jetbinname[jetbinSlNo]))
plt.savefig('RefSysRatio%s%s.png'%(Rtitle,jetbinname[jetbinSlNo]))
plt.show()
print((hist2array(h1))[1:] / (hist2array(h0))[1:])
print((hist2array(h2))[1:] / (hist2array(h0))[1:])
