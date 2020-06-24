#x = np.arange(0,2*np.pi,0.0001)
#y = np.sin(x)
#plt.plot(x,y,'--')
#plt.show()
#plt.savefig('foo.png')
#plt.savefig('foo.pdf')#, bbox_inches='tight')

#import matplotlib
#matplotlib.use('Agg')
#matplotlib.use('TkAgg')
##matplotlib.use('Qt5Agg')
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
# 1: run zRawSys.py
# 2: which runs zrun_mainsys.csh
# 3: which runs zrun_sys.csh
# reading bin contents of each bin of a histogram
#--------------------------------------------------

if len(sys.argv)==1:
    print("   === Usage example: python file.py R jetbin ")
    print("   === e.g.: python file.py 02 0")
    exit()

# variable/constant values used
fptbinsZN = 5
fptbinsZA = [0.5,0.65,0.75,0.85,0.95]
fptbinsZlo = [0.4,0.6,0.7,0.8,0.9]
fptbinsZhi = [0.6,0.7,0.8,0.9,1.0]
fptbinsZlh = [0.4,0.6,0.7,0.8,0.9,1.02]
colors = [RT.kGreen+1, RT.kRed+2, RT.kGreen+2, RT.kBlue+2, RT.kOrange+2, RT.kViolet+1, RT.kYellow+1]
#--------------------------------------------------
def drawthis(histPerJet):
    # Plotting the Dpt histos from one trial
    for i in range(len(histPerJet)):
        exec("hh = histPerJet[i].Rebin(fptbinsZN,'h%s',array.array('d',fptbinsZlh))"%(i)) #hh = histPerJet[i].Rebin(fptbinsZN,'h%s',array.array('d',fptbinsZlh))
        exec("hh.SetMarkerColor(colors[i])") #hh.SetMarkerColor(colors[i])
        exec("hh.SetLineColor(colors[i])")   #hh.GetYaxis().SetRangeUser(0.4,1000)
        exec("hh.GetXaxis().SetRangeUser(0.4,1.02)")
        c1.cd(i+1)
        if(i==0 and tries==0):
            exec("hh.Draw()")
        else:
            exec("hh.Draw('same')")
#--------------------------------
def draw2this(hist2PerJet):
    # Plotting the Dpt histos from one trial
    for i in range(len(hist2PerJet)):
#        hh = hist2PerJet[i].Rebin(fptbinsZN,"h",array.array("d",fptbinsZlh))
        hh = hist2PerJet[i]
        hh.SetMarkerColor(colors[i%len(colors)])
        hh.SetLineColor(colors[i%len(colors)])
        #hh.GetYaxis().SetRangeUser(0.4,1000)
        hh.GetXaxis().SetRangeUser(0.4,1.02)
        if(i==0):
            hh.SetTitle('R='+Radius+', jet #it{p}_{T}: '+jetbintitle[jetbinSlNo])
            hh.Draw()
        else:
            hh.Draw("same")
#--------------------------------
# FinalRMS()
def FinalRMS2(histarrays, defaulthist): # needs array of arrays, i.e. array of TH1 histograms
    squaresum, trials = 0, len(histarrays)
    sumarr = 0
    arr = []
    hdef = defaulthist.Rebin(fptbinsZN,"",array.array("d",fptbinsZlh))
    hdef.Scale(1,'width')
    defhist = hist2array(hdef)
    #defhist = hist2array(defaulthist.Rebin(fptbinsZN,"",array.array("d",fptbinsZlh)))
    for i in range(trials):
        h = histarrays[i].Rebin(fptbinsZN,"",array.array("d",fptbinsZlh))
        h.Scale(1,'width')
        arr.append(hist2array(h))
#        arr.append(hist2array(histarrays[i].Rebin(fptbinsZN,"",array.array("d",fptbinsZlh))))
        sumarr += ((arr[i]-defhist)/defhist)**2
    meansq = sumarr/trials
    return np.sqrt(meansq)*100
#--------------------------------
# CONFIG SETTINGS
#--------------------------------
jetbinSlNo=int(sys.argv[2]) # 0 1 2 3
Radius=str(sys.argv[1]) #'06'
Rtitle=Radius
jetbinname = ["5_7","7_10","10_15","15_50"]
if Rtitle == '02':
    djetbin = [5,5,6,4]
else:
    djetbin = [5,6,6,6]
dptbins = djetbin[jetbinSlNo]
jetbin = jetbinname[jetbinSlNo]
jetbintitle=['5-7 GeV', '7-10 GeV', '10-15 GeV', '15-50 GeV']
Rtitle+='_finaltry'
#--------------------------------FIGURE 1 Individual DPT bins
c1 = RT.TCanvas("c","c",800,600)
c1.Divide(3,2)
## TRIAL N
#trials = 200
trials = 180
##---- store all dpt histograms and draw them per trial
for tries in range(0,trials+1):
    # Get a list of all Dpt histos from each trial
    # and store in a list
    exec("histPerJet%d = []"%(tries))
    exec("trialfile%d = RT.TFile('/media/jackbauer/data/z_out/R_%s_Multi/signalExtraction/plots/Z0to102_jetbin_%s/JetPtSpectra_SB_eff%d.root','read')"%(tries,Rtitle,jetbin,tries))
    for i in range(dptbins):
        exec("histPerJet%d.append(trialfile%d.hjetptsub_%d.Clone('hjetptsub_%d'))"%(tries,tries,i,i))
    exec("drawthis(histPerJet%d)"%(tries))

c1.Update()
c1.SaveAs('dpt'+Rtitle+jetbinname[jetbinSlNo]+'.png')
c1.SaveAs('dpt'+Rtitle+jetbinname[jetbinSlNo]+'.pdf')
#wait=input()
#--------------------------------FIGURE 2 SUM of DPT bins
c2 = RT.TCanvas("c2","c2",800,600)
c2.SetLogy()
hist2PerJet = []
for tries in range(0,trials+1):
    exec( "hist2PerJet.append(trialfile%d.hjetptspectrumReb.Clone('h%s'))" %(tries,tries) )
draw2this(hist2PerJet)
c2.Update()
c2.SaveAs('sum'+Rtitle+jetbinname[jetbinSlNo]+'.png')
c2.SaveAs('sum'+Rtitle+jetbinname[jetbinSlNo]+'.pdf')

#--------------------------------
hist2PerJetRand = [hist2PerJet[0]]
##mylist = [];print(mylist)

for i in range(0,160):
    x = random.randint(1,trials)
    print(x)
    hist2PerJetRand.append(hist2PerJet[x])
print(len(hist2PerJet))
print(len(hist2PerJetRand))

#wait=input()
#exit()

#--------------------------------FIGURE 3 RMS
plt.figure()
#plt.hlines(FinalRMS2(hist2PerJet,hist2PerJet[0])[1:], fptbinsZlo, fptbinsZhi)
plt.hlines(FinalRMS2(hist2PerJetRand,hist2PerJet[0])[:], fptbinsZlo, fptbinsZhi)
#plt.show()
#print(len(hist2PerJet))
#print(FinalRMS2(hist2PerJet,hist2PerJet[0])[1:])
print(FinalRMS2(hist2PerJetRand,hist2PerJet[0])[:])
for a,b in zip(fptbinsZA,FinalRMS2(hist2PerJetRand,hist2PerJetRand[0])[:]):
    plt.text(a,b, '%.2f'%(abs(b)))

axes = plt.gca()
axes.set_xlim(0.4,1.0)
axes.set_ylim(0,20)

plt.ylabel('RMS in %')
plt.xlabel('$z_{||}^{ch}$')
plt.title('R=%s, jet $p_{T}$: %s'%(Radius, jetbintitle[jetbinSlNo]))

plt.grid()
plt.savefig('rawYieldRMS%s%s.pdf'%(Rtitle,jetbinname[jetbinSlNo]))
plt.savefig('rawYieldRMS%s%s.png'%(Rtitle,jetbinname[jetbinSlNo]))

print(FinalRMS2(hist2PerJet,hist2PerJet[0])[:])
print(FinalRMS2(hist2PerJetRand,hist2PerJet[0])[:])

## For one jet pt interval:
## ----------------------------
## we have individual Dpt plots CHECK
## Now, we need multiple trials
## And plot them
## Then take the average from each Dpt bin
## And plot in a single canvas
## then take RMS of the averages

##----------------
## plan
## 1. write bash to run signalExtraction_SBz.C for different variations
## and create output files in appropriate places
## 2. read those files from python
## 3. Compute and plot them
