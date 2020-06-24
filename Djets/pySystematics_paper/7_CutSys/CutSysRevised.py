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

if len(sys.argv)==1:
    print("   === Usage example: python file.py R jetbin ")
    print("   === e.g.: python file.py 02 0")
    exit()

##----------------------------------
## FUNCTIONS
def setHistoDetails(hh, color, Mstyle, Msize, Lwidth, Lstyle): 
       # setHistoDetails(hhdata_sub  ,2,20,0.8,1,1);
    hh.SetMarkerColor(color);
    hh.SetMarkerStyle(Mstyle);
    hh.SetLineColor(color);
    hh.SetLineWidth(Lwidth);
    hh.SetMarkerSize(Msize);
    hh.SetLineStyle(Lstyle);
    hh.SetTitle("");
    hh.GetXaxis().SetTitle("#it{z}_{||}");
    return hh

##----------------------------------
#--------------------------------------------------
# variable/constant values used
fptbinsZN = 5
fptbinsZA = [0.5,0.65,0.75,0.85,0.95]
fptbinsZlo = [0.4,0.6,0.7,0.8,0.9]
fptbinsZhi = [0.6,0.7,0.8,0.9,1.0]
fptbinsZlh = [0.4,0.6,0.7,0.8,0.9,1.02]
#pyColor = ['red','chocolate','olive','darkorange','goldenrod','yellowgreen','indianred','cadetblue','dodgerblue','purple','crimson']
pyColor = ['olive','blue','violet','darkorange','indianred','cadetblue','dodgerblue','purple','crimson']
#RTcolors = [RT.kGreen+1, RT.kRed+2, RT.kGreen+2, RT.kBlue+2, RT.kOrange+2, RT.kViolet+1, RT.kYellow+1,RT.kRed+1, RT.kOrange+1]
RTcolors = [ RT.kGreen+2,RT.kBlue+2,RT.kViolet+1,RT.kOrange+2,RT.kRed+2,  RT.kYellow+1,RT.kRed+1, RT.kOrange+1]
#--------------------------------------------------
#--------------------------------
# CONFIG SETTINGS
#--------------------------------
jetbinSlNo=int(sys.argv[2]) #1
R=str(sys.argv[1])#'06'
Rtitle=R+'_finaltry'
jetbinname = ["5_7","7_10","10_15","15_50"]
jetbin = jetbinname[jetbinSlNo]
jetbintitle=['5-7 GeV', '7-10 GeV', '10-15 GeV', '15-50 GeV']

#--------------------------------------------------
#'_cutfix_data_3sig'#'_cutfix_data'#'_cutfix_MC' # mc or data or whatever, when fixing the sigma or mu
suffix='_cutfix_data'#'_cutfix_data' # mc or data or whatever, when fixing the sigma or mu
cutsprefix=['','_L2','_L3','_T2','_T3']
cut_legend=['Default','L2','L3','T2','T3']
files,hists,hi=[],[],[]
for i in range(5):
    files.append(RT.TFile('/media/jackbauer/data/z_out/R_%s%s%s/FDsubtraction/outFD_1.root' %(Rtitle,cutsprefix[i],suffix), 'read'))
    exec("hists.append(files[i].hsub_c"+str(jetbinSlNo+2)+".Clone('hist'+str(i)))")
    hi.append(hists[i].Rebin(fptbinsZN,"",array.array("d",fptbinsZlh)))
    hi[i].Scale(1,'width')
    print('/media/jackbauer/data/z_out/R_%s%s%s/FDsubtraction/outFD_1.root' %(Rtitle,cutsprefix[i],suffix))

print( (hist2array(hi[1]))[:]/(hist2array(hi[0]))[:] )
print( (hist2array(hi[2]))[:]/(hist2array(hi[0]))[:] )
print( (hist2array(hi[3]))[:]/(hist2array(hi[0]))[:] )
print( (hist2array(hi[4]))[:]/(hist2array(hi[0]))[:] )
print( hist2array(hists[2]) )
print( hist2array(hists[1]) )
print( hist2array(hists[0]) )

##---------------- FIGURE 1: The Spectra
c1 = RT.TCanvas('c1','c1', 800, 600);c1.SetLogy()
l1 = RT.TLegend(0.15,0.7,0.35,0.9)

for i in [2,1,0,3,4]:#range(len(hi)):
    setHistoDetails(hi[i],RTcolors[i],20,1,2,1)
    if i == 2:
        hi[i].Draw()
    else:
        hi[i].Draw('same')
    l1.AddEntry(hi[i],cut_legend[i])
l1.Draw('same')
hi[2].SetTitle('R='+"%.1f"%(0.1*int(R))+', '+jetbintitle[jetbinSlNo])
c1.Update()

c1.SaveAs('CutsSpectra'+Rtitle+jetbinname[jetbinSlNo]+suffix+'.pdf')
c1.SaveAs('CutsSpectra'+Rtitle+jetbinname[jetbinSlNo]+suffix+'.png')
wait=input()

##---------------- FIGURE 2: The Ratio

plt.figure()

for i in [2,1,3,4]:
    plt.hlines((hist2array(hi[i]))[:]/(hist2array(hi[0]))[:],fptbinsZlo, fptbinsZhi,colors=pyColor[i],label=cut_legend[i])

plt.legend(loc='lower right')
plt.plot((0.4,1.0),(1.0,1.0),'k.-.')
plt.xlabel('$z_{||}^{ch}$')
plt.ylabel('Cut systematics: ratios')
plt.title('R=%s, jet $p_{T}$: %s'%(R, jetbintitle[jetbinSlNo]))

axes = plt.gca()
axes.set_ylim(0.4,1.4)
if R == '02':
    axes.set_ylim(0,3)

plt.savefig('CutsRatio%s%s%s.pdf'%(Rtitle,jetbinname[jetbinSlNo],suffix))
plt.savefig('CutsRatio%s%s%s.png'%(Rtitle,jetbinname[jetbinSlNo],suffix))
plt.draw()
plt.waitforbuttonpress(1);input();plt.close()
##---------------- FIGURE 3: The RMS
RMSval = []
sumval,count = 0,0
sumval += ( ( (hist2array(hi[2]))[:]-(hist2array(hi[0]))[:] ) / (hist2array(hi[0]))[:] )**2; count+=1
sumval += ( ( (hist2array(hi[1]))[:]-(hist2array(hi[0]))[:] ) / (hist2array(hi[0]))[:] )**2; count+=1
sumval += ( ( (hist2array(hi[3]))[:]-(hist2array(hi[0]))[:] ) / (hist2array(hi[0]))[:] )**2; count+=1
sumval += ( ( (hist2array(hi[4]))[:]-(hist2array(hi[0]))[:] ) / (hist2array(hi[0]))[:] )**2; count+=1  #

rms = np.array(np.sqrt(sumval/count))*100

plt.figure()
plt.hlines( rms, fptbinsZlo, fptbinsZhi)
for a,b in zip(fptbinsZA,rms):
    plt.text(a,b, '%.2f'%(abs(b-1)))

plt.xlabel('$z_{||}^{ch}$')
plt.ylabel('Cut systematics: RMS in %')
plt.title('R=%s, jet $p_{T}$: %s'%(R, jetbintitle[jetbinSlNo]))

axes = plt.gca()
axes.set_ylim(0,70)
if R == '02':
    axes.set_ylim(0,120)
if R == '02' and jetbinSlNo==1:
    sumval_2 = sumval - ( ( (hist2array(hi[2]))[:]-(hist2array(hi[0]))[:] ) / (hist2array(hi[0]))[:] )**2
    rms_2 = np.array(np.sqrt(sumval_2/(count-1)))*100
    for a,b in zip(fptbinsZA[:1],rms_2[:1]):
        plt.text(a,b, '%.2f w/o L3'%(abs(b-1)))



plt.savefig('CutsRMS%s%s%s.pdf'%(Rtitle,jetbinname[jetbinSlNo],suffix))
plt.savefig('CutsRMS%s%s%s.png'%(Rtitle,jetbinname[jetbinSlNo],suffix))

plt.draw()
plt.waitforbuttonpress(1);input();plt.close();

#plt.show()
print(rms)
exit()

##############################################################
###########      statistical uncertainties

unc0 = np.array([hist0.GetBinError(1), hist0.GetBinError(2),hist0.GetBinError(3),hist0.GetBinError(4)])/np.array([hist0.GetBinContent(1), hist0.GetBinContent(2), hist0.GetBinContent(3), hist0.GetBinContent(4)])
unc1 = np.array([hist1.GetBinError(1), hist1.GetBinError(2),hist1.GetBinError(3),hist1.GetBinError(4)])/np.array([hist1.GetBinContent(1), hist1.GetBinContent(2), hist1.GetBinContent(3), hist1.GetBinContent(4)])
unc2 = np.array([hist2.GetBinError(1), hist2.GetBinError(2),hist2.GetBinError(3),hist2.GetBinError(4)])/np.array([hist2.GetBinContent(1), hist2.GetBinContent(2), hist2.GetBinContent(3), hist2.GetBinContent(4)])
unc3 = np.array([hist3.GetBinError(1), hist3.GetBinError(2),hist3.GetBinError(3),hist3.GetBinError(4)])/np.array([hist3.GetBinContent(1), hist3.GetBinContent(2), hist3.GetBinContent(3), hist3.GetBinContent(4)])
unc4 = np.array([hist4.GetBinError(1), hist4.GetBinError(2),hist4.GetBinError(3),hist4.GetBinError(4)])/np.array([hist4.GetBinContent(1), hist4.GetBinContent(2), hist4.GetBinContent(3), hist4.GetBinContent(4)])
unc5 = np.array([hist5.GetBinError(1), hist5.GetBinError(2),hist5.GetBinError(3),hist5.GetBinError(4)])/np.array([hist5.GetBinContent(1), hist5.GetBinContent(2), hist5.GetBinContent(3), hist5.GetBinContent(4)])

plt.figure()
plt.hlines(unc0, fptbinsZlo, fptbinsZhi)
plt.savefig('DefCutsUnc%s%s.png'%(Rtitle,jetbinname[jetbinSlNo]))

plt.figure()
plt.hlines(unc1, fptbinsZlo, fptbinsZhi)
plt.savefig('LL1CutsUnc%s%s.png'%(Rtitle,jetbinname[jetbinSlNo]))

plt.figure()
plt.hlines(unc2, fptbinsZlo, fptbinsZhi)
plt.savefig('LL2CutsUnc%s%s.png'%(Rtitle,jetbinname[jetbinSlNo]))

plt.figure()
plt.hlines(unc3, fptbinsZlo, fptbinsZhi)
plt.savefig('LL3CutsUnc%s%s.png'%(Rtitle,jetbinname[jetbinSlNo]))

plt.figure()
plt.hlines(unc4, fptbinsZlo, fptbinsZhi)
plt.savefig('TT2CutsUnc%s%s.png'%(Rtitle,jetbinname[jetbinSlNo]))

plt.figure()
plt.hlines(unc5, fptbinsZlo, fptbinsZhi)
plt.savefig('T3CutsUnc%s%s.png'%(Rtitle,jetbinname[jetbinSlNo]))

#wait=input()



