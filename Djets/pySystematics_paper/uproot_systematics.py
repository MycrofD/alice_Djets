import sys
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatch
import matplotlib
#matplotlib.use("TkAgg")
import uproot
from hist import Hist
from hist.intervals import ratio_uncertainty
import boost_histogram as bh
import array
import numpy
from scipy.optimize import curve_fit

from uproot_funcsettings import funcLine, HistoStyle, rootRMS, GetDigitText 
import mplhep
#mplhep.style.use("ALICE")
plt.style.use(mplhep.style.ALICE)
def uncertainty(hist_A, hist_B, operation):
    if operation == 'divide':
        f = hist_A.values() / hist_B.values()
    if operation == 'mupltiply':
        f = hist_A.values() * hist_B.values()
    if operation == 'add':
        f = hist_A.values() + hist_B.values()
    if operation == 'subtract':
        f = hist_A.values() - hist_B.values()
    sigma_A = hist_A.errors()
    sigma_B = hist_B.errors()
    value_A = hist_A.values()
    value_B = hist_B.values()
    relative_f = numpy.sqrt( (sigma_A/value_A)**2 + (sigma_B/value_B)**2 )
    sigma_f = numpy.sqrt(f**2) * relative_f
    return {'relative_uncertainty': relative_f, 'sigma': sigma_f}
#--------------------------------------------------
if len(sys.argv)==1:
    print("   === Usage example: python file.py R ")
    print("   === e.g.: python file.py 02")
    print("   === OR ===   ")
    print("   === Usage example: python file.py R jetbin")
    print("   === R=02,03,04,06. jetbin=0,1,2,3,4")
    print("   === e.g.: python file.py 02 1")
    exit()
#--------------------------------------------------
R = sys.argv[1] #'02'
lensysin = len(sys.argv)
if(lensysin>2): 
    jetbin = sys.argv[2]
    fptbinsJlh = [0.4,0.6,0.7,0.8,0.9,1.02]
    fptbinsJN=5
    fptbinsJC=[0.5,0.65,0.75,0.85,0.96]
    fDptbins = [3,5,5,6,4]
else: 
    jetbin='';print(jetbin);
    fptbinsJlh = [2,3,4,5,6,8,10,14,20,30,50]
    fptbinsJN = 10
    fptbinsJC = [2.5,3.5,4.5,5.5,7,9,12,17,25,40]
    fDptbins = [12];jetbin='0'
fptbinsJlo = fptbinsJlh[:-1]
fptbinsJhi = fptbinsJlh[1:]
##--------------------------------------------------
Rtitle=R #for accessing the directories


#ROOTColors = [ROOT.kRed+2, ROOT.kGreen+2, ROOT.kBlue+2, ROOT.kOrange+2, ROOT.kViolet+2, ROOT.kYellow+2, ROOT.kCyan-6, ROOT.kAzure+2, ROOT.kMagenta-6,ROOT.kGreen-8,ROOT.kYellow-8]
cernbox = os.environ['cernbox']
datafileFinal = uproot.open(f"{cernbox}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR{R}_paperCuts/Default/unfolding_Bayes_5/finalSpectra/JetPtSpectrum_final.root")
datafileFD = uproot.open(f"{cernbox}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR{R}_paperCuts/Default/FDsubtraction/JetPtSpectrum_FDsub.root")
datafileUnf = uproot.open(f"{cernbox}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR{R}_paperCuts/Default/unfolding_Bayes_5/unfoldedSpectrum_unfoldedJetSpectrum.root")
datafileJES = uproot.open(f"{cernbox}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR{R}_paperCuts/JES/Default/unfolding_Bayes_5/unfoldedSpectrum_unfoldedJetSpectrum.root")
if(lensysin>2):
    datafileUnf = uproot.open(f"{cernbox}/media/jackbauer/data/z_out/R_{R}_finaltry/unfolding/Bayes/alljetz2D/unfold2DoutFileRegBayesSys5.root")
    datafileJES = uproot.open(f"{cernbox}/media/jackbauer/data/z_out/R_{R}_finaltry/JES/unfolding/Bayes/alljetz2D/unfold2DoutFileRegBayesSys5.root")

#######
flagMulti=0
flagRef=0
flagClos=0
flagSigSB=0
flagCutsys=1
flagFD=0
flagUIter=0
flagUPrior=0
flagUSvd=0
flagJES=0
flagStat=0
#######
############### -----------------------------
##histos:
# 3. Cutsys
cutindices=['','2','3','5','6']
Cutstitles=['def','cut1','cut2','cut3','cut4']
sizeCutsys=5#
if(flagCutsys):
    if(lensysin>2):
        datafileCutsys=[uproot.open(f"{cernbox}/media/jackbauer/data/z_out/R_{R}_finaltry/FDsubtraction/outFD_1.root")]
        hCuts=[datafileCutsys[0][f'hsub_c{int(jetbin)+1}']]
        for i in range(1,sizeCutsys):
            datafileCutsys.append(uproot.open(f"{cernbox}/media/jackbauer/data/z_out/R_{R}_finaltry/SQ{cutindices[i]}/FDsubtraction/outFD_1.root"))
            hCuts.append(datafileCutsys[i][f'hsub_c{int(jetbin)+1}'])
    else:
        datafileCutsys = [datafileFD]
        hCuts=[datafileFD['hData_binned_sub']];
        for i in range(1,sizeCutsys):
            datafileCutsys.append(
                    uproot.open(f"{cernbox}/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR{R}_paperCuts/SQ{cutindices[i]}/Default/FDsubtraction/JetPtSpectrum_FDsub.root")
                    )
            hCuts.append(datafileCutsys[i]['hData_binned_sub'])
    #### -----------------------------
    hhCutsratio=[]
    print(uncertainty(hCuts[1],hCuts[0], 'divide')['relative_uncertainty'])
    for i in range(sizeCutsys):
        hhCutsratio.append(Hist(hCuts[i]))
        hhCutsratio[i].values = hCuts[i].values()/hCuts[0].values()
        hhCutsratio[i].errors = uncertainty(hCuts[i],hCuts[0], 'divide')['sigma']


    # plot:
    fig, ax = plt.subplots()
    
    y    = hhCutsratio[1].values
    yerr = hhCutsratio[1].errors
    xerr = [0.1,0.05,0.05,0.05,0.06]
    x = fptbinsJC
    ax.errorbar(x, y, yerr, xerr, fmt='o', linewidth=2)
    
    ax.set(xlim=(0.4, 1.02), #xticks=numpy.arange(1, 8),
           ylim=(-0.1, 2.5), #yticks=numpy.arange(1, 8)
    )
    
    plt.draw()
    plt.waitforbuttonpress(1)
    input()
    plt.close()


    
    #print(vars(hCuts[0]))
    #hhCutsratio.extend(bh.Histogram(hCuts[i])/bh.Histogram(hCuts[0]) for i in range(sizeCutsys))
    #print(hhCutsratio[1])

#        hhCutsratio[i].Divide(hCuts[0])
#        HistoStyle(hhCutsratio[i],1,2,21,0,i);
#        hhCutsratio[i].GetYaxis().SetTitle('ratio');hhCutsratio[i].SetTitle("Topological Cuts: R=0."+str(int(R)))
#        hhCutsratio[i].SetLineColor(ROOTColors[i])
#        if(i==0):
#            hhCutsratio[i].GetYaxis().SetRangeUser(0.4,2)
#            if(lensysin>2):hhCutsratio[i].GetYaxis().SetRangeUser(0.4,3)
#            hhCutsratio[i].Draw()
#        else:
#            hhCutsratio[i].Draw('same')
#        lCuts1.AddEntry(hhCutsratio[i],Cutstitles[i],"l");
#    lCuts1.Draw("same");
#
#    #c_Cuts.SaveAs('plots/3_Cuts/3_Cuts_ratio'+R+'.pdf')
#    #c_Cuts.SaveAs('plots/3_Cuts/3_Cuts_ratio'+R+'.png')
#    c_Cuts.SaveAs('plots/'+JET_or_Z+'3_Cuts/3_Cuts_ratio'+R+whichjetbin+'.pdf')
#    c_Cuts.SaveAs('plots/'+JET_or_Z+'3_Cuts/3_Cuts_ratio'+R+whichjetbin+'.png')
#    #### -----------------------------RMS
#    c_CutsRMS = ROOT.TCanvas("cCutsRMS","cCutsRMS",900,600)
#    histCutsRMS = rootRMS(hhCutsratio)
#    histCutsRMS.GetYaxis().SetRangeUser(0,1);histCutsRMS.SetLineColor(ROOT.kBlue+2);
#    histCutsRMS.SetFillColor(ROOT.kBlue+2);histCutsRMS.SetFillStyle(3654)
#    histCutsRMS.GetYaxis().SetTitle("RMS");histCutsRMS.Draw()
#    # extra for R 02
#    #histCutsRMS2= rootRMS(hhCutsratio[0:1]+hhCutsratio[3:])
#    if(lensysin==2):
#        histCutsRMS2 = rootRMS(hhCutsratio,0.2)
#    elif(lensysin>2):
#        histCutsRMS2 = rootRMS(hhCutsratio,0.6)
#    histCutsRMS3 = histCutsRMS.Clone('histCutsRMS3')
#
#    ttCuts2 = ROOT.TText(fptbinsJlh[-2],0.8,GetDigitText(histCutsRMS2,fptbinsJC))#[-10:])
#    #ttCuts2 = ROOT.TText(5,0.6,GetDigitText(histCutsRMS3,fptbinsJC))
#    histCutsRMS3.SetFillColor(ROOT.kRed+2);histCutsRMS3.SetFillStyle(3645)
#    if (lensysin==2 and R=='02'):# or (lensysin>2 and R=='02'):
#        histCutsRMS3.SetBinContent(histCutsRMS3.GetNbinsX(),histCutsRMS2.GetBinContent(histCutsRMS2.GetXaxis().FindBin(40)))
#        histCutsRMS3.Draw('same');ttCuts2.Draw('same')
#    elif (lensysin>2 and R=='02'):
#        histCutsRMS3.SetBinContent(1,histCutsRMS2.GetBinContent(histCutsRMS2.GetXaxis().FindBin(0.5)))
#        histCutsRMS3.Draw('same');
#        if(jetbin=='2'):
#            ttCuts2=ROOT.TText(fptbinsJlh[1],0.7,"w/o cut2")
#            ttCuts2.Draw('same') 
#        elif(jetbin=='4'):
#            ttCuts2 = ROOT.TText(fptbinsJlo[1],0.7,GetDigitText(histCutsRMS2,fptbinsJC)+" w/o cut2")#[-10:])
#            #ttCuts2=ROOT.TText(fptbinsJlh[1],0.7,"w/o cut2")
#            ttCuts2.Draw('same') 
#
#    # text
#    ttCuts = ROOT.TText(fptbinsJlo[1],0.9,GetDigitText(histCutsRMS,fptbinsJC))
#    ttCuts.Draw('same')
#
#    #c_CutsRMS.SaveAs('plots/3_Cuts/3_Cuts_sysRMS'+R+'.pdf')
#    #c_CutsRMS.SaveAs('plots/3_Cuts/3_Cuts_sysRMS'+R+'.png')
#    c_CutsRMS.SaveAs('plots/'+JET_or_Z+'3_Cuts/3_Cuts_sysRMS'+R+whichjetbin+'.pdf')
#    c_CutsRMS.SaveAs('plots/'+JET_or_Z+'3_Cuts/3_Cuts_sysRMS'+R+whichjetbin+'.png')
#
#    ##FITTING plot
#    yCUTS = array.array('d',histCutsRMS)[1:-1]
#    yCUTS2= array.array('d',histCutsRMS)[1:-1]
#    xCUTS = array.array('d',fptbinsJC)
#    if(lensysin>2):
#        print(f'xCUTS{xCUTS}')
#        print(f'yCUTS{yCUTS}')
#        if(R=='02'):
#            if(whichJetInZ==2 or 4):
#                xCUTS=xCUTS[1:]
#                yCUTS=yCUTS[1:]
#                print('this works!!!!!!!!')
#    else:
#        yCUTS = yCUTS[3:-1]
#        xCUTS = xCUTS[3:-1]
#    xCUTSlh= numpy.array(array.array('d',fptbinsJlh[:]))
#    print(xCUTS, yCUTS)
#    funcCUTS = funcLine #funcJES = funcBola
#    print(f'curve_fit parameters:{funcCUTS}, {xCUTS}, {yCUTS}')
#    try:
#        popt, pcov = curve_fit(funcCUTS, numpy.array(xCUTS), numpy.array(yCUTS))
#        print(f'popt:{popt}')
#        print(f'pcov:{pcov}')
#    except RuntimeError:
#        print("Error")
#
#    fig=plt.figure(1);
#    fig.tight_layout()
#    fig.subplots_adjust(bottom=0.15)
#    plt.grid(axis='both',color='0.95');
#    plt.ylim(0,0.5);#plt.xticks(np.arange(min(xJESall),max(xJESall)+5,5.0));
#    plt.xlim(fptbinsJlh[0],fptbinsJlh[-1])
#    plt.plot(xCUTSlh, funcCUTS(xCUTSlh,*popt),'r-', label='fit');
#    plt.plot(fptbinsJlo[0:2],[array.array('d',histCutsRMS)[0],array.array('d',histCutsRMS)[0]],'k')
#    if(lensysin==2 and R=='02'):
#        plt.step(fptbinsJlh[:], array.array('d',histCutsRMS2),'k', label="test");
#    else:
#        plt.step(fptbinsJlh[:], array.array('d',histCutsRMS)[:-1],'k', label="test");
#    plt.errorbar(fptbinsJC,numpy.array(array.array('d',histCutsRMS)[1:-1]),fmt='.k');
#    if(lensysin==2):
#        plt.xlabel(r'jet $p_T$', fontsize=12)
#        locCUTS = 0.4
#        locCUTS_x = 3
#    elif(lensysin>2):
#        plt.xlabel(r'$z_{||}^{ch}$', fontsize=12)
#        locCUTS = 0.47
#        locCUTS_x = 0.42
#    plt.ylabel('RMS', fontsize=12);
#    #plt.title("Topological selections, R = 0."+str(int(R)))
#    uncCUTSvals = ""
#    valsCUTS = abs(100*(funcCUTS(numpy.array(array.array('d',fptbinsJC)),*popt))) #100*(valsJES[i]-1)
#    for i in range(len(valsCUTS)):
#        uncCUTSvals += "%.1f, "%(valsCUTS[i]) #"R%s"%(R))
#    uncCUTSvals = uncCUTSvals[:-2]+' in %'
#    text_title = 'Topological selections, $R$=0.'+str(int(R))
#    text_title +='\n'+JET_or_Z_text
#    text_title +='\n'+uncCUTSvals
#    plt.text(locCUTS_x, locCUTS, text_title,#uncCUTSvals,
#         #rotation=45,
#         horizontalalignment='left',
#         verticalalignment='top',
#         #multialignment='center'
#         fontsize=12,
#         linespacing = 1.5,
#         )
#    plt.draw()
#    plt.savefig('plots/'+JET_or_Z+'3_Cuts/3_Cut_FIT'+R+whichjetbin+'.pdf')
#    plt.savefig('plots/'+JET_or_Z+'3_Cuts/3_Cut_FIT'+R+whichjetbin+'.png')
#    plt.waitforbuttonpress(1);input()
#    plt.close()
#    #if(wait):input()
############### -----------------------------
