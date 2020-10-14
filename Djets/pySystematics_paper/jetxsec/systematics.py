import os, os.path, sys
import matplotlib.pyplot as plt
import ROOT as RT
import ROOT as ROOT
import rootpy as rp
import numpy as np
import scipy as sp
from rootpy.io import root_open
import array

from matplotlib import colors as mcolors
##----------------------------------------------------------------
from ROOT import TCanvas, TLegend, TLine
##----------------------------------------------------------------
from funcsettings import * # HistoMarkers#from stylesettings import * # HistoMarkers
#--------------------------------------------------
if len(sys.argv)==1:
    print("   === Usage example: python file.py R ")
    print("   === e.g.: python file.py 02")
    print("   === OR ===   ")
    print("   === Usage example: python file.py R jetbin")
    print("   === e.g.: python file.py 02 1")
    exit()
wait=1
#--------------------------------------------------
#fptbinsJN = 7
#fptbinsJC = [5.5,7,9,12,17,25,40]
#fptbinsJlo = [5,6,8,10,14,20,30]
#fptbinsJhi = [6,8,10,14,20,30,50]
#fptbinsJlh = [5,6,8,10,14,20,30,50]
##--------------------------------------------------
fptbinsJNreal = 10
fptbinsJlhreal = [2,3,4,5,6,8,10,14,20,30,50]
##--------------------------------------------------
R = str(sys.argv[1]) #'02'
lensysin = len(sys.argv)
if(lensysin>2): 
    jetbin = str(sys.argv[2]);print(jetbin);print(type(jetbin))
    fptbinsJlh = [0.2,0.4,0.6,0.7,0.8,0.9,1.02]
    fptbinsJN=6
    fptbinsJC=[0.3,0.5,0.65,0.75,0.85,0.96]
else: 
    jetbin='';print(jetbin);print(type(jetbin))
    fptbinsJlh = [2,3,4,5,6,8,10,14,20,30,50]
    fptbinsJN = 10
    fptbinsJC = [2.5,3.5,4.5,5.5,7,9,12,17,25,40]
fptbinsJlo = fptbinsJlh[:-1]
fptbinsJhi = fptbinsJlh[1:]
##--------------------------------------------------
Rtitle=R #for accessing the directories

RTColors = [RT.kRed+2, RT.kGreen+2, RT.kBlue+2, RT.kOrange+2, RT.kViolet+2, RT.kYellow+2, RT.kCyan-6, RT.kAzure+2, RT.kMagenta-6,RT.kGreen-8,RT.kYellow-8]

datafileFinal = RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_Bayes_5/finalSpectra/JetPtSpectrum_final.root"%(R))
datafileFD = RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/FDsubtraction/JetPtSpectrum_FDsub.root"%(R))
datafileSig = RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/signalExtraction/JetPtSpectra_SB_eff.root"%(R))

datafileUnf = RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_Bayes_5/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R))
datafileJES = RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/JESsysFinal_DzeroR%s_paperCuts/Default/unfolding_Bayes_5/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R))
#######
flagMulti=1
flagRef=0
flagSigSB=0
flagCutsys=0
flagFD=0
flagUIter=0
flagUPrior=0
flagUSvd=0
flagJES=0
#######
RT.gStyle.SetOptStat(0000)
ROOT.gStyle.SetLegendBorderSize(0)
############## -----------------------------
#histos:
# 0. Multi trial
## raw systematics files
sizeMulti=973
Multititles=[
        ]
if(flagMulti):
    datafileMulti = [
            RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/RawSysFinal_DzeroR%s_paperCuts/Default/signalExtraction_multitrial/JetPtSpectra_SB_eff.root"%(R))
            ]
    hMulti=[
            datafileMulti[0].Get('hjetptspectrumReb').Clone('hMulti_0')
            ];
    for i in range(1,sizeMulti):
        try:
            datafileMulti.append(
                RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/RawSysFinal_DzeroR%s_paperCuts/Default/signalExtraction_multitrial/JetPtSpectra_SB_eff%d.root"%(R,i))
                )
            hMulti.append(
                datafileMulti[i].Get('hjetptspectrumReb').Clone('hMulti_'+str(i))
                )
        except:
            pass
    #### -----------------------------
    c_MultiRaw = TCanvas("cMultiRaw","cMultiRaw",1200,900)
    c_MultiRaw.Divide(3,4)
    lMulti1 = TLegend(0.12,0.60,0.38,0.88);
    hhMultiRaw=[]
    colorindex=30
    #for i in range(sizeMulti):
    for i in range(len(datafileMulti)):
        try:
            hhDptSpec = []
            if(i%20==1): colorindex += 1
            for dbins in range(12):
                hhdptspecRaw = datafileMulti[i].Get('hjetptsub_'+str(dbins)).Clone('hjetptsub_'+str(dbins))
                hhdptspectra=hhdptspecRaw.Rebin(fptbinsJNreal,'hjetptsubReb_'+str(dbins),array.array('d',fptbinsJlhreal)) 
                hhDptSpec.append(hhdptspectra)
                c_MultiRaw.cd(dbins+1)
                c_MultiRaw.cd(dbins+1).SetLogy()
                hhDptSpec[dbins].SetMinimum(0.5)
                hhDptSpec[dbins].SetLineColor(colorindex)#RTColors[colorindex])
                if(i==0): hhDptSpec[dbins].Draw() #default
                else: hhDptSpec[dbins].Draw('same')
            hhMultiRaw.append(hhDptSpec)
        except:
            pass

    c_MultiRaw.SaveAs('plots/0_Multi/0_Multi_sysRaw'+R+'.pdf')
    c_MultiRaw.SaveAs('plots/0_Multi/0_Multi_sysRaw'+R+'.png')
    c_MultiRaw.SaveAs('plots/0_Multi/0_Multi_sysRaw'+R+'.svg')
    #### -----------------------------Distribution of yields
    c_MultiYieldDist = TCanvas("cMultiYieldDist","cMultiYieldDist",900,600)
    c_MultiYieldDist.Divide(4,3)
    hTrialsDist=[]
    for i in range(fptbinsJN):
        c_MultiYieldDist.cd(i+1)
        c_MultiYieldDist.cd(i+1).SetMargin(2,0,0.9,0.9)
        hTrialsDist.append(yieldDist(hMulti,fptbinsJC[i]))#hTrialsPerPt.append(hYPTValues[0])
        hTrialsDist[i].Draw()
    c_MultiYieldDist.SaveAs('plots/0_Multi/0_Multi_yieldDist'+R+'.pdf')
    c_MultiYieldDist.SaveAs('plots/0_Multi/0_Multi_yieldDist'+R+'.png')
    #### -----------------------------
    c_MultiYieldTrialsJetpt = TCanvas("cMultiYieldTrialsJetpt","cMultiYieldTrialsJetpt",900,600)
    c_MultiYieldTrialsJetpt.Divide(4,3)
    hTrialsPerPt = []

    meanYPT,line_YPT=[],[]
    sigYPTplus,line_YPTplus,sigUp=[],[],1.5
    sigYPTminu,line_YPTminu,sigDo=[],[],1.5
    for i in range(fptbinsJN):
        c_MultiYieldTrialsJetpt.cd(i+1)
        c_MultiYieldTrialsJetpt.cd(i+1).SetMargin(2,0,0.9,0.9)
        hTrialsPerPt.append(yieldPtrials(hMulti,fptbinsJC[i]))#hTrialsPerPt.append(hYPTValues[0])
        hTrialsPerPt[i].SetTitle(str(fptbinsJlo[i])+' < #it{p}_{T,ch.jet} < '+str(fptbinsJhi[i])+' GeV')
        hTrialsPerPt[i].Draw()
        
        c_MultiYieldTrialsJetpt.cd(i+1).Update()
        meanYPT.append(hMulti[0].GetBinContent(hMulti[0].GetXaxis().FindBin(fptbinsJC[i]))) 
        sigYPTplus.append(hMulti[0].GetBinContent(hMulti[0].GetXaxis().FindBin(fptbinsJC[i])) + sigUp*hMulti[0].GetBinError(hMulti[0].GetXaxis().FindBin(fptbinsJC[i])) )
        sigYPTminu.append(hMulti[0].GetBinContent(hMulti[0].GetXaxis().FindBin(fptbinsJC[i])) - sigDo*hMulti[0].GetBinError(hMulti[0].GetXaxis().FindBin(fptbinsJC[i])) )
        line_YPT.append(ROOT.TLine(0,meanYPT[i],len(datafileMulti), meanYPT[i]))
        line_YPTplus.append(ROOT.TLine(0,sigYPTplus[i],len(datafileMulti), sigYPTplus[i]))
        line_YPTminu.append(ROOT.TLine(0,sigYPTminu[i],len(datafileMulti), sigYPTminu[i]))
        line_YPT[i].SetLineColor(ROOT.kRed+2)
        line_YPTplus[i].SetLineColor(ROOT.kRed+2)
        line_YPTminu[i].SetLineColor(ROOT.kRed+2)
        c_MultiYieldTrialsJetpt.cd(i+1)
        line_YPT[i].Draw('same')
        line_YPTplus[i].Draw('same')
        line_YPTminu[i].Draw('same')

    c_MultiYieldTrialsJetpt.cd(fptbinsJN+1);
    ttMultiYPT= [];delYPTtt=0.066
    YPTvariations=['#sigma = free, #sigma_{MC}, #sigma_{MC}*1.1, #sigma_{MC}*0.9 ', 
            'bkg = 0(exp), 1(lin), 2(poly2)',
            'mass = free, m_{PDG}',
            'lower limit = 1.71, 1.72, 1.70',
            'upper limit = 2.10, 2.09, 2.11',
            'mass bin width= 2, 4, 1   ']
    for i in range(len(YPTvariations)):
        ttMultiYPT.append(ROOT.TLatex(0.5,.5-i*delYPTtt,str(i)+": "+YPTvariations[i]))
        ttMultiYPT[i].SetTextFont(82);
        ttMultiYPT[i].SetTextAlign(22);
        ttMultiYPT[i].SetTextSize(0.06);
        ttMultiYPT[i].SetTextAngle(0)
        ttMultiYPT[i].SetTextColor(ROOT.kRed+2)
        ttMultiYPT[i].Draw('same')
    
    c_MultiYieldTrialsJetpt.cd(fptbinsJN+2)
    ttMultiYPT2 = [ROOT.TLatex(.12,.8,'R = 0.'+str(int(R))+', Multi-trial'),
            ROOT.TLatex(.12,.7,'pp, 5.02 TeV'),
            ROOT.TLatex(.12,.6,'')
            ]
    for i in range(len(ttMultiYPT2)):ttMultiYPT2[i].SetTextSize(0.1);ttMultiYPT2[i].SetTextFont(42);ttMultiYPT2[i].Draw('same')

    c_MultiYieldTrialsJetpt.SaveAs('plots/0_Multi/0_Multi_YieldTrialsJetpt'+R+'.pdf')
    c_MultiYieldTrialsJetpt.SaveAs('plots/0_Multi/0_Multi_YieldTrialsJetpt'+R+'.png')
    c_MultiYieldTrialsJetpt.SaveAs('plots/0_Multi/0_Multi_YieldTrialsJetpt'+R+'.svg')
    #### -----------------------------
    c_MultiRatio = TCanvas("cMultiRatio","cMultiRatio",900,600)
    lMulti3 = TLegend(0.12,0.60,0.38,0.88);
    hhMultiratio=[]
    #for i in range(sizeMulti):
    for i in range(len(datafileMulti)):
        try:
            hhMultiratio.append(hMulti[i].Clone('hMultiratio_'+str(i)))
            hhMultiratio[i].Divide(hMulti[0])
            HistoStyle(hhMultiratio[i],1,2,21,0,RTColors[0]);
            hhMultiratio[i].GetYaxis().SetTitle('ratio');hhMultiratio[i].SetTitle("Multi-trial: R=0."+str(int(R)))
            hhMultiratio[i].SetLineColor((i%20)+29)
            hhMultiratio[i].GetYaxis.SetRangeUser(0,2)
            if(i==1): hhMultiratio[i].Draw(); 
            elif(i>1): hhMultiratio[i].Draw('same')
        except: pass
    c_MultiRatio.SaveAs('plots/0_Multi/0_Multi_sysRatio'+R+'.pdf')
    c_MultiRatio.SaveAs('plots/0_Multi/0_Multi_sysRatio'+R+'.png')
    #### -----------------------------RMS
    c_MultiRMS = TCanvas("cMultiRMS","cMultiRMS",900,600)
    lMulti4 = TLegend(0.12,0.60,0.38,0.88);
    histMultiRMS = rootRMS_mt(hhMultiratio,sigUp,sigDo);
    histMultiRMS.SetLineColor(ROOT.kRed+2);histMultiRMS.SetFillColor(ROOT.kRed+2);histMultiRMS.SetFillStyle(3354);
    histMultiRMS.GetYaxis().SetRangeUser(0,0.25);histMultiRMS.GetYaxis().SetTitle("RMS");
    histMultiRMS.Draw()
    #lMulti4.AddEntry(histMultiRMS,'');
    lMulti4.Draw('same')

    # text
    ytextplace=0.22
    ttMultiRMS= ROOT.TText(18,ytextplace,GetDigitText(histMultiRMS,fptbinsJC));ttMultiRMS.SetTextSize(0.04);
    ttMultiRMS.Draw('same')
    c_MultiRMS.SaveAs('plots/0_Multi/0_Multi_sysRMS'+R+'.pdf')
    c_MultiRMS.SaveAs('plots/0_Multi/0_Multi_sysRMS'+R+'.png')

    #if(wait):input()
############## -----------------------------
#histos:
# 1. Reflection
## 
if(flagRef):
    datafileRf0 = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/RawSysFinal_DzeroR%s_paperCuts/Default/signalExtraction_refsys/JetPtSpectra_SB_effrefSys0.root"%(R),"read")
    datafileRf1 = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/RawSysFinal_DzeroR%s_paperCuts/Default/signalExtraction_refsys/JetPtSpectra_SB_effrefSys1.root"%(R),"read")
    hRdf = datafileSig.Get('hjetptspectrumRebScaled').Clone('hRdf')
    hRf0 = datafileRf0.Get('hjetptspectrumRebScaled').Clone('hRf0')
    hRf1 = datafileRf1.Get('hjetptspectrumRebScaled').Clone('hRf1')
    #### -----------------------------
    c_Ref = TCanvas("cRef","cRef",900,600)
    leg_Ref = ROOT.TLegend(0.7,0.7,0.88,0.88)
    hRf0ratio=hRf0.Clone('hRf0ratio');hRf0ratio.Divide(hRdf);HistoStyle(hRf0ratio,1,2,21,0,ROOT.kRed+2);
    hRf1ratio=hRf1.Clone('hRf1ratio');hRf1ratio.Divide(hRdf);HistoStyle(hRf1ratio,1,2,21,0,ROOT.kGreen+2);
    for i in range(hRf0ratio.GetNbinsX()): hRf0ratio.SetBinError(i+1,0);hRf1ratio.SetBinError(i+1,0)
    hRf0ratio.GetYaxis().SetRangeUser(0.9,1.1)
    hRf0ratio.Draw();hRf1ratio.Draw('same')
    leg_Ref.AddEntry(hRf0ratio,"refl,-50%","l");leg_Ref.AddEntry(hRf1ratio,"refl,+50%","l");leg_Ref.Draw('same')
    lRef = ROOT.TLine(2, 1, 50, 1);lRef.SetLineStyle(2);lRef.Draw("same")
    ttRefsys0 = ROOT.TText(5,1.05,GetDigitTextFromRatio(hRf0ratio,fptbinsJC))
    ttRefsys1 = ROOT.TText(5,0.95,GetDigitTextFromRatio(hRf1ratio,fptbinsJC))
    ttRefsys0.Draw('same');ttRefsys1.Draw('same')
    c_Ref.SaveAs('plots/1_Ref/1_Ref_rat'+R+'.pdf')
    c_Ref.SaveAs('plots/1_Ref/1_Ref_rat'+R+'.png')

    if(wait):input()
############## -----------------------------
#histos:
# 2. Signal and SB ranges
## raw systematics files
sigSBsize=18
sigSBtitles=[
        'S:2\sigma, SB:4-9\sigma',#default, not necessary
        'S:3\sigma, SB:4-9\sigma',
        'S:2\sigma, SB:4-8\sigma',
        'S:3\sigma, SB:4-8\sigma',
        'S:2\sigma, SB:4-7\sigma',
        'S:3\sigma, SB:4-7\sigma',
        'S:2\sigma, SB:4.5-9\sigma',
        'S:3\sigma, SB:4.5-9\sigma',
        'S:2\sigma, SB:4.5-8\sigma',
        'S:3\sigma, SB:4.5-8\sigma',
        'S:2\sigma, SB:4.5-7\sigma',
        'S:3\sigma, SB:4.5-7\sigma',
        'S:2\sigma, SB:3.5-9\sigma',
        'S:3\sigma, SB:3.5-9\sigma',
        'S:2\sigma, SB:3.5-8\sigma',
        'S:3\sigma, SB:3.5-8\sigma',
        'S:2\sigma, SB:3.5-7\sigma',
        'S:3\sigma, SB:3.5-7\sigma',
        ]
if(flagSigSB):
    datafileSigSB = [
            RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/RawSysFinal_DzeroR%s_paperCuts/Default/signalExtraction_SigSBranges/JetPtSpectra_SB_eff.root"%(R))
            ]
    hSigSB=[
            datafileSigSB[0].Get('hjetptspectrumRebScaled').Clone('hSigSB_0')
            ];
    for i in range(1,sigSBsize):
        datafileSigSB.append(
            RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/RawSysFinal_DzeroR%s_paperCuts/Default/signalExtraction_SigSBranges/JetPtSpectra_SB_eff%d.root"%(R,i))
            )
        hSigSB.append(
            datafileSigSB[i].Get('hjetptspectrumRebScaled').Clone('hSigSB_'+str(i))
            )
    #### -----------------------------
    c_SigSB = TCanvas("cSigSB","cSigSB",900,600)
    lSigSB1 = TLegend(0.12,0.50,0.38,0.88);
    hhSigSBratio=[]
    excludedSigSB=[]
    #if(R=='04'):excludedSigSB=[6,7,8,9,10,11]
    if(R=='04' or R=='06'):excludedSigSB=[6,7,8,9,10,11]
    #elif(R=='06'):excludedSigSB=[7,9,11]
    #elif(R=='06'):excludedSigSB=[11]
    for i in range(sigSBsize):
    #for i in range(5):#for i in [0,1,2,6,7]:
        hhSigSBratio.append(hSigSB[i].Clone('hSigSBratio_'+str(i)))
        hhSigSBratio[i].Divide(hSigSB[0])
        HistoStyle(hhSigSBratio[i],1,2,21,0,22+2*i);
        hhSigSBratio[i].GetYaxis().SetTitle('ratio');hhSigSBratio[i].SetTitle("Signal-SB ranges: R=0."+str(int(R)))
        if(i==1):
            hhSigSBratio[i].GetYaxis().SetRangeUser(0.4,2)
            #hhSigSBratio[i].GetYaxis().SetRangeUser(-20,5)
            hhSigSBratio[i].Draw()
        elif(i>1 and i not in excludedSigSB):
            hhSigSBratio[i].Draw('same')
        if(i>=1 and i not in excludedSigSB):
            lSigSB1.AddEntry(hhSigSBratio[i],sigSBtitles[i],"l");
    lSigSB1.Draw("same");

    c_SigSB.SaveAs('plots/2_SigSB/2_SigSB_sys'+R+'.pdf')
    c_SigSB.SaveAs('plots/2_SigSB/2_SigSB_sys'+R+'.png')
    #### -----------------------------RMS
    hhSigSBratio_o=[]
    hhSigSBratio_o2=[]
    for i in range(len(hhSigSBratio)):
        if(i in excludedSigSB):continue
        else:hhSigSBratio_o.append(hhSigSBratio[i])
    for i in range(len(hhSigSBratio)):
        if(i in excludedSigSB or i%2!=0):continue
        else:hhSigSBratio_o2.append(hhSigSBratio[i])
    c_SigSBRMS = TCanvas("cSigSBRMS","cSigSBRMS",900,600)
    lSigSB2 = TLegend(0.12,0.80,0.38,0.88);
    if(R=='04' or R=='06'):lSigSB2 = TLegend(0.12,0.715,0.3,0.88);
    histSigSBRMS = rootRMS(hhSigSBratio_o)
    histSigSBRMS2= rootRMS(hhSigSBratio_o2)
    histSigSBRMS.GetYaxis().SetRangeUser(0,0.2);histSigSBRMS.SetLineColor(ROOT.kGreen+2);
    histSigSBRMS2.GetYaxis().SetRangeUser(0,0.2);histSigSBRMS2.SetLineColor(ROOT.kRed+2);
    histSigSBRMS.SetFillColor(ROOT.kGreen+2);histSigSBRMS.SetFillStyle(3654)
    histSigSBRMS2.SetFillColor(ROOT.kRed+2); histSigSBRMS2.SetFillStyle(3645)
    histSigSBRMS.GetYaxis().SetTitle("RMS");histSigSBRMS.Draw();lSigSB2.AddEntry(histSigSBRMS,'all')
    if(len(excludedSigSB)>0):histSigSBRMS2.Draw('same');lSigSB2.AddEntry(histSigSBRMS2,'only 2sig')
    lSigSB2.Draw('same')
    # text
    ttSigSB = ROOT.TText(14,0.18,GetDigitText(histSigSBRMS,fptbinsJC))
    ttSigSB2 = ROOT.TText(14,0.16,GetDigitText(histSigSBRMS2,fptbinsJC))
    ttSigSB.Draw('same')
    if(R=='04' or R=='06'):ttSigSB2.Draw('same')
    c_SigSBRMS.SaveAs('plots/2_SigSB/2_SigSB_sysRMS'+R+'.pdf')
    c_SigSBRMS.SaveAs('plots/2_SigSB/2_SigSB_sysRMS'+R+'.png')

    if(wait):input()
############## -----------------------------
#histos:
# 3. Cutsys
cutindices=['','2','3','5','6']
Cutstitles=['def','cut1','cut2','cut3','cut4']
sizeCutsys=5#
if(flagCutsys):
    datafileCutsys = [datafileFD]
    hCuts=[datafileFD.Get('hData_binned_sub').Clone('hCuts_0')];
    for i in range(1,sizeCutsys):
        datafileCutsys.append(
            RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/CutSys%sFinal_DzeroR%s_paperCuts/Default/FDsubtraction/JetPtSpectrum_FDsub.root"%(cutindices[i],R))
            )
        hCuts.append(datafileCutsys[i].Get('hData_binned_sub').Clone('hCuts_'+str(i)))
    #### -----------------------------
    c_Cuts = TCanvas("cCuts","cCuts",900,600)
    lCuts1 = TLegend(0.12,0.70,0.38,0.88);
    hhCutsratio=[]
    for i in range(sizeCutsys):
        hhCutsratio.append(hCuts[i].Clone('hCutsRatio_'+str(i)))
        hhCutsratio[i].Divide(hCuts[0])
        HistoStyle(hhCutsratio[i],1,2,21,0,i);
        hhCutsratio[i].GetYaxis().SetTitle('ratio');hhCutsratio[i].SetTitle("Topological Cuts: R=0."+str(int(R)))
        hhCutsratio[i].SetLineColor(RTColors[i])
        if(i==0):
            hhCutsratio[i].GetYaxis().SetRangeUser(0.4,2)
            hhCutsratio[i].Draw()
        elif(i>0):
            hhCutsratio[i].Draw('same')
    
        if(i>=0):
            lCuts1.AddEntry(hhCutsratio[i],Cutstitles[i],"l");
    lCuts1.Draw("same");

    c_Cuts.SaveAs('plots/3_Cuts/3_Cuts_ratio'+R+'.pdf')
    c_Cuts.SaveAs('plots/3_Cuts/3_Cuts_ratio'+R+'.png')
    #### -----------------------------RMS
    c_CutsRMS = TCanvas("cCutsRMS","cCutsRMS",900,600)
    histCutsRMS = rootRMS(hhCutsratio)
    histCutsRMS.GetYaxis().SetRangeUser(0,1);histCutsRMS.SetLineColor(ROOT.kBlue+2);
    histCutsRMS.SetFillColor(ROOT.kBlue+2);histCutsRMS.SetFillStyle(3654)
    histCutsRMS.GetYaxis().SetTitle("RMS");histCutsRMS.Draw()
    # extra for R 02
    histCutsRMS2= rootRMS(hhCutsratio[0:1]+hhCutsratio[3:])
    histCutsRMS3 = histCutsRMS.Clone('histCutsRMS3')
    histCutsRMS3.SetBinContent(histCutsRMS3.GetNbinsX(),histCutsRMS2.GetBinContent(histCutsRMS2.GetXaxis().FindBin(40)))

    ttCuts2 = ROOT.TText(33,0.8,GetDigitText(histCutsRMS2,fptbinsJC)[-10:])
    #ttCuts2 = ROOT.TText(5,0.6,GetDigitText(histCutsRMS3,fptbinsJC))
    histCutsRMS3.SetFillColor(ROOT.kRed+2);histCutsRMS3.SetFillStyle(3645)
    if (R=='02'):histCutsRMS3.Draw('same');ttCuts2.Draw('same')

    # text
    ttCuts = ROOT.TText(5,0.9,GetDigitText(histCutsRMS,fptbinsJC))
    ttCuts.Draw('same')

    c_CutsRMS.SaveAs('plots/3_Cuts/3_Cuts_sysRMS'+R+'.pdf')
    c_CutsRMS.SaveAs('plots/3_Cuts/3_Cuts_sysRMS'+R+'.png')

    if(wait):input()
############## -----------------------------
#histos:
# 4. B-FD
if(flagFD):
    hFD_ce = datafileFD.Get('hData_binned_sub').Clone('hFD_ce')
    hFD_up = datafileFD.Get('hData_binned_sub_up').Clone('hFD_up')
    hFD_do = datafileFD.Get('hData_binned_sub_down').Clone('hFD_do')
    #### -----------------------------
    c_FD = TCanvas("cFD","cFD",900,600)
    hFD_upratio = hFD_up.Clone('hFD_upratio')
    hFD_doratio = hFD_do.Clone('hFD_doratio')
    hFD_upratio.Divide(hFD_ce)
    hFD_doratio.Divide(hFD_ce)
    HistoStyle(hFD_upratio,1,2,21,0,ROOT.kRed+2);
    HistoStyle(hFD_doratio,1,2,21,0,ROOT.kGreen+2);
    hFD_upratio.GetYaxis().SetTitle('ratio');hFD_upratio.SetTitle("B feed-down: R=0."+str(int(R)))
    hFD_upratio.Draw()
    hFD_doratio.Draw('same')
    # text
    ttup = ROOT.TText(5,1.5,GetDigitTextFromRatio(hFD_upratio,fptbinsJC))
    ttdo = ROOT.TText(5,0.5,GetDigitTextFromRatio(hFD_doratio,fptbinsJC))
    ttup.Draw('same');ttdo.Draw('same')
    c_FD.SaveAs('plots/4_FD/4_FD_sys'+R+'.pdf')
    c_FD.SaveAs('plots/4_FD/4_FD_sys'+R+'.png')

    if(wait):input()
############## -----------------------------
#histos:
# 5. UnfBayes Iterations
if(flagUIter):
    datafileUIter4=RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_Bayes_4/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R))
    datafileUIter6=RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_Bayes_6/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R))
    datafileUIter = [datafileUnf,datafileUIter4,datafileUIter6]
    hIters=[datafileUIter[0].Get('unfoldedSpectrum').Clone('hIter_5'),datafileUIter[1].Get('unfoldedSpectrum').Clone('hIter_4'),datafileUIter[2].Get('unfoldedSpectrum').Clone('hIter_6')]
    hIterRat0=hIters[0].Clone('hIterRat0');hIterRat0.Divide(hIters[0])
    hIterRat1=hIters[1].Clone('hIterRat1');hIterRat1.Divide(hIters[0])
    hIterRat2=hIters[2].Clone('hIterRat2');hIterRat2.Divide(hIters[0])
    for i in range(hIterRat1.GetNbinsX()): hIterRat1.SetBinError(i+1,0);hIterRat2.SetBinError(i+1,0);
    hIterRat1.SetLineColor(ROOT.kRed+2);hIterRat2.SetLineColor(ROOT.kBlue+2)
    #### -----------------------------
    c_Iter = TCanvas("cIter","cIter",900,500)
    hIterRat1.GetYaxis().SetRangeUser(0.990,1.01);hIterRat1.GetYaxis().SetTitle('ratio to B iter=5')
    hIterRat1.SetTitle("Bayes Iterations: R=0."+str(int(R)))
    leg_Iter=ROOT.TLegend(0.7,0.7,0.85,0.85)
    leg_Iter.AddEntry(hIterRat1,'B iter=4','l')
    leg_Iter.AddEntry(hIterRat2,'B iter=6','l')
    hIterRat1.Draw();hIterRat2.Draw('same')
    lIter = ROOT.TLine(2, 1, 50, 1);lIter.SetLineStyle(2);lIter.Draw("same")
    leg_Iter.Draw('same')
    c_Iter.SaveAs('plots/5_Iter/5_Iter_ratio'+R+'.pdf')
    c_Iter.SaveAs('plots/5_Iter/5_Iter_ratio'+R+'.png')
    #### -----------------------------
    c_IterSys= TCanvas("cIterSys","cIterSys",900,500)
    histIterRMS = rootRMS([hIterRat0,hIterRat1,hIterRat2])
    histIterRMS.GetYaxis().SetRangeUser(0,0.016);histIterRMS.SetLineColor(ROOT.kBlue+2);
    histIterRMS.SetTitle("Bayes Iterations: R=0."+str(int(R)))
    histIterRMS.SetFillColor(ROOT.kBlue+2);histIterRMS.SetFillStyle(3654)
    histIterRMS.GetYaxis().SetTitle("RMS");histIterRMS.Draw()
    ttIter = ROOT.TText(5,0.010,GetDigitText(histIterRMS,fptbinsJC))
    ttIter.Draw('same')
    c_IterSys.SaveAs('plots/5_Iter/5_Iter_sys'+R+'.pdf')
    c_IterSys.SaveAs('plots/5_Iter/5_Iter_sys'+R+'.png')

    if(wait):input()
############## -----------------------------
#histos:
# 6. UnfBayes Priors
sizePrior=10
if(flagUPrior):
    dataUPrior=[datafileUnf]
    hUPrior=[dataUPrior[0].Get('unfoldedSpectrum').Clone('hPriorMC')]
    hUPriorratio=[hUPrior[0].Clone('hratPriorMC')]
    for i in range(1,sizePrior):
        dataUPrior.append(RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_Bayes_5_priorType%d/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R,i-1)))
        hUPrior.append(dataUPrior[i].Get('unfoldedSpectrum').Clone('hPrior'+str(i-1)))
        hUPriorratio.append(dataUPrior[i].Get('unfoldedSpectrum').Clone('hPriorratio'+str(i-1)))
    #### -----------------------------
    c_Prior = TCanvas("cPrior","cPrior",900,500)
    leg_Prior=ROOT.TLegend(0.15,0.55,0.3,0.87)
    hUPriorratio[0].GetYaxis().SetRangeUser(0.94,1.07)
    hUPriorratio[0].GetYaxis().SetTitle('ratio');hUPriorratio[0].SetTitle('Bayes unfolding: priors')
    for i in range(sizePrior):
        hUPriorratio[i].Divide(hUPrior[0])
        HistoStyle(hUPriorratio[i],1,2,21,0,RTColors[i])
        for j in range(hUPriorratio[i].GetNbinsX()+1):hUPriorratio[i].SetBinError(j,0)
        if(i==0):hUPriorratio[i].Draw()
        elif(i==1 or i==6):continue
        else:hUPriorratio[i].Draw('same');leg_Prior.AddEntry(hUPriorratio[i],'Prior: '+str(i-1))
    leg_Prior.Draw('same')
    c_Prior.SaveAs('plots/6_Priors/6_Priors_ratio'+R+'.pdf')
    c_Prior.SaveAs('plots/6_Priors/6_Priors_ratio'+R+'.png')
    #### -----------------------------
    c_PriorRMS = TCanvas("cPriorRMS","cPriorRMS",900,500)
    histPriorRMS=rootRMS(hUPriorratio[0:1]+hUPriorratio[2:6]+hUPriorratio[7:])
    histPriorRMS.GetYaxis().SetRangeUser(0,0.06)
    histPriorRMS.SetFillColor(ROOT.kRed+2);histPriorRMS.SetFillStyle(3354);histPriorRMS.GetYaxis().SetTitle('RMS');histPriorRMS.SetTitle('Bayes unfolding: priors sys.')
    histPriorRMS.Draw()
    ttPrior = ROOT.TText(5,0.04,GetDigitText(histPriorRMS,fptbinsJC))
    ttPrior.Draw('same')
    c_PriorRMS.SaveAs('plots/6_Priors/6_Priors_RMS'+R+'.pdf')
    c_PriorRMS.SaveAs('plots/6_Priors/6_Priors_RMS'+R+'.png')

    if(wait):input()
############## -----------------------------
#histos:
# 7. UnfSvd to Bayes
sizeSVD=4;radiusSVD=0
if R=='04' or R=='06':radiusSVD=1 #d vector for R=0.2,0.3 is 7. for 0.4,0.6 it is 8
if(flagUSvd):
    dataUSvd=[datafileUnf]
    hUSvd=[dataUSvd[0].Get('unfoldedSpectrum').Clone('hBayes')]
    hUSvdratio=[hUSvd[0].Clone('hratSvd0')]
    for i in range(sizeSVD-1): #default bayes already in
        dataUSvd.append(RT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_SVD_%d/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R,i+6+radiusSVD)))
        hUSvd.append(dataUSvd[i+1].Get('unfoldedSpectrum').Clone('hSvd_'+str(i+1)))
        hUSvdratio.append(dataUSvd[i+1].Get('unfoldedSpectrum').Clone('hratSvd'+str(i+1)))
    #### -----------------------------
    c_Svd = TCanvas("cSvd","cSvd",900,500)
    leg_Svd=ROOT.TLegend(0.7,0.7,0.85,0.85)
    hUSvdratio[0].GetYaxis().SetRangeUser(0.6,1.4)
    hUSvdratio[0].GetYaxis().SetTitle('ratio to Bayes 5')
    hUSvdratio[0].SetTitle('SVD unfolding systematic unc.')
    for i in range(sizeSVD):
        hUSvdratio[i].Divide(hUSvd[0])
        HistoStyle(hUSvdratio[i],1,2,21,0,RTColors[i]);
        for j in range(hUSvdratio[i].GetNbinsX()+1):hUSvdratio[i].SetBinError(j,0)
        if(i==0):hUSvdratio[i].Draw()
        else:hUSvdratio[i].Draw('same');leg_Svd.AddEntry(hUSvdratio[i],'SVD: reg='+str(i+5+radiusSVD))
    leg_Svd.Draw('same')
    c_Svd.SaveAs('plots/7_Svd/7_Svd_ratio'+R+'.pdf')
    c_Svd.SaveAs('plots/7_Svd/7_Svd_ratio'+R+'.png')
    #### -----------------------------
    c_SvdSys = TCanvas("cSvdSys","cSvdSys",900,500)
    leg_SvdSys=ROOT.TLegend(0.13,0.7,0.5,0.85)
    histSvdRMS = rootRMS(hUSvdratio);histSvdMean= rootMEAN(hUSvdratio)
    histSvdRMS.GetYaxis().SetRangeUser(0,0.35);histSvdMean.GetYaxis().SetRangeUser(0,0.35)
    histSvdRMS.SetFillColor(ROOT.kRed+2);histSvdRMS.SetFillStyle(3354)
    histSvdMean.SetFillColor(ROOT.kBlue+2);histSvdMean.SetFillStyle(3345);histSvdMean.SetLineColor(ROOT.kBlue+2)
    histSvdRMS.GetYaxis().SetTitle('RMS and Mean deviation')
    histSvdRMS.SetTitle('SVD unfolding systematic unc.')
    histSvdRMS.Draw()
    histSvdMean.Draw('same')
    leg_SvdSys.AddEntry(histSvdRMS,'RMS')
    leg_SvdSys.AddEntry(histSvdMean,'Mean')
    leg_SvdSys.Draw('same')
    ttSvdR = ROOT.TText(14,0.305,GetDigitText(histSvdRMS,fptbinsJC))
    ttSvdM = ROOT.TText(14,0.275,GetDigitText(histSvdMean,fptbinsJC))
    ttSvdR.Draw('same');ttSvdM.Draw('same')
    c_SvdSys.SaveAs('plots/7_Svd/7_Svd_sys'+R+'.pdf')
    c_SvdSys.SaveAs('plots/7_Svd/7_Svd_sys'+R+'.png')

    if(wait):input()
############## -----------------------------
#histos:
# 8. JES
##
if(flagJES):
    hJES_0 = datafileUnf.Get('unfoldedSpectrum').Clone('hJES_0') #default
    hJES_1 = datafileJES.Get('unfoldedSpectrum').Clone('hJES_1') #jes
    #### -----------------------------
    c_JES = TCanvas("cJES","cJES",900,600)
    hJES_0ratio = hJES_0.Clone('hJES_0ratio')
    hJES_1ratio = hJES_1.Clone('hJES_1ratio')
    hJES_1ratio.Divide(hJES_0)
    HistoStyle(hJES_1ratio,1,2,21,0,ROOT.kRed+2);
    hJES_1ratio.GetYaxis().SetTitle('ratio');hJES_1ratio.SetTitle("JES: R=0."+str(int(R)))
    hJES_1ratio.Draw()
    # text
    ttJES = ROOT.TText(5,1.5,GetDigitTextFromRatio(hJES_1ratio,fptbinsJC))
    ttJES.Draw('same');
    c_JES.SaveAs('plots/8_JES/8_JES_rat'+R+'.pdf')
    c_JES.SaveAs('plots/8_JES/8_JES_rat'+R+'.png')
    #### -----------------------------
    c_JESsys = TCanvas("cJESsys","cJESsys",900,600)
    lJES = ROOT.TLine(2, 1, 50, 1);lJES.SetLineStyle(2);lJES.Draw("same")
    hJES_sys = hJES_1ratio.Clone("hJES_sys");hJES_sys.GetYaxis().SetRangeUser(0.95,1.15)
    for i in range(hJES_sys.GetNbinsX()): hJES_sys.SetBinError(i+1,0)
    ttJESsys = ROOT.TText(5,1.05,GetDigitTextFromRatio(hJES_sys,fptbinsJC))
    hJES_sys.Draw();ttJESsys.Draw('same');lJES.Draw("same")
    c_JESsys.SaveAs('plots/8_JES/8_JES_sys'+R+'.pdf')
    c_JESsys.SaveAs('plots/8_JES/8_JES_sys'+R+'.png')

    if(wait):input()
