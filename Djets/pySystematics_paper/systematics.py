import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import ROOT
import array
import numpy
from scipy.optimize import curve_fit

#from matplotlib import colors as mcolors
##----------------------------------------------------------------
#from funcsettings import * # HistoMarkers#from stylesettings import * # HistoMarkers
from funcsettings import funcLine, HistoStyle, rootRMS, GetDigitText, yieldPtrials # HistoMarkers#from stylesettings import * # HistoMarkers
#--------------------------------------------------
if len(sys.argv)==1:
    print("   === Usage example: python file.py R ")
    print("   === e.g.: python file.py 02")
    print("   === OR ===   ")
    print("   === Usage example: python file.py R jetbin")
    print("   === R=02,03,04,06. jetbin=0,1,2,3,4")
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
R = str(sys.argv[1]) #'02'
lensysin = len(sys.argv)
if(lensysin>2): 
    jetbin = str(sys.argv[2]);print(jetbin);
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
ROOTColors = [ROOT.kRed+2, ROOT.kGreen+2, ROOT.kBlue+2, ROOT.kOrange+2, ROOT.kViolet+2, ROOT.kYellow+2, ROOT.kCyan-6, ROOT.kAzure+2, ROOT.kMagenta-6,ROOT.kGreen-8,ROOT.kYellow-8]
datafileFinal = ROOT.TFile("$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_Bayes_5/finalSpectra/JetPtSpectrum_final.root"%(R))
datafileFD = ROOT.TFile("$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/FDsubtraction/JetPtSpectrum_FDsub.root"%(R))
#datafileJES = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/JESsysFinal_DzeroR%s_paperCuts/Default/unfolding_Bayes_5/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R))

datafileUnf = ROOT.TFile("$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_Bayes_5/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R))
if(lensysin>2):
    datafileUnf= ROOT.TFile("$cernbox/media/jackbauer/data/z_out/R_"+R+"_finaltry/unfolding/Bayes/alljetz2D/unfold2DoutFileRegBayesSys5.root","read")

datafileJES = ROOT.TFile("$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/JES/Default/unfolding_Bayes_5/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R))
if(lensysin>2):
    datafileJES = ROOT.TFile("$cernbox/media/jackbauer/data/z_out/R_"+R+"_finaltry/JES/unfolding/Bayes/alljetz2D/unfold2DoutFileRegBayesSys5.root","read")


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
ROOT.gStyle.SetOptStat(0000)
ROOT.gStyle.SetLegendBorderSize(0)
############## -----------------------------
JET_or_Z=''
JET_or_Z_title=''
whichjetbin = ''
whichjetbintitle = whichjetbin
whichJetInZBins=['2_5','5_7','7_10','10_15','15_50']
if(lensysin>2):
    whichJetInZ=int(str(sys.argv[2]))
    JET_or_Z+='Z/'
    whichjetbin = '_'+whichJetInZBins[whichJetInZ]
    whichjetbintitle = whichjetbin[1:]+' GeV'
    JET_or_Z_title='#it{z}_{||,ch.jet}'
    JET_or_Z_text=whichjetbintitle.replace('_', ' < $p_\mathrm{T,ch.jet}$ < ')
else:
    JET_or_Z+='JET/'
    JET_or_Z_title='#it{p}_{T,ch.jet}'
    JET_or_Z_text=''
############## -----------------------------
if(lensysin>2):
    string_SYSTEMATICS_MULTI='$cernbox/media/jackbauer/data/z_out/R_'+R+'_finaltry/RawSys_Multi/signalExtraction/plots/Z0to102_jetbin_'+whichJetInZBins[whichJetInZ]
else:
    string_SYSTEMATICS_MULTI='$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/RawSysFinal_DzeroR'+R+'_paperCuts/Default/signalExtraction_multitrial'
#histos:
# 0. Multi trial
## raw systematics files
if(lensysin>2): sizeMulti=650
else: sizeMulti=973
Multititles=[
        ]
if(flagMulti):
    datafileMulti = [
            ROOT.TFile(string_SYSTEMATICS_MULTI+"/JetPtSpectra_SB_eff.root")
            ]
    hMulti=[
            datafileMulti[0].Get('hjetptspectrumReb').Clone('hMulti_0')
            ];
    for i in range(1,sizeMulti):
        try:
            datafileMulti.append(
                ROOT.TFile(string_SYSTEMATICS_MULTI+"/JetPtSpectra_SB_eff%d.root"%(i))
                )
            hMulti.append(
                datafileMulti[i].Get('hjetptspectrumReb').Clone('hMulti_'+str(i))
                )
        except:
            pass
    #### -----------------------------
    c_MultiRaw = ROOT.TCanvas("cMultiRaw","cMultiRaw",900,600)
    c_MultiRaw.Divide(4,3)
    if(lensysin>2): c_MultiRaw = ROOT.TCanvas("cMultiRaw","cMultiRaw",900,450);c_MultiRaw.Divide(3,2)
    lMulti1 = ROOT.TLegend(0.12,0.60,0.38,0.88);
    hhMultiRaw=[]
    colorindex=30
    #for i in range(sizeMulti):
    for i in range(len(datafileMulti)):
        try:
            hhDptSpec = []
            if(i%20==1): colorindex += 1
            for dbins in range(fDptbins[int(jetbin)]):
                hhdptspecRaw = datafileMulti[i].Get('hjetptsub_'+str(dbins)).Clone('hjetptsub_'+str(dbins))
                hhdptspectra=hhdptspecRaw.Rebin(fptbinsJN,'hjetptsubReb_'+str(dbins),array.array('d',fptbinsJlh)) 
                hhDptSpec.append(hhdptspectra)
                c_MultiRaw.cd(dbins+1)
                #if(lensysin==2): c_MultiRaw.cd(dbins+1).SetLogy()
                hhDptSpec[dbins].SetMinimum(-5)
                hhDptSpec[dbins].SetLineColor(colorindex)#ROOTColors[colorindex])
                if(i==0): hhDptSpec[dbins].Draw() #default
                else: hhDptSpec[dbins].Draw('same')
            hhMultiRaw.append(hhDptSpec)
        except:
            pass

    c_MultiRaw.SaveAs('plots/'+JET_or_Z+'/0_Multi/0_Multi_'+R+whichjetbin+'_sysRaw'+'.pdf')
    c_MultiRaw.SaveAs('plots/'+JET_or_Z+'/0_Multi/0_Multi_'+R+whichjetbin+'_sysRaw'+'.svg')
    #### -----------------------------Distribution of yields
    c_MultiYieldDist = ROOT.TCanvas("cMultiYieldDist","cMultiYieldDist",900,600)
    c_MultiYieldDist.Divide(4,3)
    if(lensysin>2): c_MultiYieldDist = ROOT.TCanvas("cMultiYieldDist","cMultiYieldDist",1000,600); c_MultiYieldDist.Divide(3,2)
    hTrialsDist=[]
    for i in range(fptbinsJN):
        c_MultiYieldDist.cd(i+1)
        c_MultiYieldDist.cd(i+1).SetMargin(2,0,0.9,0.9)
        hTrialsDist.append(yieldDist(hMulti,fptbinsJC[i]))#hTrialsPerPt.append(hYPTValues[0])
        hTrialsDist[i].SetTitle(str(fptbinsJlo[i])+' < '+JET_or_Z_title+' < '+str(fptbinsJhi[i])+' GeV')
        hTrialsDist[i].GetYaxis().SetTitle('Frequency')
        hTrialsDist[i].GetXaxis().SetTitle('Yields of D-jets')
        hTrialsDist[i].Draw()
    #c_MultiYieldDist.cd();padMYD = TPad("","",0,0,1,1.2);padMYD.SetFillStyle(4000);padMYD.Draw();padMYD.cd();latMYD = TLatex();latMYD.SetTextSize(0.04);latMYD.DrawLatexNDC(.4,0.99,"Frequency of raw yields")
    c_MultiYieldDist.SaveAs('plots/'+JET_or_Z+'/0_Multi/0_Multi_'+R+whichjetbin+'_yieldDist'+'.pdf')
    c_MultiYieldDist.SaveAs('plots/'+JET_or_Z+'/0_Multi/0_Multi_'+R+whichjetbin+'_yieldDist'+'.svg')
    #### -----------------------------
    c_MultiYieldTrialsJetpt = ROOT.TCanvas("cMultiYieldTrials","cMultiYieldTrials",900,600)
    c_MultiYieldTrialsJetpt.Divide(4,3)
    if(lensysin>2): c_MultiYieldTrialsJetpt = ROOT.TCanvas("cMultiYieldTrials","cMultiYieldTrials",1000,600); c_MultiYieldTrialsJetpt.Divide(3,2)
    hTrialsPerPt = []

    meanYPT,line_YPT=[],[]
    sigYPTplus,line_YPTplus,sigUp=[],[],1.5
    sigYPTminu,line_YPTminu,sigDo=[],[],1.5
    for i in range(fptbinsJN):
        c_MultiYieldTrialsJetpt.cd(i+1)
        c_MultiYieldTrialsJetpt.cd(i+1).SetMargin(2,0,0.9,0.9)
        hTrialsPerPt.append(yieldPtrials(hMulti,fptbinsJC[i]))#hTrialsPerPt.append(hYPTValues[0])
        hTrialsPerPt[i].SetTitle(str(fptbinsJlo[i])+' < '+JET_or_Z_title+' < '+str(fptbinsJhi[i])+' GeV')
        hTrialsPerPt[i].Draw()
        
        c_MultiYieldTrialsJetpt.cd(i+1).Update()
        meanYPT.append(hMulti[0].GetBinContent(hMulti[0].GetXaxis().FindBin(fptbinsJC[i]))) 
        sigYPTplus.append(hMulti[0].GetBinContent(hMulti[0].GetXaxis().FindBin(fptbinsJC[i])) + sigUp*hMulti[0].GetBinError(hMulti[0].GetXaxis().FindBin(fptbinsJC[i])) )
        sigYPTminu.append(hMulti[0].GetBinContent(hMulti[0].GetXaxis().FindBin(fptbinsJC[i])) - sigDo*hMulti[0].GetBinError(hMulti[0].GetXaxis().FindBin(fptbinsJC[i])) )
        line_YPT.append(ROOT.TLine(0,meanYPT[i],len(datafileMulti), meanYPT[i]))
        line_YPTplus.append(ROOT.TLine(0,sigYPTplus[i],len(datafileMulti), sigYPTplus[i]))
        line_YPTminu.append(ROOT.TLine(0,max(0,sigYPTminu[i]),len(datafileMulti), max(0,sigYPTminu[i])))
        line_YPT[i].SetLineColor(ROOT.kRed+2)
        line_YPTplus[i].SetLineColor(ROOT.kRed+2)
        line_YPTminu[i].SetLineColor(ROOT.kRed+2)
        c_MultiYieldTrialsJetpt.cd(i+1)
        line_YPT[i].Draw('same')
        line_YPTplus[i].Draw('same')
        line_YPTminu[i].Draw('same')

    c_MultiYieldTrialsJetpt.cd(fptbinsJN+1);
    ttMultiYPT2 = [ROOT.TLatex(.12,.8,'R = 0.'+str(int(R))+', Multi-trial')
            ,ROOT.TLatex(.12,.7,'pp, 5.02 TeV')
            ,ROOT.TLatex(.12,.6,''+whichjetbintitle)
            ]
    for i in range(len(ttMultiYPT2)):ttMultiYPT2[i].SetTextSize(0.1);ttMultiYPT2[i].SetTextFont(42);ttMultiYPT2[i].Draw('same')
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
    
    #c_MultiYieldTrialsJetpt.cd(fptbinsJN+2)
    #ttMultiYPT2 = [ROOT.TLatex(.12,.8,'R = 0.'+str(int(R))+', Multi-trial')
    #        ,ROOT.TLatex(.12,.7,'pp, 5.02 TeV')
    #        ,ROOT.TLatex(.12,.6,''+whichJetInZBins[whichJetInZ]+' GeV')
    #        ]
    #for i in range(len(ttMultiYPT2)):ttMultiYPT2[i].SetTextSize(0.1);ttMultiYPT2[i].SetTextFont(42);ttMultiYPT2[i].Draw('same')

    c_MultiYieldTrialsJetpt.SaveAs('plots/'+JET_or_Z+'/0_Multi/0_Multi_'+R+whichjetbin+'_YieldTrials'+'.pdf')
    c_MultiYieldTrialsJetpt.SaveAs('plots/'+JET_or_Z+'/0_Multi/0_Multi_'+R+whichjetbin+'_YieldTrials'+'.svg')
    #### -----------------------------
    c_MultiRatio = ROOT.TCanvas("cMultiRatio","cMultiRatio",900,600)
    lMulti3 = ROOT.TLegend(0.12,0.60,0.38,0.88);
    hhMultiratio=[]
    for i in range(len(datafileMulti)):
        try:
            hhMultiratio.append(hMulti[i].Clone('hMultiratio_'+str(i)))
            hhMultiratio[i].Divide(hMulti[0])
            HistoStyle(hhMultiratio[i],1,2,21,0,ROOTColors[0]);
            hhMultiratio[i].GetYaxis().SetTitle('ratio');hhMultiratio[i].SetTitle("Multi-trial: R=0."+str(int(R))+' '+whichjetbintitle)
            hhMultiratio[i].SetLineColor((i%20)+29)
            hhMultiratio[i].GetYaxis().SetRangeUser(0,2)
            if(i==1): hhMultiratio[i].Draw(); 
            elif(i>1): hhMultiratio[i].Draw('same')
        except: pass
    c_MultiRatio.SaveAs('plots/'+JET_or_Z+'/0_Multi/0_Multi_'+R+whichjetbin+'_sysRatio'+'.pdf')
    c_MultiRatio.SaveAs('plots/'+JET_or_Z+'/0_Multi/0_Multi_'+R+whichjetbin+'_sysRatio'+'.svg')
    #### -----------------------------RMS
    c_MultiRMS = ROOT.TCanvas("cMultiRMS","cMultiRMS",900,600)
    lMulti4 = ROOT.TLegend(0.12,0.60,0.38,0.88);
    histMultiRMS = rootRMS_mt(hhMultiratio,sigUp,sigDo);
    histMultiRMS.SetLineColor(ROOT.kRed+2);histMultiRMS.SetFillColor(ROOT.kRed+2);histMultiRMS.SetFillStyle(3354);
    histMultiRMS.GetYaxis().SetRangeUser(0,0.25);histMultiRMS.GetYaxis().SetTitle("RMS");
    histMultiRMS.Draw()
    #lMulti4.AddEntry(histMultiRMS,'');lMulti4.Draw('same')

    # text
    ytextplace=0.22
    ttMultiRMS= ROOT.TText(15,ytextplace,GetDigitText(histMultiRMS,fptbinsJC));ttMultiRMS.SetTextSize(0.04);
    if(lensysin>2):ttMultiRMS= ROOT.TText(0.6,0.2,GetDigitText(histMultiRMS,fptbinsJC));ttMultiRMS.SetTextSize(0.04);
    ttMultiRMS.Draw('same')
    c_MultiRMS.SaveAs('plots/'+JET_or_Z+'/0_Multi/0_Multi_'+R+whichjetbin+'_sysRMS'+'.pdf')
    c_MultiRMS.SaveAs('plots/'+JET_or_Z+'/0_Multi/0_Multi_'+R+whichjetbin+'_sysRMS'+'.svg')

    if(wait):input()
############## -----------------------------
#histos:
# 9. Closure test
## 
if(flagClos):
    size_closure=10
    dataClosure=[]
    hhClosure=[]
    #### -----------------------------
    c_close = ROOT.TCanvas("cClose","cClose",900,600)
    lClose = ROOT.TLegend(0.12,0.80,0.38,0.88);
    #### -----------------------------
    for i in range(size_closure):
        if(lensysin>2):
            dataClosure.append(ROOT.TFile("$cernbox/media/jackbauer/data/z_out/R_"+R+"_finaltry/unfolding/Bayes/alljetz2D/MCClosure_"+str(i)+".root","read"))
            hhClosure.append(dataClosure[i].Get('unfByTruth'+str(int(jetbin)+1)).Clone('hRawVTrue'+str(i)))
        else:
            dataClosure.append(ROOT.TFile("$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_Bayes_5/plots/Closure/unfoldedSpectrum_closure%d.root"%(R,i),"read"))
            hhClosure.append(dataClosure[i].Get('hRawVTrue').Clone('hRawVTrue'+str(i)))
        hhClosure[i].GetYaxis().SetTitle("ratio")
        hhClosure[i].GetYaxis().SetRangeUser(0,2)
        if(i==0):hhClosure[i].Draw()
        else:hhClosure[i].Draw("same")
    c_close.SaveAs('plots/'+JET_or_Z+'9_Closure/9_Close'+R+whichjetbin+'.pdf')
    c_close.SaveAs('plots/'+JET_or_Z+'9_Closure/9_Close'+R+whichjetbin+'.png')
    #### -----------------------------RMS
    c_CloRMS = ROOT.TCanvas("cCloRMS","cCloRMS",900,600)
    lCloRMS = ROOT.TLegend(0.12,0.80,0.38,0.88);
    #skip = 0;histCloRMS = ratioRMS(hhClosure[:skip]+hhClosure[skip+1:])
    histCloRMS = ratioRMS(hhClosure)
    histCloRMS.GetYaxis().SetRangeUser(0,0.5);histCloRMS.SetLineColor(ROOT.kGreen+2);
    histCloRMS.SetFillColor(ROOT.kGreen+2);histCloRMS.SetFillStyle(3654)
    histCloRMS.GetYaxis().SetTitle("RMS");histCloRMS.Draw();
    # text
    ttCloRMS = ROOT.TText(fptbinsJlo[1],0.4,GetDigitText(histCloRMS,fptbinsJC))
    ttCloRMS.Draw('same')
    c_CloRMS.SaveAs('plots/'+JET_or_Z+'9_Closure/9_CloRMS'+R+whichjetbin+'.pdf')
    c_CloRMS.SaveAs('plots/'+JET_or_Z+'9_Closure/9_CloRMS'+R+whichjetbin+'.png')
    #### -----------------------------Mean
    c_CloMean = ROOT.TCanvas("cCloMean","cCloMean",900,600)
    lCloMean = ROOT.TLegend(0.12,0.80,0.38,0.88);
    histCloMean = ratioMEAN(hhClosure)
    histCloMean.GetYaxis().SetRangeUser(0,0.5);histCloMean.SetLineColor(ROOT.kGreen+2);
    histCloMean.SetFillColor(ROOT.kGreen+2);histCloMean.SetFillStyle(3654)
    histCloMean.GetYaxis().SetTitle("Mean");histCloMean.Draw();
    # text
    ttCloMean = ROOT.TText(fptbinsJlo[1],0.4,GetDigitText(histCloMean,fptbinsJC))
    ttCloMean.Draw('same')
    c_CloMean.SaveAs('plots/'+JET_or_Z+'9_Closure/9_CloMean'+R+whichjetbin+'.pdf')
    c_CloMean.SaveAs('plots/'+JET_or_Z+'9_Closure/9_CloMean'+R+whichjetbin+'.png')

    #### -----------------------------
    if(wait):input()
############## -----------------------------
############## -----------------------------
if(lensysin>2):
    string_SYSTEMATICS_REF="$cernbox/media/jackbauer/data/z_out/R_"+R+"_finaltry/signalExtraction/plots/Z0to102_jetbin_"+whichJetInZBins[whichJetInZ]
    datafileSig = ROOT.TFile(string_SYSTEMATICS_REF+"/JetPtSpectra_SB_eff.root","read")
else:
    string_SYSTEMATICS_REF="$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/RawSysFinal_DzeroR"+R+"_paperCuts/Default/signalExtraction_refsys"
    datafileSig = ROOT.TFile("$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/signalExtraction/JetPtSpectra_SB_eff.root"%(R))
#histos:
# 1. Reflection
## 
if(flagRef):
    datafileRf0 = ROOT.TFile(string_SYSTEMATICS_REF+"/JetPtSpectra_SB_effrefSys0.root","read")
    datafileRf1 = ROOT.TFile(string_SYSTEMATICS_REF+"/JetPtSpectra_SB_effrefSys1.root","read")
    hRdf = datafileSig.Get('hjetptspectrumRebScaled').Clone('hRdf')
    hRf0 = datafileRf0.Get('hjetptspectrumRebScaled').Clone('hRf0')
    hRf1 = datafileRf1.Get('hjetptspectrumRebScaled').Clone('hRf1')
    #### -----------------------------
    c_Ref = ROOT.TCanvas("cRef","cRef",900,600)
    leg_Ref = ROOT.TLegend(0.7,0.7,0.88,0.88)
    hRf0ratio=hRf0.Clone('hRf0ratio');hRf0ratio.Divide(hRdf);HistoStyle(hRf0ratio,1,2,21,0,ROOT.kRed+2);
    hRf1ratio=hRf1.Clone('hRf1ratio');hRf1ratio.Divide(hRdf);HistoStyle(hRf1ratio,1,2,21,0,ROOT.kGreen+2);
    for i in range(hRf0ratio.GetNbinsX()): hRf0ratio.SetBinError(i+1,0);hRf1ratio.SetBinError(i+1,0)
    hRf0ratio.GetYaxis().SetRangeUser(0.9,1.1)
    hRf0ratio.SetTitle("R=0."+str(int(R))+": Reflections")
    hRf0ratio.GetYaxis().SetTitle("Ratio")
    hRf0ratio.Draw();hRf1ratio.Draw('same')
    leg_Ref.AddEntry(hRf0ratio,"refl,-50%","l");leg_Ref.AddEntry(hRf1ratio,"refl,+50%","l");leg_Ref.Draw('same')
    lRef = ROOT.TLine(fptbinsJlh[0], 1, fptbinsJlh[-1], 1);lRef.SetLineStyle(2);lRef.Draw("same")
    ttRefsys0 = ROOT.TText(fptbinsJC[0],1.05,GetDigitTextFromRatio(hRf0ratio,fptbinsJC))
    ttRefsys1 = ROOT.TText(fptbinsJC[0],0.95,GetDigitTextFromRatio(hRf1ratio,fptbinsJC))
    ttRefsys0.Draw('same');ttRefsys1.Draw('same')
    c_Ref.SaveAs('plots/'+JET_or_Z+'1_Ref/1_Ref_rat'+R+whichjetbin+'.pdf')
    c_Ref.SaveAs('plots/'+JET_or_Z+'1_Ref/1_Ref_rat'+R+whichjetbin+'.png')

    if(wait):input()
############## -----------------------------
if(lensysin>2):
    string_SYSTEMATICS_SIGSB='$cernbox/media/jackbauer/data/z_out/R_'+R+'_finaltry/RawSys_SBSig/signalExtraction/plots/Z0to102_jetbin_'+whichJetInZBins[whichJetInZ]
else:
    #string_SYSTEMATICS_SIGSB="/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/RawSysFinal_DzeroR"+R+"_paperCuts/Default/signalExtraction_SigSBranges/"
    string_SYSTEMATICS_SIGSB="$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR"+R+"_paperCuts/RawSys_SBSig/Default/signalExtraction/"
#histos:
# 2. Signal and SB ranges
## raw systematics files
if(lensysin<=3): # for both jetpt and zch
    sigSBtitles=[
        #'S:2\sigma, SB:4.0-9\sigma',#default, not necessary
        'S:2\sigma, SB:3.5-9\sigma',
        'S:2\sigma, SB:4.5-9\sigma',
        'S:2\sigma, SB:4.0-8\sigma',
        'S:2\sigma, SB:3.5-8\sigma',
        'S:2\sigma, SB:4.5-8\sigma',
        'S:3\sigma, SB:4.0-9\sigma',
        'S:3\sigma, SB:3.5-9\sigma',
        'S:3\sigma, SB:4.5-9\sigma',
        'S:3\sigma, SB:4.0-8\sigma',
        'S:3\sigma, SB:3.5-8\sigma',
        'S:3\sigma, SB:4.5-8\sigma',
        ]

if(flagSigSB):
    datafileSigSB = [
            ROOT.TFile(string_SYSTEMATICS_SIGSB+"/JetPtSpectra_SB_eff.root")
            ]
    hSigSB=[datafileSigSB[0].Get('hjetptspectrumRebScaled').Clone('hSigSB_0')];
    for i in range(1,len(sigSBtitles)+1):
        try:
            datafileSigSB.append( ROOT.TFile(string_SYSTEMATICS_SIGSB+"/JetPtSpectra_SB_eff%d.root"%(i)) )
            hSigSB.append( datafileSigSB[i].Get('hjetptspectrumRebScaled').Clone('hSigSB_'+str(i)) )
        except: pass
    #### -----------------------------
    c_SigSB = ROOT.TCanvas("cSigSB","cSigSB",900,600)
    lSigSB1 = ROOT.TLegend(0.12,0.50,0.38,0.88);
    hhSigSBratio=[]
    excludedSigSB=[]#6,7,8,9,10,11]
    excludedSigSB=[2,5,8,11]
    for i in range(len(hSigSB)):
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
            lSigSB1.AddEntry(hhSigSBratio[i],sigSBtitles[i-1],"l");
    lSigSB1.Draw("same");

    c_SigSB.SaveAs('plots/'+JET_or_Z+'2_SigSB/2_SigSB_sys'+R+whichjetbin+'.pdf')
    c_SigSB.SaveAs('plots/'+JET_or_Z+'2_SigSB/2_SigSB_sys'+R+whichjetbin+'.png')
    #### -----------------------------RMS
    hhSigSBratio_o=[]
    for i in range(len(hhSigSBratio)):
        if(i in excludedSigSB):continue
        else:hhSigSBratio_o.append(hhSigSBratio[i])
    c_SigSBRMS = ROOT.TCanvas("cSigSBRMS","cSigSBRMS",900,600)
    lSigSB2 = ROOT.TLegend(0.12,0.80,0.38,0.88);
    histSigSBRMS = rootRMS(hhSigSBratio_o)
    histSigSBRMS.GetYaxis().SetRangeUser(0,0.6);histSigSBRMS.SetLineColor(ROOT.kGreen+2);
    histSigSBRMS.SetFillColor(ROOT.kGreen+2);histSigSBRMS.SetFillStyle(3654)
    histSigSBRMS.GetYaxis().SetTitle("RMS");histSigSBRMS.Draw();lSigSB2.AddEntry(histSigSBRMS,'all')
    #lSigSB2.Draw('same')
    # text
    ttSigSB = ROOT.TText(8,0.18,GetDigitText(histSigSBRMS,fptbinsJC))
    if(lensysin>=3): ttSigSB = ROOT.TText(0.6,0.3,GetDigitText(histSigSBRMS,fptbinsJC))
    ttSigSB.Draw('same')
    c_SigSBRMS.SaveAs('plots/'+JET_or_Z+'2_SigSB/2_SigSB_sysRMS'+R+whichjetbin+'.pdf')
    c_SigSBRMS.SaveAs('plots/'+JET_or_Z+'2_SigSB/2_SigSB_sysRMS'+R+whichjetbin+'.png')

    if(wait):input()
############## -----------------------------
#histos:
# 3. Cutsys
cutindices=['','2','3','5','6']
Cutstitles=['def','cut1','cut2','cut3','cut4']
sizeCutsys=5#
if(flagCutsys):
    if(lensysin>2):
        datafileCutsys=[ROOT.TFile("$cernbox/media/jackbauer/data/z_out/R_"+R+"_finaltry/FDsubtraction/outFD_1.root ","read")]
        hCuts=[datafileCutsys[0].Get('hsub_c'+str(int(jetbin)+1)).Clone('hCuts_0')];
        for i in range(1,sizeCutsys):
            datafileCutsys.append(ROOT.TFile("$cernbox/media/jackbauer/data/z_out/R_"+R+"_finaltry/SQ"+cutindices[i]+"/FDsubtraction/outFD_1.root ","read"))
            hCuts.append(datafileCutsys[i].Get('hsub_c'+str(int(jetbin)+1)).Clone('hCuts_'+str(i)))
        print(hCuts)
    else:
        datafileCutsys = [datafileFD]
        hCuts=[datafileFD.Get('hData_binned_sub').Clone('hCuts_0')];
        for i in range(1,sizeCutsys):
            datafileCutsys.append(
                    #ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/CutSys%sFinal_DzeroR%s_paperCuts/Default/FDsubtraction/JetPtSpectrum_FDsub.root"%(cutindices[i],R))
                    ROOT.TFile("$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR"+R+"_paperCuts/SQ"+cutindices[i]+"/Default/FDsubtraction/JetPtSpectrum_FDsub.root")
                    )
            hCuts.append(datafileCutsys[i].Get('hData_binned_sub').Clone('hCuts_'+str(i)))
    #### -----------------------------
    c_Cuts = ROOT.TCanvas("cCuts","cCuts",900,600)
    lCuts1 = ROOT.TLegend(0.12,0.70,0.38,0.88);
    if(lensysin>2):lCuts1 = ROOT.TLegend(0.72,0.70,0.88,0.88);
    hhCutsratio=[]
    for i in range(sizeCutsys):
        hhCutsratio.append(hCuts[i].Clone('hCutsRatio_'+str(i)))
        hhCutsratio[i].Divide(hCuts[0])
        HistoStyle(hhCutsratio[i],1,2,21,0,i);
        hhCutsratio[i].GetYaxis().SetTitle('ratio');hhCutsratio[i].SetTitle("Topological Cuts: R=0."+str(int(R)))
        hhCutsratio[i].SetLineColor(ROOTColors[i])
        if(i==0):
            hhCutsratio[i].GetYaxis().SetRangeUser(0.4,2)
            if(lensysin>2):hhCutsratio[i].GetYaxis().SetRangeUser(0.4,3)
            hhCutsratio[i].Draw()
        else:
            hhCutsratio[i].Draw('same')
        lCuts1.AddEntry(hhCutsratio[i],Cutstitles[i],"l");
    lCuts1.Draw("same");
    print('errors:')
    for i in range(hhCutsratio[1].GetNbinsX()):
        print(hhCutsratio[1].GetBinContent(i+1), ':', hhCutsratio[1].GetBinError(i+1)/hhCutsratio[1].GetBinContent(i+1))

    #c_Cuts.SaveAs('plots/3_Cuts/3_Cuts_ratio'+R+'.pdf')
    #c_Cuts.SaveAs('plots/3_Cuts/3_Cuts_ratio'+R+'.png')
    c_Cuts.SaveAs('plots/'+JET_or_Z+'3_Cuts/3_Cuts_ratio'+R+whichjetbin+'.pdf')
    c_Cuts.SaveAs('plots/'+JET_or_Z+'3_Cuts/3_Cuts_ratio'+R+whichjetbin+'.png')
    #### -----------------------------RMS
    c_CutsRMS = ROOT.TCanvas("cCutsRMS","cCutsRMS",900,600)
    l_CutsRMS = ROOT.TLegend(0.12,0.70,0.38,0.88);
    histCutsRMS = rootRMS(hhCutsratio)
    histCutsRMS.GetYaxis().SetRangeUser(0,1);histCutsRMS.SetLineColor(ROOT.kBlue+2);
    histCutsRMS.SetFillColor(ROOT.kBlue+2);histCutsRMS.SetFillStyle(3654)
    histCutsRMS.GetYaxis().SetTitle("RMS      ");
    histCutsRMS.SetTitle("")
    histCutsRMS.Draw()
    # extra for R 02
    #histCutsRMS2= rootRMS(hhCutsratio[0:1]+hhCutsratio[3:])
    if(lensysin==2):
        histCutsRMS2 = rootRMS(hhCutsratio,0.2)
    elif(lensysin>2):
        histCutsRMS2 = rootRMS(hhCutsratio,0.6)
    histCutsRMS3 = histCutsRMS.Clone('histCutsRMS3')


    ttCuts2 = ROOT.TText(fptbinsJlh[-2],0.8,GetDigitText(histCutsRMS2,fptbinsJC))#[-10:])
    #ttCuts2 = ROOT.TText(5,0.6,GetDigitText(histCutsRMS3,fptbinsJC))
    histCutsRMS3.SetFillColor(ROOT.kRed+2);histCutsRMS3.SetFillStyle(3645)
    if (lensysin==2 and R=='02'):# or (lensysin>2 and R=='02'):
        histCutsRMS3.SetBinContent(histCutsRMS3.GetNbinsX(),histCutsRMS2.GetBinContent(histCutsRMS2.GetXaxis().FindBin(40)))
        histCutsRMS3.Draw('same');ttCuts2.Draw('same')
    elif (lensysin>2 and R=='02'):
        histCutsRMS3.SetBinContent(1,histCutsRMS2.GetBinContent(histCutsRMS2.GetXaxis().FindBin(0.5)))
        histCutsRMS3.Draw('same');
        if(jetbin=='2'):
            ttCuts2=ROOT.TText(fptbinsJlh[1],0.7,"w/o cut2")
            ttCuts2.Draw('same') 
        elif(jetbin=='4'):
            ttCuts2 = ROOT.TText(fptbinsJlo[1],0.7,GetDigitText(histCutsRMS2,fptbinsJC)+" w/o cut2")#[-10:])
            #ttCuts2=ROOT.TText(fptbinsJlh[1],0.7,"w/o cut2")
            ttCuts2.Draw('same') 

    # text
    ttCuts = ROOT.TText(fptbinsJlo[1],0.9,GetDigitText(histCutsRMS,fptbinsJC))
    ttCuts.Draw('same')

    #c_CutsRMS.SaveAs('plots/3_Cuts/3_Cuts_sysRMS'+R+'.pdf')
    #c_CutsRMS.SaveAs('plots/3_Cuts/3_Cuts_sysRMS'+R+'.png')
    c_CutsRMS.SaveAs('plots/'+JET_or_Z+'3_Cuts/3_Cuts_sysRMS'+R+whichjetbin+'.pdf')
    c_CutsRMS.SaveAs('plots/'+JET_or_Z+'3_Cuts/3_Cuts_sysRMS'+R+whichjetbin+'.png')

    ##FITTING plot
    yCUTS = array.array('d',histCutsRMS)[1:-1]
    yCUTS2= array.array('d',histCutsRMS)[1:-1]
    xCUTS = array.array('d',fptbinsJC)
    if(lensysin>2):
        print('xCUTS',xCUTS)
        print('yCUTS',yCUTS)
        if(R=='02'):
            if(whichJetInZ==2 or 4):
                xCUTS=xCUTS[1:]
                yCUTS=yCUTS[1:]
                print('this works!!!!!!!!')
    else:
        yCUTS = yCUTS[3:-1]
        xCUTS = xCUTS[3:-1]
    xCUTSlh= numpy.array(array.array('d',fptbinsJlh[:]))
    print(xCUTS, yCUTS)
    funcCUTS = funcLine  
    #funcCUTS = funcBola
    print('curve_fit parameters:',funcCUTS, xCUTS, yCUTS)
    try:
        popt, pcov = curve_fit(funcCUTS, numpy.array(xCUTS), numpy.array(yCUTS))
        print('popt:',popt)
        print('pcov:',pcov)
    except RuntimeError:
        print("Error")

    fig=plt.figure(1);
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.15)
    plt.grid(axis='both',color='0.95');
    plt.ylim(0,0.5);#plt.xticks(np.arange(min(xJESall),max(xJESall)+5,5.0));
    plt.xlim(fptbinsJlh[0],fptbinsJlh[-1])
    plt.plot(xCUTSlh, funcCUTS(xCUTSlh,*popt),'r-', label='fit');
    plt.plot(fptbinsJlo[0:2],[array.array('d',histCutsRMS)[0],array.array('d',histCutsRMS)[0]],'k')
    if(lensysin==2 and R=='02'):
        plt.step(fptbinsJlh[:], array.array('d',histCutsRMS2),'k', label="test");
    else:
        plt.step(fptbinsJlh[:], array.array('d',histCutsRMS)[:-1],'k', label="test");
    plt.errorbar(fptbinsJC,numpy.array(array.array('d',histCutsRMS)[1:-1]),fmt='.k');
    if(lensysin==2):
        plt.xlabel(r'jet $p_T$', fontsize=12)
        locCUTS = 0.4
        locCUTS_x = 3
    elif(lensysin>2):
        plt.xlabel(r'$z_{||}^{ch}$', fontsize=12)
        locCUTS = 0.47
        locCUTS_x = 0.42
    plt.ylabel('RMS', fontsize=12);
    #plt.title("Topological selections, R = 0."+str(int(R)))
    uncCUTSvals = ""
    valsCUTS = abs(100*(funcCUTS(numpy.array(array.array('d',fptbinsJC)),*popt))) #100*(valsJES[i]-1)
    for i in range(len(valsCUTS)):
        uncCUTSvals += "%.1f, "%(valsCUTS[i]) #"R%s"%(R))
    uncCUTSvals = uncCUTSvals[:-2]+' in %'
    text_title = 'Topological selections, $R$=0.'+str(int(R))
    text_title +='\n'+JET_or_Z_text
    text_title +='\n'+uncCUTSvals
    plt.text(locCUTS_x, locCUTS, text_title,#uncCUTSvals,
         #rotation=45,
         horizontalalignment='left',
         verticalalignment='top',
         #multialignment='center'
         fontsize=12,
         linespacing = 1.5,
         )
    plt.draw()
    plt.savefig('plots/'+JET_or_Z+'3_Cuts/3_Cut_FIT_2'+R+whichjetbin+'.pdf')
    plt.savefig('plots/'+JET_or_Z+'3_Cuts/3_Cut_FIT_2'+R+whichjetbin+'.png')
    #plt.waitforbuttonpress(1);input()
    plt.close()
    #if(wait):input()
############## -----------------------------
#histos:
# 4. B-FD
if(flagFD):
    if(lensysin>2):
        xtitle = "z_{||}"
        textpos = 0.5
        datafileFD=ROOT.TFile("$cernbox/media/jackbauer/data/z_out/R_"+R+"_finaltry/FDsubtraction/outFD_1.root","read")
        hFD_ce = datafileFD.Get('hsub_c'+str(int(jetbin)+1)).Clone('hFD_ce')
        hFD_up = datafileFD.Get('hsub_u'+str(int(jetbin)+1)).Clone('hFD_up')
        hFD_do = datafileFD.Get('hsub_d'+str(int(jetbin)+1)).Clone('hFD_do')
    else:
        textpos=5
        xtitle = "jet-p_T"
        hFD_ce = datafileFD.Get('hData_binned_sub').Clone('hFD_ce')
        hFD_up = datafileFD.Get('hData_binned_sub_up').Clone('hFD_up')
        hFD_do = datafileFD.Get('hData_binned_sub_down').Clone('hFD_do')
    #### -----------------------------
    c_FD = ROOT.TCanvas("cFD","cFD",900,600)
    hFD_upratio = hFD_up.Clone('hFD_upratio')
    hFD_doratio = hFD_do.Clone('hFD_doratio')
    hFD_upratio.Divide(hFD_ce)
    hFD_doratio.Divide(hFD_ce)
    HistoStyle(hFD_upratio,1,2,21,0,ROOT.kRed+2);
    HistoStyle(hFD_doratio,1,2,21,0,ROOT.kGreen+2);
    hFD_upratio.GetYaxis().SetTitle('ratio');hFD_upratio.SetTitle("B feed-down: R=0."+str(int(R)))
    hFD_upratio.GetYaxis().SetRangeUser(0,2);
    if(lensysin>2):hFD_upratio.GetXaxis().SetTitle(xtitle)
    hFD_upratio.Draw()
    hFD_doratio.Draw('same')
    # text
    ttup = ROOT.TText(textpos,1.5,GetDigitTextFromRatio(hFD_upratio,fptbinsJC))
    ttdo = ROOT.TText(textpos,0.5,GetDigitTextFromRatio(hFD_doratio,fptbinsJC))
    ttup.Draw('same');ttdo.Draw('same')
    #c_FD.SaveAs('plots/4_FD/4_FD_sys'+R+'.pdf')
    #c_FD.SaveAs('plots/4_FD/4_FD_sys'+R+'.png')
    c_FD.SaveAs('plots/'+JET_or_Z+'4_FD/4_FD_sys'+R+whichjetbin+'.pdf')
    c_FD.SaveAs('plots/'+JET_or_Z+'4_FD/4_FD_sys'+R+whichjetbin+'.png')

    if(wait):input()
############## -----------------------------
#histos:
# 5. UnfBayes Iterations
if(flagUIter):
    if(lensysin>2):
        datafileUIter4 = ROOT.TFile("$cernbox/media/jackbauer/data/z_out/R_"+R+"_finaltry/unfolding/Bayes/alljetz2D/unfold2DoutFileRegBayesSys4.root","read")
        datafileUIter6= ROOT.TFile("$cernbox/media/jackbauer/data/z_out/R_"+R+"_finaltry/unfolding/Bayes/alljetz2D/unfold2DoutFileRegBayesSys6.root","read")
        datafileUIter5= ROOT.TFile("$cernbox/media/jackbauer/data/z_out/R_"+R+"_finaltry/unfolding/Bayes/alljetz2D/unfold2DoutFileRegBayesSys5.root","read")
        datafileUIter = [datafileUIter5,datafileUIter4,datafileUIter6]
        hIters=[datafileUIter[0].Get('UnfProjectX_'+jetbin).Clone('hIter_5'),datafileUIter[1].Get('UnfProjectX_'+jetbin).Clone('hIter_4'),datafileUIter[2].Get('UnfProjectX_'+jetbin).Clone('hIter_6')]
    else:
        #datafileUIter4=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_Bayes_4/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R))
        #datafileUIter6=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_Bayes_6/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R))
        #datafileUIter = [datafileUnf,datafileUIter4,datafileUIter6]
        #hIters=[datafileUIter[0].Get('unfoldedSpectrum').Clone('hIter_5'),datafileUIter[1].Get('unfoldedSpectrum').Clone('hIter_4'),datafileUIter[2].Get('unfoldedSpectrum').Clone('hIter_6')]
        datafileUIter = [datafileUnf,datafileUnf,datafileUnf]
        hIters=[datafileUIter[0].Get('unfoldedSpectrum').Clone('hIter_5'),datafileUIter[1].Get('unfoldedSpectrum_minusone').Clone('hIter_4'),datafileUIter[2].Get('unfoldedSpectrum_plusone').Clone('hIter_6')]
    hIterRat0=hIters[0].Clone('hIterRat0');hIterRat0.Divide(hIters[0])
    hIterRat1=hIters[1].Clone('hIterRat1');hIterRat1.Divide(hIters[0])
    hIterRat2=hIters[2].Clone('hIterRat2');hIterRat2.Divide(hIters[0])
    for i in range(hIterRat1.GetNbinsX()): hIterRat1.SetBinError(i+1,0);hIterRat2.SetBinError(i+1,0);
    hIterRat1.SetLineColor(ROOT.kRed+2);hIterRat2.SetLineColor(ROOT.kBlue+2)
    #### -----------------------------
    c_Iter = ROOT.TCanvas("cIter","cIter",900,500)
    if(lensysin>2):
        lower_limit=0.9;upper_limit=1.1
    else:
        lower_limit=0.99;upper_limit=1.02
    hIterRat1.GetYaxis().SetRangeUser(lower_limit,upper_limit);hIterRat1.GetYaxis().SetTitle('ratio to B iter=5')
    hIterRat1.SetTitle("Bayes Iterations: R=0."+str(int(R))+whichjetbintitle)
    leg_Iter=ROOT.TLegend(0.7,0.7,0.85,0.85)
    leg_Iter.AddEntry(hIterRat1,'B iter=4','l')
    leg_Iter.AddEntry(hIterRat2,'B iter=6','l')
    hIterRat1.Draw();hIterRat2.Draw('same')
    lIter = ROOT.TLine(fptbinsJlh[0], 1, fptbinsJlh[-1], 1);lIter.SetLineStyle(2);lIter.Draw("same")
    leg_Iter.Draw('same')
    #c_Iter.SaveAs('plots/5_Iter/5_Iter_ratio'+R+'.pdf')
    #c_Iter.SaveAs('plots/5_Iter/5_Iter_ratio'+R+'.png')
    c_Iter.SaveAs('plots/'+JET_or_Z+'5_Iter/5_Iter_ratio'+R+whichjetbin+'.pdf')
    c_Iter.SaveAs('plots/'+JET_or_Z+'5_Iter/5_Iter_ratio'+R+whichjetbin+'.png')
    #### -----------------------------
    c_IterSys= ROOT.TCanvas("cIterSys","cIterSys",900,500)
    histIterRMS = rootRMS([hIterRat0,hIterRat1,hIterRat2])
    if(lensysin>2):histIterRMS.GetYaxis().SetRangeUser(0,0.2);histIterRMS.SetLineColor(ROOT.kBlue+2);
    else:histIterRMS.GetYaxis().SetRangeUser(0,0.016);histIterRMS.SetLineColor(ROOT.kBlue+2);
    histIterRMS.SetTitle("Bayes Iterations: R=0."+str(int(R))+" . "+whichjetbintitle)
    histIterRMS.SetFillColor(ROOT.kBlue+2);histIterRMS.SetFillStyle(3654)
    histIterRMS.GetYaxis().SetTitle("RMS");histIterRMS.Draw()
    ttIter = ROOT.TText(fptbinsJlh[1],0.010,GetDigitText(histIterRMS,fptbinsJC))
    ttIter.Draw('same')
    #c_IterSys.SaveAs('plots/5_Iter/5_Iter_sys'+R+'.pdf')
    #c_IterSys.SaveAs('plots/5_Iter/5_Iter_sys'+R+'.png')
    c_IterSys.SaveAs('plots/'+JET_or_Z+'5_Iter/5_Iter_rat'+R+whichjetbin+'.pdf')
    c_IterSys.SaveAs('plots/'+JET_or_Z+'5_Iter/5_Iter_rat'+R+whichjetbin+'.png')
    ##### -----------------------------
    c_It= ROOT.TCanvas("cIt","cIt",900,500)
    histIterRMSIT = histIterRMS.Clone('histIterRMSIT')
    histIterRMSIT.SetMarkerSize(0)
    for i in range(histIterRMS.GetNbinsX()):
        histIterRMSIT.SetBinContent(i+1,histIterRMS.GetBinContent(i+1)+lower_limit)
    histIterRMSIT.GetYaxis().SetRangeUser(lower_limit,upper_limit);
    histIterRMSIT.GetYaxis().SetTitle('RMS (+'+str(lower_limit)+').                       Ratio to B iter=5             ')
    #histIterRMSIT.SetTitle('ratio to B iter=5 and RMS above lower limit')
    histIterRMSIT.Draw()
    hIterRat1.Draw('same');hIterRat2.Draw('same')
    leg_IT=ROOT.TLegend(0.6,0.7,0.85,0.85)
    leg_IT.AddEntry(hIterRat1,'B iter=4, ratio to iter=5','l')
    leg_IT.AddEntry(hIterRat2,'B iter=6, ratio to iter=5','l')
    leg_IT.AddEntry(histIterRMSIT,'RMS (+ lower limit on Y)')
    leg_IT.Draw('same')
    ttIT = ROOT.TText(fptbinsJlh[1],0.013+lower_limit,GetDigitText(histIterRMS,fptbinsJC))
    ttIT.Draw('same')
    c_It.SaveAs('plots/'+JET_or_Z+'5_Iter/5_ItFinal'+R+whichjetbin+'.pdf')
    c_It.SaveAs('plots/'+JET_or_Z+'5_Iter/5_ItFinal'+R+whichjetbin+'.png')
    if(wait):input()
############## -----------------------------
#histos:
# 6. UnfBayes Priors
if(flagUPrior):
    dataUPrior=[datafileUnf]
    if(lensysin>2):
        sizePrior=9
        hUPrior=[dataUPrior[0].Get('UnfProjectX_'+jetbin).Clone('hPriorMC')]
        hUPriorratio=[hUPrior[0].Clone('hratPriorMC')]
        for i in range(1,sizePrior):
            dataUPrior.append(ROOT.TFile("$cernbox/media/jackbauer/data/z_out/R_%s_finaltry/unfolding/Bayes/alljetz2D/unfold2DoutFileAPW%d.root"%(R,i)))
            hUPrior.append(dataUPrior[i].Get('UnfProjectX_'+jetbin).Clone('hPrior'+str(i-1)))
            hUPriorratio.append(dataUPrior[i].Get('UnfProjectX_'+jetbin).Clone('hPriorratio'+str(i-1)))
    else:
        sizePrior=10
        hUPrior=[dataUPrior[0].Get('unfoldedSpectrum').Clone('hPriorMC')]
        hUPriorratio=[hUPrior[0].Clone('hratPriorMC')]
        for i in range(1,sizePrior):
            dataUPrior.append(ROOT.TFile("$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_Bayes_5_priorType%d/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R,i-1)))
            hUPrior.append(dataUPrior[i].Get('unfoldedSpectrum').Clone('hPrior'+str(i-1)))
            hUPriorratio.append(dataUPrior[i].Get('unfoldedSpectrum').Clone('hPriorratio'+str(i-1)))
    #### -----------------------------
    limYpriorDo=0.94;limYpriorUp=1.07
    c_Prior = ROOT.TCanvas("cPrior","cPrior",900,500)
    leg_Prior=ROOT.TLegend(0.15,0.55,0.3,0.87)
    hUPriorratio[0].GetYaxis().SetRangeUser(limYpriorDo,limYpriorUp)
    hUPriorratio[0].GetYaxis().SetTitle('ratio');hUPriorratio[0].SetTitle('Bayes unfolding: priors')
    for i in range(sizePrior):
        hUPriorratio[i].Divide(hUPrior[0])
        HistoStyle(hUPriorratio[i],1,2,21,0,ROOTColors[i])
        for j in range(hUPriorratio[i].GetNbinsX()+1):hUPriorratio[i].SetBinError(j,0)

        if(lensysin>2):
            if(i==0):hUPriorratio[i].Draw()
            else:hUPriorratio[i].Draw('same');leg_Prior.AddEntry(hUPriorratio[i],'Prior: '+str(i))
        else:
            if(i==0):hUPriorratio[i].Draw()
            elif(i==1 or i==6):continue
            else:hUPriorratio[i].Draw('same');leg_Prior.AddEntry(hUPriorratio[i],'Prior: '+str(i-1))
    leg_Prior.Draw('same')
    c_Prior.SaveAs('plots/'+JET_or_Z+'6_Priors/6_Priors_ratio'+R+whichjetbin+'.pdf')
    c_Prior.SaveAs('plots/'+JET_or_Z+'6_Priors/6_Priors_ratio'+R+whichjetbin+'.png')
    #### -----------------------------
    c_PriorRMS = ROOT.TCanvas("cPriorRMS","cPriorRMS",900,500)
    if(lensysin>2):histPriorRMS=rootRMS(hUPriorratio[:])
    else:histPriorRMS=rootRMS(hUPriorratio[0:1]+hUPriorratio[2:6]+hUPriorratio[7:])
    #else:histPriorRMS=rootRMS(hUPriorratio[:])
    histPriorRMS.GetYaxis().SetRangeUser(0,0.06)
    histPriorRMS.SetFillColor(ROOT.kRed+2);histPriorRMS.SetFillStyle(3354);histPriorRMS.GetYaxis().SetTitle('RMS');histPriorRMS.SetTitle('Bayes unfolding: priors sys.')
    histPriorRMS.Draw()
    ttPrior = ROOT.TText(fptbinsJlo[1],0.04,GetDigitText(histPriorRMS,fptbinsJC))
    ttPrior.Draw('same')
    c_PriorRMS.SaveAs('plots/'+JET_or_Z+'6_Priors/6_Priors_RMS'+R+whichjetbin+'.pdf')
    c_PriorRMS.SaveAs('plots/'+JET_or_Z+'6_Priors/6_Priors_RMS'+R+whichjetbin+'.png')
    #### -----------------------------
    c_PriorRMSFinal = ROOT.TCanvas("cPriorRMSFinal","cPriorRMSFinal",900,500)
    leg_PriorFinal=ROOT.TLegend(0.15,0.55,0.3,0.87)
    histPriorRMSIT = histPriorRMS.Clone('histPriorsRMSIT')
    for i in range(histPriorRMSIT.GetNbinsX()):
        histPriorRMSIT.SetBinContent(i+1,histPriorRMS.GetBinContent(i+1)+limYpriorDo)
    histPriorRMSIT.SetFillColor(ROOT.kRed+2);histPriorRMSIT.SetFillStyle(3354);
    histPriorRMSIT.GetYaxis().SetTitle('RMS (+0.94).                Ratio (from all priors)         ');
    histPriorRMSIT.SetTitle('Bayes unfolding: priors sys.')
    histPriorRMSIT.GetYaxis().SetRangeUser(limYpriorDo,limYpriorUp)
    histPriorRMSIT.Draw()
    for i in range(sizePrior):
        if(lensysin>2):
            if(i==0):hUPriorratio[i].Draw('same')
            else:hUPriorratio[i].Draw('same');leg_PriorFinal.AddEntry(hUPriorratio[i],'Prior: '+str(i))
        else:
            if(i==0):hUPriorratio[i].Draw('same')
            elif(i==1 or i==6):continue
            else:hUPriorratio[i].Draw('same');leg_PriorFinal.AddEntry(hUPriorratio[i],'Prior: '+str(i-1))
    leg_PriorFinal.Draw('same')
    ttPrior2 = ROOT.TText(fptbinsJlo[1],limYpriorDo+0.02,GetDigitText(histPriorRMS,fptbinsJC))
    ttPrior2.Draw('same')
#    histPriorRMSIT.GetYaxis().SetRangeUser(limYprior,0.06)
#    ttPrior = ROOT.TText(fptbinsJlo[1],0.04,GetDigitText(histPriorRMS,fptbinsJC))
#    ttPrior.Draw('same')
    c_PriorRMSFinal.SaveAs('plots/'+JET_or_Z+'6_Priors/6_Priors_RMS_Final'+R+whichjetbin+'.pdf')
#    c_PriorRMS.SaveAs('plots/'+JET_or_Z+'6_Priors/6_Priors_RMS_Final'+R+whichjetbin+'.png')


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
        dataUSvd.append(ROOT.TFile("$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR%s_paperCuts/Default/unfolding_SVD_%d/unfoldedSpectrum_unfoldedJetSpectrum.root"%(R,i+6+radiusSVD)))
        hUSvd.append(dataUSvd[i+1].Get('unfoldedSpectrum').Clone('hSvd_'+str(i+1)))
        hUSvdratio.append(dataUSvd[i+1].Get('unfoldedSpectrum').Clone('hratSvd'+str(i+1)))
    #### -----------------------------
    c_Svd = ROOT.TCanvas("cSvd","cSvd",900,500)
    leg_Svd=ROOT.TLegend(0.7,0.7,0.85,0.85)
    hUSvdratio[0].GetYaxis().SetRangeUser(0.6,1.4)
    hUSvdratio[0].GetYaxis().SetTitle('ratio to Bayes 5')
    hUSvdratio[0].SetTitle('SVD unfolding systematic unc.')
    for i in range(sizeSVD):
        hUSvdratio[i].Divide(hUSvd[0])
        HistoStyle(hUSvdratio[i],1,2,21,0,ROOTColors[i]);
        for j in range(hUSvdratio[i].GetNbinsX()+1):hUSvdratio[i].SetBinError(j,0)
        if(i==0):hUSvdratio[i].Draw()
        else:hUSvdratio[i].Draw('same');leg_Svd.AddEntry(hUSvdratio[i],'SVD: reg='+str(i+5+radiusSVD))
    leg_Svd.Draw('same')
    c_Svd.SaveAs('plots/7_Svd/7_Svd_ratio'+R+'.pdf')
    c_Svd.SaveAs('plots/7_Svd/7_Svd_ratio'+R+'.png')
    #### -----------------------------
    c_SvdSys = ROOT.TCanvas("cSvdSys","cSvdSys",900,500)
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
    if(lensysin>2):
        hJES_0 = datafileUnf.Get('UnfProjectX_'+jetbin).Clone('hJES_0')
        hJES_1 = datafileJES.Get('UnfProjectX_'+jetbin).Clone('hJES_1')
    else:
        hJES_0 = datafileUnf.Get('unfoldedSpectrum').Clone('hJES_0') #default
        hJES_1 = datafileJES.Get('unfoldedSpectrum').Clone('hJES_1') #jes
    #### -----------------------------
    c_JES = ROOT.TCanvas("cJES","cJES",900,600)
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
    ## print ratio values for correlated JES
    print("correlated JES start")
    corrJES = []
    for i in range(len(fptbinsJC)):
        corJ = hJES_1ratio.GetBinContent(i+1)
        corrJES.append(corJ)
    print(array.array('d',corrJES))
    print("correlated JES end")
    #### -----------------------------
    c_JESsys = ROOT.TCanvas("cJESsys","cJESsys",900,600)
    lJES = ROOT.TLine(2, 1, 50, 1);lJES.SetLineStyle(2);lJES.Draw("same")
    hJES_sys = hJES_1ratio.Clone("hJES_sys");hJES_sys.GetYaxis().SetRangeUser(0.95,1.15)
    for i in range(hJES_sys.GetNbinsX()): hJES_sys.SetBinError(i+1,0)
    ttJESsys = ROOT.TText(5,1.05,GetDigitTextFromRatio(hJES_sys,fptbinsJC))
    hJES_sys.Draw();ttJESsys.Draw('same');lJES.Draw("same")
    c_JESsys.SaveAs('plots/8_JES/8_JES_sys'+R+'.pdf')
    c_JESsys.SaveAs('plots/8_JES/8_JES_sys'+R+'.png')

    yJES = array.array('d',hJES_1ratio)[:]
    xJES = array.array('d',fptbinsJC)[:]
    xJESlh= array.array('d',fptbinsJlh[:])
    print(xJES, yJES)
    funcJES = funcLine #funcJES = funcBola
    try:
        popt, pcov = curve_fit(funcJES, xJES, yJES)
    except RuntimeError:
        print("Error")

    yerrJES = [];
    for i in range(hJES_1ratio.GetNbinsX()): 
        yerrJES.append(hJES_1ratio.GetBinError(i+1));
    yerrJES=array.array('d',yerrJES)
    fig=plt.figure(1);plt.grid(axis='both',color='0.95');plt.ylim(0.6,1.8);#plt.xticks(np.arange(min(xJESall),max(xJESall)+5,5.0));
    plt.xlim(fptbinsJlh[0],fptbinsJlh[-1])
    plt.plot(xJESlh, funcJES(xJESlh,*popt),'r-', label='fit');plt.step(xJESlh[1:], yJES,'k', label="test");plt.plot(fptbinsJlh[0:2],[yJES[0],yJES[0]],'k')
    plt.errorbar(xJES,yJES,yerr=yerrJES,fmt='.k');
    if(lensysin==2):
        plt.xlabel(r'jet $p_T$')
        locJES = [2,1.4]
    elif(lensysin>2):
        plt.xlabel(r'z$_||$')
        locJES = [0.5,1.4]
    plt.ylabel('ratio');
    plt.title("JES unc: R=0."+str(int(R)))
    uncJESvals = ""
    valsJES = abs(100*(funcJES(xJES,*popt)-1)) #100*(valsJES[i]-1)
    for i in range(len(valsJES)):
        uncJESvals += "%.1f, "%(valsJES[i]) #"R%s"%(R))
    uncJESvals = uncJESvals[:-1]+' in %'
    plt.text(locJES[0], locJES[1], uncJESvals,
         #rotation=45,
         horizontalalignment='left',
         verticalalignment='top',
         #multialignment='center'
         )
    plt.draw()
    plt.savefig('plots/'+JET_or_Z+'8_JES/8_JES_rat'+R+whichjetbin+'.pdf')
    plt.savefig('plots/'+JET_or_Z+'8_JES/8_JES_rat'+R+whichjetbin+'.png')
    plt.waitforbuttonpress(1);input()
    plt.close()

    if(wait):input()

############## -----------------------------
#histos:
# . Statistical uncertainties
if(flagStat):
    if(lensysin>2):
        ### post unfolding
        fileString="$cernbox/media/jackbauer/data/z_out/R_"+R+"_finaltry/unfolding/Bayes/alljetz2D/unfold2DoutFileAPW.root"
        statSignal=ROOT.TFile(fileString)
        hSignal = statSignal.Get('UnfProjectX_'+jetbin).Clone('hStat')
        ### pre unfolding
        fileStringFD="$cernbox/media/jackbauer/data/z_out/R_"+R+"_finaltry/FDsubtraction/outFD_1.root"
        statSignalFD=ROOT.TFile(fileStringFD)
        hSignalFD = statSignalFD.Get('hsub_c'+jetbin).Clone('hStatFD')

    else:
        fileString="$cernbox/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR"+R+"_paperCuts/Default/unfolding_Bayes_5/unfoldedSpectrum_unfoldedJetSpectrum.root"
        statSignal=ROOT.TFile(fileString)
        hSignal = statSignal.Get('unfoldedSpectrum').Clone('hStat')

    statErr,statErrFD,staterror,staterrorFD=[],[],0,0
    hStat = ROOT.TH1D("hStat", "hStat", fptbinsJN, array.array('d',fptbinsJlh))
    hStatFD=ROOT.TH1D("hStatFD","hStatFD",fptbinsJN,array.array('d',fptbinsJlh))
    fptbinsJN=5
    for i in range(hSignal.GetNbinsX()):
        staterror = round(100*hSignal.GetBinError(i+1)/hSignal.GetBinContent(i+1),1)
        statErr.append(staterror)
        hStat.SetBinContent(i+1,staterror);hStat.SetBinError(i+1,0)
        if(lensysin>2):
            staterrorFD=round(100*hSignalFD.GetBinError(i+1)/hSignalFD.GetBinContent(i+1),1)
            statErrFD.append(staterrorFD)
            hStatFD.SetBinContent(i+1,staterrorFD);hStatFD.SetBinError(i+1,0)
    c_Stat = ROOT.TCanvas("cStat","cStat",900,600)
    hStat.SetLineColor(ROOT.kRed+2)
    hStat.SetMinimum(0)
    hStat.SetMaximum(100)
    hStat.Draw();
    lstat = ROOT.TLegend(0.12,0.60,0.38,0.88);
    lstat.AddEntry(hStat,"unfolded")
    if(lensysin>2): 
        hStatFD.SetLineColor(ROOT.kBlue+2)
        hStatFD.Draw("same")
        lstat.AddEntry(hStatFD,"FD sub")
    lstat.Draw("same")
    
    c_Stat.SaveAs('plots/'+JET_or_Z+'10_Stat/10_Stat_'+R+whichjetbin+'.pdf')
    c_Stat.SaveAs('plots/'+JET_or_Z+'10_Stat/10_Stat_'+R+whichjetbin+'.png')

    print(statErr)
    if(lensysin>2):print(statErrFD)
