import os, os.path, sys
import ROOT
ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetBatch()
import style_settings
import array
import yaml
#import numpy as np

ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetPadLeftMargin(0.19)
ROOT.gStyle.SetPadRightMargin(0.03)

os.system('cp ../../finalSpectra/5TeVunc.yaml 5TeVunc.yaml')

### Sanity check
if len(sys.argv)!=3:
    print("python file.py R(2,3,4,6) sys.yaml")
### SETTINGS
flag_paper = 0 # turn off this flag if you want to get all four R together a, b, c, d
Rpar=int(sys.argv[1]) # R=2,3,4,6
incIndex = {2:1,3:2,4:3,6:4} # index for radius of inclusive jets
enumIndex= {2:"(a)",3:"(b)",4:"(c)",6:"(d)"}
if flag_paper:
    enumIndex = {2:"(a)",3:"",4:"(b)",6:"(c)"}
energy = "5.02"
plotmin, plotmax = 5,50
plotYmin, plotYmax = 0,0.17;
fptbinsJlh = [5,6,8,10,14,20,30,50]
fptbinsJN = len(fptbinsJlh)-1
Colors = [ROOT.kRed+1,ROOT.kBlue+2,ROOT.kGreen+2,ROOT.kViolet+2]
Markers = [20,21,22,23,24,4]
x_title = "#it{p}_{T,ch. jet} (GeV/#it{c})" 
y_title = "#frac{d^{2}#it{#sigma}^{D^{0}-jets}}{d#it{p}_{T,ch. jet}d#it{#eta}_{ch. jet}} / #frac{d^{2}#it{#sigma}^{inclusive jets}}{d#it{p}_{T,ch. jet}d#it{#eta}_{ch. jet}}"
#y_title = "#sigma(D^{0}-jets) / #sigma(inclusive jets)"
#y_title = "D^{0}-jet fraction "
#y_title = "D^{0} tagged charged particle jet fraction  "
## Settings others
textsize=0.035
markersize=1.5
ptval, value, ptvalunc, valueerrup, valueerrdown = [], [], [], [], []
valDA, valDAerrup, valDAerrdo = [], [], []
for i in range(fptbinsJN):
    ## powheg 
    ptval.append((fptbinsJlh[i]+fptbinsJlh[i+1])/2.0)
    ptvalunc.append((fptbinsJlh[i+1]-fptbinsJlh[i])/2.0)
###### yaml systematics
yaml_config_file = sys.argv[2]
# /home/jackbauer/ALICE_HeavyFlavour/work/Djets/kvapilj/alice_Djets/Djets/finalSpectra/5TeVunc.yaml
with open(yaml_config_file, 'r') as stream:
    unc_loaded = yaml.safe_load(stream)
systUncD_up = unc_loaded[0]["0"+str(Rpar)]['systUncD_up']
systUncD_do = unc_loaded[0]["0"+str(Rpar)]['systUncD_down']
systUncD_JES = unc_loaded[0]["0"+str(Rpar)]['systUncD_JES']
systUncD_CUTS = unc_loaded[0]["0"+str(Rpar)]['systUncD_CUTS']
print(type(systUncD_up))
## numpy does not work: so alternative.
djetsUp=[0.0,0.0,0.0,0.0,0.0,0.0,0.0]
djetsDo=[0.0,0.0,0.0,0.0,0.0,0.0,0.0]
BR=0.008
LUMI=0.021
for i in range(fptbinsJN):
    djetsUp[i] = ROOT.TMath.Sqrt(systUncD_up[i]**2+systUncD_JES[i]**2+systUncD_CUTS[i]**2+BR**2+LUMI**2)
    djetsDo[i] = ROOT.TMath.Sqrt(systUncD_do[i]**2+systUncD_JES[i]**2+systUncD_CUTS[i]**2+BR**2+LUMI**2)
print(djetsUp)
print(djetsDo)
####### GATHERING DATA
fileDJ = ROOT.TFile("/eos/user/a/amohanty/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Default/unfolding_Bayes_5/final/JetPtSpectrum_final_fullGlobal_addedCUTandJES.root","read")
fileIJ = ROOT.TFile("/eos/user/a/amohanty/home/jackbauer/ALICE_HeavyFlavour/work/Djets/alice_Djets/InclusiveJetsPP502TeV/HEPData-ins1733689-v1-root.root","read")
####### Creating DATA histograms
# D-jets
hDJ = fileDJ.Get("hData_binned");
# Inclusive jets
dirIJ = fileIJ.Get("Jets in pp 5.02 TeV") #dirIJ = fileIJ.Get("Jets with UE subtraction in pp 5.02 TeV")
hIJ  = dirIJ.Get("Hist1D_y"+str(incIndex[Rpar]))
eIJ1 = dirIJ.Get("Hist1D_y"+str(incIndex[Rpar])+"_e1")
eIJ2 = dirIJ.Get("Hist1D_y"+str(incIndex[Rpar])+"_e2")

hDJ.Sumw2()
hIJ.Sumw2()
hIJsyst = hIJ.Clone("hIJsyst")
#### UNSCALE DATA by bin width
print("===== unscaling by bin-width")
## inclusive jets: setting bin error for the histograms
for i in range(hIJ.GetNbinsX()):
    hIJ.SetBinError(i+1,eIJ1.GetBinContent(i+1))
    hIJsyst.SetBinError(i+1,eIJ2.GetBinContent(i+1))
    print("IncJets:",hIJ.GetBinWidth(i+1), hIJ.GetBinCenter(i+1))
    hIJ.SetBinContent(i+1,hIJ.GetBinContent(i+1)*hIJ.GetBinWidth(i+1))
    hIJsyst.SetBinContent(i+1,hIJsyst.GetBinContent(i+1)*hIJsyst.GetBinWidth(i+1))
    hIJ.SetBinError(i+1,hIJ.GetBinError(i+1)*hIJ.GetBinWidth(i+1))
    hIJsyst.SetBinError(i+1,hIJsyst.GetBinError(i+1)*hIJsyst.GetBinWidth(i+1))
## D-jets: setting bin error for the histos according to their bin width 
for i in range(hDJ.GetNbinsX()):
    hDJ.SetBinContent(i+1,hDJ.GetBinContent(i+1)*hDJ.GetBinWidth(i+1))
    hDJ.SetBinError(i+1,hDJ.GetBinError(i+1)*hDJ.GetBinWidth(i+1))
    print("D-Jets:",hDJ.GetBinWidth(i+1), hDJ.GetBinCenter(i+1), hDJ.GetBinContent(i+1))
print("===== unscaling by bin-width completed")
#### REBIN DATA
hIJsystReb = hIJsyst.Rebin(fptbinsJN,'hIJsystReb',array.array('d',fptbinsJlh));
hIJReb = hIJ.Rebin(fptbinsJN,'hIJReb',array.array('d',fptbinsJlh));
hDJReb = hDJ.Rebin(fptbinsJN,'hDJReb',array.array('d',fptbinsJlh));

hIJsystReb.Sumw2()
hIJReb.Sumw2()
hDJReb.Sumw2()

hR_Data = hDJReb.Clone("hR_Data")
hR_Data.Divide(hIJReb)
### assigning systematic uncertainties of D-jets as statistical uncertainties
hR_DatasystUp = hDJReb.Clone("hR_DataSyst")
hR_DatasystDo = hDJReb.Clone("hR_DataSyst")
for i in range(fptbinsJN):
    print("i",i,"i:",hR_DatasystUp.GetXaxis().FindBin(ptval[i]))
    hR_DatasystUp.SetBinError(
            (hR_DatasystUp.GetXaxis().FindBin(ptval[i])),
            djetsUp[i]*hR_DatasystUp.GetBinContent(hR_DatasystUp.GetXaxis().FindBin(ptval[i])))
    hR_DatasystDo.SetBinError(
            (hR_DatasystDo.GetXaxis().FindBin(ptval[i])),
            djetsDo[i]*hR_DatasystDo.GetBinContent(hR_DatasystDo.GetXaxis().FindBin(ptval[i])))

hR_DatasystUp.Divide(hIJsystReb)
hR_DatasystDo.Divide(hIJsystReb)
#######################################################################################################################################################
###############     Theory
###----1. Pythia 8 Monash
fileDJmonash = ROOT.TFile("/eos/user/a/amohanty/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Default/unfolding_Bayes_5/final/JetPtSpectrum_final_fullGlobal_addedCUTandJES.root","read")
fileIJmonash = ROOT.TFile("/eos/user/a/amohanty/home/jackbauer/ALICE_HeavyFlavour/work/Djets/alice_Djets/InclusiveJetsPP502TeV/theoryMonashR0"+str(Rpar)+".root","read")

hDJmb   = fileDJmonash.Get("hsimPythia8_central")
hIJmb   = fileIJmonash.Get("fSpecPythia8_R0"+str(Rpar)+"_PYTHIA8_Monash_2013_reb")

###----3. Powheg+Pythia8
hDJpow8 = fileDJmonash.Get("hsimPowhegPythia8_central")
hDJpow8Err = [fileDJmonash.Get("hsimPowhegPythia8_up"),fileDJmonash.Get("hsimPowhegPythia8_down")]

fileIJpowheg = ROOT.TFile("/eos/user/a/amohanty/home/jackbauer/ALICE_HeavyFlavour/work/Djets/alice_Djets/InclusiveJetsPP502TeV/theoryPowhegR0"+str(Rpar)+".root")
hIJpow8 = fileIJpowheg.Get("R0"+str(Rpar)+"PwgScaleErr")
###### Histograms UNSCALE AND REBIN:
##--D-jets min bias pythia 8
for i in range(hDJmb.GetNbinsX()):
    hDJmb.SetBinContent(i+1,hDJmb.GetBinContent(i+1)*hDJmb.GetBinWidth(i+1))
    hDJmb.SetBinError(i+1,hDJmb.GetBinError(i+1)*hDJmb.GetBinWidth(i+1))

hDJmbReb = hDJmb.Rebin(fptbinsJN,'hDJmbReb',array.array('d',fptbinsJlh))
##--inclusive jets  min bias pythia 8
for i in range(hIJmb.GetNbinsX()):
    hIJmb.SetBinContent(i+1,hIJmb.GetBinContent(i+1)*hIJmb.GetBinWidth(i+1))
    hIJmb.SetBinError(i+1,hIJmb.GetBinError(i+1)*hIJmb.GetBinWidth(i+1))

hIJmbReb = hIJmb.Rebin(fptbinsJN,'hIJmbReb',array.array('d',fptbinsJlh))
##--D-jets  powheg pythia 8
for i in range(hDJpow8.GetNbinsX()):
    hDJpow8.SetBinContent(i+1,hDJpow8.GetBinContent(i+1)*hDJpow8.GetBinWidth(i+1))
    hDJpow8Err[0].SetBinContent(i+1,hDJpow8Err[0].GetBinContent(i+1)*hDJpow8Err[0].GetBinWidth(i+1))
    hDJpow8Err[1].SetBinContent(i+1,hDJpow8Err[1].GetBinContent(i+1)*hDJpow8Err[1].GetBinWidth(i+1))
    hDJpow8.SetBinError(i+1,0)
    hDJpow8Err[0].SetBinError(i+1,0)
    hDJpow8Err[1].SetBinError(i+1,0)

hDJpowReb = hDJpow8.Rebin(fptbinsJN,'hDJpow8',array.array('d',fptbinsJlh))
hDJpowRebErr= [hDJpow8Err[0].Rebin(fptbinsJN,'hDJpow8RebErrUp',array.array('d',fptbinsJlh)), hDJpow8Err[1].Rebin(fptbinsJN,'hDJpow8RebErrDo',array.array('d',fptbinsJlh))]
##--inclusive jets  powheg pythia 8
for i in range(hIJpow8.GetNbinsX()):
    hIJpow8.SetBinContent(i+1,hIJpow8.GetBinContent(i+1)*hIJpow8.GetBinWidth(i+1))
    hIJpow8.SetBinError(i+1,hIJpow8.GetBinError(i+1)*hIJpow8.GetBinWidth(i+1))

hIJpowReb = hIJpow8.Rebin(fptbinsJN,'hIJpow8',array.array('d',fptbinsJlh))

###########--- DIVISION
##--- minimum bias/pythia 8
hR_mb = hDJmbReb.Clone("hR_mb") #hR_mb0= hDJmb.Clone("hR_mb0")
hR_mb.Divide(hIJmbReb) #hR_mb0.Divide(hIJmbReb)
##--- powheg pythia 8
hR_pow8 = hDJpowReb.Clone("hR_pow8")
hR_pow8Err = [hDJpowRebErr[0].Clone("hR_pow8ErrUp"), hDJpowRebErr[1].Clone("hR_pow8ErrDo")]
hR_pow8.Divide(hIJpowReb)
for i in range(hR_pow8.GetNbinsX()):
    print(hDJpowRebErr[0].GetBinContent(i+1)-hDJpowReb.GetBinContent(i+1))
    hR_pow8.SetBinError(i+1,0)
    hR_pow8Err[0].SetBinContent(i+1,
            (
            (hR_pow8Err[0].GetBinContent(i+1)-hDJpowReb.GetBinContent(i+1))
            /(hDJpowReb.GetBinContent(i+1))
            +hIJpowReb.GetBinError(i+1)/hIJpowReb.GetBinContent(i+1)
            )*0.5*hR_pow8.GetBinContent(i+1) + hR_pow8.GetBinContent(i+1)
                )

    hR_pow8Err[1].SetBinContent(i+1,
            -(
            (-hR_pow8Err[1].GetBinContent(i+1)+hDJpowReb.GetBinContent(i+1))
            /(hDJpowReb.GetBinContent(i+1))
            +hIJpowReb.GetBinError(i+1)/hIJpowReb.GetBinContent(i+1)
            )*0.5*hR_pow8.GetBinContent(i+1) + hR_pow8.GetBinContent(i+1)
                )
### style
hR_pow8Err[0].SetLineColor(ROOT.kOrange+20)
hR_pow8Err[1].SetLineColor(ROOT.kOrange+20)
hR_pow8.SetLineColor(ROOT.kBlue+2)
hR_mb.SetLineColor(ROOT.kOrange+2)
hR_mb.SetMarkerColor(ROOT.kOrange+2)
hR_mb.SetMarkerSize(0)
hR_mb.SetLineWidth(3)

hR_Data.SetLineWidth(3)
hR_Data.SetLineColor(ROOT.kRed+2)
hR_Data.SetMarkerSize(markersize)
hR_Data.SetMarkerStyle(Markers[0])
hR_Data.SetMarkerColor(ROOT.kRed+2)

hR_Data.GetYaxis().SetTitle(" cross section ratio for R=0."+str(Rpar)+", #frac{cross section_{D-jets}}{cross section_{inclusive jets}}")

hR_DatasystUp.SetLineColor(ROOT.kBlue+2)
hR_DatasystDo.SetLineColor(ROOT.kBlue+2)
####################################################
## theory and data uncertainty
##-------------------
for i in range(fptbinsJN):
    ## powheg 
    value.append(hR_pow8.GetBinContent(hR_pow8.GetXaxis().FindBin(ptval[i])))
    valueerrup.append(hR_pow8Err[0].GetBinContent(hR_pow8Err[0].GetXaxis().FindBin(ptval[i]))-hR_pow8.GetBinContent(hR_pow8Err[0].GetXaxis().FindBin(ptval[i])))
    valueerrdown.append(hR_pow8.GetBinContent(hR_pow8Err[0].GetXaxis().FindBin(ptval[i]))-hR_pow8Err[1].GetBinContent(hR_pow8Err[1].GetXaxis().FindBin(ptval[i])))
    ## data sys
    #hR_DatasystUp.Divide(hIJsystReb)
    #hR_DatasystDo.Divide(hIJsystReb)
    valDA.append(hR_Data.GetBinContent(hR_Data.GetXaxis().FindBin(ptval[i])))
    valDAerrup.append(hR_DatasystUp.GetBinError(hR_DatasystUp.GetXaxis().FindBin(ptval[i])))
    valDAerrdo.append(hR_DatasystDo.GetBinError(hR_DatasystDo.GetXaxis().FindBin(ptval[i])))

### POWHEG SYS
grsys = ROOT.TGraphAsymmErrors(fptbinsJN,array.array('d',ptval),array.array('d',value),array.array('d',ptvalunc),array.array('d',ptvalunc),array.array('d',valueerrdown),array.array('d',valueerrup))
grsys.SetFillColor(Colors[1])
grsys.SetLineColor(Colors[1])
grsys.SetLineWidth(2)
grsys.SetFillStyle(0)
grsys.SetMarkerStyle(Markers[5])#71)#107)
grsys.SetMarkerColor(Colors[1])
grsys.SetMarkerSize(markersize);
grsys.SetTitle('')
#grsys.GetXaxis().SetTitle(x_title);
#grsys.GetYaxis().SetTitle(y_title);
#grsys.GetXaxis().SetLabelSize(0.04);
#grsys.GetXaxis().SetTitleSize(0.04);
#grsys.GetXaxis().SetTitleOffset(1.1);
#grsys.GetYaxis().SetTitleOffset(1.6);
#grsys.GetYaxis().SetLabelSize(0.04);
#grsys.GetYaxis().SetTitleSize(0.04);
#grsys.GetXaxis().SetRangeUser(plotmin,plotmax);
#grsys.GetYaxis().SetRangeUser(plotYmin,plotYmax);
### DATA SYS
grdat = ROOT.TGraphAsymmErrors(fptbinsJN,array.array('d',ptval),array.array('d',valDA),array.array('d',ptvalunc),array.array('d',ptvalunc),array.array('d',valDAerrdo),array.array('d',valDAerrup))
grdat.SetFillColor(ROOT.TColor.GetColor("#cccccc"))
grdat.SetLineColor(ROOT.TColor.GetColor("#cccccc"))#ROOT.TColor.GetColor("#000099")) 
#grdat.SetLineWidth(2)
grdat.SetFillStyle(1001)
grdat.SetMarkerStyle(Markers[0])
grdat.SetMarkerColor(ROOT.kRed+2)
grdat.SetMarkerSize(markersize)
grdat.SetTitle('')
grdat.GetXaxis().SetTitle(x_title);
grdat.GetYaxis().SetTitle(y_title);
grdat.GetXaxis().SetLabelSize(0.04);
grdat.GetXaxis().SetTitleSize(0.04);
grdat.GetXaxis().SetTitleOffset(1.1);
grdat.GetYaxis().SetTitleOffset(1.9);
grdat.GetYaxis().SetLabelSize(0.04);
grdat.GetYaxis().SetTitleSize(0.04);
grdat.GetXaxis().SetRangeUser(plotmin,plotmax);
grdat.GetYaxis().SetRangeUser(plotYmin,plotYmax);

############################## FORMAL STUFF
ptL, ptR, pt_row = [], [], []
for i in range(4): # number of jet radii to be plotted
    if(i==1):ptL.append(ROOT.TPaveText(0.22,0.80,0.65,0.95,"NB NDC"))
    #elif(i==2):ptL.append(ROOT.TPaveText(0.22,0.80,0.65,0.95,"NB NDC"))
    else:ptL.append(ROOT.TPaveText(0.22,0.87,0.65,0.95,"NB NDC"))

    ptL[i].SetBorderSize(0)
    ptL[i].SetFillStyle(0)
    ptL[i].SetTextAlign(13)
    ptL[i].SetTextFont(43)
    ptL[i].SetTextSize(20)

    ptR.append(ROOT.TPaveText(0.75,0.87,0.85,0.95,"NB NDC"))
    ptR[i].SetBorderSize(0)
    ptR[i].SetFillStyle(0)
    ptR[i].SetTextAlign(13)
    ptR[i].SetTextFont(43)
    ptR[i].SetTextSize(20)
    
    ptR[i].AddText(enumIndex[Rpar]+" #it{R} = 0.%d"%(int(Rpar)))

ptL[0].AddText("ALICE, pp, #sqrt{#it{s}} = 5.02 TeV")
if flag_paper:
    ptL[2].AddText("pp, #sqrt{#it{s}} = 5.02 TeV")
    ptL[2].AddText("Charged jets, Anti-#it{k}_{T}, |#it{#eta}_{ch. jet}| < 0.9 - #it{R}")
    ptL[2].AddText("2 < #it{p}_{T,D^{0}} < 36 GeV/#it{c}") 
else:
    ptL[1].AddText("pp, #sqrt{#it{s}} = 5.02 TeV")
    ptL[1].AddText("Charged jets, Anti-#it{k}_{T}, |#it{#eta}_{ch. jet}| < 0.9 - #it{R}")
    ptL[2].AddText("2 < #it{p}_{T,D^{0}} < 36 GeV/#it{c}") 


############################## CANVAS
leg = ROOT.TLegend(0.22,0.7,0.65,0.85);
leg.SetTextSize(textsize);
leg.SetFillColor(0);
leg.SetFillStyle(0);

cDvI = ROOT.TCanvas("cDvI","cDvI",620,510);
grdat.GetYaxis().SetDecimals() #POWHEG SYSTEMATICS
grdat.Draw("2ap") # DATA SYSTEMATICS
grsys.Draw("p2same")
hR_Data.Draw("same p e0 x0")
hR_mb.Draw("same")

leg.AddEntry(grdat,"Data","pf");
leg.AddEntry(hR_mb,"PYTHIA 8 HardQCD Monash 2013","l");
leg.AddEntry(grsys,"POWHEG hvq + PYTHIA 8","pf");
if(Rpar==2):leg.Draw("same");
ptL[incIndex[Rpar]-1].Draw()
ptR[incIndex[Rpar]-1].Draw()

cDvI.Update()

cDvI.SaveAs("plots/DJetsVIncJets_R0"+str(Rpar)+".png")
cDvI.SaveAs("plots/DJetsVIncJets_R0"+str(Rpar)+".pdf")
ofile = ROOT.TFile("DvInc_R0"+str(Rpar)+".root","recreate")
hR_Data.Write("hR_Data_R0"+str(Rpar))
ofile.Close()
