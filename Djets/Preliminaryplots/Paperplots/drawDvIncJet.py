import os, os.path, sys
import ROOT
ROOT.TH1.AddDirectory(False)
import style_settings
import array


ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetPadLeftMargin(0.135)
ROOT.gStyle.SetPadRightMargin(0.03)
  
### SETTINGS
Rpar=int(sys.argv[1]) #R=2,3,4,6
incIndex = {2:1,3:2,4:3,6:4}
energy = "5.02"
plotmin, plotmax = 2,36;
plotYmin, plotYmax = 0,0.4;

fptbinsJlh = [5,6,8,10,14,20,30,50]
fptbinsJN = len(fptbinsJlh)-1



###############     data
fileDJ3 = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_D3zeroR0"+str(Rpar)+"_paperCuts/Default/unfolding_Bayes_5/final/JetPtSpectrum_final_fullGlobal_addedCUTandJES.root","read")
fileDJ = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Default/unfolding_Bayes_5/final/JetPtSpectrum_final_fullGlobal_addedCUTandJES.root","read")
fileIJ = ROOT.TFile("/home/jackbauer/ALICE_HeavyFlavour/work/Djets/alice_Djets/InclusiveJetsPP502TeV/HEPData-ins1733689-v1-root.root","read")

## histograms
# D-jets
hDJ = fileDJ.Get("hData_binned");
hDJ3 = fileDJ3.Get("hData_binned")
# Inclusive jets
dirIJ = fileIJ.Get("Jets in pp 5.02 TeV") #dirIJ = fileIJ.Get("Jets with UE subtraction in pp 5.02 TeV")
hIJ  = dirIJ.Get("Hist1D_y"+str(incIndex[Rpar]))
eIJ1 = dirIJ.Get("Hist1D_y"+str(incIndex[Rpar])+"_e1")
eIJ2 = dirIJ.Get("Hist1D_y"+str(incIndex[Rpar])+"_e2")

## inclusive jets: setting bin error for the histograms
for i in range(hIJ.GetNbinsX()):
    hIJ.SetBinError(i+1,eIJ1.GetBinContent(i+1))
    print(hIJ.GetBinWidth(i+1), hIJ.GetBinCenter(i+1))
    hIJ.SetBinContent(i+1,hIJ.GetBinContent(i+1)*hIJ.GetBinWidth(i+1))
    hIJ.SetBinError(i+1,hIJ.GetBinError(i+1)*hIJ.GetBinWidth(i+1))

## D-jets: setting bin error for the histos according to their bin width 
for i in range(hIJ.GetNbinsX()):
    hDJ.SetBinContent(i+1,hDJ.GetBinContent(i+1)*hDJ.GetBinWidth(i+1))
    hDJ.SetBinError(i+1,hDJ.GetBinError(i+1)*hDJ.GetBinWidth(i+1))
#### UNSCALE by bin width
print("===== unscaling by bin-width")

#exit()
#### REBIN
hIJReb = hIJ.Rebin(fptbinsJN,'hIJReb',array.array('d',fptbinsJlh));
hDJReb = hDJ.Rebin(fptbinsJN,'hDJReb',array.array('d',fptbinsJlh));
hDJReb3 = hDJ3.Rebin(fptbinsJN,'hDJReb3',array.array('d',fptbinsJlh))

hR_Data = hDJReb.Clone("hR_Data")
hR_Data3 = hDJReb3.Clone("hR_Data3")
hR_Data.Divide(hIJReb)
hR_Data3.Divide(hIJReb)
###############     Theory
fileDJsoft = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Simulations/Prompt/JetPt_AnalysisResults_FastSim_pythia8_charm_soft2color_1598639645_Dpt2_36_Dzero.root","read")
fileDJmb = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Simulations/Prompt/JetPt_AnalysisResults_FastSim_pythia8_charm_1594137862_Dpt2_36_Dzero.root","read")
fileIJsoft = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Simulations/Prompt/JetPt_AnalysisResults_FastSim_pythia8_soft2color_1611748149_Dpt2_36_Dzero.root","read")
fileIJmb = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Simulations/Prompt/JetPt_AnalysisResults_FastSim_pythia8_mb_1611743722_Dpt2_36_Dzero.root","read")
fileIJmonash = ROOT.TFile("/home/jackbauer/ALICE_HeavyFlavour/work/Djets/alice_Djets/InclusiveJetsPP502TeV/theoryMonashR04.root","read")

hDJsoft = fileDJsoft.Get("hPt")
hDJmb   = fileDJmb.Get("hPt")
hIJsoft = fileIJsoft.Get("hPt")
hIJmb   = fileIJmb.Get("hPt")
hIJmonash=fileIJmonash.Get("fSpecPythia8_R04_PYTHIA8_Monash_2013_reb")

hDJsoftReb = hDJsoft.Rebin(fptbinsJN,'hDJsoftReb',array.array('d',fptbinsJlh))
hIJsoftReb = hIJsoft.Rebin(fptbinsJN,'hIJsoftReb',array.array('d',fptbinsJlh))
hR_soft = hDJsoftReb.Clone("hR_soft")
hR_soft.Divide(hIJsoftReb)

hDJmbReb = hDJmb.Rebin(fptbinsJN,'hDJmbReb',array.array('d',fptbinsJlh))
hIJmbReb = hIJmb.Rebin(fptbinsJN,'hIJmbReb',array.array('d',fptbinsJlh))
hR_mb = hDJmbReb.Clone("hR_mb")
hR_mb.Divide(hIJmbReb)
hR_mbNew = hDJmbReb.Clone("hR_mbNew")
### style
hR_mb.SetLineColor(ROOT.kRed+2)
hR_mb.SetLineWidth(3)
hR_soft.SetLineColor(ROOT.kGreen+2)
hR_soft.SetLineWidth(3)
hR_Data.SetLineWidth(3)
hR_Data3.SetLineWidth(3);hR_Data3.SetLineColor(ROOT.kBlue+2);hR_Data3.SetMarkerColor(ROOT.kBlue+2)

hR_Data.GetYaxis().SetTitle("R=0.4, #frac{D-jets}{Inclusive jets}")

leg = ROOT.TLegend(0.5, 0.6, 0.9, 0.9)
### TCANVAS
cRatio = ROOT.TCanvas("cRatio","cRatio",800,500)
hR_Data.GetYaxis().SetRangeUser(0.01,0.14)
hR_Data.Draw()
#hR_Data3.Draw('same')
#hR_soft.Scale(1/10)
#hR_mb.Scale(1/10)
hR_soft.Draw("same")
hR_mb.Draw("same")
hR_mbNew.Draw("same")
hIJmonash.Draw("same")
hDJmbReb.Draw("same")

leg.AddEntry(hR_Data,"Data Dpt>2 GeV")
#leg.AddEntry(hR_Data3,"Data Dpt>3 GeV")
leg.AddEntry(hR_mb,"Pythia 8 MB (scaled by 1/10)")
leg.AddEntry(hR_soft,"Pythia 8 soft, mode 2 (scaled by 1/10)")
leg.Draw("same")

cRatio.SaveAs("plots/DvIncJets.pdf")
cRatio.SaveAs("plots/DvIncJets.png")

wait=input()
