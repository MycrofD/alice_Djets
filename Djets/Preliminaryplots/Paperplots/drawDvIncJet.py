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
fileDJ = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Default/unfolding_Bayes_5/final/JetPtSpectrum_final_fullGlobal_addedCUTandJES.root","read")
fileIJ = ROOT.TFile("/home/jackbauer/ALICE_HeavyFlavour/work/Djets/alice_Djets/InclusiveJetsPP502TeV/HEPData-ins1733689-v1-root.root","read")

## histograms
# D-jets
hDJ = fileDJ.Get("hData_binned");
# Inclusive jets
dirIJ = fileIJ.Get("Jets in pp 5.02 TeV") #dirIJ = fileIJ.Get("Jets with UE subtraction in pp 5.02 TeV")
hIJ  = dirIJ.Get("Hist1D_y"+str(incIndex[Rpar]))
eIJ1 = dirIJ.Get("Hist1D_y"+str(incIndex[Rpar])+"_e1")
eIJ2 = dirIJ.Get("Hist1D_y"+str(incIndex[Rpar])+"_e2")

hDJ.Sumw2()
hIJ.Sumw2()
#### UNSCALE by bin width
print("===== unscaling by bin-width")
## inclusive jets: setting bin error for the histograms
for i in range(hIJ.GetNbinsX()):
    hIJ.SetBinError(i+1,eIJ1.GetBinContent(i+1))
    print("IncJets:",hIJ.GetBinWidth(i+1), hIJ.GetBinCenter(i+1))
    hIJ.SetBinContent(i+1,hIJ.GetBinContent(i+1)*hIJ.GetBinWidth(i+1))
    hIJ.SetBinError(i+1,hIJ.GetBinError(i+1)*hIJ.GetBinWidth(i+1))

## D-jets: setting bin error for the histos according to their bin width 
for i in range(hIJ.GetNbinsX()):
    hDJ.SetBinContent(i+1,hDJ.GetBinContent(i+1)*hDJ.GetBinWidth(i+1))
    hDJ.SetBinError(i+1,hDJ.GetBinError(i+1)*hDJ.GetBinWidth(i+1))
    print("D-Jets:",hDJ.GetBinWidth(i+1), hDJ.GetBinCenter(i+1), hDJ.GetBinContent(i+1))
print("===== unscaling by bin-width completed")

#exit()
#### REBIN
hIJReb = hIJ.Rebin(fptbinsJN,'hIJReb',array.array('d',fptbinsJlh));
hDJReb = hDJ.Rebin(fptbinsJN,'hDJReb',array.array('d',fptbinsJlh));

hIJReb.Sumw2()
hDJReb.Sumw2()

hIJReb = hIJ.Rebin(fptbinsJN,'hIJReb',array.array('d',fptbinsJlh));
hR_Data = hDJReb.Clone("hR_Data")
hR_Data.Divide(hIJReb)
###############     Theory
fileIJmonash = ROOT.TFile("/home/jackbauer/ALICE_HeavyFlavour/work/Djets/alice_Djets/InclusiveJetsPP502TeV/theoryMonashR0"+str(Rpar)+".root","read")
fileDJmonash = ROOT.TFile("~/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Default/unfolding_Bayes_5/final/JetPtSpectrum_final_fullGlobal_addedCUTandJES.root","read")

hDJmb   = fileDJmonash.Get("hsimPythia8_central")
hIJmb   = fileIJmonash.Get("fSpecPythia8_R0"+str(Rpar)+"_PYTHIA8_Monash_2013_reb")

## inclusive jets
for i in range(hIJmb.GetNbinsX()):
    hIJmb.SetBinContent(i+1,hIJmb.GetBinContent(i+1)*hIJmb.GetBinWidth(i+1))
    hIJmb.SetBinError(i+1,hIJmb.GetBinError(i+1)*hIJmb.GetBinWidth(i+1))

## D-jets
for i in range(hDJmb.GetNbinsX()):
    hDJmb.SetBinContent(i+1,hDJmb.GetBinContent(i+1)*hDJmb.GetBinWidth(i+1))
    hDJmb.SetBinError(i+1,hDJmb.GetBinError(i+1)*hDJmb.GetBinWidth(i+1))

hDJmbReb = hDJmb.Rebin(fptbinsJN,'hDJmbReb',array.array('d',fptbinsJlh))
hIJmbReb = hIJmb.Rebin(fptbinsJN,'hIJmbReb',array.array('d',fptbinsJlh))

for i in range(hIJmb.GetNbinsX()+1):
    if(hIJmb.GetBinCenter(i)>5):
        print(i,":",hIJmb.GetBinCenter(i),hIJmb.GetBinContent(i))

for i in range(hIJmbReb.GetNbinsX()+1):
    print(i,":",hIJmbReb.GetBinCenter(i),hIJmbReb.GetBinContent(i))


hR_mb0= hDJmb.Clone("hR_mb0")
hR_mb = hDJmbReb.Clone("hR_mb")
hR_mb.Divide(hIJmbReb)
hR_mb0.Divide(hIJmbReb)
### style
hR_mb.SetLineColor(ROOT.kRed+2)
hR_mb.SetLineWidth(3)
hR_Data.SetLineWidth(3)

hR_Data.GetYaxis().SetTitle("#frac{cross section_{D-jets}}{cross section_{inclusive jets}} cross section ratio for R=0."+str(Rpar)+", ")

leg = ROOT.TLegend(0.5, 0.6, 0.9, 0.9)
### TCANVAS
cRatio = ROOT.TCanvas("cRatio","cRatio",800,500)
hR_Data.GetYaxis().SetRangeUser(0.01,0.14)
hR_Data.Draw()
hR_mb.Draw("same")
hR_mb0.Draw("same")

leg.AddEntry(hR_Data,"Data Dpt>2 GeV")
leg.AddEntry(hR_mb,"Pythia 8 MB")
leg.Draw("same")

cRatio.SaveAs("plots/DvIncJets.pdf")
cRatio.SaveAs("plots/DvIncJets.png")

wait=input()
