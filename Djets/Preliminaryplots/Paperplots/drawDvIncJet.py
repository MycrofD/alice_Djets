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
#######################################################################################################################################################
###############     Theory
###----1. Pythia 8 Monash
fileDJmonash = ROOT.TFile("~/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Default/unfolding_Bayes_5/final/JetPtSpectrum_final_fullGlobal_addedCUTandJES.root","read")
fileIJmonash = ROOT.TFile("/home/jackbauer/ALICE_HeavyFlavour/work/Djets/alice_Djets/InclusiveJetsPP502TeV/theoryMonashR0"+str(Rpar)+".root","read")

hDJmb   = fileDJmonash.Get("hsimPythia8_central")
hIJmb   = fileIJmonash.Get("fSpecPythia8_R0"+str(Rpar)+"_PYTHIA8_Monash_2013_reb")

###----2. Pythia 8 Soft mode 2
#fileDJsoft = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Simulations/Prompt/JetPt_AnalysisResults_FastSim_pythia8_charm_soft2color_1598639645_Dpt2_36_Dzero.root","read")
#fileIJsoft = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Simulations/Prompt/JetPt_AnalysisResults_FastSim_pythia8_soft2color_1611748149_Dpt2_36_Dzero.root","read")
#fileIJmbRaw = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Simulations/Prompt/JetPt_AnalysisResults_FastSim_pythia8_mb_1611743722_Dpt2_36_Dzero.root","read")

#hDJsm2   = fileDJsoft.Get("hPt")
#hIJsm2   = fileIJsoft.Get("hPt")
#hIJmbRaw = fileIJmbRaw.Get("hPt")

###----3. Powheg+Pythia8
hDJpow8 = [fileDJmonash.Get("hsimPowhegPythia8_central"),fileDJmonash.Get("hsimPowhegPythia8_central"),fileDJmonash.Get("hsimPowhegPythia8_central")]
hDJpow8Err = [fileDJmonash.Get("hsimPowhegPythia8_up"),fileDJmonash.Get("hsimPowhegPythia8_down")]

fileIJpowheg = ROOT.TFile("/home/jackbauer/ALICE_HeavyFlavour/work/Djets/alice_Djets/InclusiveJetsPP502TeV/theoryPowhegR0"+str(Rpar)+".root")
hIJpow8 = fileIJpowheg.Get("R0"+str(Rpar)+"PwgScaleErr")
###### Histograms:
## D-jets min bias pythia 8
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
for i in range(hDJpow8[0].GetNbinsX()):
    hDJpow8[0].SetBinContent(i+1,hDJpow8[0].GetBinContent(i+1)*hDJpow8[0].GetBinWidth(i+1))
    hDJpow8[1].SetBinContent(i+1,hDJpow8[1].GetBinContent(i+1)*hDJpow8[1].GetBinWidth(i+1))
    hDJpow8[2].SetBinContent(i+1,hDJpow8[2].GetBinContent(i+1)*hDJpow8[2].GetBinWidth(i+1))
    hDJpow8[0].SetBinError(i+1,0)
    hDJpow8[1].SetBinError(i+1,hDJpow8Err[0].GetBinContent(i+1)*hDJpow8[0].GetBinWidth(i+1))
    hDJpow8[2].SetBinError(i+1,hDJpow8Err[1].GetBinContent(i+1)*hDJpow8[1].GetBinWidth(i+1))

hDJpowReb = hDJpow8[0].Rebin(fptbinsJN,'hDJpow8',array.array('d',fptbinsJlh))
hDJpowRebErr= [hDJpow8[1].Rebin(fptbinsJN,'hDJpow8RebErrUp',array.array('d',fptbinsJlh)), hDJpow8[2].Rebin(fptbinsJN,'hDJpow8RebErrDo',array.array('d',fptbinsJlh))]
##--inclusive jets  powheg pythia 8
for i in range(hIJpow8.GetNbinsX()):
    hIJpow8.SetBinContent(i+1,hIJpow8.GetBinContent(i+1)*hIJpow8.GetBinWidth(i+1))
    hIJpow8.SetBinError(i+1,hIJpow8.GetBinError(i+1)*hIJpow8.GetBinWidth(i+1))

hIJpowReb = hIJpow8.Rebin(fptbinsJN,'hIJpow8',array.array('d',fptbinsJlh))
### printing for debugging purpose
#for i in range(hIJmb.GetNbinsX()+1):
#    if(hIJmb.GetBinCenter(i)>5):
#        print(i,":",hIJmb.GetBinCenter(i),hIJmb.GetBinContent(i))
#
#for i in range(hIJmbReb.GetNbinsX()+1):
#    print(i,":",hIJmbReb.GetBinCenter(i),hIJmbReb.GetBinContent(i))

###-- soft mode 2 pythia 8
#hIJmbRebRaw = hIJmbRaw.Rebin(fptbinsJN,'hIJmbRebRaw',array.array('d',fptbinsJlh))
#hDJsmReb = hDJsm2.Rebin(fptbinsJN,'hDJsmReb',array.array('d',fptbinsJlh))
#hIJsmReb = hIJsm2.Rebin(fptbinsJN,'hIJsmReb',array.array('d',fptbinsJlh))

###########--- DIVISION
##--- minimum bias/pythia 8
hR_mb = hDJmbReb.Clone("hR_mb") #hR_mb0= hDJmb.Clone("hR_mb0")
hR_mb.Divide(hIJmbReb) #hR_mb0.Divide(hIJmbReb)
###--- soft mode 2/pythia 8
#hR_sm2 = hDJsm2.Clone("hR_sm2")
#hR_sm  = hDJsmReb.Clone("hR_sm")
#hR_sm2.Divide(hIJsm2)
#hR_sm.Divide(hIJsmReb)
##--- powheg pythia 8
hR_pow8 = hDJpowReb.Clone("hR_pow8")
hR_pow8Err = [hDJpowRebErr[0].Clone("hR_pow8ErrUp"), hDJpowRebErr[1].Clone("hR_pow8ErrDo")]
hR_pow8.Divide(hIJpowReb)
hR_pow8Err[0].Divide(hIJpowReb)
hR_pow8Err[1].Divide(hIJpowReb)
### style
hR_mb.SetLineColor(ROOT.kRed+2)
hR_mb.SetLineWidth(3)
hR_Data.SetLineWidth(3)

hR_Data.GetYaxis().SetTitle(" cross section ratio for R=0."+str(Rpar)+", #frac{cross section_{D-jets}}{cross section_{inclusive jets}}")

leg = ROOT.TLegend(0.5, 0.6, 0.9, 0.9)
################################################### TCANVAS
#hIJmbRebScale = hIJmbReb.Clone("hIJmbRebScale")
#hIJmbRebRaw.Divide(hIJmbRebScale)
##hIJmbRebRaw = hIJmbRebRaw.Scale(1,"width")
cRatio2 = ROOT.TCanvas("cRatio2","cRatio2",800,500)
#hIJpowReb.Draw()
hDJpowReb.Draw()
#hIJmbRebRaw.Draw()
##hIJmbReb.Draw("same")

for i in range(hDJpow8[0].GetNbinsX()):
    hR_pow8Err[0].SetBinContent(i+1,hR_pow8Err[0].GetBinContent(i+1)+hR_pow8Err[0].GetBinError(i+1))
    hR_pow8Err[0].SetBinError(i+1,0)
    hR_pow8Err[1].SetBinContent(i+1,hR_pow8Err[1].GetBinContent(i+1)-hR_pow8Err[1].GetBinError(i+1))
    hR_pow8Err[1].SetBinError(i+1,0)
################################################### TCANVAS
cRatio = ROOT.TCanvas("cRatio","cRatio",800,500)
hR_Data.GetYaxis().SetRangeUser(0.0,0.44)
hR_Data.Draw()
hR_mb.Draw("same")
hR_pow8.Draw("same");hR_pow8Err[0].Draw("same");hR_pow8Err[1].Draw("same")
#

leg.AddEntry(hR_Data,"Data")
leg.AddEntry(hR_mb,"Pythia 8 Monash 2013")
leg.AddEntry(hR_pow8,"POWHEG+Pythia 8")
leg.Draw("same")

cRatio.SaveAs("plots/DvIncJets.pdf")
cRatio.SaveAs("plots/DvIncJets.png")

wait=input()
