import os, os.path, sys
import ROOT
ROOT.TH1.AddDirectory(False)
import style_settings
import array

Colors = [ROOT.kRed+1,ROOT.kBlue+2,ROOT.kGreen+2,ROOT.kViolet+2]
Markers = [20,21,22,23]
plotmin, plotmax = 2,36
plotYmin, plotYmax = 0.00,0.4

ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetPadLeftMargin(0.135)
ROOT.gStyle.SetPadRightMargin(0.03)

### SETTINGS
R=int(sys.argv[1]) #R=2,3,4,6
PnP=int(sys.argv[2]) #0 for non-prompt, 1 for prompt
if(PnP):
    dirtext,labeltext="prompt","Prompt"
else:
    dirtext,labeltext="nonPrompt","Non-prompt"
energy = "5.02"
ptjetbins = [5,7,10,15,50]
ptbinsDNlist = [5,6,6,6]
jetbinnames=['5_7','7_10','10_15','15_50']
jetptLegends=["5 < #it{p}_{T,jet}^{ch} < 7 GeV/#it{c}",
        "7 < #it{p}_{T,jet}^{ch} < 10 GeV/#it{c}",
        "10 < #it{p}_{T,jet}^{ch} < 15 GeV/#it{c}",
        "15 < #it{p}_{T,jet}^{ch} < 50 GeV/#it{c}"]
fptbinsJN = [5,6,6,6]
fptbinsJlh=[
        [2,3,4,5,6,7],
        [3,4,5,6,7,8,10],
        [5,6,7,8,10,12,15],
        [5,8,10,12,16,24,36]
        ]

hEmpty = ROOT.TH1D("hE","hE",100,plotmin,plotmax)
## HIST SETTINGS
hEmpty.SetTitle('')
hEmpty.GetXaxis().SetTitle("#it{p}_{T,D^{0}} (GeV/#it{c})");
hEmpty.GetYaxis().SetTitle("Acceptance #times Efficiency");
hEmpty.GetXaxis().SetLabelSize(0.04);
hEmpty.GetXaxis().SetTitleSize(0.05);
hEmpty.GetXaxis().SetTitleOffset(1.);
hEmpty.GetYaxis().SetTitleOffset(1.3);
hEmpty.GetYaxis().SetLabelSize(0.045);
hEmpty.GetYaxis().SetTitleSize(0.05);
hEmpty.GetXaxis().SetRangeUser(plotmin,plotmax);
hEmpty.GetYaxis().SetRangeUser(plotYmin,plotYmax);
### READING INPUT FILE
inFiles, hEffPrompts = [], []
leg = ROOT.TLegend(0.6,0.17,0.7,0.40);
leg.SetTextSize(0.04);
for i in range(len(jetbinnames)):
    jetbinname = jetbinnames[i]
    inFiles.append(
            ROOT.TFile('/media/jackbauer/data/z_out/R_0'+str(R)+'_finaltry/efficiency/DjetEff_'+dirtext+'_jetpt'+str(jetbinname)+'.root','read')
            )
    hh = inFiles[i].Get('hEff_reb').Clone("h"+str(i))
    hh.GetXaxis().SetRangeUser(plotmin,plotmax);
    hEffPrompts.append(hh)
    ## HIST SETTINGS
    hEffPrompts[i].SetTitle('')
    hEffPrompts[i].SetMarkerColor(Colors[i]);
    hEffPrompts[i].SetLineColor(Colors[i]);
    hEffPrompts[i].SetMarkerStyle(Markers[i]);
    hEffPrompts[i].SetMarkerSize(1.3);
    hEffPrompts[i].GetXaxis().SetTitle("#it{p}_{T,D^{0}} (GeV/#it{c})");
    hEffPrompts[i].GetYaxis().SetTitle("Acceptance #times Efficiency");
    hEffPrompts[i].GetXaxis().SetLabelSize(0.04);
    hEffPrompts[i].GetXaxis().SetTitleSize(0.05);
    hEffPrompts[i].GetXaxis().SetTitleOffset(1.);
    hEffPrompts[i].GetYaxis().SetTitleOffset(1.3);
    hEffPrompts[i].GetYaxis().SetLabelSize(0.045);
    hEffPrompts[i].GetYaxis().SetTitleSize(0.05);
    hEffPrompts[i].GetXaxis().SetRangeUser(fptbinsJlh[i][0],fptbinsJlh[i][-1]);
    hEffPrompts[i].Rebin(fptbinsJN[i],'h_'+str(i),array.array('d',fptbinsJlh[i]))
    #hEffPrompts[i].SetMaximum(hEffPrompts[i].GetMaximum()*2.5);
    #hEffPrompts[i].SetMaximum(0.45);
    leg.AddEntry(hEffPrompts[i],jetptLegends[i],"p");

### LEGEND STUFF
textsize = 0.04
shift = -0.0;
pvALICE = ROOT.TPaveText(0.15,0.85,0.8,0.9,"brNDC");
pvALICE.SetFillStyle(0);
pvALICE.SetBorderSize(0);
pvALICE.SetTextFont(42);
pvALICE.SetTextSize(textsize);
pvALICE.SetTextAlign(11);
pvALICE.AddText("ALICE, PYTHIA 6, pp, #sqrt{#it{s}} = "+energy+" TeV");

#shift += 0.05;
#pvEn= ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
#pvEn.SetFillStyle(0);
#pvEn.SetBorderSize(0);
#pvEn.SetTextFont(42);
#pvEn.SetTextSize(textsize)
#pvEn.SetTextAlign(11);
#pvEn.AddText("PYTHIA 6, pp, #sqrt{#it{s}} = "+energy+" TeV");

shift += 0.05;
pvD = ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvD.SetFillStyle(0);
pvD.SetBorderSize(0);
pvD.SetTextFont(42);
pvD.SetTextSize(textsize)
pvD.SetTextAlign(11);
pvD.AddText(labeltext+" D^{0} #rightarrow K^{-}#pi^{+} and charge conj.");

shift += 0.05;
pvJet = ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvJet.SetFillStyle(0);
pvJet.SetBorderSize(0);
pvJet.SetTextFont(42);
pvJet.SetTextSize(textsize)
pvJet.SetTextAlign(11);
pvJet.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0."+str(R));


shift += 0.05;
pvEta = ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvEta.SetFillStyle(0);
pvEta.SetBorderSize(0);
pvEta.SetTextFont(42);
pvEta.SetTextSize(textsize)
pvEta.SetTextAlign(11);
pvEta.AddText("|#it{#eta}_{lab}^{jet}| < 0."+str(9-R));

print(hEffPrompts[1].GetBinContent(1))
### CANVAS
cEff = ROOT.TCanvas("cEff","cEff",1000,800);
#cEff.SetLogy();
hEmpty.Draw()
for i in range(len(jetbinnames)):
    hEffPrompts[i].Draw("same");

pvALICE.Draw("same");
#pvEn.Draw("same");
pvJet.Draw("same");
pvD.Draw("same");
pvEta.Draw("same");
leg.Draw("same");

cEff.SaveAs('plots/z_'+dirtext+'Efficiency_R0'+str(R)+'_'+str(int(float(energy)))+'.pdf')
cEff.SaveAs('plots/z_'+dirtext+'Efficiency_R0'+str(R)+'_'+str(int(float(energy)))+'.png')
cEff.SaveAs('plots/z_'+dirtext+'Efficiency_R0'+str(R)+'_'+str(int(float(energy)))+'.eps')

input()
