import os, os.path, sys
import ROOT
ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetBatch()
import style_settings
import array


promptColor = ROOT.kRed+1;
nonpromptColor = ROOT.kBlue+1;

ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetPadLeftMargin(0.135)
ROOT.gStyle.SetPadRightMargin(0.03)

### SETTINGS
Rpar=int(sys.argv[1]) #R=2,3,4,6
energy = "5.02"
plotmin, plotmax = 2,36;
plotYmin, plotYmax = 0,0.45;
#markerstyle[20,21]
#markerstyle=[20,72]
markerstyle=[71,21]
markersize=2

inFilePrompt = ROOT.TFile("/eos/user/a/amohanty/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Default/efficiency/DjetEff_prompt_jetpt5_50.root","read");
inFileFD = ROOT.TFile("/eos/user/a/amohanty/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(Rpar)+"_paperCuts/Default/efficiency/DjetEff_nonPrompt_jetpt5_50.root","read");

hEffPrompt = inFilePrompt.Get("hEff_reb")
hEffNonPrompt = inFileFD.Get("hEff_reb")

hEffPrompt.SetTitle("");
hEffPrompt.SetMarkerColor(promptColor);
hEffPrompt.SetLineColor(promptColor);
hEffPrompt.SetMarkerStyle(markerstyle[0]);
hEffPrompt.SetMarkerSize(markersize);
hEffPrompt.GetXaxis().SetTitle("#it{p}_{T,D^{0}} (GeV/#it{c})");
#hEffPrompt.GetYaxis().SetTitle("Acceptance #times Efficiency");
hEffPrompt.GetYaxis().SetTitle("Acc #times #scale[1.25]{#varepsilon} ");
hEffPrompt.GetXaxis().SetLabelSize(0.045);
hEffPrompt.GetXaxis().SetTitleSize(0.05);
hEffPrompt.GetXaxis().SetTitleOffset(1.);
hEffPrompt.GetYaxis().SetTitleOffset(1.3);
hEffPrompt.GetYaxis().SetLabelSize(0.045);
hEffPrompt.GetYaxis().SetTitleSize(0.05);
hEffPrompt.GetXaxis().SetRangeUser(plotmin,plotmax);
#hEffPrompt.SetMaximum(hEffPrompt.GetMaximum()*3.5);#for logy version
#hEffPrompt.SetMaximum(hEffPrompt.GetMaximum()*1.5);#for linear y version
hEffPrompt.SetMaximum(plotYmax);
hEffPrompt.SetMinimum(plotYmin);

hEffNonPrompt.SetTitle("");
hEffNonPrompt.SetMarkerColor(nonpromptColor);
hEffNonPrompt.SetLineColor(nonpromptColor);
hEffNonPrompt.SetMarkerStyle(markerstyle[1]);
hEffNonPrompt.SetMarkerSize(markersize);
hEffNonPrompt.GetXaxis().SetTitle("#it{p}_{T}^{D^{0}} (GeV/#it{c})");
hEffNonPrompt.GetXaxis().SetLabelSize(0.045);
hEffNonPrompt.GetXaxis().SetTitleSize(0.05);
hEffNonPrompt.GetXaxis().SetTitleOffset(1.);
hEffNonPrompt.GetYaxis().SetLabelSize(0.045);
hEffNonPrompt.GetYaxis().SetTitleSize(0.05);
hEffNonPrompt.GetXaxis().SetRangeUser(plotmin,plotmax);
#hEffNonPrompt.SetMaximum(hEffNonPrompt.GetMaximum()*2);
hEffNonPrompt.SetMaximum(plotYmax);
hEffNonPrompt.SetMinimum(plotYmin);

leg = ROOT.TLegend(0.5,0.25,0.85,0.40);
leg.SetTextSize(0.04);
leg.AddEntry(hEffPrompt,"Prompt D^{0}","p");
leg.AddEntry(hEffNonPrompt,"Non-prompt D^{0}","p");

######
shift=0
textsize=0.04
pvALICE = ROOT.TPaveText(0.15,0.85,0.8,0.9,"brNDC");
pvALICE.SetFillStyle(0);
pvALICE.SetBorderSize(0);
pvALICE.SetTextFont(42);
pvALICE.SetTextSize(textsize);
pvALICE.SetTextAlign(11);
pvALICE.AddText("ALICE, pp, #sqrt{#it{s}} = 5.02 TeV");

#shift+=0.1
#pvEn= ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
#pvEn.SetFillStyle(0);
#pvEn.SetBorderSize(0);
#pvEn.SetTextFont(42);
#pvEn.SetTextSize(0.045);
#pvEn.SetTextAlign(11);
#pvEn.AddText("PYTHIA 6, pp, #sqrt{#it{s}} = 5.02 TeV");

shift+=0.05
pvJet= ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvJet.SetFillStyle(0);
pvJet.SetBorderSize(0);
pvJet.SetTextFont(42);
pvJet.SetTextSize(textsize);
pvJet.SetTextAlign(11);
pvJet.AddText("charged jets, anti-#scale[0.5]{ }#it{k}_{T}, #it{R} = 0."+str(Rpar));

shift+=0.05
pvD= ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvD.SetFillStyle(0);
pvD.SetBorderSize(0);
pvD.SetTextFont(42);
pvD.SetTextSize(textsize);
pvD.SetTextAlign(11);
pvD.AddText("with D^{0} #rightarrow K^{#font[122]{-}}#pi^{+} and charge conj.");

shift+=0.05
pvEta= ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvEta.SetFillStyle(0);
pvEta.SetBorderSize(0);
pvEta.SetTextFont(42);
pvEta.SetTextSize(textsize);
pvEta.SetTextAlign(11);
pvEta.AddText("|#it{#eta}_{ch jet}| < 0."+str(9-Rpar));

cEff = ROOT.TCanvas("cEff","cEff",1000,800);
hEffPrompt.Draw();
hEffNonPrompt.Draw("same");

pvALICE.Draw("same");
#pvEn.Draw("same");
pvJet.Draw("same");
pvD.Draw("same");
pvEta.Draw("same");
leg.Draw("same");


cEff.SaveAs("plots/jet_DjetEff_Sim_"+str(int(float(energy)))+".png");
cEff.SaveAs("plots/jet_DjetEff_Sim_"+str(int(float(energy)))+".pdf");
cEff.SaveAs("plots/jet_DjetEff_Sim_"+str(int(float(energy)))+".eps");

#input()
