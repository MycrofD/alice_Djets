import ROOT as RT
import style_settings


promptColor = RT.kRed+1
nonpromptColor = RT.kBlue+1

RT.gStyle.SetOptStat(000)
RT.gStyle.SetLegendFont(42)
RT.gStyle.SetPadLeftMargin(0.135)
RT.gStyle.SetPadRightMargin(0.03)

### SETTINGS
R=4
jetbin=1
ptjetbins = [5,7,10,15,50]
ptbinsDNlist = [5,6,6,6]
jetbinnames=['5_7','7_10','10_15','15_50']
jetbinname = jetbinnames[jetbin-1]

### READING INPUT FILE
inFilePrompt = RT.TFile('/media/jackbauer/data/z_out/R_0'+str(R)+'_finaltry/efficiency/DjetEff_prompt_jetpt'+str(jetbinname)+'.root','read')
inFileFD = RT.TFile('/media/jackbauer/data/z_out/R_0'+str(R)+'_finaltry/efficiency/DjetEff_nonPrompt_jetpt'+str(jetbinname)+'.root','read')

hEffPrompt=inFilePrompt.Get('hEff_reb')
hEffNonPrompt=inFileFD.Get('hEff_reb')

## PROMPT HIST SETTINGS
hEffPrompt.SetTitle('')
hEffPrompt.SetMarkerColor(promptColor);
hEffPrompt.SetLineColor(promptColor);
hEffPrompt.SetMarkerStyle(20);
hEffPrompt.SetMarkerSize(1.2);
hEffPrompt.GetXaxis().SetTitle("#it{p}_{T,D^{0}} (GeV/#it{c})");
hEffPrompt.GetYaxis().SetTitle("Acceptance #times Efficiency");
hEffPrompt.GetXaxis().SetLabelSize(0.04);
hEffPrompt.GetXaxis().SetTitleSize(0.05);
hEffPrompt.GetXaxis().SetTitleOffset(1.);
hEffPrompt.GetYaxis().SetTitleOffset(1.3);
hEffPrompt.GetYaxis().SetLabelSize(0.045);
hEffPrompt.GetYaxis().SetTitleSize(0.05);
#hEffPrompt.GetXaxis().SetRangeUser(plotmin,plotmax);
hEffPrompt.SetMaximum(hEffPrompt.GetMaximum()*2.5);


## NON PROMPT HIST SETTINGS
hEffNonPrompt.SetTitle('');
hEffNonPrompt.SetMarkerColor(nonpromptColor);
hEffNonPrompt.SetLineColor(nonpromptColor);
hEffNonPrompt.SetMarkerStyle(21);
hEffNonPrompt.SetMarkerSize(1.2);
hEffNonPrompt.GetXaxis().SetTitle("#it{p}_{T}^{D} (GeV/#it{c})");
hEffNonPrompt.GetYaxis().SetTitle("Acceptance #times Efficiency");
hEffNonPrompt.GetXaxis().SetLabelSize(0.04);
hEffNonPrompt.GetXaxis().SetTitleSize(0.05);
hEffNonPrompt.GetXaxis().SetTitleOffset(1.);
hEffNonPrompt.GetYaxis().SetLabelSize(0.045);
hEffNonPrompt.GetYaxis().SetTitleSize(0.05);
#hEffNonPrompt.GetXaxis().SetRangeUser(plotmin,plotmax);
hEffNonPrompt.SetMaximum(hEffNonPrompt.GetMaximum()*2);

leg = RT.TLegend(0.15,0.60,0.45,0.70);
leg.SetTextSize(0.045);
leg.AddEntry(hEffPrompt,"Prompt D^{0}","p");
leg.AddEntry(hEffNonPrompt,"Feed-down D^{0}","p");

pvALICE = RT.TPaveText(0.15,0.85,0.8,0.9,"brNDC");
pvALICE.SetFillStyle(0);
pvALICE.SetBorderSize(0);
pvALICE.SetTextFont(42);
pvALICE.SetTextSize(0.045);
pvALICE.SetTextAlign(11);
pvALICE.AddText("ALICE Preliminary");

pvEn= RT.TPaveText(0.15,0.80,0.8,0.85,"brNDC");
pvEn.SetFillStyle(0);
pvEn.SetBorderSize(0);
pvEn.SetTextFont(42);
pvEn.SetTextSize(0.045);
pvEn.SetTextAlign(11);
pvEn.AddText("pp, #sqrt{#it{s}} = 5.02 TeV");


shift = -0.15;
pvD = RT.TPaveText(0.46,0.66-shift,0.9,0.70-shift,"brNDC");
pvD.SetFillStyle(0);
pvD.SetBorderSize(0);
pvD.SetTextFont(42);
pvD.SetTextSize(0.045);
pvD.SetTextAlign(11);
pvD.AddText("D^{0} #rightarrow K^{-}#pi^{+} and charge conj.");

pvJet = RT.TPaveText(0.46,0.6-shift,0.9,0.64-shift,"brNDC");
pvJet.SetFillStyle(0);
pvJet.SetBorderSize(0);
pvJet.SetTextFont(42);
pvJet.SetTextSize(0.045);
pvJet.SetTextAlign(11);
pvJet.AddText("in charged jets, anti-#it{k}_{T}, #it{R} = 0."+str(R));

pvjetpt = RT.TPaveText(0.465,0.54-shift,0.8,0.59-shift,"brNDC");
pvjetpt.SetFillStyle(0);
pvjetpt.SetBorderSize(0);
pvjetpt.SetTextFont(42);
pvjetpt.SetTextSize(0.045);
pvjetpt.SetTextAlign(11);
pvjetpt.AddText("%.0f < #it{p}_{T, ch jet} < %.0f GeV/#it{c}"%(ptjetbins[jetbin-1], ptjetbins[jetbin]));

pvEta = RT.TPaveText(0.465,0.47-shift,0.8,0.52-shift,"brNDC");
pvEta.SetFillStyle(0);
pvEta.SetBorderSize(0);
pvEta.SetTextFont(42);
pvEta.SetTextSize(0.045);
pvEta.SetTextAlign(11);
pvEta.AddText("|#it{#eta}_{lab}^{jet}| < 0."+str(9-R));

### CANVAS
cEff = RT.TCanvas("cEff","cEff",1000,800);
#cEff.SetLogy();
hEffPrompt.Draw();
hEffNonPrompt.Draw("same");

pvALICE.Draw("same");
pvEn.Draw("same");
pvJet.Draw("same");
pvD.Draw("same");
pvjetpt.Draw("same");
pvEta.Draw("same");
leg.Draw("same");

cEff.SaveAs('Efficiency_R0'+str(R)+'_'+str(jetbin)+'.pdf')
cEff.SaveAs('Efficiency_R0'+str(R)+'_'+str(jetbin)+'.png')
cEff.SaveAs('Efficiency_R0'+str(R)+'_'+str(jetbin)+'.eps')

input()
