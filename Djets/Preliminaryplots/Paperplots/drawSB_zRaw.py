import os, os.path, sys
import ROOT
from style_settings import *
import array

def setHistoDetails(h, color, mstyle, size, width):
    h.SetMarkerStyle(mstyle)
    h.SetMarkerColor(color)
    h.SetMarkerSize(size)
    h.SetLineColor(color)
    h.SetLineWidth(width)
    h.SetTitle(0)
    h.GetXaxis().SetTitle('#it{z}_{||}^{ch}')

def SaveCanvas(c, name):
    #c.SaveAs(name+".png");
    c.SaveAs(name+".pdf");
    c.SaveAs(name+".eps");
    #c.SaveAs(name+".root");
    #c.SaveAs(name+".C");
############
if len(sys.argv) != 3:
    print("""
            == usage example: python file.py R jetbin
            R: 2,3,4,6
            jetbin: 1,2,3,4 (5-7, 7-10, 10-15, 15-50)
            """)


### SETTINGS
energy = "5.02"
Rpar=int(sys.argv[1])  # R=2,3,4,5
jetbin =int(sys.argv[2])   # 1,2,3,4:   ## 5-7-10-15-50
ptjetbins = [5, 7, 10, 15, 50]
bin1, bin2, bin3 = 1,2,3
islog = 0
ptbinsDNlist = [5,6,6,6];
jetbinnames=['5_7','7_10','10_15','15_50']
ptDbins_list = [[2,3,4,5,6,7], [3,4,5,6,7,8,10], [5,6,7,8,10,12,15], [5,6,7,8,10,12,16,24,36]]
jetbinname=jetbinnames[jetbin-1]
ptbinsDN = ptbinsDNlist[jetbin-1]
ptDbins = ptDbins_list[jetbin-1]

ptbinsJetN = 4;
ptJetbins = [0.4,0.6,0.8,0.9,1.0]

massColor = ROOT.kBlack; signalColor = ROOT.kRed+1; SBColor = ROOT.kGreen+3; subColor = ROOT.kBlue+1;
markersize = 2.;
markerstyle= [24,25,27]

ltextsize = 0.06;
ROOT.gStyle.SetOptStat(000);ROOT.gStyle.SetLegendFont(42);ROOT.gStyle.SetPadLeftMargin(0.15);ROOT.gStyle.SetPadRightMargin(0.02);ROOT.gStyle.SetPadTopMargin(0.02);

### READING INPUT FILE
inFile = ROOT.TFile('/media/jackbauer/data/z_out/R_0'+str(Rpar)+'_finaltry/signalExtraction/plots/Z0to102_jetbin_'+str(jetbinname)+'/JetPtSpectra_SB_eff.root','read')
hjetpt=[] #ptbinsDN;
hjetpt_s=[] #ptbinsDN;
hjetptsub=[] #ptbinsDN
hjetptcorr=[] #ptbinsDN

for i in range(ptbinsDN):
    hjet = inFile.Get("hjetpt_"+str(i)).Clone("hjetpt_"+str(i))
    hjetReb=hjet.Rebin(ptbinsJetN,"hjetptReb_"+str(i),array.array('d',ptJetbins))
    hjetpt.append(hjetReb)
    hjetpt[i].SetTitle("")
    hjetpt[i].SetMarkerColor(signalColor);
    hjetpt[i].SetLineColor(signalColor);
    hjetpt[i].SetLineWidth(2);
    hjetpt[i].SetMarkerStyle(markerstyle[0]);
    hjetpt[i].SetMarkerSize(markersize);

    hjetpt[i].GetXaxis().SetTitle('#it{z}_{||}^{ch}')
    hjetpt[i].GetYaxis().SetTitle("Entries");
    hjetpt[i].GetXaxis().SetLabelSize(0.045);
    hjetpt[i].GetXaxis().SetTitleSize(0.055);
    hjetpt[i].GetYaxis().SetLabelSize(0.045);
    hjetpt[i].GetYaxis().SetTitleSize(0.055);
    hjetpt[i].GetYaxis().SetTitleOffset(1.4);

    hjetpt[i].SetMaximum(hjetpt[i].GetMaximum()*2.2);
    #hjetpt[i].SetMaximum(1600);
    if(islog):
        hjetpt[i].SetMaximum(hjetpt[i].GetMaximum()*10);
        hjetpt[i].SetMinimum(1);
    hjetpt[i].SetMinimum(-hjetpt[i].GetMaximum()*0.10);

    hjets = inFile.Get("hjetpt_sb_"+str(i)).Clone("hjetpt_sb_"+str(i))
    hjetsReb = hjets.Rebin(ptbinsJetN,"hjetpt_s_"+str(i),array.array('d',ptJetbins))
    hjetpt_s.append(hjetsReb)
    hjetpt_s[i].SetTitle("");
    hjetpt_s[i].SetMarkerColor(SBColor);
    hjetpt_s[i].SetLineColor(SBColor);
    hjetpt_s[i].SetLineWidth(2);
    hjetpt_s[i].SetMarkerStyle(markerstyle[1]);
    hjetpt_s[i].SetMarkerSize(markersize);

    hjetsub = inFile.Get("hjetptsub_"+str(i)).Clone("hjetptsub_"+str(i))
    hjetsubReb = hjetsub.Rebin(ptbinsJetN,"hjetptsub_"+str(i),array.array('d',ptJetbins))
    hjetptsub.append(hjetsubReb)
    hjetptsub[i].SetTitle("")
    hjetptsub[i].SetMarkerColor(subColor);
    hjetptsub[i].SetLineColor(subColor);
    hjetptsub[i].SetLineWidth(2);
    hjetptsub[i].SetMarkerStyle(markerstyle[2]);
    hjetptsub[i].SetMarkerSize(markersize+1.8);


if(islog):
    hjetpt[bin3].SetMaximum(hjetpt[bin3].GetMaximum()*1000);
legBands1 = ROOT.TLegend(0.45,0.88,0.85,0.92)#0.15,0.76,0.7,0.86)
legBands1.SetTextSize(ltextsize)
legBands1.SetFillStyle(0);
legBands1.SetTextAlign(13);
#legBands1.AddEntry(hjetpt[0],"#splitline{Signal region}{|#it{M}(K#pi)-#it{M}_{D^{0}}|<2#sigma}","p");
legBands1.AddEntry(hjetpt[0],"Signal region","p");

legBands2 = ROOT.TLegend(0.45,0.88,0.85,0.92)#0.15,0.65,0.7,0.75)
legBands2.SetTextSize(ltextsize)
legBands2.SetFillStyle(0)
legBands2.SetTextAlign(13)
#legBands2.AddEntry(hjetpt_s[0],"#splitline{Side bands (SB)}{4#sigma<|#it{M}(K#pi)-#it{M}_{D^{0}}|<9#sigma}","p");
legBands2.AddEntry(hjetpt_s[0],"Side bands (SB)","p");

legBands3 = ROOT.TLegend(0.45,0.8,0.85,0.85)#0.15,0.55,0.7,0.6)
legBands3.SetTextSize(ltextsize)
legBands3.SetFillStyle(0)
legBands3.SetTextAlign(13)
legBands3.AddEntry(hjetptsub[0],"Signal - SB","p")

pvALICE = ROOT.TPaveText(0.187,0.88,0.6,0.92,"brNDC")
pvALICE.SetFillStyle(0);
pvALICE.SetBorderSize(0);
pvALICE.SetTextFont(42);
pvALICE.SetTextSize(ltextsize);
pvALICE.SetTextAlign(11);
pvALICE.AddText("ALICE")

pvEn = ROOT.TPaveText(0.25,0.88,0.8,0.92,"brNDC")
pvEn.SetFillStyle(0);
pvEn.SetBorderSize(0);
pvEn.SetTextFont(42);
pvEn.SetTextSize(ltextsize);
pvEn.SetTextAlign(11);
pvEn.AddText("pp, #sqrt{#it{s}} = "+energy+" TeV");

pvD = ROOT.TPaveText(0.25,0.88,0.55,0.92,"brNDC")#0.187,0.82,0.55,0.87,"brNDC")
pvD.SetFillStyle(0)
pvD.SetBorderSize(0)
pvD.SetTextFont(42)
pvD.SetTextSize(ltextsize)
pvD.SetTextAlign(11)
pvD.AddText("D^{0} #rightarrow K^{-}#pi^{+} and charge conj.")

pvJet = ROOT.TPaveText(0.25,0.80,0.55,0.84,"brNDC")#0.189,0.77,0.55,0.82,"brNDC")
pvJet.SetFillStyle(0);
pvJet.SetBorderSize(0);
pvJet.SetTextFont(42);
pvJet.SetTextSize(ltextsize);
pvJet.SetTextAlign(11);
pvJet.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0."+str(int(Rpar)));

pvEta = ROOT.TPaveText(0.6,0.63,0.85,0.67,"brNDC")#0.187,0.72,0.6,0.77,"brNDC");
pvEta.SetFillStyle(0);
pvEta.SetBorderSize(0);
pvEta.SetTextFont(42);
pvEta.SetTextSize(ltextsize);
pvEta.SetTextAlign(11);
pvEta.AddText("|#it{#eta}_{lab}^{jet}| < 0."+str(int(9-Rpar)));

pvpt1 = ROOT.TPaveText(0.58,0.66,0.85,0.71,"brNDC")
pvpt1.SetFillStyle(0);
pvpt1.SetBorderSize(0);
pvpt1.SetTextFont(42);
pvpt1.SetTextSize(ltextsize);
pvpt1.AddText("%.0f < #it{p}_{T,D^{0}} < %.0f GeV/#it{c}"%(ptDbins[bin1],ptDbins[bin1+1]))

pvpt2 = ROOT.TPaveText(0.58,0.66,0.85,0.71,"brNDC");
pvpt2.SetFillStyle(0);
pvpt2.SetBorderSize(0);
pvpt2.SetTextFont(42);
pvpt2.SetTextSize(ltextsize);
pvpt2.AddText("%.0f < #it{p}_{T,D^{0}} < %.0f GeV/#it{c}"%(ptDbins[bin2],ptDbins[bin2+1]))

pvpt3 = ROOT.TPaveText(0.31,0.88,0.65,0.92,"brNDC")
pvpt3.SetFillStyle(0);
pvpt3.SetBorderSize(0);
pvpt3.SetTextFont(42);
pvpt3.SetTextSize(ltextsize);
pvpt3.AddText("%.0f < #it{p}_{T,D^{0}} < %.0f GeV/#it{c}"%(ptDbins[bin3],ptDbins[bin3+1]))


#cMass = ROOT.TCanvas("cMass","cMass",2160,1008)
cMass = ROOT.TCanvas("cMass","cMass",3000,587)
cMass.Divide(3,1);
cMass.cd(1);
ROOT.gPad.SetLogy(islog);

# the DrawFrame method takes the parameters (x1,y1,x2,y2).
hjetpt[bin1].Draw();
hjetpt[bin1].Draw("same");
hjetpt[bin1].Draw("same");
hjetpt_s[bin1].Draw("same");
hjetpt_s[bin1].Draw("same");
hjetpt_s[bin1].Draw("same");
hjetptsub[bin1].Draw("same");
hjetptsub[bin1].Draw("same");
hjetptsub[bin1].Draw("same");

#change this line and leave out the "A" for axis.
#pvpt1.Draw("same");
#pvALICE.Draw("same");
pvD.Draw("same");
pvJet.Draw("same");
pvEta.Draw("same");

cMass.cd(2);
ROOT.gPad.SetLogy(islog);
hjetpt[bin2].Draw();
hjetpt[bin2].Draw("same");
hjetpt_s[bin2].Draw("same");
hjetpt_s[bin2].Draw("same");
hjetptsub[bin2].Draw("same");
hjetptsub[bin2].Draw("same");

legBands1.Draw("same");
#pvpt2.Draw("same");
#pvEn.Draw("same");

cMass.cd(3);
ROOT.gPad.SetLogy(islog);
hjetpt[bin3].Draw();
hjetpt[bin3].Draw("same");
hjetpt_s[bin3].Draw("same");
hjetpt_s[bin3].Draw("same");
hjetptsub[bin3].Draw("same");
hjetptsub[bin3].Draw("same");

#pvpt3.Draw("same");

legBands2.Draw("same");
legBands3.Draw("same");

if(islog):
    SaveCanvas(cMass,"plots/RawZ_R0"+str(Rpar)+"_"+str(jetbin)+"_"+str(int(float(energy)))+"TeV_Perf_log_2");
else:
    SaveCanvas(cMass,"plots/RawZ_R0"+str(Rpar)+"_"+str(jetbin)+"_"+str(int(float(energy)))+"TeV_Perf_2");
