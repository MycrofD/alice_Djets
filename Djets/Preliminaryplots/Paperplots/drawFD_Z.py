import os, os.path, sys
import ROOT
ROOT.TH1.AddDirectory(False)
import style_settings
import array

Colors = [ROOT.kRed+1,ROOT.kBlue+2,ROOT.kGreen+2,ROOT.kViolet+2]
Markers = [20,21,22,23]
plotmin, plotmax = 0.4,1.0
plotYmin, plotYmax = 0.00001,1.0

ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetPadLeftMargin(0.135)
ROOT.gStyle.SetPadRightMargin(0.03)

############
if len(sys.argv) != 2:
    print("""
            == usage example: python file.py R
            R: 2,3,4,6
            """)
### SETTINGS
R=int(sys.argv[1]) #R=2,3,4,6
energy = "5.02"
jetbinnames=['5_7','7_10','10_15','15_50']
jetptLegends=["5 < #it{p}_{T,jet}^{ch} < 7 GeV/#it{c}",
        "7 < #it{p}_{T,jet}^{ch} < 10 GeV/#it{c}",
        "10 < #it{p}_{T,jet}^{ch} < 15 GeV/#it{c}",
        "15 < #it{p}_{T,jet}^{ch} < 50 GeV/#it{c}"]
fptbinsJN = [5,5,5,5]
fptbinsJlh=[
        [0.4,0.6,0.7,0.8,0.9,1.0],
        [0.4,0.6,0.7,0.8,0.9,1.0],
        [0.4,0.6,0.7,0.8,0.9,1.0],
        [0.4,0.6,0.7,0.8,0.9,1.0],
        ]

toPlot = [0,3]

hEmpty = ROOT.TH1D("hE","hE",100,plotmin,plotmax)
## HIST SETTINGS
hEmpty.SetTitle('')
hEmpty.GetXaxis().SetTitle("#it{z}_{||}^{ch}");
hEmpty.GetYaxis().SetTitle("B Feed-Down Fraction");
hEmpty.GetXaxis().SetLabelSize(0.03);
hEmpty.GetXaxis().SetTitleSize(0.03);
hEmpty.GetXaxis().SetTitleOffset(1.);
hEmpty.GetYaxis().SetTitleOffset(1.3);
hEmpty.GetYaxis().SetLabelSize(0.03);
hEmpty.GetYaxis().SetTitleSize(0.03);
hEmpty.GetXaxis().SetRangeUser(plotmin,plotmax);
hEmpty.GetYaxis().SetRangeUser(plotYmin,plotYmax);
### READING INPUT FILE
inFiles, hFD, hFD_up, hFD_do = [], [], [], []
leg = ROOT.TLegend(0.65,0.50,0.85,0.80);
leg.SetTextSize(0.03);
for i in range(len(jetbinnames)):
    jetbinname = jetbinnames[i]
    inFiles.append(
            ROOT.TFile('/media/jackbauer/data/z_out/R_0'+str(R)+'_finaltry/FDsubtraction/Jetbin_'+str(jetbinname)+'/plots/JetPtSpectrum_FDsub.root','read')
            )
    hh = inFiles[i].Get('hFD_ratio').Clone("h"+str(i))
    hhUp = inFiles[i].Get('hFD_ratio_up').Clone("hup"+str(i))
    hhDo = inFiles[i].Get('hFD_ratio_down').Clone("hdo"+str(i))
    hFD.append(hh)
    hFD_up.append(hhUp)
    hFD_do.append(hhDo)
    ## HIST SETTINGS
    hFD[i].SetTitle('')
    hFD[i].SetMarkerColor(Colors[i]);
    hFD[i].SetLineColor(Colors[i]);
    hFD[i].SetMarkerStyle(Markers[i]);
    hFD[i].SetMarkerSize(1.2);
    hFD[i].GetXaxis().SetLabelSize(0.04);
    hFD[i].GetXaxis().SetTitleSize(0.05);
    hFD[i].GetXaxis().SetTitleOffset(1.);
    hFD[i].GetYaxis().SetTitleOffset(1.3);
    hFD[i].GetYaxis().SetLabelSize(0.045);
    hFD[i].GetYaxis().SetTitleSize(0.05);
    hFD[i].GetXaxis().SetRangeUser(fptbinsJlh[i][0],fptbinsJlh[i][-1]);
    hFD[i]=hFD[i].Rebin(fptbinsJN[i],'h_'+str(i),array.array('d',fptbinsJlh[i]))
    hFD[i].Scale(0.1,"width");
    #hFD[i].SetMaximum(hFD[i].GetMaximum()*2.5);
    if(i in toPlot):
        leg.AddEntry(hFD[i],jetptLegends[i],"p");

    hFD_up[i].SetLineStyle(1);
    hFD_up[i].SetMarkerColor(Colors[i]);
    hFD_up[i].SetLineColor(Colors[i]);
    hFD_up[i]=hFD_up[i].Rebin(fptbinsJN[i],'h_'+str(i),array.array('d',fptbinsJlh[i]))
    hFD_up[i].Scale(0.1,"width");
    hFD_do[i].SetLineStyle(1);
    hFD_do[i].SetMarkerColor(Colors[i]);
    hFD_do[i].SetLineColor(Colors[i]);
    hFD_do[i]=hFD_do[i].Rebin(fptbinsJN[i],'h_'+str(i),array.array('d',fptbinsJlh[i]))
    hFD_do[i].Scale(0.1,"width");
shift = -0.0;
pvALICE = ROOT.TPaveText(0.15,0.85,0.8,0.9,"brNDC");
pvALICE.SetFillStyle(0);
pvALICE.SetBorderSize(0);
pvALICE.SetTextFont(42);
pvALICE.SetTextSize(0.03);
pvALICE.SetTextAlign(11);
pvALICE.AddText("ALICE, pp, #sqrt{#it{s}} = "+energy+" TeV");

shift += 0.05;
pvJet = ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvJet.SetFillStyle(0);
pvJet.SetBorderSize(0);
pvJet.SetTextFont(42);
pvJet.SetTextSize(0.03);
pvJet.SetTextAlign(11);
pvJet.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0."+str(R)+", |#it{#eta}_{lab}^{jet}| < 0."+str(9-R));

shift += 0.05;
pvD = ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvD.SetFillStyle(0);
pvD.SetBorderSize(0);
pvD.SetTextFont(42);
pvD.SetTextSize(0.03);
pvD.SetTextAlign(11);
pvD.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and charge conj.");


#shift += 0.05;
#pvEta = ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
#pvEta.SetFillStyle(0);
#pvEta.SetBorderSize(0);
#pvEta.SetTextFont(42);
#pvEta.SetTextSize(0.03);
#pvEta.SetTextAlign(11);
#pvEta.AddText("|#it{#eta}_{lab}^{jet}| < 0."+str(9-R));

### CANVAS
cEff = ROOT.TCanvas("cEff","cEff",1000,800);
#cEff.SetLogy();
hEmpty.Draw()
#for i in range(len(jetbinnames)):
for i in toPlot:
    hFD[i].Draw("same");
    hFD_up[i].Draw("same");
    hFD_do[i].Draw("same");

pvALICE.Draw("same");
pvJet.Draw("same");
pvD.Draw("same");
#pvEta.Draw("same");
leg.Draw("same");

cEff.SaveAs('plots/zFDratio_R0'+str(R)+'_'+str(int(float(energy)))+'.pdf')
cEff.SaveAs('plots/zFDratio_R0'+str(R)+'_'+str(int(float(energy)))+'.png')
cEff.SaveAs('plots/zFDratio_R0'+str(R)+'_'+str(int(float(energy)))+'.eps')

input()
