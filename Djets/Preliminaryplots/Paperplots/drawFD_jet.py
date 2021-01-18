import os, os.path, sys
import ROOT
ROOT.TH1.AddDirectory(False)
import style_settings
import array

Colors = [ROOT.kRed+1,ROOT.kBlue+2,ROOT.kGreen+2,ROOT.kViolet+2]
Markers = [20,21,22,23]
plotmin, plotmax = 5,50
plotYmin, plotYmax = 0.00001,1.0

ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetPadLeftMargin(0.135)
ROOT.gStyle.SetPadRightMargin(0.03)
# Graph settings
#------------------
#def TGRAPH(bins, histcenters, syscent, exlh, sysdo, sysup, graphcolor, fillstyle, histedges, histsetting, ytitle, xtitle):
#    #      bins, histXcenters, histYcenters, xerr, xerr, yerrdown, yerrup, color, fillstyle, edges of hist, base histogram
#    graph = ROOT.TGraphAsymmErrors(bins, array.array('d',histcenters), array.array('d',syscent), array.array('d',exlh), array.array('d',exlh), array.array('d',sysdo), array.array('d',sysup))
#    graph.SetFillColor(graphcolor)
#    graph.SetFillStyle(fillstyle)
#    graphhist = ROOT.TH1F("","",fptbinsZN,array.array('d',histedges))
#    graphhist.SetMinimum(histsetting.GetMaximum()*(-0.04))
#    graphhist.SetMaximum(histsetting.GetMaximum()*1.2)
#
#    graphhist.GetYaxis().SetTitle(ytitle)
#    graphhist.GetXaxis().SetTitle(xtitle)
##    graphhist.SetTitleSize(0.0,"XYZ")
##    graphhist.SetTitleOffset(1.6,"XYZ")
#    graphhist.SetTitle("")
#    graphhist.SetMinimum(0)
#
#    graph.SetHistogram(graphhist)
#    return [graph, graphhist]
############
if len(sys.argv) != 2:
    print("""
            == usage example: python file.py R
            R: 2,3,4,6
            """)
    exit()
### SETTINGS
R=int(sys.argv[1]) #R=2,3,4,6
energy = "5.02"
fptbinsJN = 7
fptbinsJlh=[5,6,8,10,14,20,30,50]
jetLegFD="Raw B Feed-Down Fraction"
jetLegFDsys="Sys. Unc. (POWHEG+PYTHIA6)"

hEmpty = ROOT.TH1D("hE","hE",100,plotmin,plotmax)
## HIST SETTINGS
hEmpty.SetTitle('')
hEmpty.GetXaxis().SetTitle("#it{z}_{||}^{ch}");
hEmpty.GetYaxis().SetTitle("B Feed-Down Fraction");
#hEmpty.GetXaxis().SetLabelSize(0.04);
#hEmpty.GetXaxis().SetTitleSize(0.04);
#hEmpty.GetXaxis().SetTitleOffset(1.);
#hEmpty.GetYaxis().SetTitleOffset(1.3);
#hEmpty.GetYaxis().SetLabelSize(0.04);
#hEmpty.GetYaxis().SetTitleSize(0.04);
hEmpty.GetXaxis().SetRangeUser(plotmin,plotmax);
hEmpty.GetYaxis().SetRangeUser(plotYmin,plotYmax);
### READING INPUT FILE
leg = ROOT.TLegend(0.15,0.6,0.65,0.75);
#pvALICE = ROOT.TPaveText(0.15,0.85,0.8,0.9,"brNDC");
leg.SetTextSize(0.035);

inFile=ROOT.TFile('/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0'+str(R)+'_paperCuts/Default/FDsubtraction/JetPtSpectrum_FDsub.root','read')
hh = inFile.Get('hFD_ratio').Clone("h")
hhUp = inFile.Get('hFD_ratio_up').Clone("hup")
hhDo = inFile.Get('hFD_ratio_down').Clone("hdo")
hFD=hh
hFD_up=hhUp
hFD_do=hhDo
## HIST SETTINGS
hFD.SetTitle('')
hFD.SetMarkerColor(Colors[0]);
hFD.SetLineColor(Colors[0]);
hFD.SetMarkerStyle(Markers[0]);
hFD.SetMarkerSize(1.2);
#hFD.GetXaxis().SetLabelSize(0.04);
#hFD.GetXaxis().SetTitleSize(0.05);
#hFD.GetXaxis().SetTitleOffset(1.);
#hFD.GetYaxis().SetTitleOffset(1.3);
#hFD.GetYaxis().SetLabelSize(0.045);
#hFD.GetYaxis().SetTitleSize(0.05);
#hFD.GetXaxis().SetRangeUser(fptbinsJlh[0],fptbinsJlh[-1]);
hFD=hFD.Rebin(fptbinsJN,'h_',array.array('d',fptbinsJlh))
#hFD.Scale(1,"width");
#hFD.SetMaximum(hFD.GetMaximum()*2.5);
leg.AddEntry(hFD,jetLegFD,"p");

hFD_up.SetLineStyle(1);
hFD_up.SetMarkerColor(Colors[0]);
hFD_up.SetLineColor(Colors[0]);
hFD_up=hFD_up.Rebin(fptbinsJN,'hup_',array.array('d',fptbinsJlh))
#hFD_up.Scale(0.1,"width");
hFD_do.SetLineStyle(1);
hFD_do.SetMarkerColor(Colors[0]);
hFD_do.SetLineColor(Colors[0]);
hFD_do=hFD_do.Rebin(fptbinsJN,'hdo_',array.array('d',fptbinsJlh))
#hFD_do.Scale(0.1,"width");
#leg.AddEntry(hFD_up,jetLegFDsys,"p");

## theory uncertainty
##-------------------
ptval, value, ptvalunc, valueerrup, valueerrdown = [], [], [], [], []
for i in range(fptbinsJN):
    ptval.append((fptbinsJlh[i]+fptbinsJlh[i+1])/2.0)
    ptvalunc.append((fptbinsJlh[i+1]-fptbinsJlh[i])/2.0)
    value.append(hFD.GetBinContent(hFD.GetXaxis().FindBin(ptval[i])))
    valueerrup.append(hFD_up.GetBinContent(hFD_up.GetXaxis().FindBin(ptval[i]))-hFD.GetBinContent(hFD_up.GetXaxis().FindBin(ptval[i])))
    valueerrdown.append(hFD.GetBinContent(hFD_up.GetXaxis().FindBin(ptval[i]))-hFD_do.GetBinContent(hFD_do.GetXaxis().FindBin(ptval[i])))

grsys = ROOT.TGraphAsymmErrors(fptbinsJN,array.array('d',ptval),array.array('d',value),array.array('d',ptvalunc),array.array('d',ptvalunc),array.array('d',valueerrdown),array.array('d',valueerrup))
grsys.SetFillColor(Colors[0])
grsys.SetLineColor(Colors[0])
grsys.SetFillStyle(0)
grsys.SetTitle('')
grsys.GetXaxis().SetTitle("#it{z}_{||}^{ch}");
grsys.GetYaxis().SetTitle("B Feed-Down Fraction");
grsys.GetXaxis().SetLabelSize(0.04);
grsys.GetXaxis().SetTitleSize(0.04);
grsys.GetXaxis().SetTitleOffset(1.);
grsys.GetYaxis().SetTitleOffset(1.3);
grsys.GetYaxis().SetLabelSize(0.04);
grsys.GetYaxis().SetTitleSize(0.04);
grsys.GetXaxis().SetRangeUser(plotmin,plotmax);
grsys.GetYaxis().SetRangeUser(plotYmin,plotYmax);

#grsys.SetHistogram(hFD)
grsys.SetMinimum(0)
grsys.SetMaximum(1)
leg.AddEntry(grsys,jetLegFDsys,"f");

## legend stuff
##----------------
shift = -0.0;
pvALICE = ROOT.TPaveText(0.15,0.85,0.8,0.9,"brNDC");
pvALICE.SetFillStyle(0);
pvALICE.SetBorderSize(0);
pvALICE.SetTextFont(42);
pvALICE.SetTextSize(0.035);
pvALICE.SetTextAlign(11);
pvALICE.AddText("ALICE, pp, #sqrt{#it{s}} = "+energy+" TeV");

shift += 0.05;
pvJet = ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvJet.SetFillStyle(0);
pvJet.SetBorderSize(0);
pvJet.SetTextFont(42);
pvJet.SetTextSize(0.035);
pvJet.SetTextAlign(11);
pvJet.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0."+str(R)+", |#it{#eta}_{lab}^{jet}| < 0."+str(9-R));

shift += 0.05;
pvD = ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvD.SetFillStyle(0);
pvD.SetBorderSize(0);
pvD.SetTextFont(42);
pvD.SetTextSize(0.035);
pvD.SetTextAlign(11);
pvD.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and charge conj.");

### CANVAS
cEff = ROOT.TCanvas("cEff","cEff",1000,800);
#cEff.SetLogy();
hEmpty.Draw()
grsys.Draw("a2same")
hFD.Draw("same");
#hFD_up.Draw("same");
#hFD_do.Draw("same");

pvALICE.Draw("same");
pvJet.Draw("same");
pvD.Draw("same");
#pvEta.Draw("same");
leg.Draw("same");

cEff.SaveAs('jetFDratio_R0'+str(R)+'_'+str(int(float(energy)))+'.pdf')
cEff.SaveAs('jetFDratio_R0'+str(R)+'_'+str(int(float(energy)))+'.png')
cEff.SaveAs('jetFDratio_R0'+str(R)+'_'+str(int(float(energy)))+'.eps')

input()
