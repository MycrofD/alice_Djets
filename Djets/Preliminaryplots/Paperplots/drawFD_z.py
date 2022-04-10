# trying using TGraphAsymErrors
import os, os.path, sys
import ROOT
ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetBatch()
import style_settings
import array

Colors = [ROOT.kRed+1,ROOT.kBlue+2,ROOT.kGreen+2,ROOT.kViolet+2]
Styles = [3003, 3365, 3005, 3004]
#Markers = [20,23,22,21]
Markers = [71,72,75,85]
plotmin, plotmax = 0.4,1.0
plotYmin, plotYmax = 0.00001,1.0

ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetPadLeftMargin(0.135)
ROOT.gStyle.SetPadRightMargin(0.03)

############
if len(sys.argv) < 2:
    print("""
            == usage example for 0 and 3 from [0,1,2,3] jetpt bins: python file.py R
            R: 2,3,4,6
            == usage example for individual jetpt bins: python file.py R jetpt_index
            R: 2,3,4,6
            jetpt_index: 0,1,2,3
            """)
### SETTINGS
R=int(sys.argv[1]) #R=2,3,4,6
energy = "5.02"
jetbinnames=['5_7','7_10','10_15','15_50']
jetptLegends=["  5 < #it{p}_{T,ch jet} <   7 GeV/ #it{c}",
        "  7 < #it{p}_{T,ch jet} < 10 GeV/ #it{c}",
        "10 < #it{p}_{T,ch jet} < 15 GeV/ #it{c}",
        "15 < #it{p}_{T,ch jet} < 50 GeV/ #it{c}"]
fptbinsJN = 5
fptbinsJlh=[0.4,0.6,0.7,0.8,0.9,1.0]
ptD_min = [2,3,5,5]
if R == 2:
    ptD_min = [2,4,5,10]

toPlot = [0,3] #paper
if len(sys.argv) == 3:
    jetpt_index=int(sys.argv[2])
    toPlot = [jetpt_index]
x_title = "#it{z}_{||}^{ch}"
y_title = "b-hadron feed-down fraction"
jetLegFD="b-hadron feed-down fraction"
jetLegFDsys="Sys. unc. (POWHEG+PYTHIA6)"
textsize=0.035

### READING INPUT FILE
inFiles, hFD, hFD_up, hFD_do = [], [], [], []
leg = ROOT.TLegend(0.52,0.57,0.85,0.67);
leg.SetTextSize(textsize);
for i in range(len(jetbinnames)):
    inFiles.append(
            ROOT.TFile('/eos/user/a/amohanty/media/jackbauer/data/z_out/R_0'+str(R)+'_finaltry/FDsubtraction/Jetbin_'+str(jetbinnames[i])+'/plots/JetPtSpectrum_FDsub.root','read')
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
    hFD[i].SetMarkerSize(2);
    hFD[i]=hFD[i].Rebin(fptbinsJN,'h_'+str(i),array.array('d',fptbinsJlh))
    hFD[i].Scale(0.1,"width");
    if(i in toPlot):
        leg.AddEntry(hFD[i],jetptLegends[i],"p");

    hFD_up[i].SetLineStyle(1);
    hFD_up[i].SetMarkerColor(Colors[i]);
    hFD_up[i].SetLineColor(Colors[i]);
    hFD_up[i]=hFD_up[i].Rebin(fptbinsJN,'h_'+str(i),array.array('d',fptbinsJlh))
    hFD_up[i].Scale(0.1,"width");
    hFD_do[i].SetLineStyle(1);
    hFD_do[i].SetMarkerColor(Colors[i]);
    hFD_do[i].SetLineColor(Colors[i]);
    hFD_do[i]=hFD_do[i].Rebin(fptbinsJN,'h_'+str(i),array.array('d',fptbinsJlh))
    hFD_do[i].Scale(0.1,"width");


## theory uncertainty
##-------------------
grsys, ptval, ptvalunc =[], [], []
for i in range(fptbinsJN):
    ptval.append((fptbinsJlh[i]+fptbinsJlh[i+1])/2.0)
    ptvalunc.append((fptbinsJlh[i+1]-fptbinsJlh[i])/2.0)
for j in range(len(jetbinnames)):
    value,  valueerrup, valueerrdown = [], [], []
    for i in range(fptbinsJN):
        value.append(hFD[j].GetBinContent(hFD[j].GetXaxis().FindBin(ptval[i])))
        valueerrup.append(hFD_up[j].GetBinContent(hFD_up[j].GetXaxis().FindBin(ptval[i]))-hFD[j].GetBinContent(hFD_up[j].GetXaxis().FindBin(ptval[i])))
        valueerrdown.append(hFD[j].GetBinContent(hFD_up[j].GetXaxis().FindBin(ptval[i]))-hFD_do[j].GetBinContent(hFD_do[j].GetXaxis().FindBin(ptval[i])))

    # theory unc for all jetbins
    graphsys = ROOT.TGraphAsymmErrors(fptbinsJN,array.array('d',ptval),array.array('d',value),array.array('d',ptvalunc),array.array('d',ptvalunc),array.array('d',valueerrdown),array.array('d',valueerrup))
    grsys.append(graphsys)
    grsys[j].SetFillColor(Colors[j])
    grsys[j].SetLineColor(Colors[j])
    #grsys[j].SetLineWidth(1504)
    #grsys[j].SetFillStyle(0)
    grsys[j].SetFillStyle(Styles[j])


shift = -0.0;
pvALICE = ROOT.TPaveText(0.15,0.85,0.8,0.9,"brNDC");
pvALICE.SetFillStyle(0);
pvALICE.SetBorderSize(0);
pvALICE.SetTextFont(42);
pvALICE.SetTextSize(0.04);
pvALICE.SetTextAlign(11);
pvALICE.AddText("ALICE, pp, #sqrt{#it{s}} = "+energy+" TeV");

shift += 0.05;
pvJet = ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvJet.SetFillStyle(0);
pvJet.SetBorderSize(0);
pvJet.SetTextFont(42);
pvJet.SetTextSize(textsize);
pvJet.SetTextAlign(11);
pvJet.AddText("charged jets, anti-#scale[0.5]{ }#it{k}_{T}, #it{R} = 0."+str(R)+", | #it{#eta}_{ch jet}| < 0."+str(9-R));

shift += 0.05;
pvD = ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvD.SetFillStyle(0);
pvD.SetBorderSize(0);
pvD.SetTextFont(42);
pvD.SetTextSize(textsize);
pvD.SetTextAlign(11);
pvD.AddText("with D^{0} #rightarrow K^{#font[122]{-}}#pi^{+} and charge conj.");


shift += 0.05;
pvPtD = ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
pvPtD.SetFillStyle(0);
pvPtD.SetBorderSize(0);
pvPtD.SetTextFont(42);
pvPtD.SetTextSize(textsize);
pvPtD.SetTextAlign(11);
text_pvPtD = "#it{p}_{T, D^{0}} > 2 GeV/ #it{c}"
if len(sys.argv) == 3:
    text_pvPtD = "#it{p}_{T, D^{0}} > "+str(ptD_min[jetpt_index])+" GeV/ #it{c}"
pvPtD.AddText(text_pvPtD)
######
#pvHist = ROOT.TPaveText(0.58,0.65,0.85,0.75,"brNDC");
pvHist = ROOT.TPaveText(0.52,0.65,0.85,0.75,"brNDC");
pvHist.SetFillStyle(0);
pvHist.SetBorderSize(0);
pvHist.SetTextFont(42);
pvHist.SetTextSize(textsize);
pvHist.SetTextAlign(11);
pvHist.AddText(jetLegFD);

#pvSyst = ROOT.TPaveText(0.58,0.47,0.85,0.57,"brNDC");
pvSyst = ROOT.TPaveText(0.52,0.45,0.85,0.57,"brNDC");
pvSyst.SetFillStyle(0);
pvSyst.SetBorderSize(0);
pvSyst.SetTextFont(42);
pvSyst.SetTextSize(textsize);
pvSyst.SetTextAlign(11);
pvSyst.AddText(jetLegFDsys);
########MULTIGRAPH
mg = ROOT.TMultiGraph()
mg.SetTitle(";"+x_title+";"+y_title)
#mg.Add(grsys[toPlot[0]])
#mg.Add(grsys[toPlot[1]])
for i in range(len(toPlot)):
    mg.Add(grsys[toPlot[i]])
leg2 = ROOT.TLegend(0.52,0.35,0.88,0.5);
#leg2.AddEntry(grsys[toPlot[0]],jetptLegends[toPlot[0]],"f")
#leg2.AddEntry(grsys[toPlot[1]],jetptLegends[toPlot[1]],"f")
for i in range(len(toPlot)):
    leg2.AddEntry(grsys[toPlot[i]],jetptLegends[toPlot[i]],"f")
mg.GetXaxis().SetLabelSize(0.04);
mg.GetXaxis().SetTitleSize(0.045);
mg.GetXaxis().SetTitleOffset(1.);
mg.GetYaxis().SetTitleOffset(1.3);
mg.GetYaxis().SetLabelSize(0.04);
mg.GetYaxis().SetTitleSize(0.045);
mg.SetMinimum(0)
mg.SetMaximum(1)
### CANVAS
cEff = ROOT.TCanvas("cEff","cEff",1000,800);
mg.Draw("a2")
for i in toPlot:
    hFD[i].Draw("same");
    #hFD_up[i].Draw("same");
    #hFD_do[i].Draw("same");

pvALICE.Draw("same");
pvJet.Draw("same");
pvD.Draw("same");
pvPtD.Draw("same");
pvHist.Draw("same")
pvSyst.Draw("same")
leg.Draw("same");
leg2.Draw("same");

suffix=''
if len(sys.argv) == 3:
    suffix = '_'+str(jetpt_index)

cEff.SaveAs('plots/zFDratio_R0'+str(R)+'_'+str(int(float(energy)))+suffix+'.pdf')
cEff.SaveAs('plots/zFDratio_R0'+str(R)+'_'+str(int(float(energy)))+suffix+'.png')
cEff.SaveAs('plots/zFDratio_R0'+str(R)+'_'+str(int(float(energy)))+suffix+'.eps')
