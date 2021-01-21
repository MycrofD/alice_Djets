import os, os.path, sys
import ROOT 
import style_settings
import array
ROOT.TH1.AddDirectory(False)
############
############
def Rebin2D(name, h, nx, binx, ny, biny): #rebin in 2d variable size
    #TH2D * Rebin2D(const char* name, TH2D *h, int nx, const double *binx, int ny, const double *biny, bool crop)
    xaxis = h.GetXaxis();
    yaxis = h.GetYaxis();

    hre = ROOT.TH2D(name,name,nx,binx,ny,biny);
    hre.Sumw2();
    for i in range(1,xaxis.GetNbins()+1):
        for j in range(1,yaxis.GetNbins()+1):
            hre.Fill(xaxis.GetBinCenter(i),yaxis.GetBinCenter(j),h.GetBinContent(i,j));


    for j in range(hre.GetNbinsY()+2):
        hre.SetBinContent(0,j,0);
        hre.SetBinError(0,j,0);
        hre.SetBinContent(hre.GetNbinsX()+1,j,0);
        hre.SetBinError(hre.GetNbinsX()+1,j,0);

    return hre
############
if len(sys.argv) != 2:
    print("""
            == usage example: python file.py R
            R: 2,3,4,6
            """)
    exit()
### SETTINGS
plotmin,plotmax=2,50
fptbinsJN = 10
fptbinsJlh=[2,3,4,5,6,8,10,14,20,30,50]
R=int(sys.argv[1]) #R=2,3,4,6
energy = "5.02"
jetorz="jet"
jetbin=""

ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetPadBottomMargin(0.1)
ROOT.gStyle.SetPadTopMargin(1)
ROOT.gStyle.SetPadRightMargin(1.)
##########################
## the histograms
inFileJ = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(R)+"_paperCuts/Default/unfolding_Bayes_5/unfoldedSpectrum_unfoldedJetSpectrum.root","read")
inFileJresol = ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(R)+"_paperCuts/Default/ResponseMatrix/DetMatrix_prompt.root","read")


hResp = inFileJ.Get("fMatrixProd")
hResp.GetXaxis().SetLabelSize(0.04);
hResp.GetXaxis().SetTitleSize(0.04);
hResp.GetXaxis().SetTitleOffset(1.);
hResp.GetYaxis().SetTitleOffset(1.3);
hResp.GetYaxis().SetLabelSize(0.04);
hResp.GetYaxis().SetTitleSize(0.04);
#hResp.GetXaxis().SetRangeUser(2,100);
#hResp.GetYaxis().SetRangeUser(2,100);
#print(hResp.Integral())
#hResp.Scale(1./hResp.Integral())
## HIST SETTINGS
#hEmpty = ROOT.TH2D("hE","hE",100,plotmin,plotmax,100,plotmin,plotmax)
hEmpty = ROOT.TH2D("hE","hE",fptbinsJN,array.array('d',fptbinsJlh),fptbinsJN,array.array('d',fptbinsJlh))
hEmpty.SetTitle('')
#hEmpty.GetXaxis().SetTitle("#it{z}_{||}^{ch}");
#hEmpty.GetYaxis().SetTitle("B Feed-Down Fraction");
hEmpty.GetXaxis().SetLabelSize(0.04);
hEmpty.GetXaxis().SetTitleSize(0.04);
hEmpty.GetXaxis().SetTitleOffset(1.);
hEmpty.GetYaxis().SetTitleOffset(1.3);
hEmpty.GetYaxis().SetLabelSize(0.04);
hEmpty.GetYaxis().SetTitleSize(0.04);
#hEmpty.GetXaxis().SetRangeUser(plotmin,plotmax);
#hEmpty.GetYaxis().SetRangeUser(plotmin,plotmax);

##########################
y1,y2=0.8,0.85
x1,x2=0.12,0.8
pvALICE = ROOT.TPaveText(x1,y1,x2,y2,"brNDC");
pvALICE.SetFillStyle(0);
pvALICE.SetBorderSize(0);
pvALICE.SetTextFont(42);
pvALICE.SetTextSize(0.04);
pvALICE.SetTextAlign(11);
pvALICE.AddText("ALICE, pp, #sqrt{#it{s}} = 5.02 TeV");

shift = 0.05;
pvD = ROOT.TPaveText(x1,y1-shift,x2,y2-shift,"brNDC");
pvD.SetFillStyle(0);
pvD.SetBorderSize(0);
pvD.SetTextFont(42);
pvD.SetTextSize(0.04);
pvD.SetTextAlign(11);
pvD.AddText("D^{0} #rightarrow K^{-}#pi^{+} and charge conj.");

shift += 0.05;
pvJet = ROOT.TPaveText(x1,y1-shift,x2,y2-shift,"brNDC");
pvJet.SetFillStyle(0);
pvJet.SetBorderSize(0);
pvJet.SetTextFont(42);
pvJet.SetTextSize(0.04);
pvJet.SetTextAlign(11);
pvJet.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0."+str(R));

#shift += 0.05;
#pvjetpt = ROOT.TPaveText(0.15,0.85-shift,0.8,0.9-shift,"brNDC");
#pvjetpt.SetFillStyle(0);
#pvjetpt.SetBorderSize(0);
#pvjetpt.SetTextFont(42);
#pvjetpt.SetTextSize(0.04);
#pvjetpt.SetTextAlign(11);
#pvjetpt.AddText("")#%.0f < #it{p}_{T, ch jet} < %.0f GeV/#it{c}"%(ptjetbins[jetbin-1], ptjetbins[jetbin]));

shift += 0.05;
pvEta = ROOT.TPaveText(x1,y1-shift,x2,y2-shift,"brNDC");
pvEta.SetFillStyle(0);
pvEta.SetBorderSize(0);
pvEta.SetTextFont(42);
pvEta.SetTextSize(0.04);
pvEta.SetTextAlign(11);
pvEta.AddText("|#it{#eta}_{lab}^{jet}| < 0."+str(9-R));

### CANVAS
cResp = ROOT.TCanvas("cResp","cResp",950,800);
cResp.SetLogz()
#hEmpty.Draw('')
hResp.Draw('colz')

pvALICE.Draw("same");
pvJet.Draw("same");
pvD.Draw("same");
#pvjetpt.Draw("same");
pvEta.Draw("same");

cResp.SaveAs('plots/Resp_R0'+str(R)+jetorz+'_'+str(jetbin)+'.pdf')
cResp.SaveAs('plots/Resp_R0'+str(R)+jetorz+'_'+str(jetbin)+'.png')
cResp.SaveAs('plots/Resp_R0'+str(R)+jetorz+'_'+str(jetbin)+'.eps')
#################################################################################
#################################################################################
#################################################################################
y1,y2=0.8,0.85
x1,x2=0.12,0.8
textsize=0.032
pvALICE2 = ROOT.TPaveText(x1,y1,x2,y2,"brNDC");
pvALICE2.SetFillStyle(0);
pvALICE2.SetBorderSize(0);
pvALICE2.SetTextFont(42);
pvALICE2.SetTextSize(textsize);
pvALICE2.SetTextAlign(11);
pvALICE2.AddText("ALICE PYTHIA6");

shift = 0.05;
pvpp = ROOT.TPaveText(x1,y1-shift,x2,y2-shift,"brNDC");
pvpp.SetFillStyle(0);
pvpp.SetBorderSize(0);
pvpp.SetTextFont(42);
pvpp.SetTextSize(textsize);
pvpp.SetTextAlign(11);
pvpp.AddText("pp, #sqrt{#it{s}} = 5.02 TeV");

shift += 0.05;
pvD2 = ROOT.TPaveText(x1,y1-shift,x2,y2-shift,"brNDC");
pvD2.SetFillStyle(0);
pvD2.SetBorderSize(0);
pvD2.SetTextFont(42);
pvD2.SetTextSize(textsize);
pvD2.SetTextAlign(11);
pvD2.AddText("Prompt D^{0} #rightarrow K^{-}#pi^{+}");

shift += 0.05;
pvD3 = ROOT.TPaveText(x1,y1-shift,x2,y2-shift,"brNDC");
pvD3.SetFillStyle(0);
pvD3.SetBorderSize(0);
pvD3.SetTextFont(42);
pvD3.SetTextSize(textsize);
pvD3.SetTextAlign(11);
pvD3.AddText("and charge conj.");

shift += 0.05;
pvD4 = ROOT.TPaveText(x1,y1-shift,x2,y2-shift,"brNDC");
pvD4.SetFillStyle(0);
pvD4.SetBorderSize(0);
pvD4.SetTextFont(42);
pvD4.SetTextSize(textsize);
pvD4.SetTextAlign(11);
pvD4.AddText("#it{p}_{T, D^{0}} > 2 GeV/#it{c}");

pvjet1 = ROOT.TPaveText(x1+0.5,y1,x2+0.5,y2,"brNDC");
pvjet1.SetFillStyle(0);
pvjet1.SetBorderSize(0);
pvjet1.SetTextFont(42);
pvjet1.SetTextSize(textsize);
pvjet1.SetTextAlign(11);
pvjet1.AddText("Charged Jets");
#pvjet1.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0."+str(R));

pvjet2 = ROOT.TPaveText(x1+0.5,y1-0.05,x2+0.5,y2-0.05,"brNDC");
pvjet2.SetFillStyle(0);
pvjet2.SetBorderSize(0);
pvjet2.SetTextFont(42);
pvjet2.SetTextSize(textsize);
pvjet2.SetTextAlign(11);
pvjet2.AddText("Anti-#it{k}_{T}, #it{R} = 0."+str(R));

pvEta1 = ROOT.TPaveText(x1+0.5,y1-2*0.05,x2+0.5,y2-2*0.05,"brNDC");
pvEta1.SetFillStyle(0);
pvEta1.SetBorderSize(0);
pvEta1.SetTextFont(42);
pvEta1.SetTextSize(textsize);
pvEta1.SetTextAlign(11);
pvEta1.AddText("|#it{#eta}_{lab}^{jet}| < 0."+str(9-R));

hDel2D = inFileJresol.Get("hjetRecGenResol")
hDel2D_noeff= inFileJresol.Get("hjetRGResol")

def GetDeltaProb(hDel2D,binLimLow,binLimHig,markerstyle,markersize,color,legend):
    hDel = ROOT.TH1D("hDel","hDel",50,-1.0,1.0)
    hDel.GetXaxis().SetLabelSize(0.04);
    hDel.GetXaxis().SetTitleSize(0.04);
    hDel.GetXaxis().SetTitleOffset(1.);
    hDel.GetYaxis().SetTitleOffset(1.3);
    hDel.GetYaxis().SetLabelSize(0.04);
    hDel.GetYaxis().SetTitleSize(0.04);
    hDel.SetMarkerStyle(markerstyle);
    hDel.SetMarkerSize(markersize);
    hDel.SetMarkerColor(color);
    hDel.SetLineColor(color);
    hDel.GetYaxis().SetTitle("Probability Density");
    hDel.GetXaxis().SetTitle("#Delta_{#it{p}_{T}}");
    legend.AddEntry(hDel,str(binLimLow-1)+" < #it{p}_{T, gen jet}^{ch} < "+str(binLimHig)+" GeV/#it{c}","p");

    htemp = hDel2D.ProjectionX("",binLimLow,binLimHig,"")
    hDel.Add(hDel,htemp,1,1)

    if(hDel.Integral()!=0):
        hDel.Scale(1./hDel.Integral())
    hDel.Scale(1,"width")
    return hDel

leg = ROOT.TLegend(0.15,0.23,0.4,0.40);
leg.SetTextSize(textsize);
hbin1 = GetDeltaProb(hDel2D,6,6,24,1.2,ROOT.kRed+2,leg)
hbin2 = GetDeltaProb(hDel2D,9,10,25,1.2,ROOT.kBlue+2,leg)
hbin3 = GetDeltaProb(hDel2D,21,30,27,1.8,ROOT.kGreen+2,leg)
### CANVAS
cDel = ROOT.TCanvas("cDel","cDel",950,800);
cDel.SetLogz()
hbin1.GetYaxis().SetRangeUser(0,30)
hbin1.Draw()
hbin2.Draw('same')
hbin3.Draw('same')
pvALICE2.Draw('same')
pvpp.Draw('same')
pvD2.Draw('same')
pvD3.Draw('same')
pvD4.Draw('same')
pvjet1.Draw('same')
pvjet2.Draw('same')
pvEta1.Draw('same')
leg.Draw('same')
cDel.SaveAs('plots/Delta_R0'+str(R)+jetorz+'_.pdf')
#########################################################
leg2 = ROOT.TLegend(0.15,0.23,0.4,0.40);
leg2.SetTextSize(textsize);
hbin12 = GetDeltaProb(hDel2D_noeff,7,8,24,1.2,ROOT.kRed+2,leg2)
hbin22 = GetDeltaProb(hDel2D_noeff,9,10,25,1.2,ROOT.kBlue+2,leg2)
hbin32 = GetDeltaProb(hDel2D_noeff,21,30,27,1.8,ROOT.kGreen+2,leg2)
#### CANVAS
cDel2 = ROOT.TCanvas("cDel2","cDel2",950,800);
cDel2.SetLogz()
hbin12.GetYaxis().SetRangeUser(0,30)
hbin12.Draw()
hbin22.Draw('same')
hbin32.Draw('same')
pvALICE2.Draw('same')
pvpp.Draw('same')
pvD2.Draw('same')
pvD3.Draw('same')
pvD4.Draw('same')
pvjet1.Draw('same')
pvjet2.Draw('same')
pvEta1.Draw('same')
leg2.Draw('same')
cDel2.SaveAs('plots/Delta2_R0'+str(R)+jetorz+'_.pdf')

input()
