import os, os.path, sys
import ROOT 
import style_settings
import array
ROOT.TH1.AddDirectory(False)
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
R=int(sys.argv[1]) #R=2,3,4,6
energy = "5.02"
jetorz="z"

#ROOT.gStyle.SetPalette(57,ROOT.nullptr,1.0)
ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadTopMargin(1)
ROOT.gStyle.SetPadRightMargin(1.)
##########################
## the histograms
inFileZ = ROOT.TFile("/media/jackbauer/data/z_out/R_0"+str(R)+"_finaltry/unfolding/Bayes/alljetz2D/responseobject.root","read")

hResp = inFileZ.Get("responseNow")
hResp.GetXaxis().SetTitle("#it{z}_{||}^{ch} in #it{p}_{T,jet}^{ch} (RECO)")
hResp.GetXaxis().SetLabelSize(0.0);
hResp.GetXaxis().SetTitleSize(0.04);
hResp.GetXaxis().SetTitleOffset(1.7);
hResp.GetYaxis().SetTitle("#it{z}_{||}^{ch} in #it{p}_{T,jet}^{ch} (GEN)")
hResp.GetYaxis().SetTitleOffset(1.5);
hResp.GetYaxis().SetLabelSize(0.0);
hResp.GetYaxis().SetTitleSize(0.04);
print(hResp.Integral())
#hResp.Scale(1./hResp.Integral())

## HIST SETTINGS
####################################
y1,y2=0.83,0.88
x1,x2=0.12,0.8
textsize=0.035
pvALICE = ROOT.TPaveText(x1,y1,x2,y2,"brNDC");
pvALICE.SetFillStyle(0);
pvALICE.SetBorderSize(0);
pvALICE.SetTextFont(42);
pvALICE.SetTextSize(textsize);
pvALICE.SetTextAlign(11);
pvALICE.AddText("ALICE, pp, #sqrt{#it{s}} = 5.02 TeV");

shift = 0.05;
pvD = ROOT.TPaveText(x1,y1-shift,x2,y2-shift,"brNDC");
pvD.SetFillStyle(0);
pvD.SetBorderSize(0);
pvD.SetTextFont(42);
pvD.SetTextSize(textsize);
pvD.SetTextAlign(11);
pvD.AddText("D^{0} #rightarrow K^{-}#pi^{+} and charge conj.");

shift += 0.05;
pvJet = ROOT.TPaveText(x1,y1-shift,x2,y2-shift,"brNDC");
pvJet.SetFillStyle(0);
pvJet.SetBorderSize(0);
pvJet.SetTextFont(42);
pvJet.SetTextSize(textsize);
pvJet.SetTextAlign(11);
pvJet.AddText("Charged Jets, Anti-#it{k}_{T}");

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
pvEta.SetTextSize(textsize);
pvEta.SetTextAlign(11);
pvEta.AddText("#it{R} = 0."+str(R)+", |#it{#eta}_{lab}^{jet}| < 0."+str(9-R));

### CANVAS
cResp = ROOT.TCanvas("cResp","cResp",950,800);
cResp.SetLogz()
hResp.Draw('colz')
##########################
axisV1 = ROOT.TGaxis(0, 0, 0, 1, 0.4, 0.6, 1, "-L");axisV1.SetLabelFont(43);axisV1.SetLabelSize(17);axisV1.SetLabelOffset(0.03);axisV1.Draw();
axisV2 = ROOT.TGaxis(0, 2, 0, 5, 0.7, 1.0, 3, "-L");axisV2.SetLabelFont(43);axisV2.SetLabelSize(17);axisV2.SetLabelOffset(0.03);axisV2.Draw();
axisV3 = ROOT.TGaxis(0, 6, 0,10, 0.6, 1.0, 4, "-L");axisV3.SetLabelFont(43);axisV3.SetLabelSize(17);axisV3.SetLabelOffset(0.03);axisV3.Draw();
axisV4 = ROOT.TGaxis(0,11, 0,15, 0.6, 1.0, 4, "-L");axisV4.SetLabelFont(43);axisV4.SetLabelSize(17);axisV4.SetLabelOffset(0.03);axisV4.Draw();
axisV5 = ROOT.TGaxis(0,16, 0,20, 0.6, 1.0, 4, "-L");axisV5.SetLabelFont(43);axisV5.SetLabelSize(17);axisV5.SetLabelOffset(0.03);axisV5.Draw();
axisV6 = ROOT.TGaxis(0,21, 0,25, 0.6, 1.0, 4, "-L");axisV6.SetLabelFont(43);axisV6.SetLabelSize(17);axisV6.SetLabelOffset(0.03);axisV6.Draw();
axisV7 = ROOT.TGaxis(0, 0, 0, 5,   2,   5, 1, "-L");axisV7.SetLabelFont(43);axisV7.SetLabelSize(26);axisV7.SetLabelOffset(0.06);axisV7.Draw();
axisV8 = ROOT.TGaxis(0,10, 0,15,   7,  10, 1, "-L");axisV8.SetLabelFont(43);axisV8.SetLabelSize(26);axisV8.SetLabelOffset(0.06);axisV8.Draw();
axisV9 = ROOT.TGaxis(0,20, 0,25,  15,  50, 1, "-L");axisV9.SetLabelFont(43);axisV9.SetLabelSize(26);axisV9.SetLabelOffset(0.06);axisV9.Draw();

axisH1 = ROOT.TGaxis( 0, 0, 1, 0, 0.4,0.6, 1, "-+L");axisH1.SetLabelFont(43);axisH1.SetLabelSize(17);axisH1.SetLabelOffset(0.005);axisH1.Draw();
axisH2 = ROOT.TGaxis( 2, 0, 5, 0, 0.7,1.0, 3, "-+L");axisH2.SetLabelFont(43);axisH2.SetLabelSize(17);axisH2.SetLabelOffset(0.005);axisH2.Draw();
axisH3 = ROOT.TGaxis( 6, 0,10, 0, 0.6,1.0, 4, "-+L");axisH3.SetLabelFont(43);axisH3.SetLabelSize(17);axisH3.SetLabelOffset(0.005);axisH3.Draw();
axisH4 = ROOT.TGaxis(11, 0,15, 0, 0.6,1.0, 4, "-+L");axisH4.SetLabelFont(43);axisH4.SetLabelSize(17);axisH4.SetLabelOffset(0.005);axisH4.Draw();
axisH5 = ROOT.TGaxis(16, 0,20, 0, 0.6,1.0, 4, "-+L");axisH5.SetLabelFont(43);axisH5.SetLabelSize(17);axisH5.SetLabelOffset(0.005);axisH5.Draw();
axisH6 = ROOT.TGaxis(21, 0,25, 0, 0.6,1.0, 4, "-+L");axisH6.SetLabelFont(43);axisH6.SetLabelSize(17);axisH6.SetLabelOffset(0.005);axisH6.Draw();
axisH7 = ROOT.TGaxis( 0, 0, 5, 0,   2,  5, 1, "-+L");axisH7.SetLabelFont(43);axisH7.SetLabelSize(26);axisH7.SetLabelOffset(0.04); axisH7.Draw();
axisH8 = ROOT.TGaxis(10, 0,15, 0,   7, 10, 1, "-+L");axisH8.SetLabelFont(43);axisH8.SetLabelSize(26);axisH8.SetLabelOffset(0.04); axisH8.Draw();
axisH9 = ROOT.TGaxis(20, 0,25, 0,  15, 50, 1, "-+L");axisH9.SetLabelFont(43);axisH9.SetLabelSize(26);axisH9.SetLabelOffset(0.04); axisH9.Draw();
##########################
pvALICE.Draw("same");
pvJet.Draw("same");
pvD.Draw("same");
#pvjetpt.Draw("same");
pvEta.Draw("same");

cResp.SaveAs('plots/Resp_R0'+str(R)+jetorz+'_.pdf')
cResp.SaveAs('plots/Resp_R0'+str(R)+jetorz+'_.png')
cResp.SaveAs('plots/Resp_R0'+str(R)+jetorz+'_.eps')

####################################################
####################################################
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

pvD4 = ROOT.TPaveText(x1+0.02,y1-0.4,x2+0.02,y2-0.4,"brNDC");
pvD4.SetFillStyle(0);
pvD4.SetBorderSize(0);
pvD4.SetTextFont(42);
pvD4.SetTextSize(textsize);
pvD4.SetTextAlign(11);
pvD4.AddText("#it{p}_{T, D^{0}} > 2 GeV/#it{c}");

pvD4J = ROOT.TPaveText(x1+0.02,y1-0.35,x2+0.02,y2-0.35,"brNDC");
pvD4J.SetFillStyle(0);
pvD4J.SetBorderSize(0);
pvD4J.SetTextFont(42);
pvD4J.SetTextSize(textsize);
pvD4J.SetTextAlign(11);
pvD4J.AddText("5 < #it{p}_{T, gen jet}^{ch} < 7 GeV/#it{c}");

pvD5 = ROOT.TPaveText(x1+0.42,y1-0.4,x2+0.42,y2-0.4,"brNDC");
pvD5.SetFillStyle(0);
pvD5.SetBorderSize(0);
pvD5.SetTextFont(42);
pvD5.SetTextSize(textsize);
pvD5.SetTextAlign(11);
pvD5.AddText("#it{p}_{T, D^{0}} > 5 GeV/#it{c}");

pvD5J = ROOT.TPaveText(x1+0.42,y1-0.35,x2+0.42,y2-0.35,"brNDC");
pvD5J.SetFillStyle(0);
pvD5J.SetBorderSize(0);
pvD5J.SetTextFont(42);
pvD5J.SetTextSize(textsize);
pvD5J.SetTextAlign(11);
pvD5J.AddText("10 < #it{p}_{T, gen jet}^{ch} < 15 GeV/#it{c}");

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

hDel2D = inFileZ.Get("resolNow")
hDel2D_noeff = inFileZ.Get("resol_noeffNow")

def GetDeltaProb(hDel2D,binLimLow,binLimHig,markerstyle,markersize,color,legend,jetlow,jethig):
    hDelTemp = ROOT.TH1D("hDelTemp","hDelTemp",250,0,250)
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
    hDel.GetXaxis().SetTitle("#Delta_{#it{z}}");
    legend.AddEntry(hDel,str(jetlow)+" < #it{z}_{||, gen}^{ch} < "+str(jethig),"p");

    htemp = hDel2D.ProjectionX("",binLimLow+1,binLimHig,"")
    hDelTemp.Add(hDelTemp,htemp,1,1)

    for i in range(1,hDel.GetNbinsX()+1):
        bincontent = 0
        for j in range(int(hDel2D.GetNbinsX()/hDel.GetNbinsX())):
            bincontent += hDelTemp.GetBinContent(i+j*hDel.GetNbinsX())
        hDel.SetBinContent(i,bincontent)

    if(hDel.Integral()!=0):
        hDel.Scale(1./hDel.Integral())
    hDel.Scale(1,"width")

    return hDel


leg1 = ROOT.TLegend(0.15,0.23,0.4,0.40);
leg1.SetTextSize(textsize);
leg2 = ROOT.TLegend(0.55,0.23,0.8,0.40);
leg2.SetTextSize(textsize);
hbin1 = GetDeltaProb(hDel2D,5,6,24,1.2,ROOT.kRed+2,leg1,0.4,0.6)
hbin2 = GetDeltaProb(hDel2D,6,8,25,1.2,ROOT.kBlue+2,leg1,0.7,0.8)
hbin3 = GetDeltaProb(hDel2D,8,10,27,1.8,ROOT.kGreen+2,leg1,0.9,1.0)
hbin4 = GetDeltaProb(hDel2D,15,16,20,1.2,ROOT.kRed+2,leg2,0.4,0.6)
hbin5 = GetDeltaProb(hDel2D,16,17,21,1.2,ROOT.kBlue+2,leg2,0.6,0.7)
hbin6 = GetDeltaProb(hDel2D,18,19,33,1.8,ROOT.kGreen+2,leg2,0.8,0.9)
####################
######## CANVAS
cDel = ROOT.TCanvas("cDel","cDel",950,800);
cDel.SetLogy();hbin1.GetYaxis().SetRangeUser(0.0001,10000)
hbin1.Draw()
hbin2.Draw('same')
hbin3.Draw('same')
pvD4.Draw('same')
pvD4J.Draw('same')
leg1.Draw('same')
####=========
#hbin4.Draw('same')
#hbin5.Draw('same')
#hbin6.Draw('same')
#pvD5.Draw('same')
#pvD5J.Draw('same')
#leg2.Draw('same')
####=========
pvALICE2.Draw('same')
pvpp.Draw('same')
pvD2.Draw('same')
pvD3.Draw('same')
pvjet1.Draw('same')
pvjet2.Draw('same')
pvEta1.Draw('same')
cDel.SaveAs('plots/Delta_R0'+str(R)+jetorz+'_.pdf')
input()
