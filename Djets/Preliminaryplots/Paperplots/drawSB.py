import os, os.path, sys
import ROOT
from style_settings import *
def setHistoDetails(h, color, mstyle, size, width):
    h.SetMarkerStyle(mstyle)
    h.SetMarkerColor(color)
    h.SetMarkerSize(size)
    h.SetLineColor(color)
    h.SetLineWidth(width)
    h.SetTitle(0)
    h.GetXaxis().SetTitle('#it{z}_{||}^{ch}')

############
if len(sys.argv) != 3:
    print("""
            == usage example: python file.py R jetbin
            R: 2,3,4,6
            jetbin: 1,2,3,4 (5-7, 7-10, 10-15, 15-50)
            """)

### SETTINGS
energy = "5.02"
R=int(sys.argv[1])  # R=2,3,4,5
jetbin =int(sys.argv[2])   # 1,2,3,4:   ## 5-7-10-15-50
ptjetbins = [5, 7, 10, 15, 50]
bin1, bin2, bin3 = 1,2,3
ptbinsDNlist = [5,6,6,6]  ## 5-7-10-15-50
jetbinnames=['5_7','7_10','10_15','15_50']
ptDbins_list = [[2,3,4,5,6,7], [3,4,5,6,7,8,10], [5,6,7,8,10,12,15], [5,6,7,8,10,12,16,24,36]]
jetbinname=jetbinnames[jetbin-1]
ptbinsDN = ptbinsDNlist[jetbin-1]
ptDbins = ptDbins_list[jetbin-1]


massColor = ROOT.kBlack; signalColor = ROOT.kRed+1; SBColor = ROOT.kGreen+3

ROOT.gStyle.SetOptStat(000);ROOT.gStyle.SetLegendFont(42);ROOT.gStyle.SetPadLeftMargin(0.15);ROOT.gStyle.SetPadRightMargin(0.02);ROOT.gStyle.SetPadTopMargin(0.02);

### READING INPUT FILE
inFile = ROOT.TFile('/media/jackbauer/data/z_out/R_0'+str(R)+'_finaltry/signalExtraction/plots/Z0to102_jetbin_'+str(jetbinname)+'/JetPtSpectra_SB_eff.root','read')

hmean = inFile.Get('hmean')
hsigma= inFile.Get('hsigma')
hsign = inFile.Get('hsign')
hsb   = inFile.Get('hsb')

hmass, hmass_l, hmass_u, hmass_c, fullfit, bkgfit = [], [], [], [], [], []
for i in range(ptbinsDN):
    hmass.append(inFile.Get('hmass_'+str(i)))
    hmass[i].SetTitle('')
    hmass[i].SetMarkerColor(massColor)
    hmass[i].SetLineColor(massColor)
    hmass[i].SetMarkerStyle(20)
    hmass[i].SetMarkerSize(1.2)
    hmass[i].GetXaxis().SetRangeUser(1.72,2.0)
    hmass[i].GetXaxis().SetTitle('#it{M}(K#pi) (GeV/#it{c}^{2})')
    hmass[i].GetYaxis().SetTitle('Entries/%.0f MeV/#it{c}^{2}'%(hmass[i].GetBinWidth(1)*1000))
    hmass[i].GetXaxis().SetLabelSize(0.045)
    hmass[i].GetXaxis().SetTitleSize(0.055)
    hmass[i].GetYaxis().SetLabelSize(0.045)
    hmass[i].GetYaxis().SetTitleSize(0.055)
    hmass[i].SetMaximum(hmass[i].GetMaximum()*1.3)
    hmass[i].SetMinimum(1)

    hmass_l.append(inFile.Get('hmass_l_'+str(i)))
    hmass_l[i].SetTitle('')
    hmass_l[i].SetMarkerColor(massColor)
    hmass_l[i].SetLineColor(SBColor)
    hmass_l[i].SetFillColor(SBColor)
    hmass_l[i].SetFillStyle(3004)
    hmass_l[i].SetLineWidth(1)

    hmass_u.append(inFile.Get('hmass_u_'+str(i)))
    hmass_u[i].SetTitle('')
    hmass_u[i].SetMarkerColor(massColor)
    hmass_u[i].SetLineColor(SBColor)
    hmass_u[i].SetFillColor(SBColor)
    hmass_u[i].SetFillStyle(3004)
    hmass_u[i].SetLineWidth(1)

    hmass_c.append(inFile.Get('hmass_c_'+str(i)))
    hmass_c[i].SetTitle('')
    hmass_c[i].SetMarkerColor(massColor)
    hmass_c[i].SetLineColor(signalColor)
    hmass_c[i].SetFillColor(signalColor)
    hmass_c[i].SetFillStyle(3004)
    hmass_c[i].SetLineWidth(1)

    fullfit.append(inFile.Get('fullfit_'+str(i)))
    fullfit[i].SetLineWidth(2)
    fullfit[i].SetLineColor(4)

    bkgfit.append(inFile.Get('bkgFitWRef_'+str(i)))
    bkgfit[i].SetLineWidth(2)
    bkgfit[i].SetLineStyle(1)
    bkgfit[i].SetLineColor(2)


hmass[bin1].SetMaximum(hmass[bin1].GetMaximum()*1.7)
hmass[bin2].SetMaximum(hmass[bin2].GetMaximum()*1.4)
hmass[bin3].SetMaximum(hmass[bin3].GetMaximum()*1.45)

hmass[bin1].GetYaxis().SetTitleOffset(1.4)
hmass[bin2].GetYaxis().SetTitleOffset(1.4)
hmass[bin3].GetYaxis().SetTitleOffset(1.4)

hmass[bin1].GetXaxis().SetTitleOffset(1.1)
hmass[bin2].GetXaxis().SetTitleOffset(1.1)
hmass[bin3].GetXaxis().SetTitleOffset(1.1)

legBands = ROOT.TLegend(0.18,0.86,0.8,0.95)
legBands = ROOT.TLegend(0.181,0.86,0.8,0.95)
legBands.SetTextSize(0.04)
legBands.AddEntry(hmass_c[0],"Signal region","f")
legBands.AddEntry(hmass_l[0],"Side bands","f")

lines = ROOT.TLegend(0.65,0.86,0.90,0.95)
lines.SetTextSize(0.04)
lines.AddEntry(fullfit[bin3],"Signal + bkg","l")
lines.AddEntry(bkgfit[bin3],"Background","l")

pvALICE = ROOT.TPaveText(0.18,0.890,0.6,0.93,"brNDC")
pvALICE.SetFillStyle(0)
pvALICE.SetBorderSize(0)
pvALICE.SetTextFont(42)
pvALICE.SetTextSize(0.048)
pvALICE.SetTextAlign(11)
pvALICE.AddText("ALICE")

pvEn = ROOT.TPaveText(0.168,0.88,0.8,0.92,"brNDC")
pvEn.SetFillStyle(0)
pvEn.SetBorderSize(0)
pvEn.SetTextFont(42)
pvEn.SetTextSize(0.048)
pvEn.SetTextAlign(11)
pvEn.AddText("pp, #sqrt{#it{s}} = "+energy+" TeV")

pvD = ROOT.TPaveText(0.18,0.84,0.55,0.875,"brNDC")
pvD.SetFillStyle(0)
pvD.SetBorderSize(0)
pvD.SetTextFont(42)
pvD.SetTextSize(0.047)
pvD.SetTextAlign(11)
pvD.AddText("D^{0} #rightarrow K^{-}#pi^{+} and charge conj.")

pvJet = ROOT.TPaveText(0.181,0.79,0.55,0.83,"brNDC")
pvJet.SetFillStyle(0)
pvJet.SetBorderSize(0)
pvJet.SetTextFont(42)
pvJet.SetTextSize(0.046)
pvJet.SetTextAlign(11)
pvJet.AddText("in charged jets, anti-#it{k}_{T}, #it{R} = 0."+str(R))

pvEta = ROOT.TPaveText(0.1881,0.725,0.4,0.78,"brNDC")
pvEta.SetFillStyle(0)
pvEta.SetBorderSize(0)
pvEta.SetTextFont(42)
pvEta.SetTextSize(0.046)
pvEta.SetTextAlign(11)
pvEta.AddText("|#it{#eta}_{lab}^{jet}| < 0."+str(9-R))

pvpt1 = ROOT.TPaveText(0.62,0.74,0.9,0.775,"brNDC");
pvpt1.SetFillStyle(0);
pvpt1.SetBorderSize(0);
pvpt1.SetTextFont(42);
pvpt1.SetTextSize(0.046);
pvpt1.AddText("%.0f < #it{p}_{T,D^{0}} < %.0f GeV/#it{c}"%(ptDbins[bin1],ptDbins[bin1+1]));

pvpt2 = ROOT.TPaveText(0.62,0.86,0.85,0.86,"brNDC");
pvpt2.SetFillStyle(0);
pvpt2.SetBorderSize(0);
pvpt2.SetTextFont(42);
pvpt2.SetTextSize(0.046);
pvpt2.AddText("%.0f < #it{p}_{T,D^{0}} < %.0f GeV/#it{c}"%(ptDbins[bin2],ptDbins[bin2+1]));

pvpt3 = ROOT.TPaveText(0.62,0.82,0.85,0.81,"brNDC");
pvpt3.SetFillStyle(0);
pvpt3.SetBorderSize(0);
pvpt3.SetTextFont(42);
pvpt3.SetTextSize(0.046);
pvpt3.AddText("%.0f < #it{p}_{T,D^{0}} < %.0f GeV/#it{c}"%(ptDbins[bin3],ptDbins[bin3+1]));

pvsig1 = ROOT.TPaveText(0.57,0.60,0.9,0.65,"brNDC");
pvsig1.SetFillStyle(0);
pvsig1.SetBorderSize(0);
pvsig1.SetTextFont(42);
pvsig1.SetTextSize(0.04);
pvsig1.AddText("S (2#sigma) = %.1f #pm %.1f"%(hsign.GetBinContent(hsign.FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2.)),hsign.GetBinError(hsign.FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2.))))

pvjetpt1 = ROOT.TPaveText(0.60,0.59,0.89,0.63,"brNDC");
pvjetpt1.SetFillStyle(0);
pvjetpt1.SetBorderSize(0);
pvjetpt1.SetTextFont(42);
pvjetpt1.SetTextSize(0.046);
pvjetpt1.AddText("%.0f < #it{p}_{T, ch.jet} < %.0f GeV/#it{c}"%(ptjetbins[jetbin-1],ptjetbins[jetbin]));

pvsb1 = ROOT.TPaveText(0.17,0.45,0.5,0.49,"brNDC");
pvsb1.SetFillStyle(0);
pvsb1.SetBorderSize(0);
pvsb1.SetTextFont(42);
pvsb1.SetTextSize(0.046);
pvsb1.AddText("S/B (2#sigma) = %.2f"%(hsb.GetBinContent(hsb.FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2. ))));

pvmean1 = ROOT.TPaveText(0.178,0.67,0.62,0.72,"brNDC");
pvmean1.SetFillStyle(0);
pvmean1.SetBorderSize(0);
pvmean1.SetTextFont(42);
pvmean1.SetTextSize(0.046);
pvmean1.SetTextAlign(11);
pvmean1.AddText("#mu = %.2f #pm %.2f GeV/#it{c}^{2}"%(hmean.GetBinContent(hmean.FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2. )),hmean.GetBinError(hmean.FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2. ))));

pvsigma1 = ROOT.TPaveText(0.179,0.62,0.62,0.68,"brNDC");
pvsigma1.SetFillStyle(0);
pvsigma1.SetBorderSize(0);
pvsigma1.SetTextFont(42);
pvsigma1.SetTextSize(0.046);
pvsigma1.SetTextAlign(11);
pvsigma1.AddText("#sigma = %.1f #pm %.1f MeV/#it{c}^{2}"%(hsigma.GetBinContent(hsigma.FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2. )),hsigma.GetBinError(hsigma.FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2. ))));

pvsig2 = ROOT.TPaveText(0.57,0.73,0.9,0.78,"brNDC");
pvsig2.SetFillStyle(0);
pvsig2.SetBorderSize(0);
pvsig2.SetTextFont(42);
pvsig2.SetTextSize(0.04);
pvsig2.AddText("S (2#sigma) = %.1f #pm %.1f"%(hsign.GetBinContent(hsign.FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2.)),hsign.GetBinError(hsign.FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2.))));

pvsb2 = ROOT.TPaveText(0.170,0.45,0.5,0.49,"brNDC");
pvsb2.SetFillStyle(0);
pvsb2.SetBorderSize(0);
pvsb2.SetTextFont(42);
pvsb2.SetTextSize(0.046);
pvsb2.AddText("S/B (2#sigma) = %.2f"%(hsb.GetBinContent(hsb.FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2. ))));

pvmean2 = ROOT.TPaveText(0.178,0.785,0.62,0.81,"brNDC");
pvmean2.SetFillStyle(0);
pvmean2.SetBorderSize(0);
pvmean2.SetTextFont(42);
pvmean2.SetTextSize(0.046);
pvmean2.SetTextAlign(11);

pvmean2.AddText("#mu = %.2f #pm %.2f GeV/#it{c}^{2}"%(hmean.GetBinContent(hmean.FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2. )),hmean.GetBinError(hmean.FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2. ))));

pvsigma2 = ROOT.TPaveText(0.1797,0.735,0.62,0.77,"brNDC");
pvsigma2.SetFillStyle(0);
pvsigma2.SetBorderSize(0);
pvsigma2.SetTextFont(42);
pvsigma2.SetTextSize(0.046);
pvsigma2.SetTextAlign(11);
pvsigma2.AddText("#sigma = %.1f #pm %.1f MeV/#it{c}^{2}"%(hsigma.GetBinContent(hsigma.FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2. )),hsigma.GetBinError(hsigma.FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2. ))));

pvsig3 = ROOT.TPaveText(0.17,0.82,0.9,0.87,"brNDC");
pvsig3.SetFillStyle(0);
pvsig3.SetBorderSize(0);
pvsig3.SetTextFont(42);
pvsig3.SetTextSize(0.04);
pvsig3.AddText("S (2#sigma) = %.1f #pm %.1f"%(hsign.GetBinContent(hsign.FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2.)),hsign.GetBinError(hsign.FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2.))));

pvsb3 = ROOT.TPaveText(0.170,0.45,0.5,0.49,"brNDC");
pvsb3.SetFillStyle(0);
pvsb3.SetBorderSize(0);
pvsb3.SetTextFont(42);
pvsb3.SetTextSize(0.046);
pvsb3.AddText("S/B (2#sigma) = %.2f"%(hsb.GetBinContent(hsb.FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2. ))));

pvmean3 = ROOT.TPaveText(0.17,0.73,0.62,0.76,"brNDC");
pvmean3.SetFillStyle(0);
pvmean3.SetBorderSize(0);
pvmean3.SetTextFont(42);
pvmean3.SetTextSize(0.046);
pvmean3.SetTextAlign(11);
pvmean3.AddText("#mu = %.2f #pm %.2f GeV/#it{c}^{2}"%(hmean.GetBinContent(hmean.FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2. )),hmean.GetBinError(hmean.FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2. ))));

pvsigma3 = ROOT.TPaveText(0.17,0.68,0.62,0.72,"brNDC");
pvsigma3.SetFillStyle(0);
pvsigma3.SetBorderSize(0);
pvsigma3.SetTextFont(42);
pvsigma3.SetTextSize(0.046);
pvsigma3.SetTextAlign(11);
pvsigma3.AddText("#sigma = %.1f #pm %.1f MeV/#it{c}^{2}"%(hsigma.GetBinContent(hsigma.FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2. )),hsigma.GetBinError(hsigma.FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2. ))));

###   CANVAS
cMass = ROOT.TCanvas('cMass','cMass',3000,900)
cMass.Divide(3,1);
cMass.cd(1)
hmass[bin1].Draw();
hmass_l[bin1].Draw("hsame");
hmass_u[bin1].Draw("hsame");
hmass_c[bin1].Draw("hsame");
hmass[bin1].Draw("same");
bkgfit[bin1].Draw("same");
fullfit[bin1].Draw("same");

pvpt1.Draw("same");
pvALICE.Draw("same");
pvD.Draw("same");
pvJet.Draw("same");
pvmean1.Draw("same");
pvsigma1.Draw("same");
pvsb1.Draw("same");
pvEta.Draw("same");
pvjetpt1.Draw('same')


cMass.cd(2)
hmass[bin2].Draw();
hmass_l[bin2].Draw("hsame");
hmass_u[bin2].Draw("hsame");
hmass_c[bin2].Draw("hsame");
hmass[bin2].Draw("same");
bkgfit[bin2].Draw("same");
fullfit[bin2].Draw("same");

pvpt2.Draw("same");
pvEn.Draw("same");
pvmean2.Draw("same");
pvsigma2.Draw("same");
pvsb2.Draw("same");

cMass.cd(3)
hmass[bin3].Draw();
hmass_l[bin3].Draw("hsame");
hmass_u[bin3].Draw("hsame");
hmass_c[bin3].Draw("hsame");
hmass[bin3].Draw("same");
bkgfit[bin3].Draw("same");
fullfit[bin3].Draw("same");

pvpt3.Draw("same");
pvmean3.Draw("same");
pvsigma3.Draw("same");
pvsb3.Draw("same");
legBands.Draw("same");
lines.Draw("same");

cMass.Print('DjetZinMass'+str(R)+'_'+str(jetbin)+'_'+str(int(float(energy)))+'tev.pdf')
cMass.Print('DjetZinMass'+str(R)+'_'+str(jetbin)+'_'+str(int(float(energy)))+'tev.png')
cMass.Print('DjetZinMass'+str(R)+'_'+str(jetbin)+'_'+str(int(float(energy)))+'tev.eps')
