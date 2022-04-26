import os, os.path, sys
import ROOT 
import style_settings
import array
ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetBatch()
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
textsize=0.032
textsize1=0.04
textsize2=0.042

ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadTopMargin(1)
ROOT.gStyle.SetPadRightMargin(0.15)
##########################
## the histograms
inFileJ = ROOT.TFile("/eos/user/a/amohanty/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(R)+"_paperCuts/Default/unfolding_Bayes_5/unfoldedSpectrum_unfoldedJetSpectrum.root","read")
inFileJresol = ROOT.TFile("/eos/user/a/amohanty/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0"+str(R)+"_paperCuts/Default/ResponseMatrix/DetMatrix_prompt.root","read")

hResp = inFileJ.Get("fMatrixProd")
hResp.GetXaxis().SetTitle("#it{p}_{T,ch jet}^{det} (GeV/#it{c})")
hResp.GetXaxis().SetLabelSize(textsize1);
hResp.GetXaxis().SetTitleSize(textsize2);
hResp.GetXaxis().SetTitleOffset(1.3);
hResp.GetYaxis().SetTitle("#it{p}_{T,ch jet}^{part} (GeV/#it{c})")
hResp.GetYaxis().SetTitleOffset(1.3);
hResp.GetYaxis().SetLabelSize(textsize1);
hResp.GetYaxis().SetTitleSize(textsize2);
## HIST SETTINGS
#hEmpty = ROOT.TH2D("hE","hE",100,plotmin,plotmax,100,plotmin,plotmax)
hEmpty = ROOT.TH2D("hE","hE",fptbinsJN,array.array('d',fptbinsJlh),fptbinsJN,array.array('d',fptbinsJlh))
hEmpty.SetTitle('')
#hEmpty.GetXaxis().SetTitle("#it{z}_{||}^{ch}");
#hEmpty.GetYaxis().SetTitle("B Feed-Down Fraction");
hEmpty.GetXaxis().SetLabelSize(0.04);
#hEmpty.GetXaxis().SetTitleSize(0.041);
hEmpty.GetXaxis().SetTitleOffset(1.3);
hEmpty.GetYaxis().SetTitleOffset(1.3);
hEmpty.GetYaxis().SetLabelSize(0.04);
#hEmpty.GetYaxis().SetTitleSize(0.041);
hEmpty.GetXaxis().SetRangeUser(plotmin,plotmax);
hEmpty.GetYaxis().SetRangeUser(plotmin,plotmax);

##########################
y1,y2=0.3,0.35
x1,x2=0.46,0.83
#pvALICE = ROOT.TPaveText(x1,y1,x2,y2,"brNDC");
pvALICE = ROOT.TPaveText(0.15,0.8,0.52,0.85,"brNDC");
pvALICE.SetFillStyle(0);
pvALICE.SetBorderSize(0);
pvALICE.SetTextFont(42);
pvALICE.SetTextSize(textsize1);
pvALICE.SetTextAlign(11);
pvALICE.AddText("ALICE, pp, #sqrt{#it{s}} = 5.02 TeV");

shift = 0.05;
pvJet = ROOT.TPaveText(x1,y1-shift,x2,y2-shift,"brNDC");
pvJet.SetFillStyle(0);
pvJet.SetBorderSize(0);
pvJet.SetTextFont(42);
pvJet.SetTextSize(textsize);
pvJet.SetTextAlign(11);
pvJet.AddText("charged jets, anti-#scale[0.5]{ }#it{k}_{T}, #it{R} = 0."+str(R));

shift += 0.05;
pvD = ROOT.TPaveText(x1,y1-shift,x2,y2-shift,"brNDC");
pvD.SetFillStyle(0);
pvD.SetBorderSize(0);
pvD.SetTextFont(42);
pvD.SetTextSize(textsize);
pvD.SetTextAlign(11);
pvD.AddText("with D^{0} #rightarrow K^{#font[122]{-}}#pi^{+} and charge conj.");

shift += 0.05;
pvEta = ROOT.TPaveText(x1,y1-shift,x2,y2-shift,"brNDC");
pvEta.SetFillStyle(0);
pvEta.SetBorderSize(0);
pvEta.SetTextFont(42);
pvEta.SetTextSize(textsize);
pvEta.SetTextAlign(11);
pvEta.AddText("#it{p}_{T, D^{0}} > 2 GeV/ #it{c}, |#it{#eta}_{ch jet}| < 0."+str(9-R));

### CANVAS
cResp = ROOT.TCanvas("cResp","cResp",950,800);
cResp.SetLogz()
hResp.Scale(1./hResp.Integral())
hResp.SetMinimum(0.00001)
hResp.SetMaximum(1)
hResp.GetXaxis().SetRangeUser(plotmin, plotmax)
hResp.GetYaxis().SetRangeUser(plotmin, plotmax)
#hResp.SetMaximum(10000)
hResp.Draw('colz')

pvALICE.Draw("same");
pvJet.Draw("same");
pvD.Draw("same");
#pvjetpt.Draw("same");
pvEta.Draw("same");

ROOT.gPad.Update()
palette = hResp.GetListOfFunctions().FindObject("palette")

# the following lines moe the paletter. Choose the values you need for the position.
palette.SetX1NDC(0.88);
palette.SetX2NDC(0.93);
palette.SetY1NDC(0.1);
palette.SetY2NDC(0.9);
ROOT.gPad.Modified();
ROOT.gPad.Update();



cResp.SaveAs('plots/Resp_R0'+str(R)+jetorz+'_prob.pdf')
cResp.SaveAs('plots/Resp_R0'+str(R)+jetorz+'_prob.png')
cResp.SaveAs('plots/Resp_R0'+str(R)+jetorz+'_prob.eps')
#################################################################################
#################################################################################
#################################################################################
