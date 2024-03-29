import os, os.path, sys
import ROOT 
import style_settings
import array
ROOT.TH1.AddDirectory(False)
############
### SETTINGS
plotmin,plotmax=0,100
fptbinsJN = 10
fptbinsJlh=[2,3,4,5,6,8,10,14,20,30,50]
R=4 #int(sys.argv[1]) #R=2,3,4,6
##########################
## the histograms
inFileJ = ROOT.TFile("unfoldedSpectrum_unfoldedJetSpectrum.root","read")

hResp = inFileJ.Get("fMatrixProd")
hResp.GetXaxis().SetCanExtend(True)
#hResp.GetXaxis().SetCanExtend(True);
hResp.GetYaxis().SetCanExtend(True);

hResp.ExtendAxis(200,hResp.GetXaxis())
hResp.ExtendAxis(200,hResp.GetYaxis())
hResp.GetXaxis().SetRangeUser(plotmin,plotmax);
hResp.GetYaxis().SetRangeUser(plotmin,plotmax);
## HIST SETTINGS
#hEmpty = ROOT.TH2D("hE","hE",100,plotmin,plotmax,100,plotmin,plotmax)
##hEmpty = ROOT.TH2D("hE","hE",fptbinsJN,array.array('d',fptbinsJlh),fptbinsJN,array.array('d',fptbinsJlh))
#hEmpty.SetTitle('')
#hEmpty.GetXaxis().SetRangeUser(plotmin,plotmax);
#hEmpty.GetYaxis().SetRangeUser(plotmin,plotmax);
##########################
### CANVAS
cResp = ROOT.TCanvas("cResp","cResp",950,800);
cResp.SetLogz()
#hEmpty.Draw()
hResp.Draw('colz')

cResp.SaveAs('minim_Resp_R04_.png')

input()
