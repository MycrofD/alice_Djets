import ROOT
ROOT.TH1.AddDirectory(False)
from style import *
from ROOT import TCanvas, gStyle,TPad,TH1D
import array
import numpy as np


def PrepareCanvas(xAxisBins,xAxis,zBin,pdf,plotRanges):
    xAxis = array.array('d',xAxis)
    style()
    #print(*args)
    #prepare main canvas
    FinalSpectrum = TCanvas("FinalSpectrum","FinalSpectrum",0,45,700,700)
    gStyle.SetOptStat(0)

    gStyle.SetOptTitle(0);
    FinalSpectrum.SetHighLightColor(2);
    FinalSpectrum.Range(0,0,1,1);
    FinalSpectrum.SetFillColor(0);
    FinalSpectrum.SetBorderMode(0);
    FinalSpectrum.SetBorderSize(2);
    FinalSpectrum.SetFrameBorderMode(0);

    #Set primitives in upper pad
    pad_top = TPad("pad_top","pad_top",0,0.35,1,1);
    pad_top.Draw();
    pad_top.cd();
    pad_top.Range(-1.986821e-07,-4.69897,33.33333,0.3499945);
    pad_top.SetFillColor(0);
    pad_top.SetBorderMode(0);
    pad_top.SetBorderSize(2);
    if(zBin==0 or pdf==False):pad_top.SetLogy();
    pad_top.SetTickx(1);
    pad_top.SetTicky(1);
    pad_top.SetLeftMargin(0.18);
    pad_top.SetBottomMargin(0);
    pad_top.SetFrameBorderMode(0);
    pad_top.SetFrameBorderMode(0);

    FinalSpectrum.cd();

    #Set primitives in bottom pad
    pad_bottom = TPad("pad_bottom", "pad_bottom",0,0,1,0.35);
    pad_bottom.Draw();
    pad_bottom.cd();
    pad_bottom.Range(-1.986821e-07,-0.9209589,33.33333,2.49);
    pad_bottom.SetFillColor(0);
    pad_bottom.SetBorderMode(0);
    pad_bottom.SetBorderSize(2);
    pad_bottom.SetGridy();
    pad_bottom.SetTickx(1);
    pad_bottom.SetTicky(1);
    pad_bottom.SetLeftMargin(0.18);
    pad_bottom.SetTopMargin(0);
    pad_bottom.SetBottomMargin(0.27);
    pad_bottom.SetFrameBorderMode(0);
    pad_bottom.SetFrameBorderMode(0);

    pad_top.cd();
    #hEmpty_up = TH1D("hEmpty_up","Central Values",xAxisBins, array.array('d',xAxis));
    hEmpty_up = TH1D("hEmpty_up","Central Values",xAxisBins, xAxis);
    hEmpty_up.SetMinimum(plotRanges[2]);
    hEmpty_up.SetMaximum(plotRanges[3]);
    hEmpty_up.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
    hEmpty_up.GetXaxis().SetLabelFont(43);
    hEmpty_up.GetXaxis().SetLabelSize(0.035);
    hEmpty_up.GetXaxis().SetTitleSize(0.035);
    hEmpty_up.GetXaxis().SetTitleFont(42);
    if(zBin==0):hEmpty_up.GetYaxis().SetTitle("#frac{d^{2}#it{#sigma}}{d#it{p}_{T}d#it{#eta}} [mb (GeV/#it{c})^{#minus1}]");
    elif(pdf==False):hEmpty_up.GetYaxis().SetTitle("#frac{d^{2}#it{#sigma}}{d#it{z}_{#parallel}d#it{#eta}} (mb)");
    else:hEmpty_up.GetYaxis().SetTitle("Probability density");
    hEmpty_up.GetYaxis().SetLabelFont(43);
    hEmpty_up.GetYaxis().SetLabelSize(22);
    hEmpty_up.GetYaxis().SetTitleSize(26);
    hEmpty_up.GetYaxis().SetLabelOffset(0.015);
    hEmpty_up.GetYaxis().SetTitleOffset(2.);
    hEmpty_up.GetYaxis().SetTitleFont(43);
    hEmpty_up.GetYaxis().SetDecimals();
    hEmpty_up.Draw("axis");

    pad_bottom.cd();
    hEmpty_bottom = TH1D("hEmpty_bottom","Central Values",xAxisBins,xAxis);
    hEmpty_bottom.SetMinimum(plotRanges[0]);
    hEmpty_bottom.SetMaximum(plotRanges[1]);
    if(zBin ==0):hEmpty_bottom.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
    else:hEmpty_bottom.GetXaxis().SetTitle("z_{#parallel}");
    hEmpty_bottom.GetXaxis().SetLabelFont(43);
    hEmpty_bottom.GetXaxis().SetLabelSize(22);
    hEmpty_bottom.GetXaxis().SetTitleSize(26);
    hEmpty_bottom.GetXaxis().SetTitleOffset(3.);
    hEmpty_bottom.GetXaxis().SetTitleFont(43);
    hEmpty_bottom.GetYaxis().SetTitle("theory / data");
    hEmpty_bottom.GetYaxis().SetNdivisions(509);
    hEmpty_bottom.GetYaxis().CenterTitle();
    hEmpty_bottom.GetYaxis().SetDecimals();
    hEmpty_bottom.GetYaxis().SetLabelOffset(0.015);
    hEmpty_bottom.GetXaxis().SetLabelOffset(0.02);
    hEmpty_bottom.GetYaxis().SetLabelFont(43);
    hEmpty_bottom.GetYaxis().SetLabelSize(22);
    hEmpty_bottom.GetYaxis().SetTitleSize(26);
    hEmpty_bottom.GetYaxis().SetTitleOffset(2.);
    hEmpty_bottom.GetYaxis().SetTitleFont(43);
    hEmpty_bottom.Draw("axis");
    """
    """
    return (FinalSpectrum,pad_top,pad_bottom,hEmpty_up,hEmpty_bottom)

def GetData(dataFile,histBase,dataScaling,nBins,xBins,systUncD_down,systUncD_up,pdf,dy):
    #def GetData():
    print("get DATA--------")
    #jetPtFile = ROOT.TFile(dataFile,"read");htmp = jetPtFile.Get(histBase)
    htmp = GetInputHist(dataFile, histBase)
    print("getting ",histBase, "from ",dataFile)
    if(htmp==None):print("no histo found")
    hData_binned = htmp.Rebin(nBins,"hData_binned",xBins)
    #data_int = hData_binned.Integral();
    if(pdf==True):
        hData_binned.Scale(1./hData_binned.Integral());
        hData_binned.Scale(1,"width");
    else:
        hData_binned.Scale(1,"width");
        hData_binned.Scale(dataScaling);
        hData_binned.Scale(1./dy);

    hData_binned.SetTitle("");
    hData_binned.SetMaximum(hData_binned.GetMaximum()*2);
    #hData_binned.GetYaxis().SetTitle("d^{2}#sigma/dp_{T}d#it{#eta} (mb)");
    if(pdf==False):hData_binned.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#it{#eta}} mb (GeV/#it{c})^{-1}")
    else:hData_binned.GetYaxis().SetTitle("Probability density");

    sysuncAbs_down = np.zeros(nBins)
    sysuncAbs_up = np.zeros(nBins)
    statunc = np.zeros(nBins)
    value = np.zeros(nBins)
    xval = np.zeros(nBins)
    xvalwidth = np.zeros(nBins)
    for j in range(nBins):
        xval[j] = (xBins[j]+xBins[j+1]) / 2.
        xvalwidth[j] = (xBins[j+1]-xBins[j]) / 2.
        value[j] = hData_binned.GetBinContent(hData_binned.GetXaxis().FindBin(xval[j]))
        error = hData_binned.GetBinError(hData_binned.GetXaxis().FindBin(xval[j]))
        #print(j," ",xval[j]," ",xvalwidth[j]," ",value[j]," ",error)
        sysuncAbs_down[j] = value[j] * systUncD_down[j];
        sysuncAbs_up[j] = value[j] * systUncD_up[j];
        if(value[j]>1e-9):statunc[j] = (error/ value[j])*100;
    
    grsysl = ROOT.TGraphAsymmErrors(nBins,xval,value,xvalwidth,xvalwidth,sysuncAbs_down,sysuncAbs_up)
    grsysl.SetName("haeData_binned_syst")
    hData_binned_ratio = hData_binned.Clone("hData_binned_ratio")
    hData_binned_ratio.SetName("hData_binned_ratio")
    sysuncRatio_up = np.zeros(nBins)
    sysuncRatio_down = np.zeros(nBins)
    valRatio = np.zeros(nBins)
    for j in range(nBins):
        xval[j] = (xBins[j]+xBins[j+1]) / 2.;
        xvalwidth[j] = (xBins[j+1]-xBins[j]) / 2.;
        valPred = hData_binned.GetBinContent(hData_binned.GetXaxis().FindBin(xval[j]))
        if(valPred<1e-9):valPred=1000
        if(valPred>1e-9):valRatio[j] = 1.0;
        err = 0;
        if(valPred>1e-9):err=hData_binned.GetBinError(hData_binned.GetXaxis().FindBin(xval[j])) / valPred;
        sysuncRatio_down[j] = valRatio[j] * systUncD_down[j];
        sysuncRatio_up[j] = valRatio[j] * systUncD_up[j];
        hData_binned_ratio.SetBinContent(hData_binned_ratio.GetXaxis().FindBin(xval[j]),valRatio[j]);
        hData_binned_ratio.SetBinError(hData_binned_ratio.GetXaxis().FindBin(xval[j]),err);
    grsysRatio = ROOT.TGraphAsymmErrors(nBins,xval,valRatio,xvalwidth,xvalwidth,sysuncRatio_down,sysuncRatio_up)
    grsysRatio.SetName("haeData_binned_syst_ratio");
    """
    """
    return (grsysl,hData_binned, grsysRatio,hData_binned_ratio) #return (hDataSys, hData_binned, hDataSysRatio, hData_binned_ratio)
#----------
def GetSim(simname, type_, nFiles, names, simDir, simScaling, nBins, xBins, DpTbins, jetpTbins, zBin, pdf,dy):
    hPrompt, hPrompt_binned = [],[] 
    tmp = "h"
    tmp+=simname
    for nr in range(nFiles):
        file_ = simDir
        file_ += "/JetPt_"
        file_ += names[nr]
        if(type_ == 0):
            file_ += "_Dpt"+str(int(DpTbins[0]))+"_"+str(int(DpTbins[1]));
        if(type_ ==1 or type_ ==2):
            file_ += "_Dpt"; file_ += DpTbins[0][zBin-1]; file_ += "_"; file_ += DpTbins[1][zBin-1];
            file_ += "_Jetpt"; file_ += jetpTbins[zBin-1]; file_ += "_"; file_ += jetpTbins[zBin];
        file_ += "_Dzero.root"
        htmp = ROOT.TH1D()
        if(type_ ==0):
            htmp = GetInputHist(file_,"hPt")
            if(htmp==False): print("no histogram in file: ",file_)
            htmp.GetYaxis().SetTitle("d#sigma/dp_{T} (mb)");
        if(type_ ==1 or type_ ==2):
            htmp = GetInputHist(file_,"hz")
            if(htmp==False): print("no histogram in file: ",file_)
            htmp.GetYaxis().SetTitle("d#sigma/dz (mb)");
        hPrompt.append(htmp.Clone("hPrompt_%d"%(nr)))
        hPrompt_binned.append(htmp.Rebin(nBins,"hPrompt_binned_%d"%(nr),xBins))
        hPrompt_binned[nr].SetName(tmp+"_var"+("%d"%(nr)))
        if(pdf==True):
            hPrompt_binned[nr].Scale(1./hPrompt_binned[nr].Integral())
            hPrompt_binned[nr].Scale(1,"width")
    hPrompt_central_binned = hPrompt_binned[0].Clone(tmp+"_central")
    setHistoDetails(hPrompt_central_binned,4,24,1,2,2)
    # get up unc
    hPrompt_up = GetUpSys(hPrompt_binned,nFiles)
    hPrompt_up.SetName(tmp+"_up")
    setHistoDetails(hPrompt_up,4,24,0,2,2)
    # get down unc
    hPrompt_down = GetDownSys(hPrompt_binned,nFiles)
    hPrompt_down.SetName(tmp+"_down")
    setHistoDetails(hPrompt_down,4,24,0,2,2)

    if(pdf==False):
        hPrompt_central_binned.Scale(simScaling)
        hPrompt_central_binned.Scale(1,"width")
        hPrompt_central_binned.Scale(1./dy) #2*jetEta;
        hPrompt_up.Scale(simScaling)
        hPrompt_up.Scale(1,"width")
        hPrompt_up.Scale(1./dy)
        hPrompt_down.Scale(simScaling)
        hPrompt_down.Scale(1,"width")
        hPrompt_down.Scale(1./dy)
        for nr in range(nFiles):
            hPrompt_binned[nr].Scale(simScaling)
            hPrompt_binned[nr].Scale(1,"width")
            hPrompt_binned[nr].Scale(1./dy)
    xval = np.zeros(nBins)
    xvalwidth = np.zeros(nBins)
    valuetheory = np.zeros(nBins)
    valuetheoryerrup = np.zeros(nBins)
    valuetheoryerrdown = np.zeros(nBins)
    for j in range(nBins):
        xval[j] = (xBins[j]+xBins[j+1]) / 2.
        xvalwidth[j] = (xBins[j+1]-xBins[j]) / 2.
        valuetheory[j] = hPrompt_central_binned.GetBinContent(hPrompt_central_binned.GetXaxis().FindBin(xval[j]))
        valuetheoryerrup[j] = hPrompt_up.GetBinContent(hPrompt_up.GetXaxis().FindBin(xval[j])) - valuetheory[j]
        valuetheoryerrdown[j] = valuetheory[j] - hPrompt_down.GetBinContent(hPrompt_up.GetXaxis().FindBin(xval[j]))

    grsystheory = ROOT.TGraphAsymmErrors(nBins,xval,valuetheory,xvalwidth,xvalwidth,valuetheoryerrdown,valuetheoryerrup)
    tmp = "hae";
    tmp+=simname;
    grsystheory.SetName(tmp)
    return (grsystheory, hPrompt_central_binned, hPrompt_up,hPrompt_down,hPrompt_binned)
#----------
def GetDataSimRatio(simname, data_cent, sim_cent, sim_up, sim_down, nBins, xBins):
    name = sim_cent.GetName()
    name_up = sim_up.GetName()
    name_down = sim_down.GetName()
    hPrompt_central_binned_ratio = sim_cent.Clone(name+"_ratio")
    hPrompt_central_binned_ratio.Divide(data_cent)
    hPrompt_down_ratio = sim_down.Clone(name_down+"_ratio")
    hPrompt_up_ratio = sim_up.Clone(name_up+"_ratio")
    hPrompt_up_ratio.Divide(data_cent)
    hPrompt_down_ratio.Divide(data_cent)
    ptvaltheoryratio = np.zeros(nBins)
    ptvalunctheoryratio = np.zeros(nBins)
    valuetheoryratio = np.zeros(nBins)
    valuetheoryerrupratio = np.zeros(nBins)
    valuetheoryerrdownratio = np.zeros(nBins)
    for j in range(nBins):
        ptvaltheoryratio[j] = (xBins[j]+xBins[j+1]) / 2.
        ptvalunctheoryratio[j] = (xBins[j+1]-xBins[j]) / 2.
        valuetheoryratio[j] = hPrompt_central_binned_ratio.GetBinContent(hPrompt_central_binned_ratio.GetXaxis().FindBin(ptvaltheoryratio[j]))
        valuetheoryerrupratio[j] = hPrompt_up_ratio.GetBinContent(hPrompt_up_ratio.GetXaxis().FindBin(ptvaltheoryratio[j])) - valuetheoryratio[j]
        valuetheoryerrdownratio[j] = valuetheoryratio[j] - hPrompt_down_ratio.GetBinContent(hPrompt_down_ratio.GetXaxis().FindBin(ptvaltheoryratio[j]))
    grsystheoryratio = ROOT.TGraphAsymmErrors(nBins,ptvaltheoryratio,valuetheoryratio,ptvalunctheoryratio,ptvalunctheoryratio,valuetheoryerrdownratio,valuetheoryerrupratio)
    tmp = "hae";
    tmp+=simname;
    grsystheoryratio.SetName(tmp+"_ratio")
    return (grsystheoryratio, hPrompt_central_binned_ratio, hPrompt_up_ratio,hPrompt_down_ratio)
#----------
def GetUpSys(hh,nFiles):
    max_ , maxerr = 0, 0
    name = hh[0].GetName()
    hh_up = hh[0].Clone(name +"_up")
    for iBin in range(1,hh[0].GetNbinsX()+1):
        max_ = hh[0].GetBinContent(iBin)
        for iFile in range(1, nFiles):
            if(hh[iFile].GetBinContent(iBin) > max_):
                max_ = hh[iFile].GetBinContent(iBin)
                maxerr = hh[iFile].GetBinError(iBin)
                
        hh_up.SetBinContent(iBin,max_)
        hh_up.SetBinError(iBin,0)

    return hh_up
#----------
def GetDownSys(hh,nFiles):
    max_ , maxerr = 0, 0
    name = hh[0].GetName()
    hh_down = hh[0].Clone(name +"_down")
    for iBin in range(1,hh[0].GetNbinsX()+1):
        max_ = hh[0].GetBinContent(iBin)
        for iFile in range(1, nFiles):
            if(hh[iFile].GetBinContent(iBin) < max_):
                max_ = hh[iFile].GetBinContent(iBin)
                maxerr = hh[iFile].GetBinError(iBin)
                
        hh_down.SetBinContent(iBin,max_)
        hh_down.SetBinError(iBin,0)

    return hh_down
#----------
def PlaceOnPadSim(pad, histo, ci, style, linestyle):
    pad.cd()
    histo.SetFillColor(1);
    histo.SetFillStyle(0);
    histo.SetLineColor(ci);
    histo.SetMarkerColor(ci);
    histo.SetMarkerStyle(style);
    histo.SetLineStyle(linestyle);
    #histo.SetFillStyle(3005);
    #histo.SetMarkerSize(markersize); #add up
    histo.SetLineWidth(2);
    odraaw = "";
    if(linestyle ==1):odraaw = "2p"
    else:
        histo.SetLineWidth(3)
        histo.SetMarkerSize(0)
        odraaw = "E"
    histo.Draw(odraaw)
#----------
def PlaceOnPadData( pad, histo1, histo2, markersize):
    pad.cd()
    ci = ROOT.TColor.GetColor("#990000")
    ci = ROOT.kBlack
    histo1.SetLineColor(ci);
    histo1.SetMarkerColor(ci);
    histo1.SetMarkerStyle(20);
    histo1.SetMarkerSize(markersize);
    ci = ROOT.TColor.GetColor("#cccccc");
    histo1.SetFillColor(ci);
    histo1.SetLineColor(ci);
    histo1.Draw("2p");
    #data central w stat. unc.
    #ci = static_cast<Color_t>(TColor::GetColor("#990000"));
    ci = ROOT.kBlack;
    histo2.SetLineColor(ci);
    histo2.SetMarkerColor(ci);
    histo2.SetMarkerStyle(20);
    histo2.SetMarkerSize(markersize);
    histo2.Draw("same p  e0 x0");
#----------
def TerminateCanvas(pad1, pad2, histo1,histo2):
    pad1.cd()
    histo1.Draw("sameaxis")
    pad2.cd()
    histo2.Draw("sameaxis")
    #histo2.Draw("sameaxig")

#----------
def GetInputHist(inFile, histName):
    jetPtFile = ROOT.TFile(inFile,"read")
    hh = jetPtFile.Get(histName)
    return hh

def setHistoDetails(hh, color, Mstyle, Msize, Lwidth, Lstyle):
    hh.SetMarkerColor(color)
    hh.SetMarkerStyle(Mstyle)
    hh.SetLineColor(color)
    hh.SetLineWidth(Lwidth)
    hh.SetMarkerSize(Msize)
    hh.SetLineStyle(Lstyle)
    hh.SetTitle("")
    hh.GetXaxis().SetTitle("p_{T}^{ch,jet} (GeV/c)")

