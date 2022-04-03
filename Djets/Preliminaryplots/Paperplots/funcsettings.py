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
    #return (FinalSpectrum,pad_top,pad_bottom,hEmpty_up,hEmpty_bottom)
    return (FinalSpectrum,pad_top,hEmpty_up)

