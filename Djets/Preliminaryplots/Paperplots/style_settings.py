import ROOT as RT
import ROOT


myStyle = RT.TStyle('myStyle','myStyle')
myStyle.SetPalette(57,RT.nullptr,1.0)
myStyle.SetOptStat(0);
myStyle.SetOptFit(0);
myStyle.SetOptTitle(0);
myStyle.SetOptDate(0);

myStyle.SetCanvasColor(0);# canvas...
myStyle.SetCanvasBorderMode(0);
myStyle.SetCanvasBorderSize(0);

myStyle.SetPadBottomMargin(0.14); #margins...
myStyle.SetPadTopMargin(0.04);
myStyle.SetPadLeftMargin(0.08);
myStyle.SetPadRightMargin(0.04);

myStyle.SetPadColor(0);
myStyle.SetPadGridX(0); # grids, tickmarks
myStyle.SetPadGridY(0);
myStyle.SetPadTickX(1);
myStyle.SetPadTickY(1);
myStyle.SetPadBorderSize(1);
myStyle.SetPadBorderMode(0);

myStyle.SetFrameBorderSize(1);
myStyle.SetFrameBorderMode(0);
myStyle.SetFrameLineColor(0);

myStyle.SetLegendBorderSize(0);

myStyle.SetTitleOffset(0.8,"x");
myStyle.SetTitleOffset(1.2,"y");
myStyle.SetTitleSize(0.1,"xyz");
myStyle.SetLabelSize(0.1,"xyz");


RT.gROOT.SetStyle("myStyle");
print("Styles are Set!\n")

