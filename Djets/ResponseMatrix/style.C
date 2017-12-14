// rootlogon.C

void style()
{// Add my own options here:
  TStyle* myStyle = new TStyle("myStyle","myStyle");
  myStyle->SetPalette(55,0); // avoid horrible default color scheme
  myStyle->SetOptStat(0);
  myStyle->SetOptFit(0);
  //myStyle->SetOptTitle(0);
  myStyle->SetOptDate(0);
  myStyle->SetStatColor(10);
  myStyle->SetStatFontSize(0.05);
  myStyle->SetTextSize(0.05);
  
  myStyle->SetTitleAlign(13);
  myStyle->SetTitleFont(62,"xyz"); // font option
  myStyle->SetLabelFont(62,"xyz");
  myStyle->SetLabelSize(0.04,"xyz"); // size of axis value font
  myStyle->SetTitleSize(0.04,"xyz"); // size of axis title font
  myStyle->SetTitleOffset(1.1,"x");
  myStyle->SetTitleOffset(1.2,"y");
  myStyle->SetTitleOffset(1.0,"z");
  
  myStyle->SetPadBottomMargin(0.13); //margins...
  myStyle->SetPadTopMargin(0.08);
  myStyle->SetPadLeftMargin(0.12);
  myStyle->SetPadRightMargin(0.1);
  
  myStyle->SetTitleFillColor(0);
  myStyle->SetLegendBorderSize(0);
  //myStyle->SetNdivisions(10, "x");
  //myStyle->SetNdivisions(10, "y");
  myStyle->SetLineWidth(2);
  // default canvas options
  //myStyle->SetCanvasDefW(800);
  //myStyle->SetCanvasDefH(600);
  //myStyle->SetCanvasColor(10);
  myStyle->SetCanvasColor(0);// canvas...
  myStyle->SetCanvasBorderMode(0);
  //myStyle->SetCanvasBorderMode(-1);
  myStyle->SetCanvasBorderSize(0);
  //myStyle->SetCanvasBorderSize(1);
  myStyle->SetPadColor(0);
  myStyle->SetPadBorderSize(1);
  myStyle->SetPadBorderMode(-1);
  myStyle->SetPadGridX(0); // grids, tickmarks
  myStyle->SetPadGridY(0);
  myStyle->SetPadTickX(1);
  myStyle->SetPadTickY(1);
  myStyle->SetFrameBorderSize(1);
  myStyle->SetFrameBorderMode(-1);
  //myStyle->SetFrameBorderMode(0);
  myStyle->SetFrameFillColor(0);
  myStyle->SetFrameLineColor(0); 
  myStyle->SetFrameLineWidth(1.1);


  //myStyle->SetPaperSize(20,24); // US letter size
  gROOT->SetStyle("myStyle");
  cout << "Styles are Set!" << endl;
  return; 
}

