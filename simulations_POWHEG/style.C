// rootlogon.C

void style()
{// Add my own options here:
  TStyle* myStyle = new TStyle("myStyle","myStyle");
  myStyle->SetPalette(55,0); // avoid horrible default color scheme
  myStyle->SetOptStat(0);
  myStyle->SetOptFit(0);
  //myStyle->SetOptTitle(0);
  myStyle->SetOptDate(0);
 
  //myStyle->SetTextSize(0.05);
  
  //myStyle->SetLegendFont(42);
  //myStyle->SetLegendTextSize(0.05);

   myStyle->SetCanvasColor(0);// canvas...
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetCanvasBorderSize(0);
  
  myStyle->SetPadBottomMargin(0.1); //margins...
  myStyle->SetPadTopMargin(0.04);
  myStyle->SetPadLeftMargin(0.08);
  myStyle->SetPadRightMargin(0.04);
  
  myStyle->SetPadColor(0);
  myStyle->SetPadGridX(0); // grids, tickmarks
  myStyle->SetPadGridY(0);
  myStyle->SetPadTickX(1);
  myStyle->SetPadTickY(1);
   myStyle->SetPadBorderSize(1);
  myStyle->SetPadBorderMode(0);
  
   myStyle->SetFrameBorderSize(1);
  myStyle->SetFrameBorderMode(0);
  myStyle->SetFrameLineColor(0); 
  
   myStyle->SetLegendBorderSize(0);
   
  myStyle->SetTitleOffset(0.8,"x");
  myStyle->SetTitleOffset(1.,"y");
  myStyle->SetTitleSize(0.08,"xyz");
   myStyle->SetLabelSize(0.08,"xyz");
  
  /*
  myStyle->SetTitleAlign(13);
  myStyle->SetTitleFont(62,"xyz"); // font option
  myStyle->SetLabelFont(62,"xyz");
  myStyle->SetLabelSize(0.04,"xyz"); // size of axis value font
  myStyle->SetTitleSize(0.04,"xyz"); // size of axis title font
  myStyle->SetTitleOffset(1.1,"x");
  myStyle->SetTitleOffset(1.2,"y");
  myStyle->SetTitleOffset(1.0,"z");
  
  
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
*/

  //myStyle->SetPaperSize(20,24); // US letter size
  gROOT->SetStyle("myStyle");
  cout << "Styles are Set!" << endl;
  return; 
}

