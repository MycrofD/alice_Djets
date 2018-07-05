

{
//=========Macro generated from canvas: FinalSpectrum/FinalSpectrum
//=========  (Tue Sep  5 14:25:17 2017) by ROOT version5.34/30
   TCanvas *FinalSpectrum = new TCanvas("FinalSpectrum", "FinalSpectrum",0,70,700,700);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   FinalSpectrum->SetHighLightColor(2);
   FinalSpectrum->Range(0,0,1,1);
   FinalSpectrum->SetFillColor(0);
   FinalSpectrum->SetBorderMode(0);
   FinalSpectrum->SetBorderSize(2);
   FinalSpectrum->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: FinalSpectrum_1
   TPad *FinalSpectrum_1 = new TPad("FinalSpectrum_1", "FinalSpectrum_1",0,0.35,1,1);
   FinalSpectrum_1->Draw();
   FinalSpectrum_1->cd();
   FinalSpectrum_1->Range(-1.986821e-07,-2.958607,33.33333,1.96998);
   FinalSpectrum_1->SetFillColor(0);
   FinalSpectrum_1->SetBorderMode(0);
   FinalSpectrum_1->SetBorderSize(2);
   FinalSpectrum_1->SetLogy();
   FinalSpectrum_1->SetTickx(1);
   FinalSpectrum_1->SetTicky(1);
   FinalSpectrum_1->SetLeftMargin(0.15);
   FinalSpectrum_1->SetBottomMargin(0);
   FinalSpectrum_1->SetFrameBorderMode(0);
   FinalSpectrum_1->SetFrameBorderMode(0);
   Double_t xAxis1[7] = {5, 6, 8, 10, 14, 20, 30}; 
   
   TH1D *CentralPointsStatisticalUncertainty__1__1 = new TH1D("CentralPointsStatisticalUncertainty__1__1","Central Values",6, xAxis1);
   CentralPointsStatisticalUncertainty__1__1->SetMinimum(0.0011);
   CentralPointsStatisticalUncertainty__1__1->SetMaximum(30);
   CentralPointsStatisticalUncertainty__1__1->SetEntries(8);
   CentralPointsStatisticalUncertainty__1__1->SetDirectory(0);
   CentralPointsStatisticalUncertainty__1__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__1__1->SetLineColor(ci);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__1__1->SetMarkerColor(ci);
   CentralPointsStatisticalUncertainty__1__1->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__1__1->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__1__1->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#it{#eta}} [mb (GeV/#it{c})^{-1}]");
   CentralPointsStatisticalUncertainty__1__1->GetYaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__1__1->GetYaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__1__1->GetYaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__1__1->GetYaxis()->SetTitleOffset(1.6);
   CentralPointsStatisticalUncertainty__1__1->GetYaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__1__1->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__1__1->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__1__1->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__1__1->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__1__1->Draw("axis");
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6);
   grae->SetName("CentralPointsSystematicUncertainty_copy");
   grae->SetTitle("Bayes, iter=4, prior=ResponseTruth Systematics");
   grae->SetFillColor(17);
   grae->SetLineColor(17);
   grae->SetMarkerColor(17);
   grae->SetPoint(0,5.5,0.9776732);
   grae->SetPointError(0,0.5,0.5,0.146651,0.146651);
   grae->SetPoint(1,7,0.4808065);
   grae->SetPointError(1,1,1,0.06250484,0.06250484);
   grae->SetPoint(2,9,0.2106707);
   grae->SetPointError(2,1,1,0.0294939,0.0294939);
   grae->SetPoint(3,12,0.09992121);
   grae->SetPointError(3,2,2,0.01199054,0.01199054);
   grae->SetPoint(4,17,0.02028462);
   grae->SetPointError(4,3,3,0.002839846,0.002839846);
   grae->SetPoint(5,25,0.004155956);
   grae->SetPointError(5,5,5,0.000664953,0.000664953);
   
   TH1F *Graph_Graph_central_syst_unc11 = new TH1F("Graph_Graph_central_syst_unc11","Bayes, iter=4, prior=ResponseTruth Systematics",100,2.5,32.5);
   Graph_Graph_central_syst_unc11->SetMinimum(4.779682e-05);
   Graph_Graph_central_syst_unc11->SetMaximum(0.02142993);
   Graph_Graph_central_syst_unc11->SetDirectory(0);
   Graph_Graph_central_syst_unc11->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_central_syst_unc11->SetLineColor(ci);
   Graph_Graph_central_syst_unc11->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   Graph_Graph_central_syst_unc11->GetXaxis()->SetLabelFont(42);
   Graph_Graph_central_syst_unc11->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_central_syst_unc11->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_central_syst_unc11->GetXaxis()->SetTitleFont(42);
   Graph_Graph_central_syst_unc11->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#it{#eta}} [mb (GeV/#it{c})^{-1}]");
   Graph_Graph_central_syst_unc11->GetYaxis()->SetLabelFont(42);
   Graph_Graph_central_syst_unc11->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_central_syst_unc11->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_central_syst_unc11->GetYaxis()->SetTitleFont(42);
   Graph_Graph_central_syst_unc11->GetZaxis()->SetLabelFont(42);
   Graph_Graph_central_syst_unc11->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_central_syst_unc11->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_central_syst_unc11->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_central_syst_unc11);
   
   grae->Draw("2");
   Double_t xAxis2[7] = {5, 6, 8, 10, 14, 20, 30}; 
   
   TH1D *CentralPointsStatisticalUncertainty__2__2 = new TH1D("CentralPointsStatisticalUncertainty__2__2","Central Values",6, xAxis2);
   CentralPointsStatisticalUncertainty__2__2->SetBinContent(0,1.834125);
   CentralPointsStatisticalUncertainty__2__2->SetBinContent(1,0.9776732);
   CentralPointsStatisticalUncertainty__2__2->SetBinContent(2,0.4808065);
   CentralPointsStatisticalUncertainty__2__2->SetBinContent(3,0.2106707);
   CentralPointsStatisticalUncertainty__2__2->SetBinContent(4,0.09992121);
   CentralPointsStatisticalUncertainty__2__2->SetBinContent(5,0.02028462);
   CentralPointsStatisticalUncertainty__2__2->SetBinContent(6,0.004155956);
   CentralPointsStatisticalUncertainty__2__2->SetBinContent(7,0.0002553179);
   CentralPointsStatisticalUncertainty__2__2->SetBinError(0,0.06249125);
   CentralPointsStatisticalUncertainty__2__2->SetBinError(1,0.05513667);
   CentralPointsStatisticalUncertainty__2__2->SetBinError(2,0.03286703);
   CentralPointsStatisticalUncertainty__2__2->SetBinError(3,0.02035879);
   CentralPointsStatisticalUncertainty__2__2->SetBinError(4,0.01021899);
   CentralPointsStatisticalUncertainty__2__2->SetBinError(5,0.00376896);
   CentralPointsStatisticalUncertainty__2__2->SetBinError(6,0.001104688);
   CentralPointsStatisticalUncertainty__2__2->SetBinError(7,0.0002005154);
   CentralPointsStatisticalUncertainty__2__2->SetMinimum(2e-05);
   CentralPointsStatisticalUncertainty__2__2->SetMaximum(4);
   CentralPointsStatisticalUncertainty__2__2->SetEntries(8);
   CentralPointsStatisticalUncertainty__2__2->SetDirectory(0);
   CentralPointsStatisticalUncertainty__2__2->SetStats(0);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__2__2->SetLineColor(ci);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__2__2->SetMarkerColor(ci);
   CentralPointsStatisticalUncertainty__2__2->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__2__2->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__2__2->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__2__2->GetXaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__2__2->GetXaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__2__2->GetXaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__2__2->GetXaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__2__2->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#it{#eta}} [mb (GeV/#it{c})^{-1}]");
   CentralPointsStatisticalUncertainty__2__2->GetYaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__2__2->GetYaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__2__2->GetYaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__2__2->GetYaxis()->SetTitleOffset(1.4);
   CentralPointsStatisticalUncertainty__2__2->GetYaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__2__2->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__2__2->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__2__2->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__2__2->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__2__2->Draw("same p e0 x0");
   Double_t xAxis3[7] = {5, 6, 8, 10, 14, 20, 30}; 
   
   TH1D *GeneratorLevel_JetPtSpectrum__3__3 = new TH1D("GeneratorLevel_JetPtSpectrum__3__3","D0_MCTruth_Charged_R040_JetPtDPtSpectrum",6, xAxis3);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(0,1.932359);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(1,0.6118639);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(2,0.3309636);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(3,0.1487388);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(4,0.05567059);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(5,0.01439281);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(6,0.00310096);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(7,0.0009589948);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(0,0.005520777);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(1,0.003106589);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(2,0.001615591);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(3,0.001083062);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(4,0.000468532);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(5,0.0001945152);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(6,6.993658e-05);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(7,3.889237e-05);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(7,3.889237e-05);
  /* GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(0,1.586436);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(1,0.5048265);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(2,0.2791936);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(3,0.1260778);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(4,0.04664339);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(5,0.01253488);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(6,0.002788148);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(7,0.00081422);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(0,0.004003268);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(1,0.002258262);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(2,0.00118752);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(3,0.0007980088);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(4,0.0003432167);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(5,0.000145274);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(6,5.307149e-05);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(7,2.867968e-05);*/
   GeneratorLevel_JetPtSpectrum__3__3->SetEntries(316731);
   GeneratorLevel_JetPtSpectrum__3__3->SetDirectory(0);
   GeneratorLevel_JetPtSpectrum__3__3->SetStats(0);

   ci = TColor::GetColor("#000099");
   GeneratorLevel_JetPtSpectrum__3__3->SetLineColor(ci);

   ci = TColor::GetColor("#000099");
   GeneratorLevel_JetPtSpectrum__3__3->SetMarkerColor(ci);
   GeneratorLevel_JetPtSpectrum__3__3->SetMarkerStyle(24);
   GeneratorLevel_JetPtSpectrum__3__3->SetMarkerSize(1.2);
   GeneratorLevel_JetPtSpectrum__3__3->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   GeneratorLevel_JetPtSpectrum__3__3->GetXaxis()->SetLabelFont(42);
   GeneratorLevel_JetPtSpectrum__3__3->GetXaxis()->SetLabelSize(0.035);
   GeneratorLevel_JetPtSpectrum__3__3->GetXaxis()->SetTitleSize(0.035);
   GeneratorLevel_JetPtSpectrum__3__3->GetXaxis()->SetTitleFont(42);
   GeneratorLevel_JetPtSpectrum__3__3->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}} #times #Delta#it{p}_{T} (mb)");
   GeneratorLevel_JetPtSpectrum__3__3->GetYaxis()->SetLabelFont(42);
   GeneratorLevel_JetPtSpectrum__3__3->GetYaxis()->SetLabelSize(0.035);
   GeneratorLevel_JetPtSpectrum__3__3->GetYaxis()->SetTitleSize(0.035);
   GeneratorLevel_JetPtSpectrum__3__3->GetYaxis()->SetTitleFont(42);
   GeneratorLevel_JetPtSpectrum__3__3->GetZaxis()->SetLabelFont(42);
   GeneratorLevel_JetPtSpectrum__3__3->GetZaxis()->SetLabelSize(0.035);
   GeneratorLevel_JetPtSpectrum__3__3->GetZaxis()->SetTitleSize(0.035);
   GeneratorLevel_JetPtSpectrum__3__3->GetZaxis()->SetTitleFont(42);
   GeneratorLevel_JetPtSpectrum__3__3->Draw("same p e0 x0");
   
   TGraphAsymmErrors *graeT = new TGraphAsymmErrors(6);
   graeT->SetName("theorySyst_copy");
   graeT->SetTitle("Graph");
   graeT->SetFillColor(1);
   graeT->SetFillStyle(0);

   ci = TColor::GetColor("#000099");
   graeT->SetLineColor(ci);
   graeT->SetLineWidth(2);
   
   graeT->SetPoint(0,5.5,0.6118639);
   graeT->SetPointError(0,0.5,0.5,0.4024011,0.7254089);
   graeT->SetPoint(1,7,0.3309636);
   graeT->SetPointError(1,1,1,0.200555,0.3737457);
   graeT->SetPoint(2,9,0.1487388);
   graeT->SetPointError(2,1,1,0.08383411,0.1545596);
   graeT->SetPoint(3,12,0.05567059);
   graeT->SetPointError(3,2,2,0.02818321,0.05231898);
   graeT->SetPoint(4,17,0.01439281);
   graeT->SetPointError(4,3,3,0.006296011,0.01428417);
   graeT->SetPoint(5,25,0.00310096);
   graeT->SetPointError(5,5,5,0.001199361,0.003019976);
   /*
   graeT->SetPoint(0,5.5,0.5048265);
   graeT->SetPointError(0,0.5,0.5,0.3307398,0.6140599);
   graeT->SetPoint(1,7,0.2791936);
   graeT->SetPointError(1,1,1,0.1709185,0.3091237);
   graeT->SetPoint(2,9,0.1260778);
   graeT->SetPointError(2,1,1,0.07066835,0.1294046);
   graeT->SetPoint(3,12,0.04664339);
   graeT->SetPointError(3,2,2,0.02322777,0.04466855);
   graeT->SetPoint(4,17,0.01253488);
   graeT->SetPointError(4,3,3,0.005750655,0.01229934);
   graeT->SetPoint(5,25,0.002788148);
   graeT->SetPointError(5,5,5,0.001190245,0.002366597);*/
   
   TH1F *Graph_Graph_theorySyst_copy22 = new TH1F("Graph_Graph_theorySyst_copy22","Graph",100,2.5,32.5);
   Graph_Graph_theorySyst_copy22->SetMinimum(2.864129e-05);
   Graph_Graph_theorySyst_copy22->SetMaximum(0.01825038);
   Graph_Graph_theorySyst_copy22->SetDirectory(0);
   Graph_Graph_theorySyst_copy22->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_theorySyst_copy22->SetLineColor(ci);
   Graph_Graph_theorySyst_copy22->GetXaxis()->SetLabelFont(42);
   Graph_Graph_theorySyst_copy22->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_theorySyst_copy22->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_theorySyst_copy22->GetXaxis()->SetTitleFont(42);
   Graph_Graph_theorySyst_copy22->GetYaxis()->SetLabelFont(42);
   Graph_Graph_theorySyst_copy22->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_theorySyst_copy22->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_theorySyst_copy22->GetYaxis()->SetTitleFont(42);
   Graph_Graph_theorySyst_copy22->GetZaxis()->SetLabelFont(42);
   Graph_Graph_theorySyst_copy22->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_theorySyst_copy22->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_theorySyst_copy22->GetZaxis()->SetTitleFont(42);
   graeT->SetHistogram(Graph_Graph_theorySyst_copy22);
   
   graeT->Draw("2");
   
   TLegend *leg = new TLegend(0.45,0.35,0.8,0.61,NULL,"NB NDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(43);
   leg->SetTextSize(23);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("NULL","Data","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#990000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(0.9);
   entry->SetTextFont(43);
   entry=leg->AddEntry("CentralPointsSystematicUncertainty_copy","Syst. Unc. (data)","f");
   entry->SetFillColor(17);
   entry->SetFillStyle(1001);
   entry->SetLineColor(17);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(43);
   entry=leg->AddEntry("NULL","POWHEG+PYTHIA6 #times A","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#000099");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(43);
   entry=leg->AddEntry("theorySyst_copy","Syst. Unc. (theory)","f");
   entry->SetFillColor(1);

   ci = TColor::GetColor("#000099");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(43);
   leg->Draw();
   
   TPaveText *pt = new TPaveText(0.16,0.64,0.55,0.9,"NB NDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(13);
   pt->SetTextFont(43);
   pt->SetTextSize(22);
   TText *text = pt->AddText("ALICE Preliminary");
   text = pt->AddText("p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
   text = pt->AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4, |#it{#eta}_{jet}| < 0.5");
   text = pt->AddText("with D^{*+}, #it{p}_{T,D^{*+}} > 3 GeV/#it{c}");
   pt->Draw();
   Double_t xAxis4[7] = {5, 6, 8, 10, 14, 20, 30}; 
   
   TH1D *CentralPointsStatisticalUncertainty__4__4 = new TH1D("CentralPointsStatisticalUncertainty__4__4","Central Values",6, xAxis4);
   CentralPointsStatisticalUncertainty__4__4->SetMinimum(2e-05);
   CentralPointsStatisticalUncertainty__4__4->SetMaximum(0.7);
   CentralPointsStatisticalUncertainty__4__4->SetEntries(8);
   CentralPointsStatisticalUncertainty__4__4->SetDirectory(0);
   CentralPointsStatisticalUncertainty__4__4->SetStats(0);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__4__4->SetLineColor(ci);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__4__4->SetMarkerColor(ci);
   CentralPointsStatisticalUncertainty__4__4->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__4__4->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__4__4->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__4__4->GetXaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__4__4->GetXaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__4__4->GetXaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__4__4->GetXaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__4__4->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]");
   CentralPointsStatisticalUncertainty__4__4->GetYaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__4__4->GetYaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__4__4->GetYaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__4__4->GetYaxis()->SetTitleOffset(1.6);
   CentralPointsStatisticalUncertainty__4__4->GetYaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__4__4->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__4__4->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__4__4->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__4__4->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__4__4->Draw("sameaxis");
   FinalSpectrum_1->Modified();
   FinalSpectrum->cd();
  
// ------------>Primitives in pad: FinalSpectrum_2
   FinalSpectrum_2 = new TPad("FinalSpectrum_2", "FinalSpectrum_2",0,0,1,0.35);
   FinalSpectrum_2->Draw();
   FinalSpectrum_2->cd();
   FinalSpectrum_2->Range(-1.986821e-07,-0.9209589,33.33333,2.49);
   FinalSpectrum_2->SetFillColor(0);
   FinalSpectrum_2->SetBorderMode(0);
   FinalSpectrum_2->SetBorderSize(2);
   FinalSpectrum_2->SetGridy();
   FinalSpectrum_2->SetTickx(1);
   FinalSpectrum_2->SetTicky(1);
   FinalSpectrum_2->SetLeftMargin(0.15);
   FinalSpectrum_2->SetTopMargin(0);
   FinalSpectrum_2->SetBottomMargin(0.27);
   FinalSpectrum_2->SetFrameBorderMode(0);
   FinalSpectrum_2->SetFrameBorderMode(0);
   Double_t xAxis5[7] = {5, 6, 8, 10, 14, 20, 30}; 
   
   TH1D *CentralPointsStatisticalUncertainty__5__5 = new TH1D("CentralPointsStatisticalUncertainty__5__5","Central Values",6, xAxis5);
   CentralPointsStatisticalUncertainty__5__5->SetMinimum(0);
   CentralPointsStatisticalUncertainty__5__5->SetMaximum(2.49);
   CentralPointsStatisticalUncertainty__5__5->SetEntries(8);
   CentralPointsStatisticalUncertainty__5__5->SetDirectory(0);
   CentralPointsStatisticalUncertainty__5__5->SetStats(0);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__5__5->SetLineColor(ci);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__5__5->SetMarkerColor(ci);
   CentralPointsStatisticalUncertainty__5__5->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__5__5->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__5__5->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__5__5->GetXaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__5__5->GetXaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__5__5->GetXaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__5__5->GetXaxis()->SetTitleOffset(2.9);
   CentralPointsStatisticalUncertainty__5__5->GetXaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__5__5->GetYaxis()->SetTitle("data / theory");
   CentralPointsStatisticalUncertainty__5__5->GetYaxis()->SetNdivisions(509);
   CentralPointsStatisticalUncertainty__5__5->GetYaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__5__5->GetYaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__5__5->GetYaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__5__5->GetYaxis()->SetTitleOffset(1.4);
   CentralPointsStatisticalUncertainty__5__5->GetYaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__5__5->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__5__5->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__5__5->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__5__5->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__5__5->Draw("axis");
   
   TGraphAsymmErrors *graeR = new TGraphAsymmErrors(6);
   graeR->SetName("ratioSyst");
   graeR->SetTitle("Bayes, iter=4, prior=ResponseTruth Systematics");
   graeR->SetFillColor(17);
   graeR->SetLineColor(17);
   graeR->SetMarkerColor(17);
   graeR->SetPoint(0,5.5,1.597861);
   graeR->SetPointError(0,0.5,0.5,0.2396791,0.2396791);
   graeR->SetPoint(1,7,1.452747);
   graeR->SetPointError(1,1,1,0.1888571,0.1888571);
   graeR->SetPoint(2,9,1.41638);
   graeR->SetPointError(2,1,1,0.1982932,0.1982932);
   graeR->SetPoint(3,12,1.794865);
   graeR->SetPointError(3,2,2,0.2153838,0.2153838);
   graeR->SetPoint(4,17,1.409358);
   graeR->SetPointError(4,3,3,0.1973101,0.1973101);
   graeR->SetPoint(5,25,1.340216);
   graeR->SetPointError(5,5,5,0.2144346,0.2144346);
   /*
   graeR->SetPoint(0,5.5,1.936652);
   graeR->SetPointError(0,0.5,0.5,0.2904978,0.2904978);
   graeR->SetPoint(1,7,1.722126);
   graeR->SetPointError(1,1,1,0.2238763,0.2238763);
   graeR->SetPoint(2,9,1.670958);
   graeR->SetPointError(2,1,1,0.2339341,0.2339341);
   graeR->SetPoint(3,12,2.142237);
   graeR->SetPointError(3,2,2,0.2570685,0.2570685);
   graeR->SetPoint(4,17,1.618254);
   graeR->SetPointError(4,3,3,0.2265555,0.2265555);
   graeR->SetPoint(5,25,1.49058);
   graeR->SetPointError(5,5,5,0.2384927,0.2384927);
   */
   TH1F *Graph_Graph_ratioSyst33 = new TH1F("Graph_Graph_ratioSyst33","Bayes, iter=4, prior=ResponseTruth Systematics",100,2.5,32.5);
   Graph_Graph_ratioSyst33->SetMinimum(0.5839274);
   Graph_Graph_ratioSyst33->SetMaximum(2.28012);
   Graph_Graph_ratioSyst33->SetDirectory(0);
   Graph_Graph_ratioSyst33->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_ratioSyst33->SetLineColor(ci);
   Graph_Graph_ratioSyst33->GetXaxis()->SetLabelFont(42);
   Graph_Graph_ratioSyst33->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_ratioSyst33->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_ratioSyst33->GetXaxis()->SetTitleFont(42);
   Graph_Graph_ratioSyst33->GetYaxis()->SetLabelFont(42);
   Graph_Graph_ratioSyst33->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_ratioSyst33->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_ratioSyst33->GetYaxis()->SetTitleFont(42);
   Graph_Graph_ratioSyst33->GetZaxis()->SetLabelFont(42);
   Graph_Graph_ratioSyst33->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_ratioSyst33->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_ratioSyst33->GetZaxis()->SetTitleFont(42);
   graeR->SetHistogram(Graph_Graph_ratioSyst33);
   
   graeR->Draw("2");
   Double_t xAxis6[7] = {5, 6, 8, 10, 14, 20, 30}; 
   
   TH1D *CentralPointsStatisticalUncertainty__6__6 = new TH1D("CentralPointsStatisticalUncertainty__6__6","Central Values",6, xAxis6);
 
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(0,1.834125);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(1,1.597861);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(2,1.452747);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(3,1.41638);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(4,1.794865);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(5,1.409358);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(6,1.340216);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(7,0.0002553179);

   CentralPointsStatisticalUncertainty__6__6->SetBinError(0,0.06249125);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(1,0.09011263);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(2,0.09930708);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(3,0.1368761);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(4,0.1835618);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(5,0.2618641);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(6,0.3562407);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(7,0.0002005154);
   
   /*
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(0,1.83412);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(1,1.936652);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(2,1.722126);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(3,1.670958);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(4,2.142237);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(5,1.618254);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(6,1.49058);
   CentralPointsStatisticalUncertainty__6__6->SetBinContent(7,0.0002553179);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(0,0.06249125);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(1,0.109219);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(2,0.1177213);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(3,0.161478);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(4,0.2190877);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(5,0.3006778);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(6,0.3962086);
   CentralPointsStatisticalUncertainty__6__6->SetBinError(7,0.0002005154);
   CentralPointsStatisticalUncertainty__6__6->SetMinimum(1.122659e-05);
   CentralPointsStatisticalUncertainty__6__6->SetMaximum(0.06069752);*/
   
   CentralPointsStatisticalUncertainty__6__6->SetEntries(14);
   CentralPointsStatisticalUncertainty__6__6->SetDirectory(0);
   CentralPointsStatisticalUncertainty__6__6->SetStats(0);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__6__6->SetLineColor(ci);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__6__6->SetMarkerColor(ci);
   CentralPointsStatisticalUncertainty__6__6->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__6__6->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__6__6->SetMinimum(0);
   CentralPointsStatisticalUncertainty__6__6->SetMaximum(2.49);
   CentralPointsStatisticalUncertainty__6__6->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__6__6->GetXaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__6__6->GetXaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__6__6->GetXaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__6__6->GetXaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__6__6->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]");
   CentralPointsStatisticalUncertainty__6__6->GetYaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__6__6->GetYaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__6__6->GetYaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__6__6->GetYaxis()->SetTitleOffset(1.4);
   CentralPointsStatisticalUncertainty__6__6->GetYaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__6__6->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__6__6->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__6__6->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__6__6->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__6__6->Draw("same p e0 x0");
   
   TGraphAsymmErrors *graeRT = new TGraphAsymmErrors(6);
   graeRT->SetName("ratioTheorySyst");
   graeRT->SetTitle("Graph");
   graeRT->SetFillColor(1);
   graeRT->SetFillStyle(0);

   ci = TColor::GetColor("#000099");
   graeRT->SetLineColor(ci);
   graeRT->SetLineWidth(2);
   graeRT->SetPoint(0,5.5,1);
   graeRT->SetPointError(0,0.5,0.5,0.6576644,1.185572);
   graeRT->SetPoint(1,7,1);
   graeRT->SetPointError(1,1,1,0.605973,1.129265);
   graeRT->SetPoint(2,9,1);
   graeRT->SetPointError(2,1,1,0.563633,1.039134);
   graeRT->SetPoint(3,12,1);
   graeRT->SetPointError(3,2,2,0.5062495,0.9397957);
   graeRT->SetPoint(4,17,1);
   graeRT->SetPointError(4,3,3,0.4374415,0.9924519);
   graeRT->SetPoint(5,25,1);
   graeRT->SetPointError(5,5,5,0.386771,0.9738843);
   /*
   graeRT->SetPoint(0,5.5,1);
   graeRT->SetPointError(0,0.5,0.5,0.6551554,1.216378);
   graeRT->SetPoint(1,7,1);
   graeRT->SetPointError(1,1,1,0.6121863,1.107202);
   graeRT->SetPoint(2,9,1);
   graeRT->SetPointError(2,1,1,0.5605137,1.026387);
   graeRT->SetPoint(3,12,1);
   graeRT->SetPointError(3,2,2,0.4979863,0.957661);
   graeRT->SetPoint(4,17,1);
   graeRT->SetPointError(4,3,3,0.4587722,0.9812092);
   graeRT->SetPoint(5,25,1);
   graeRT->SetPointError(5,5,5,0.4268945,0.848806);
   */
   
   TH1F *Graph_Graph_ratioTheorySyst44 = new TH1F("Graph_Graph_ratioTheorySyst44","Graph",100,2.5,32.5);
   Graph_Graph_ratioTheorySyst44->SetMinimum(0.2058133);
   Graph_Graph_ratioTheorySyst44->SetMaximum(1.968171);
   Graph_Graph_ratioTheorySyst44->SetDirectory(0);
   Graph_Graph_ratioTheorySyst44->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_ratioTheorySyst44->SetLineColor(ci);
   Graph_Graph_ratioTheorySyst44->GetXaxis()->SetLabelFont(42);
   Graph_Graph_ratioTheorySyst44->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_ratioTheorySyst44->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_ratioTheorySyst44->GetXaxis()->SetTitleFont(42);
   Graph_Graph_ratioTheorySyst44->GetYaxis()->SetLabelFont(42);
   Graph_Graph_ratioTheorySyst44->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_ratioTheorySyst44->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_ratioTheorySyst44->GetYaxis()->SetTitleFont(42);
   Graph_Graph_ratioTheorySyst44->GetZaxis()->SetLabelFont(42);
   Graph_Graph_ratioTheorySyst44->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_ratioTheorySyst44->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_ratioTheorySyst44->GetZaxis()->SetTitleFont(42);
   graeRT->SetHistogram(Graph_Graph_ratioTheorySyst44);
   
   graeRT->Draw("2");
   Double_t xAxis7[7] = {5, 6, 8, 10, 14, 20, 30}; 
   
   TH1D *CentralPointsStatisticalUncertainty__7__7 = new TH1D("CentralPointsStatisticalUncertainty__7__7","Central Values",6, xAxis7);
   CentralPointsStatisticalUncertainty__7__7->SetMinimum(0);
   CentralPointsStatisticalUncertainty__7__7->SetMaximum(2.49);
   CentralPointsStatisticalUncertainty__7__7->SetEntries(8);
   CentralPointsStatisticalUncertainty__7__7->SetDirectory(0);
   CentralPointsStatisticalUncertainty__7__7->SetStats(0);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__7__7->SetLineColor(ci);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__7__7->SetMarkerColor(ci);
   CentralPointsStatisticalUncertainty__7__7->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__7__7->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__7__7->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__7__7->GetXaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__7__7->GetXaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__7__7->GetXaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__7__7->GetXaxis()->SetTitleOffset(2.9);
   CentralPointsStatisticalUncertainty__7__7->GetXaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__7__7->GetYaxis()->SetTitle("data / theory");
   CentralPointsStatisticalUncertainty__7__7->GetYaxis()->SetNdivisions(509);
   CentralPointsStatisticalUncertainty__7__7->GetYaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__7__7->GetYaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__7__7->GetYaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__7__7->GetYaxis()->SetTitleOffset(1.4);
   CentralPointsStatisticalUncertainty__7__7->GetYaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__7__7->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__7__7->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__7__7->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__7__7->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__7__7->Draw("sameaxig");
   FinalSpectrum_2->Modified();
   FinalSpectrum->cd();
   FinalSpectrum->Modified();
   FinalSpectrum->cd();
   FinalSpectrum->SetSelected(FinalSpectrum);
   
   CentralPointsStatisticalUncertainty__1__1->SetName("axis");
   CentralPointsStatisticalUncertainty__2__2->SetName("dataPoints");
   GeneratorLevel_JetPtSpectrum__3__3->SetName("theory");
   CentralPointsStatisticalUncertainty__5__5->SetName("axisRatio");
   CentralPointsStatisticalUncertainty__6__6->SetName("theoryRatio");
   grae->SetName("systUnc");
   graeT->SetName("systUncTheory");
   graeR->SetName("systUncRatio");
   graeRT->SetName("systUncRatioTheory");
   
    TFile *ofile = new TFile("DstarJet_pPb_Preliminary.root","RECREATE");
    CentralPointsStatisticalUncertainty__1__1->Write(); // axis
    CentralPointsStatisticalUncertainty__2__2->Write(); // data
    GeneratorLevel_JetPtSpectrum__3__3->Write(); // theory
    CentralPointsStatisticalUncertainty__5__5->Write(); // ratio axis
    CentralPointsStatisticalUncertainty__6__6->Write(); //ratio to theory
    grae->Write();
    graeT->Write();
    graeR->Write();
    graeRT->Write();
    ofile->Close();
  
   string name = "DstarJet_pPbxsection_v2_newSim";
   FinalSpectrum->SaveAs(Form("%s.png",name.c_str()));
   FinalSpectrum->SaveAs(Form("%s.pdf",name.c_str()));
   FinalSpectrum->SaveAs(Form("%s.eps",name.c_str()));
   FinalSpectrum->SaveAs(Form("%s.root",name.c_str()));
   
    
}

