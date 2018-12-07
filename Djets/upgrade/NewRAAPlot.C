{
      
    // systematics
    
    /*
    jets:
    track res
    1% above 60 GeV
    tracking efficiency: 1% below 20 GeV, 5% 20-50 GeV, 10% > 50 GeV; more precisly
    30-40 1.5
	40-50 2.7
	50-60 5.9
	60-70 8.1
	70-90 9.7
	delta pT: 5.1 4.1 3.8 2.9 2.2
	unf method: 4%
	reg parameter (SVD): 0.4 1.1 3.3 3.6 5.4
	prior: 4.2 5.2 5.2 4.0 1.7
*/
   
    Double_t unfSys[14] = { 
		8, // 5-6
		5, // 6-8
		5, // 8-10
		5, // 10-12
		5, // 12-14
		5, // 14-16
		5, // 16-20
		5, // 20-25
		5, // 25-30
		5, // 30-40
		5, // 40-50
		5, // 50-60
		5, // 60-80
		8 // 80-100 
		};
		
    Double_t JESSys[14] = { 
		1, // 5-6
		1, // 6-8
		2, // 8-10
		2, // 10-12
		2, // 12-14
		3, // 14-16
		3, // 16-20
		3, // 20-25
		3, // 25-30
		4, // 30-40
		4, // 40-50
		5, // 50-60
		7, // 60-80
		9 // 80-100 
		}; 
	//Double_t JetTracking = 1;
	
	Double_t DTracking = 3.5;
    Double_t DyieldExtSys = 1;
    Double_t DcutVariation = 3.5;
   
   Double_t BFDRaaSys[14] = { 
		2.69, // 5-6
		2.66, // 6-8
		2.79, // 8-10
		2.87, // 10-12
		2.95, // 12-14
		3.12, // 14-16
		3.09, // 16-20
		3.01, // 20-25
		2.21, // 25-30
		1.43, // 30-40
		1.17, // 40-50
		1.08, // 50-60
		1.09, // 60-80
		0.93 // 80-100 
		}; 
		
		Double_t BFDPowhegSys[14] = { 
		3, // 5-6
		3, // 6-8
		3, // 8-10
		3, // 10-12
		3, // 12-14
		3, // 14-16
		3, // 16-20
		3, // 20-25
		3, // 25-30
		3, // 30-40
		3, // 40-50
		3, // 50-60
		3, // 60-80
		3 // 80-100 
		}; 
   
   /*
     Double_t BFDPowhegSys[14] = { 
		8.2, // 5-6
		9.5, // 6-8
		11.0, // 8-10
		11.9, // 10-12
		12.8, // 12-14
		13.4, // 14-16
		14.2, // 16-20
		14.9, // 20-25
		13.6, // 25-30
		11.4, // 30-40
		11.6, // 40-50
		11.7, // 50-60
		12.2, // 60-80
		20.0 // 80-100 
		}; 
		
	Double_t BFDRaaSys[14] = { 
		2.9, // 5-6
		4.5, // 6-8
		6.6, // 8-10
		7.5, // 10-12
		8.6, // 12-14
		10.1, // 14-16
		11.6, // 16-20
		12.9, // 20-25
		11.2, // 25-30
		8.4, // 30-40
		7.6, // 40-50
		7.8, // 50-60
		8.5, // 60-80
		8.6 // 80-100 
		}; 
		*/
    
    TFile *D0File = new TFile("Average_DmesonRaa_010_relStatUncAndUncorrSystWPID.root","READ");
    TH1F *D0pt = (TH1F*)D0File->Get("hDmesonAverageRAB");
    TGraphAsymmErrors *D0SU = (TGraphAsymmErrors*)D0File->Get("gRAB_DmesonAverage_GlobalSystematics");
     
    TFile *D0JetFile = new TFile("D0jets.root","READ");
    TH1F *D0JetptUnc = (TH1F*)D0JetFile->Get("hjetptspectrumRebUnc");
    TH1F *D0Jet = (TH1F*)D0JetFile->Get("hjetptspectrumRebUnc");
    
    const Int_t D0bins = 9;
    //Double_t D0ptbins[D0bins+1] = {1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 24., 36., 50.};
    Double_t D0ptbins[D0bins+1] = {5., 6., 7., 8., 10., 12., 16., 24., 36., 50.};
    
    const Int_t D0Jetbins = 14;
    Double_t D0Jetptbins[D0Jetbins+1] = { 5,6,8,10,12,14,16,20,25,30,40,50,60,80,100 };
    Double_t sysunc[D0Jetbins];

    TH1D *D0Spectrum = new TH1D("D0Spectrum","D0Spectrum",D0bins,D0ptbins);
        
	TH1D *D0JetSpectrum = new TH1D("D0JetSpectrum","D0JetSpectrum",D0Jetbins,D0Jetptbins);
    TH1D *D0JetSysUnc = new TH1D("D0JetSysUnc","D0JetSysUnc",D0Jetbins,D0Jetptbins);

    for(Int_t i=1; i<=D0bins; i++)
    {
        D0Spectrum->SetBinContent(i,D0pt->GetBinContent(D0pt->FindBin(D0ptbins[i-1])));
        D0Spectrum->SetBinError(i,D0pt->GetBinError(D0pt->FindBin(D0ptbins[i-1])));
        D0SU->SetPointEXhigh(i-1,D0Spectrum->GetBinWidth(i)/4.);
        D0SU->SetPointEXlow(i-1,D0Spectrum->GetBinWidth(i)/4.);
    }
    
    
    TF1 *fD0Fit = new TF1("fD0Fit","[0]+[1]*x",8,50);
    D0Spectrum->Fit("fD0Fit","ERM0");
    fD0Fit->SetLineColor(2);
    
    TF1 *fD0Fit2 = new TF1("fD0Fit2","[0]+[1]*x",8,100);
	fD0Fit2->SetParameter(0,fD0Fit->GetParameter(0));
	fD0Fit2->SetParameter(1,fD0Fit->GetParameter(1));
    fD0Fit2->SetLineColor(2);
    
     for(Int_t i=1; i<=D0Jetbins; i++)
    {
		Double_t content;
		Double_t bin = D0JetSpectrum->GetBinCenter(D0JetSpectrum->FindBin(D0Jetptbins[i-1]));
		if(bin<8) content = D0pt->GetBinContent(D0pt->FindBin(D0Jetptbins[i-1]));
		else content = fD0Fit2->Eval(bin);
		D0JetSpectrum->SetBinContent(i,content);
		//cout << D0JetptUnc->GetBinContent(D0JetptUnc->FindBin(D0Jetptbins[i-1]))  << endl;
		Double_t unc = content * ( D0JetptUnc->GetBinContent(D0JetptUnc->FindBin(D0Jetptbins[i-1])) ) ;
        D0JetSpectrum->SetBinError(i,unc);
        
        sysunc[i-1] = unfSys[i-1]/100*unfSys[i-1]/100 + JESSys[i-1]/100*JESSys[i-1]/100 + DTracking/100*DTracking/100 + DyieldExtSys/100*DyieldExtSys/100 + DcutVariation/100*DcutVariation/100 + BFDPowhegSys[i-1]/100*BFDPowhegSys[i-1]/100 + BFDRaaSys[i-1]/100*BFDRaaSys[i-1]/100;
        sysunc[i-1] = TMath::Sqrt(sysunc[i-1]);
         cout << sysunc[i-1]*100 << endl;
        sysunc[i-1] *= content;
       
        
        //D0SU->SetPointEXhigh(i-1,D0JetSpectrum->GetBinWidth(i)/4.);
        //D0SU->SetPointEXlow(i-1,D0JetSpectrum->GetBinWidth(i)/4.);
    }
    
    
    
//=========Macro generated from canvas: FinalSpectrum/FinalSpectrum
//=========  (Tue Sep  5 14:25:17 2017) by ROOT version5.34/30

    enum EColor colorDjet = Rtypes::kGray+3;
    //Int_t colorDjet = TColor::GetColor("#990000");
    enum EColor colorD0 = Rtypes::kBlue;
    enum EColor colorChJet = Rtypes::kRed;

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
   TPad *FinalSpectrum_1 = new TPad("FinalSpectrum_1", "FinalSpectrum_1",0,0,1,1);
   FinalSpectrum_1->Draw();
   FinalSpectrum_1->cd();
   FinalSpectrum_1->Range(-1.986821e-07,-2.958607,33.33333,1.96998);
   FinalSpectrum_1->SetFillColor(0);
   FinalSpectrum_1->SetBorderMode(0);
   FinalSpectrum_1->SetBorderSize(2);
   FinalSpectrum_1->SetLogx();
   FinalSpectrum_1->SetTickx(1);
   FinalSpectrum_1->SetTicky(1);
   FinalSpectrum_1->SetLeftMargin(0.12);
   FinalSpectrum_1->SetRightMargin(0.05);
   FinalSpectrum_1->SetTopMargin(0.05);
   //FinalSpectrum_1->SetBottomMargin(0);
   FinalSpectrum_1->SetFrameBorderMode(0);
   FinalSpectrum_1->SetFrameBorderMode(0);
   Double_t xAxis1[6] = {4, 5, 10, 15, 20, 110}; 
   
   
   
   TH1D *CentralPointsStatisticalUncertainty__1__1 = new TH1D("CentralPointsStatisticalUncertainty__1__1","Central Values",5, xAxis1);
   CentralPointsStatisticalUncertainty__1__1->SetMinimum(0);
   CentralPointsStatisticalUncertainty__1__1->SetMaximum(1.65);
   CentralPointsStatisticalUncertainty__1__1->SetEntries(3);
   CentralPointsStatisticalUncertainty__1__1->SetDirectory(0);
   CentralPointsStatisticalUncertainty__1__1->SetStats(0);
   //CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetNdivisions(25);
   

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#990000");
   

   
   CentralPointsStatisticalUncertainty__1__1->SetLineColor(colorDjet);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__1__1->SetMarkerColor(colorDjet);
   CentralPointsStatisticalUncertainty__1__1->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__1__1->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetTitle("#it{p}_{T,D} and #it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   //CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__1__1->GetYaxis()->SetTitle("#it{R}_{AA}");
   CentralPointsStatisticalUncertainty__1__1->GetYaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__1__1->GetYaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__1__1->GetYaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__1__1->GetYaxis()->SetTitleOffset(1.4);
   CentralPointsStatisticalUncertainty__1__1->GetXaxis()->SetTitleOffset(1.3);
   CentralPointsStatisticalUncertainty__1__1->GetYaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__1__1->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__1__1->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__1__1->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__1__1->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__1__1->Draw("axis");
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(3);
   grae->SetName("CentralPointsSystematicUncertainty_copy");
   grae->SetTitle("Bayes, iter=4, prior=ResponseTruth Systematics");
   grae->SetFillColor(17);
   grae->SetLineColor(17);
   grae->SetMarkerColor(17);
   grae->SetPoint(2,7.5,0.209306);
   grae->SetPointError(2,0.5,0.5,0.0703684,0.0703684);
   grae->SetPoint(3,12.5,0.192116);
   grae->SetPointError(3,0.5,0.5,0.0655228,0.0655228);
   grae->SetPoint(4,17.5,0.129345);
   grae->SetPointError(4,0.5,0.5,0.0501131,0.0501131);
   
   
   TH1F *Graph_Graph_central_syst_unc11 = new TH1F("Graph_Graph_central_syst_unc11","Bayes, iter=4, prior=ResponseTruth Systematics",100,2.5,32.5);
   Graph_Graph_central_syst_unc11->SetMinimum(4.779682e-05);
   Graph_Graph_central_syst_unc11->SetMaximum(0.02142993);
   Graph_Graph_central_syst_unc11->SetDirectory(0);
   Graph_Graph_central_syst_unc11->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_central_syst_unc11->SetLineColor(kGreen+2);
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
   
   //grae->Draw("2");
   Double_t xAxis2[4] = {5, 10, 15, 20}; 
   
   ci = TColor::GetColor("#000099");
   //D0Spectrum->Draw("same"); // D0 central points
    D0Spectrum->SetMarkerStyle(kOpenCircle);
    D0Spectrum->SetMarkerColor(colorD0);
    D0Spectrum->SetLineColor(colorD0);
    D0Spectrum->SetLineWidth(2);
    
    
   
    
    TGraphAsymmErrors *gD0 = new TGraphAsymmErrors(7);
   gD0->SetName("gD0");
   gD0->SetTitle("gD0");
   gD0->SetFillColor(1);
   gD0->SetFillStyle(0);

   ci = TColor::GetColor("#000099");
   gD0->SetLineColor(colorD0);
   gD0->SetLineWidth(2);
   
   //fD0Fit2->Draw("same");
   
    gD0->SetPoint(0,5.5,D0Spectrum->GetBinContent(1));
   gD0->SetPointError(0,0.25,0.25,D0SU->GetErrorYlow(4),D0SU->GetErrorYhigh(4));
   gD0->SetPoint(1,6.5,D0Spectrum->GetBinContent(2));
   gD0->SetPointError(1,0.25,0.25,D0SU->GetErrorYlow(5),D0SU->GetErrorYhigh(5));
   gD0->SetPoint(2,7.5,D0Spectrum->GetBinContent(3));
   gD0->SetPointError(2,0.25,0.25,D0SU->GetErrorYlow(6),D0SU->GetErrorYhigh(6));
   
   gD0->SetPoint(3,9.,D0Spectrum->GetBinContent(4));
   gD0->SetPointError(3,0.5,0.5,D0SU->GetErrorYlow(7),D0SU->GetErrorYhigh(7));
   gD0->SetPoint(4,11.,D0Spectrum->GetBinContent(5));
   gD0->SetPointError(4,0.5,0.5,D0SU->GetErrorYlow(8),D0SU->GetErrorYhigh(8));
   gD0->SetPoint(5,14.,D0Spectrum->GetBinContent(6));
   gD0->SetPointError(5,1.,1.,D0SU->GetErrorYlow(9),D0SU->GetErrorYhigh(9));
   
   gD0->SetPoint(6,20.,D0Spectrum->GetBinContent(7));
   gD0->SetPointError(6,2.,2.,D0SU->GetErrorYlow(10),D0SU->GetErrorYhigh(10));
   
   gD0->SetPoint(7,30.,D0Spectrum->GetBinContent(8));
   gD0->SetPointError(7,3.,3.,D0SU->GetErrorYlow(11),D0SU->GetErrorYhigh(11));
   
   gD0->SetPoint(8,43.,D0Spectrum->GetBinContent(9));
   gD0->SetPointError(8,3.5,3.5,D0SU->GetErrorYlow(12),D0SU->GetErrorYhigh(12));
  
   //gD0->Draw("2"); // D0 sys un, wider boxes
   
   D0SU->SetFillColor(1);
   D0SU->SetFillStyle(0);

   ci = TColor::GetColor("#000099");
   D0SU->SetLineColor(colorD0);
   D0SU->SetLineWidth(2);
   
   // D0SU->Draw("2"); // D0 syst unc, narrower boxes
   
   TH1D *CentralPointsStatisticalUncertainty__2__2 = new TH1D("CentralPointsStatisticalUncertainty__2__2","Central Values",3, xAxis2);
   CentralPointsStatisticalUncertainty__2__2->SetBinContent(1,0.209306);
   CentralPointsStatisticalUncertainty__2__2->SetBinContent(2,0.192116);
   CentralPointsStatisticalUncertainty__2__2->SetBinContent(3,0.129345);
   
   CentralPointsStatisticalUncertainty__2__2->SetBinError(1,0.0149163);
   CentralPointsStatisticalUncertainty__2__2->SetBinError(2,0.0329361);
   CentralPointsStatisticalUncertainty__2__2->SetBinError(3,0.0414657);
   
   CentralPointsStatisticalUncertainty__2__2->SetMinimum(2e-05);
   CentralPointsStatisticalUncertainty__2__2->SetMaximum(4);
   CentralPointsStatisticalUncertainty__2__2->SetEntries(3);
   CentralPointsStatisticalUncertainty__2__2->SetDirectory(0);
   CentralPointsStatisticalUncertainty__2__2->SetStats(0);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__2__2->SetLineColor(colorDjet);
   CentralPointsStatisticalUncertainty__2__2->SetLineWidth(2);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__2__2->SetMarkerColor(colorDjet);
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
   CentralPointsStatisticalUncertainty__2__2->Draw("same");
   Double_t xAxis3[4] = {5, 10, 15, 20}; 
   
   TH1D *GeneratorLevel_JetPtSpectrum__3__3 = new TH1D("GeneratorLevel_JetPtSpectrum__3__3","D0_MCTruth_Charged_R040_JetPtDPtSpectrum",3, xAxis3);

   //GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(0,0.0297059);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(1,0.0941487);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(2,0.0108335);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinContent(3,0.00237272);

   //GeneratorLevel_JetPtSpectrum__3__3->SetBinError(0,0.000035151);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(1,0.0018234);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(2,0.000517178);
   GeneratorLevel_JetPtSpectrum__3__3->SetBinError(3,0.000222742);
   
   GeneratorLevel_JetPtSpectrum__3__3->SetEntries(316731);
   GeneratorLevel_JetPtSpectrum__3__3->SetDirectory(0);
   GeneratorLevel_JetPtSpectrum__3__3->SetStats(0);

   ci = TColor::GetColor("#000099");
   GeneratorLevel_JetPtSpectrum__3__3->SetLineColor(kGreen+2);

   ci = TColor::GetColor("#000099");
   GeneratorLevel_JetPtSpectrum__3__3->SetMarkerColor(kGreen+2);
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
   //GeneratorLevel_JetPtSpectrum__3__3->Draw("same p e0 x0");
   
   TGraphAsymmErrors *graeT = new TGraphAsymmErrors(3);
   graeT->SetName("theorySyst_copy");
   graeT->SetTitle("Graph");
   graeT->SetFillColor(1);
   graeT->SetFillStyle(0);

   
   graeT->SetLineColor(colorDjet);
   graeT->SetLineWidth(2);
   
   Double_t CombError[3] = {0,0,0};
   CombError[0] = TMath::Sqrt(TMath::Power(0.0703684,2)+TMath::Power(0.0251167,2));
   CombError[1] = TMath::Sqrt(TMath::Power(0.0655228,2)+TMath::Power(0.0230539,2));
   CombError[2] = TMath::Sqrt(TMath::Power(0.0501131,2)+TMath::Power(0.0181083,2));
   
    graeT->SetPoint(0,7.5,0.209306);
   graeT->SetPointError(0,0.5,0.5,CombError[0],CombError[0]);
   graeT->SetPoint(1,12.5,0.192116);
   graeT->SetPointError(1,0.5,0.5,CombError[1],CombError[1]);
   graeT->SetPoint(2,17.5,0.129345);
   graeT->SetPointError(2,0.5,0.5,CombError[2],CombError[2]);
   
   TH1F *Graph_Graph_theorySyst_copy22 = new TH1F("Graph_Graph_theorySyst_copy22","Graph",100,2.5,32.5);
   Graph_Graph_theorySyst_copy22->SetMinimum(2.864129e-05);
   Graph_Graph_theorySyst_copy22->SetMaximum(0.01825038);
   Graph_Graph_theorySyst_copy22->SetDirectory(0);
   Graph_Graph_theorySyst_copy22->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_theorySyst_copy22->SetLineColor(kGreen+2);
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
   
   graeT->Draw("2"); // D-jet sys uncertanties
   
   
   // D-jets projections
   
    D0JetSpectrum->Draw("same"); // D0 central points
    D0JetSpectrum->SetMarkerStyle(kOpenCircle);
   // D0JetSpectrum->SetMarkerColor(colorD0);
    //D0JetSpectrum->SetLineColor(colorD0);
     D0JetSpectrum->SetMarkerColor(kRed+1);
    D0JetSpectrum->SetLineColor(kRed+1);
    D0JetSpectrum->SetLineWidth(2);
    
    
    TGraphAsymmErrors *graeProj = new TGraphAsymmErrors(D0Jetbins);
   graeProj->SetName("Projections sys unc");
   graeProj->SetTitle("Graph");
   graeProj->SetFillColor(1);
   graeProj->SetFillStyle(0);

   
   graeProj->SetLineColor(kRed+1);
   graeProj->SetLineWidth(2);
   
 
   for(int i=0; i<D0Jetbins; i++) {
		double width;
		if(!i)width = 0.5;
		else if(i<8) width = 0.6;
		else if(i<10) width = 1;
		else if(i<12) width = 1.5;
		else  width = 2;
		graeProj->SetPoint(i,D0JetSpectrum->GetBinCenter(i+1),D0JetSpectrum->GetBinContent(i+1));
		graeProj->SetPointError(i,width,width,sysunc[i],sysunc[i]);
	}
   
   
   
   TH1F *Graph_Graph_Proj = new TH1F("Graph_Graph_Proj","Graph",100,2.5,32.5);
   Graph_Graph_Proj->SetMinimum(2.864129e-05);
   Graph_Graph_Proj->SetMaximum(0.01825038);
   Graph_Graph_Proj->SetDirectory(0);
   Graph_Graph_Proj->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Proj->SetLineColor(kGreen+2);
   Graph_Graph_Proj->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Proj->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Proj->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Proj->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Proj->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Proj->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Proj->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Proj->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Proj->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Proj->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Proj->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_theorySyst_copy22->GetZaxis()->SetTitleFont(42);
   graeProj->SetHistogram(Graph_Graph_Proj);
   
   graeProj->Draw("2"); // D-jet sys uncertanties
   
   
   //Ch Jets
   
   TFile *ChFile = new TFile("fileRAA.root","READ");
    TGraphAsymmErrors *ChJets = (TGraphAsymmErrors*)ChFile->Get("gRAA_C0010_R030");
    ChJets->SetName("ChJets");
    ChJets->SetTitle("ChJets");
    TGraphAsymmErrors *ChJetsSysUnc = (TGraphAsymmErrors*)ChFile->Get("gRAAc_C0010_R030");
   
    const Int_t nBinsChJet = 5;
    Double_t ptBinsChJet[nBinsChJet+1] = {50., 60., 70., 80., 90., 100.};
    
    TH1D *ChJetsCentral = new TH1D("ChJetsCentral","ChJetsCentral",nBinsChJet,ptBinsChJet);
   
   
    ChJetsCentral->SetMarkerStyle(kFullSquare);
    ChJetsCentral->SetMarkerSize(1);
    ChJetsCentral->SetMarkerColor(colorChJet);
    ChJetsCentral->SetLineColor(colorChJet);
    ChJetsCentral->SetLineWidth(2);
   
   Double_t *ChJetY = ChJets->GetY();
   Double_t *ChJetX = ChJets->GetX();
   
   for(Int_t i=1; i<=nBinsChJet; i++)
   {
       ChJetsCentral->SetBinContent(i,ChJetY[i-1]);
       ChJetsCentral->SetBinError(i,ChJets->GetErrorYhigh(i-1));
       ChJetsSysUnc->SetPointEXlow(i-1,2.5);
       ChJetsSysUnc->SetPointEXhigh(i-1,2.5);
   }
   
    
   //ChJetsCentral->Draw("same");
   ChJetsSysUnc->SetFillColor(1);
   ChJetsSysUnc->SetFillStyle(0);
   ChJetsSysUnc->SetLineColor(colorChJet);
   ChJetsSysUnc->SetLineWidth(2);
   //ChJetsSysUnc->Draw("2same");
   
   TLegend *leg = new TLegend(0.1,0.5,0.5,0.6,NULL,"NB NDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(43);
   leg->SetTextSize(20);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("NULL","ALICE Preliminary, 0-20\% Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#990000");
   entry->SetMarkerColor(colorDjet);
   entry->SetLineColor(colorDjet);
   entry->SetLineWidth(2);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(0.9);
   entry->SetTextFont(43);
   
   entry=leg->AddEntry("NULL","D^{0}-tagged jets (#it{p}_{T,D^{0}} > 3 GeV/#it{c}), p-Pb data reference","");
	//entry=leg->AddEntry("NULL","(p-Pb data reference)","");
   
   //entry=leg->AddEntry(gD0,"D^{0} Sys. Unc.","f");
   //entry=leg->AddEntry(ChJetsCentral,"Ch. Jets (#it{p}_{T}^{lead} > 5 GeV/#it{c}), 0-10\%","p");
   //entry=leg->AddEntry(ChJetsCentral,"POWHEG+PYTHIA8 reference","");
   //entry=leg->AddEntry(D0Spectrum,"D^{0}, 0-10\%","lp");
   //entry=leg->AddEntry(ChJetsSysUnc,"Charged Jets Sys. Unc.","f");
   
   leg->Draw();
   //0.14296  0.492098
   TLegend *legD0 = new TLegend(0.15,0.65,0.45,0.7,NULL,"NB NDC");
   legD0->SetBorderSize(0);
   legD0->SetTextFont(43);
   legD0->SetTextSize(23);
   legD0->SetLineColor(1);
   legD0->SetLineStyle(1);
   legD0->SetLineWidth(1);
   legD0->SetFillColor(0);
   legD0->SetFillStyle(0);
   entry=legD0->AddEntry(D0Spectrum,"Average D^{0}, D^{+}, D^{*+}, 0-10\%, arxiv:1804.09083","p");
   //legD0->Draw();
   
   TLine *line = new TLine(4.,1.,110.,1.);
   line->SetLineStyle(2);
   line->Draw();
   
   TPaveText *pt = new TPaveText(0.15,0.73,0.65,0.93,"NB NDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(13);
   pt->SetTextFont(43);
   pt->SetTextSize(22);
   TText *text = pt->AddText("ALICE Upgrade Simulation");
   text = pt->AddText("0-10% Pb-Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV, #it{L} = 10 nb^{-1}");
   text = pt->AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.3, |#it{#eta}_{jet}| < 0.6");
   text = pt->AddText("with D^{0}, #it{p}_{T,D^{0}} > 3 GeV/#it{c}");
   pt->Draw();
   
   
    
}
