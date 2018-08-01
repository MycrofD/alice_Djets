//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------

void BkgMatrixRanCones(TString datafile = "fOutDzeroEvt.root", TString histName = "hDeltaPt_ptleadbin5_excluding",
TString postfix = "Djet5Excl", double binmax = 15, TString outDir = "bkgRM03", int Rpar = 3 )
{

	gStyle->SetOptStat(0000); //Mean and RMS shown
	//gStyle->SetPadRightMargin(0.1);

	gStyle->SetLegendFont(42);
	//gStyle->SetLegendTextSize(0.05);
	gStyle->SetPadLeftMargin(0.13);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.15);


	gSystem->Exec(Form("mkdir %s",outDir.Data()));

    double rhoMin = 0.1, rhoMax = 100;
    double DPtMin = -20, DPtMax = 30;

    TPaveText *pv1 = new TPaveText(0.15,0.7,0.5,0.75,"brNDC");
    pv1->SetFillStyle(0);
    pv1->SetBorderSize(0);
    pv1->AddText("Random Cones from data");

    TPaveText *pv2 = new TPaveText(0.15,0.8,0.25,0.9,"brNDC");
    pv2->SetFillStyle(0);
    pv2->SetBorderSize(0);
    pv2->AddText(Form("R=0.%d",Rpar));

    // -------- Data -------
	  TFile *File = new TFile(Form("RandomCones/%s",datafile.Data()),"read");
	  TH1F *hDeltaPt = (TH1F*)File->Get(histName.Data());
		//TH1F *hDeltaPtDraw = (TH1F*)File->Get(histName.Data());
		TH1F *hDeltaPtDraw = (TH1F*)hDeltaPt->Clone("hDeltaPtDraw");
    hDeltaPt->SetTitle();
    hDeltaPt->Sumw2();
    hDeltaPt->GetXaxis()->SetRangeUser(-10,40);
  	hDeltaPt->Scale(1./hDeltaPt->Integral());

    TCanvas *cDPt = new TCanvas("cDPt","cDPt",800,600);
    cDPt->SetLogy();
    cDPt->cd();
    hDeltaPt->Draw();


   // ------ Data -------
   // TH2F *hBkgM = new TH2F("hBkgM","hBkgM",50,0,50,50,0,50);
   double matrixBinMin = 0, matrixBinMaxX = 60, matrixBinMaxY = 60;
   int binsX = 60, binxY = 60;
   TH2F *hBkgM = new TH2F("hBkgM","hBkgM",binsX,matrixBinMin,matrixBinMaxX,binxY,matrixBinMin,matrixBinMaxY);
   hBkgM->Sumw2();
   hBkgM->SetTitle();
   hBkgM->GetXaxis()->SetTitle("p_{T,rec.jet}^{ch} (GeV/#it{c})");
   hBkgM->GetYaxis()->SetTitle("p_{T,gen.jet}^{ch} (GeV/#it{c})");

   double binmin = -5;
	 int nb = (binmax-binmin)*10;
   TH1F *hDeltaPt2 = new TH1F("hDeltaPt2","hDeltaPt2",nb,binmin,binmax);
	 hDeltaPt2->SetTitle();

   double bin = binmin;
   for(int k=1; k<nb; k++){
     hDeltaPt2->SetBinContent(k,hDeltaPt->GetBinContent(hDeltaPt->GetXaxis()->FindBin(bin)));
     bin+=0.1;
   }

  // hDeltaPt2->SetLineColor(kMagenta);
  // hDeltaPt2->Draw();

   TRandom3 *rpt = new TRandom3(0);
   const int sampling = 100000;
   double ptgen = 0.5, ptrec = 0.5;
   double delta = 0;

   for(int i=0; i<2*binsX; i++){
     for (int j=0; j<sampling; j++){
       delta = hDeltaPt2->GetRandom();
       ptrec = ptgen + delta;
       if(ptrec<0 || ptrec>matrixBinMaxX) continue;
       hBkgM->Fill(ptrec,ptgen);
           // delta = hDeltaPt->GetRandom();
         // cout << '\n' << delta;

          //  delta = hDeltaPt->GetRandom();
          // delta = rpt->Gaus(0.5,0.29);
         //  delta = rpt->Gaus(0.94,0.15);
           //delta = rpt->Gaus(0.5,0.29);
           //delta = fP->GetRandom();
           //delta = rpt->Landau(fP->GetParameter(0),fP->GetParameter(1));

      }
       //ptgen += 0.05;
       ptgen += 0.5;
    }
    hBkgM->Scale(1./(sampling));
    hBkgM->Scale(1.,"width");
    hBkgM->GetZaxis()->SetRangeUser(0,1);
    hBkgM->SetName("hBkgM");

    TH1F *hproj = (TH1F*)hBkgM->ProjectionY();
    TCanvas *c2 = new TCanvas;
    hproj->Draw();

    TCanvas *cBkgM = new TCanvas("cBkgM","cBkgM",800,600);
    cBkgM->cd();
    cBkgM->SetLogz();
    hBkgM->Draw("colz");
    pv1->Draw("same");
    pv2->Draw("same");



		hDeltaPtDraw->Rebin(2);
		//hDeltaPtDraw->Scale(1.,"width");
		hDeltaPtDraw->Scale(1./hDeltaPtDraw->Integral());
		hDeltaPtDraw->Scale(1.,"width");
		hDeltaPtDraw->SetMarkerColor(1);
		hDeltaPtDraw->SetLineColor(1);
		hDeltaPtDraw->SetMarkerStyle(20);
		hDeltaPtDraw->SetMarkerSize(1.2);

		hDeltaPtDraw->GetXaxis()->SetTitle("#delta#it{p}_{T,ch} (GeV/#it{c})");
		hDeltaPtDraw->GetYaxis()->SetTitle("Probability density (GeV/#it{c})^{-1}");
		hDeltaPtDraw->GetXaxis()->SetLabelSize(0.04);
		hDeltaPtDraw->GetXaxis()->SetTitleSize(0.05);
		hDeltaPtDraw->GetXaxis()->SetTitleOffset(1.);
		hDeltaPtDraw->GetYaxis()->SetLabelSize(0.045);
		hDeltaPtDraw->GetYaxis()->SetTitleSize(0.05);
		hDeltaPtDraw->GetYaxis()->SetTitleOffset(1.1);
		hDeltaPtDraw->GetXaxis()->SetRangeUser(-7,25);
		hDeltaPtDraw->SetMaximum(hDeltaPt2->GetMaximum()*100);


		TH1F *hDeltaPtDrawFill = (TH1F*)hDeltaPtDraw->Clone("hDeltaPtDrawFill");
		hDeltaPtDrawFill->SetFillStyle(1001);
		hDeltaPtDrawFill->SetFillColor(kGray);

TPaveText *pvALICE = new TPaveText(0.15,0.85,0.8,0.9,"brNDC");
pvALICE->SetFillStyle(0);
pvALICE->SetBorderSize(0);
pvALICE->SetTextFont(42);
pvALICE->SetTextSize(0.045);
pvALICE->SetTextAlign(11);
pvALICE->AddText("ALICE Preliminary");

TPaveText *pvEn= new TPaveText(0.15,0.80,0.8,0.85,"brNDC");
pvEn->SetFillStyle(0);
pvEn->SetBorderSize(0);
pvEn->SetTextFont(42);
pvEn->SetTextSize(0.045);
pvEn->SetTextAlign(11);
pvEn->AddText("p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");

TPaveText *pvCent= new TPaveText(0.15,0.74,0.8,0.79,"brNDC");
pvCent->SetFillStyle(0);
pvCent->SetBorderSize(0);
pvCent->SetTextFont(42);
pvCent->SetTextSize(0.045);
pvCent->SetTextAlign(11);
pvCent->AddText("Minimum bias");

double shift = 0.;

     TPaveText *pvEta = new TPaveText(0.6,0.72,0.8,0.82,"brNDC");
    pvEta->SetFillStyle(0);
    pvEta->SetBorderSize(0);
    pvEta->SetTextFont(42);
    pvEta->SetTextSize(0.045);
    pvEta->SetTextAlign(11);
    pvEta->AddText("#splitline{|#it{#eta}_{lab}| < 0.9}{#it{p}_{T, track} > 0.15 GeV/#it{c}}");


    TPaveText *pvJet = new TPaveText(0.6,0.6,0.8,0.7,"brNDC");
    pvJet->SetFillStyle(0);
    pvJet->SetBorderSize(0);
    pvJet->SetTextFont(42);
    pvJet->SetTextSize(0.045);
    pvJet->SetTextAlign(11);
    pvJet->AddText("#splitline{#it{R} = 0.3}{Random Cones}");


TPaveText *pvD = new TPaveText(0.6,0.48,0.8,0.59,"brNDC");
pvD->SetFillStyle(0);
pvD->SetBorderSize(0);
pvD->SetTextFont(42);
pvD->SetTextSize(0.045);
pvD->SetTextAlign(11);
pvD->AddText("#splitline{with D^{0} #rightarrow K^{-}#pi^{+}}{and charge conj.}");


TCanvas *cDelta = new TCanvas("cEff","cEff",1200,900);
//cEff->SetBatch();
cDelta->SetLogy();
//cMass->Divide(3,1);
//cMass->cd(1);
hDeltaPtDrawFill->Draw("hist");
hDeltaPtDraw->Draw("ep same");
pvALICE->Draw("same");
pvEn->Draw("same");
pvJet->Draw("same");
pvD->Draw("same");
pvEta->Draw("same");

cout << "=== integral: " << hDeltaPtDraw->Integral() << "\t bin width: " << hDeltaPtDraw->GetBinWidth(1) << endl;

    // ---- SAVE ----
		cDelta->SaveAs(Form("%s/DeltaPt2_%s.png",outDir.Data(),postfix.Data()));
		cDelta->SaveAs(Form("%s/DeltaPt2_%s.pdf",outDir.Data(),postfix.Data()));
    cDelta->SaveAs(Form("%s/DeltaPt2_%s.eps",outDir.Data(),postfix.Data()));

    cBkgM->SaveAs(Form("%s/BkgMatrix_%s.png",outDir.Data(),postfix.Data()));
    c2->SaveAs(Form("%s/BkgMatrixGenProj_%s.png",outDir.Data(),postfix.Data()));


    TFile *ofile = new TFile(Form("%s/RandCones_BkgM_%s.root",outDir.Data(),postfix.Data()),"RECREATE");

    hDeltaPt2->Write();
    hBkgM->Write();
    hproj->Write();
    ofile->Close();

}
