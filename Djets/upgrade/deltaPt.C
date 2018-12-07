#include <string>
#include <sstream>
#include <iostream>


void deltaPt(int Dptmin = 0, int Dptmax = 100)
{
    
  
     gStyle->SetOptStat(0000); //Mean and RMS shown
    

    double zmin = -2., zmax = 2;
    float jetmin = 0, jetmax = 100;
    double plotmin = 0, plotmax = 50;


		TString datfile = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_Hijing_RC_Bis";
        TString histName = "histosD0UpgradeN";
        TH1D *hDeltaPt;
        // loop over D mesons
        for(int j=1; j<7; j++) {

			TFile *File = new TFile(Form("%s%d.root",datfile.Data(),j),"read");
			if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
			TDirectoryFile* dir = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
			for(int i=0;i<1; i++){
				TList *histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
				TH1D *hdelta = (TH1D*)histList->FindObject("fDeltaPT_excl_lead");
				if(i==0 && j==1) hDeltaPt = (TH1D*)hdelta->Clone("hDeltaPt");
				else hDeltaPt->Add(hdelta);
		
			}
		}
        
  
		TString datfile = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_Hijing_RC_Tres";

        // loop over D mesons
        for(int j=1; j<2; j++) {

			TFile *File = new TFile(Form("%s%d.root",datfile.Data(),j),"read");
			if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
			TDirectoryFile* dir = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
			for(int i=0;i<1; i++){
				TList *histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
				TH1D *hdelta = (TH1D*)histList->FindObject("fDeltaPT_excl_lead");
				//if(i==0 && j==1) hDeltaPt = (TH1D*)hdelta->Clone("hDeltaPt");
				//else 
				hDeltaPt->Add(hdelta);
		
			}
		}
  
  
	 // gaussian fit, just a x-check
      TF1 *fitG = new TF1("fitG","gaus",-20,20);
      hDeltaPt->Fit("fitG","0RM");
      fitG->SetLineColor(2);
  
	 TF1 *fit = new TF1("fit","gaus",-40,40);
	 fit->SetParameters(fitG->GetParameter(0),fitG->GetParameter(1),fitG->GetParameter(2));
	 fit->SetLineColor(2);
  
	TCanvas *cDelta = new TCanvas("cDelta","cDelta",1200,800);
    hDeltaPt->Draw();
    fit->Draw("same");
    
    cDelta->SaveAs("deltaPt_Hijing.png");
  
	 double matrixBinMinX = -20,  matrixBinMinY = 0, matrixBinMaxX = 100, matrixBinMaxY = 100;
   int binsX = 120, binsY = 100;
   TH2D *hBkgM = new TH2D("hBkgM","hBkgM",binsX,matrixBinMinX,matrixBinMaxX,binsY,matrixBinMinY,matrixBinMaxY);
   hBkgM->Sumw2();
   hBkgM->SetTitle();
   hBkgM->GetXaxis()->SetTitle("p_{T,rec.jet}^{ch} (GeV/#it{c})");
   hBkgM->GetYaxis()->SetTitle("p_{T,gen.jet}^{ch} (GeV/#it{c})");
   
   TH2D *hBkgM_trans = new TH2D("hBkgM_trans","hBkgM_trans",binsY,matrixBinMinY,matrixBinMaxY,binsX,matrixBinMinX,matrixBinMaxX);
   hBkgM_trans->Sumw2();
   hBkgM_trans->SetTitle();
   hBkgM_trans->GetXaxis()->SetTitle("p_{T,gen.jet}^{ch} (GeV/#it{c})");
   hBkgM_trans->GetYaxis()->SetTitle("p_{T,rec.jet}^{ch} (GeV/#it{c})");

   TRandom3 *rpt = new TRandom3(0);
   const int sampling = 100000;
   double ptgen = 3, ptrec = 0.5;
   double delta = 0;

   for(int i=0; i<2*binsX; i++){
     for (int j=0; j<sampling; j++){
       delta = rpt->Gaus(0.664,7.6912);
       //cout << "deltapt: " << delta << endl;
       ptrec = ptgen + delta;
       if(ptrec>matrixBinMaxX) continue;
       hBkgM->Fill(ptrec,ptgen);
       hBkgM_trans->Fill(ptgen,ptrec);

      }
       //ptgen += 0.05;
       ptgen += 0.5;
    }
    hBkgM->Scale(1./(sampling));
    hBkgM->Scale(1.,"width");
    hBkgM->GetZaxis()->SetRangeUser(0,1);
    hBkgM->SetName("hBkgM");
    
     hBkgM_trans->Scale(1./(sampling));
    hBkgM_trans->Scale(1.,"width");
    hBkgM_trans->GetZaxis()->SetRangeUser(0,1);
    hBkgM_trans->SetName("hBkgM_trans");

    TH1F *hproj = (TH1F*)hBkgM->ProjectionY();
    TCanvas *c2 = new TCanvas;
    hproj->Draw();

    TCanvas *cBkgM = new TCanvas("cBkgM","cBkgM",1200,800);
    cBkgM->cd();
    cBkgM->SetLogz();
    hBkgM->Draw("colz2");
    
    cBkgM->SaveAs("BkgFlucMtx.png");
    
     TCanvas *cBkgM_trans = new TCanvas("cBkgM_trans","cBkgM_trans",1200,800);
    cBkgM_trans->cd();
    cBkgM_trans->SetLogz();
    hBkgM_trans->Draw("colz2");
    
    cBkgM_trans->SaveAs("BkgFlucMtx_trans.png");

  
TFile *fOut = new TFile("BkgFluctuationMtx.root","RECREATE");
hBkgM->Write();
hDeltaPt->Write();
hBkgM_trans->Write();

fOut->Close();
   
   	return;

}
