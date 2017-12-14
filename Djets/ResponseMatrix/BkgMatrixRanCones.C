//#include "style.C"

void BkgMatrixRanCones(char *outDir = "outMatrix", int Rpar = 4 )
//void BkgMatrixRanCones(char *outDir = "outMatrix_realptshape", int Rpar = 4 )
{


    char *datafile = "fOutDjetProperties.root";
      
	gStyle->SetOptStat(0000); //Mean and RMS shown
	style();

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
    
    //TString postfix = "Incljet5Excl_noBkg";
    TString postfix = "Djet10Excl_noBkg";
        
    // -------- Data -------
    char dataFile1[200];
	sprintf(dataFile1,datafile);
	TFile *File = new TFile(dataFile1,"read");
	TH1F *hDeltaPt = (TH1F*)File->Get("hDeltaPt_ptleadbin5_exlcuding");
	//TH1F *hDeltaPt = (TH1F*)File->Get("hDeltaPt_ptleadbin10NoBkg_exlcuding");
   
    hDeltaPt->SetTitle();
    //hDeltaPt->Rebin(5);
    hDeltaPt->Sumw2();
    hDeltaPt->GetXaxis()->SetRangeUser(-10,60);
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
    

    const int nb=400;
    double binmin = -5;
    double binmax = 35;
    TH1F *hDeltaPt2 = new TH1F("hDeltaPt2","hDeltaPt2",nb,binmin,binmax);
   
    double bin = binmin;
    for(int k=1; k<nb; k++){
            hDeltaPt2->SetBinContent(k,hDeltaPt->GetBinContent(hDeltaPt->GetXaxis()->FindBin(bin)));
            bin+=0.1;
    }
   
   hDeltaPt2->SetLineColor(kMagenta);
   hDeltaPt2->Draw();
   
       

   
   /* double ptgen = 0, ptrec = 0;
    double delta = 0;
    for(int i=1; i<1000; i++){
        for (int j=0; j<10000; j++){
            delta = hDeltaPt2->GetRandom();
            ptrec = ptgen + delta;
            hBkgM->Fill(ptrec,ptgen);
        }
        ptgen += 0.1;     
    }*/
    
   /*
    char *MCfile = "/home/basia/Work/alice/analysis/pPb_run2/outMC/AnalysisResults_MCRMEff_set5.root";
    char mcfile[100];
	sprintf(mcfile,MCfile);
	TFile *File1 = new TFile(mcfile,"read");
	TDirectoryFile* dir1=(TDirectoryFile*)File1->Get("DmesonsForJetCorrelations");
    
	TList *histList1 =  (TList*)dir1->Get("histosDStarMBN0MCrec");
	THnSparseF *sparse1 = (THnSparseF*)histList1->FindObject("ResponseMatrix");
    //sparse1->GetAxis(2)->SetRangeUser(3,36); 
    
    TH1F *hMCPt = (TH1F*)sparse1->Projection(5);
    
    TCanvas *cc = new TCanvas;
    hMCPt->Draw();
    
   */
   
   
    TRandom3 *rpt = new TRandom3(0);
 
    const int sampling = 100000;
    double ptgen = 0.5, ptrec = 0.5;
    double delta = 0;
    //for(int i=0; i<500; i++){
    //for(int i=0; i<2500; i++){
    for(int i=0; i<2*binsX; i++){
        for (int j=0; j<sampling; j++){
            //ptgen = hMCPt->GetRandom();
            delta = hDeltaPt2->GetRandom();
           // delta = hDeltaPt->GetRandom();
         // cout << '\n' << delta;
          
          //  delta = hDeltaPt->GetRandom();
          // delta = rpt->Gaus(0.5,0.29);
         //  delta = rpt->Gaus(0.94,0.15);
           //delta = rpt->Gaus(0.5,0.29);
           //delta = fP->GetRandom();
           //delta = rpt->Landau(fP->GetParameter(0),fP->GetParameter(1));
            ptrec = ptgen + delta;
            if(ptrec<0 || ptrec>matrixBinMaxX) continue;
            hBkgM->Fill(ptrec,ptgen);
        }
      
       //ptgen += 0.05;     
       ptgen += 0.5;  
          
    
    }
      hBkgM->Scale(1./(sampling));
      hBkgM->Scale(1.,"width");
      hBkgM->GetZaxis()->SetRangeUser(0,1);
      
      TH1F *hproj = (TH1F*)hBkgM->ProjectionY();
      TCanvas *c2 = new TCanvas;
      hproj->Draw();
    

    hBkgM->SetName("hBkgM");

    
    TCanvas *cBkgM = new TCanvas("cBkgM","cBkgM",800,600);
    cBkgM->cd();
    cBkgM->SetLogz();
    hBkgM->Draw("colz");
  
    
    pv1->Draw("same");
    pv2->Draw("same");
    
       //cout << "-- bin widths: " << hBkgM->GetXaxis()->GetBinWidth(1) << "\t" << hBkgM->GetXaxis()->GetBinWidth(1) << endl;
       //cout << "--entries: " << hBkgM << hBkgM->GetEntries() << endl;

   
    
    // ---- SAVE ----
    
   
    cDPt->SaveAs(Form("%s/DeltaPt_%s.png",outDir,postfix.Data()));
    cBkgM->SaveAs(Form("%s/BkgMatrix_%s.png",outDir,postfix.Data()));
    c2->SaveAs(Form("%s/BkgMatrixGenProj_%s.png",outDir,postfix.Data()));
    
    TFile *ofile = new TFile(Form("%s/RandCones_BkgM_%s.root",outDir,postfix.Data()),"RECREATE");

    hDeltaPt->Write();

    hBkgM->Write();
    hproj->Write();
    ofile->Close();
    
}


