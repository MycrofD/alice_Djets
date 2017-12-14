#include "style.C"
#include <string>
#include <sstream>
#include <iostream>

const Int_t nJetBins = 9;
const double ptJetbins[nJetBins+1] = { 3,4,5,6,8,10,14,20,30,50 };

void DetRM(bool isPrompt = 1, int Dptmin = 0, int Dptmax = 100, int Rpar = 4, bool postfix = 0, TString listName = "FD" )
{
    char *outDir = "outResDetMatrixPythia";
    if(!isPrompt) outDir = "outResDetMatrixPythia_FD";
    TString datafile1 = "../outMC/AnalysisResults_fast_D0MCPythia_SMQcorr2.root";
 
    style();
    gStyle->SetOptStat(0000); //Mean and RMS shown

    double zmin = -2., zmax = 2;
    float jetmin = 0, jetmax = 60;
    double plotmin = 0, plotmax = 60;

    TPaveText *pv1 = new TPaveText(0.65,0.2,0.85,0.3,"brNDC");
    pv1->SetFillStyle(0);
    pv1->SetBorderSize(0);
    pv1->AddText(Form("%d < p_{T}^{D*} < %d GeV/c",Dptmin, Dptmax));

    TPaveText *pv2 = new TPaveText(0.15,0.8,0.25,0.9,"brNDC");
    pv2->SetFillStyle(0);
    pv2->SetBorderSize(0);
    pv2->AddText(Form("R=0.%d",Rpar));
    
    TString datafile = datafile1; 
	TFile *File = new TFile(datafile,"read");
	TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");

	TH1F *hMCpt;
	TH1F *hMCpt_reco;
    
    TH2F *hPtJet[5];
    TH1F *hPtG[5];
    TH1F *hPtR[5];
    
	TList *histList[5];
	THnSparseF *sparseMC[5];
	THnSparseF *sparsereco[5];
	
    TH2F *hPtJet2d;
    TH1F *hPtJetGen;
    TH1F *hPtJetRec;
    
	for(int i=0; i<3; i++){
	
		if(postfix) histList[i] =  (TList*)dir->Get(Form("histosD0MBN%d%sMCrec",i,listName.Data()));
        else {
			 if(isPrompt) histList[i] =  (TList*)dir->Get(Form("histosD0MBN%dMCrec",i));
			 else histList[i] =  (TList*)dir->Get(Form("histosD0MBN%dFDMCrec",i));
		}
		
		sparseMC[i] = (THnSparseF*)histList[i]->FindObject("ResponseMatrix"); 
	
        sparseMC[i]->GetAxis(2)->SetRangeUser(Dptmin,Dptmax); 
        sparseMC[i]->GetAxis(1)->SetRangeUser(jetmin,jetmax); 
        sparseMC[i]->GetAxis(6)->SetRangeUser(jetmin,jetmax); 
        hPtJet[i] = (TH2F*)sparseMC[i]->Projection(6,1,"E");
        hPtJet[i]->Sumw2();
        hPtJet[i]->SetName(Form("hPtJet_%d",i));
    
        hPtG[i] = (TH1F*)sparseMC[i]->Projection(6);
        hPtG[i]->SetName(Form("hPtG_%d",i));
        hPtR[i] = (TH1F*)sparseMC[i]->Projection(1);
        hPtR[i]->SetName(Form("hPtR_%d",i));
        hPtG[i]->Sumw2();
        hPtR[i]->Sumw2();
	
		if (!i){
			hPtJet2d = (TH2F*)hPtJet[0]->Clone("hPtJet2d");
			hPtJetGen = (TH1F*)hPtG[0]->Clone("hPtJetGen");
			hPtJetRec = (TH1F*)hPtR[0]->Clone("hPtJetRec");
            
		}
		else {
			hPtJet2d->Add(hPtJet[i]);
			hPtJetGen->Add(hPtG[i]);
			hPtJetRec->Add(hPtR[i]);
		}
	}
  
  
    hPtJet2d->SetTitle();
    hPtJet2d->SetName("hPtJet2d");
    hPtJet2d->GetXaxis()->SetTitle("p_{T,ch jet}^{rec.} (GeV/#it{c})");
    hPtJet2d->GetYaxis()->SetTitle("p_{T,ch jet}^{gen.} (GeV/#it{c})");
    
    hPtJetGen->SetName("hPtJetGen");
    hPtJetRec->SetName("hPtJetRec");
    hPtJetGen->SetLineColor(kBlue+2);
    hPtJetRec->SetLineColor(kRed+2);
    double recMean = hPtJetRec->GetMean();
    double genMean = hPtJetGen->GetMean();
    
    TCanvas *cjetPt = new TCanvas("cjetPt","cjetPt",800,600);
    cjetPt->SetLogy();
    hPtJetGen->Draw();
    hPtJetRec->Draw("same");
    cjetPt->SaveAs(Form("%s/pTdist_Dpt%d_%d.png",outDir, (int)Dptmin, (int)Dptmax));
   
    TCanvas *cjetPt2d = new TCanvas("cjetPt2d","cjetPt2d",800,600);
    cjetPt2d->SetLogz();
    hPtJet2d->Draw("colz");
    pv2->Draw("same");
    
    
    cjetPt2d->SaveAs(Form("%s/DetMatrix_Dpt%d_%d.png",outDir, (int)Dptmin, (int)Dptmax));
   
    TFile *ofile = new TFile(Form("%s/DetMatrix_Dpt%d_%d.root",outDir, (int)Dptmin, (int)Dptmax),"RECREATE");
    hPtJetGen->Write();
    hPtJetRec->Write();
    hPtJet2d->Write();
    ofile->Close();
   
   return;
   
    TH1D *proj[nJetBins];
    for(int i=0; i<nJetBins; i++){
            proj[i] = (TH1D*)hPtJet2d->ProjectionX(Form("proj_%d",i),hPtJet2d->GetYaxis()->FindBin(ptJetbins[i]), hPtJet2d->GetYaxis()->FindBin(ptJetbins[i+1]) -1);
            proj[i]->Scale(1./proj[i]->Integral());
            //proj[i-1]->SetMarkerStyle(20);
            proj[i]->SetMarkerColor(2);
            proj[i]->SetLineColor(2);
          
    }

    


  /*
    TH1F *hPtZGen = (TH1F*)sparse1_1->Projection(4);
    TH1F *hPtZReco = (TH1F*)sparse1_1->Projection(0);
     
    hPtZGen->SetName("hPtZGen");
    hPtZReco->SetName("hPtZReco");
    
    hPtZGen->SetLineColor(kBlue+2);
    hPtZReco->SetLineColor(kRed+2);
    hPtZGen->GetXaxis()->SetRangeUser(0.1,1.1);
    hPtZReco->GetXaxis()->SetRangeUser(0.1,1.1);
    hPtZGen->SetTitle("z gen");
    hPtZReco->SetTitle("z rec");
    
    
     TCanvas *cZ = new TCanvas("cZ","cZ",800,500);
     cZ->Divide(2,1);
    //cZ->SetLogy();
    cZ->cd(1);
    hPtZGen->Draw();
    cZ->cd(2);
    //hPtJetRec->Draw("same");
    hPtZReco->Draw();
    //pv1->Draw("same");
   // pv2->Draw("same");
     cZ->SaveAs(Form("%s/Zdist_Dpt%d_%d.png",outDir, (int)Dptmin, (int)Dptmax)); */
    
    //TH1F *hPtJetGen_tmp = (TH1F*)hPtJet2d->ProjectionY();
    //TH1F *hPtJetGen = (TH1F*)hPtJetGen_tmp->Clone("hPtJetGen");
    //hPtJetGen->SetLineColor(kBlue+2);
    //TH1F *hPtJetRec_tmp = (TH1F*)hPtJet2d->ProjectionX();
    //TH1F *hPtJetRec = (TH1F*)hPtJetRec_tmp->Clone("hPtJetRec");



/*   
      gStyle->SetOptStat(000);
      
       
TLegend *leg = new TLegend(0.5,0.6,0.85,0.7);
leg->SetBorderSize(0);
leg->AddEntry(proj[0],"prompt D*", "l");
//leg->AddEntry(projB[0],"non-prompt D*", "l");
TCanvas *cProj = new TCanvas("cProj","cProj",1600,1200);
cProj->Divide(3,4);





TFile *outP = new TFile(Form("%s/DetMatrixProjections_Dpt%d_%d.root",outDir, (int)Dptmin, (int)Dptmax),"RECREATE");

for(int i=0; i<nJetBins; i++){
cProj->cd(i+1);
gPad->SetLogy();
proj[i]->GetXaxis()->SetRangeUser(0,35);
proj[i]->Draw("h");
proj[i]->Write();
//projB[i]->Draw("h same");
 TPaveText *pv3 = new TPaveText(0.55,0.8,0.9,0.9,"brNDC");
    pv3->SetFillStyle(0);
    pv3->SetBorderSize(0);
    pv3->AddText(Form("%.0f < p_{T,ch jet}^{gen} < %.0f GeV/c",ptJetbins[i], ptJetbins[i+1]));
pv3->Draw("same");
if(i==1) leg->Draw("same");
}
 cProj->cd(12);  
 
   outP->Close();
   
   cProj->SaveAs(Form("%s/DetMatrixProjections_Dpt%d_%d.png",outDir, (int)Dptmin, (int)Dptmax));
   cProj->SaveAs(Form("%s/DetMatrixProjections_Dpt%d_%d.pdf",outDir, (int)Dptmin, (int)Dptmax));
*/

/*

TString outDir1 = "outResDetMatrix_Dremoved_Hijing";
TString outDir2 = "outResDetMatrix_Dremoved_Hijing_nonPrompt";
TFile *fproj1 = new TFile(Form("%s/DetMatrixProjections_Dpt%d_%d.root",outDir1.Data(), (int)Dptmin, (int)Dptmax),"READ");
TFile *fproj2 = new TFile(Form("%s/DetMatrixProjections_Dpt%d_%d.root",outDir2.Data(), (int)Dptmin, (int)Dptmax),"READ");

TCanvas *cProj2 = new TCanvas("cProj2","cProj2",1600,1200);
cProj2->Divide(3,4);

 TH1D *proj1[nJetBins];
 TH1D *proj2[nJetBins];

for(int i=0; i<nJetBins; i++){
cProj2->cd(i+1);
gPad->SetLogy();
proj1[i]=(TH1D*)fproj1->Get(Form("proj_%d",i));
proj2[i]=(TH1D*)fproj2->Get(Form("proj_%d",i));
proj1[i]->SetLineColor(2);
proj2[i]->SetLineColor(4);
proj1[i]->GetXaxis()->SetRangeUser(0,35);
proj2[i]->GetXaxis()->SetRangeUser(0,35);
proj1[i]->Draw("h");
proj2[i]->Draw("h same");

 TPaveText *pv3 = new TPaveText(0.55,0.8,0.9,0.9,"brNDC");
    pv3->SetFillStyle(0);
    pv3->SetBorderSize(0);
    pv3->AddText(Form("%.0f < p_{T,ch jet}^{gen} < %.0f GeV/c",ptJetbins[i], ptJetbins[i+1]));
pv3->Draw("same");
if(i==1) leg->Draw("same");
}

cProj2->cd(12);  

 cProj2->SaveAs(Form("%s/DetMatrixProjectionsComparison_Dpt%d_%d.pdf",outDir, (int)Dptmin, (int)Dptmax));

*/

}
