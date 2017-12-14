#include <string>
#include <sstream>
#include <iostream>

using namespace std;
	
const int ptbinsDN = 12;
double ptDbins[ptbinsDN+1] = { 1,2,3,4,5,6,7,8,10,12,16,24,36 };
float ptmin = 1, ptmax = 36;

       
void DjetEfficiency(bool isPrompt = 1, TString outDir = "SQMCorrcuts", TString effFile = "../outMC/AnalysisResults_fast_D0MCHijing_SMQcorr2.root", float jetptmin = 2, float jetptmax = 50, bool postfix = 0, TString listName = "FD", bool isSys = 0, int Rpar = 4 )
{


 gStyle->SetOptStat(0000); //Mean and RMS shown
    
stringstream sst;
sst.clear(); sst.str("");
	
	
	 // get analysis output file
    TString datafile = effFile; 

	TFile *File = new TFile(datafile,"read");
	TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");

	TH1F *hMCpt;
	TH1F *hMCpt_reco;
    
	TH1F *hhMC[5];
	TH1F * hhreco[5];
	TList *histList[5];
	THnSparseF *sparseMC[5];
	THnSparseF *sparsereco[5];
	
	for(int i=0; i<3; i++){
	
		if(postfix) histList[i] =  (TList*)dir->Get(Form("histosD0MBN%d%sMCrec",i,listName.Data()));
        else {
			 if(isPrompt) histList[i] =  (TList*)dir->Get(Form("histosD0MBN%dMCrec",i));
			 else histList[i] =  (TList*)dir->Get(Form("histosD0MBN%dFDMCrec",i));
		}
		
		sparseMC[i] = (THnSparseF*)histList[i]->FindObject("ResponseMatrix"); 
	
		sparseMC[i]->GetAxis(6)->SetRangeUser(jetptmin,jetptmax); // jet pT gen
		//sparseMC[i]->GetAxis(6)->SetRangeUser(3,36); // jet pT gen
		
		hhMC[i] = (TH1F*)sparseMC[i]->Projection(7); // Dpt gen
		hhMC[i]->SetName(Form("hhMC_%d",i));
		
		sparsereco[i] = (THnSparseF*)histList[i]->FindObject("ResponseMatrix"); 
		//sparsereco[i]->GetAxis(5)->SetRangeUser(0,36);
		
		sparsereco[i]->GetAxis(6)->SetRangeUser(jetptmin,jetptmax); // jet pT gen
		sparsereco[i]->GetAxis(1)->SetRangeUser(0,jetptmax); // jet pT gen
		
		hhreco[i] = (TH1F*)sparsereco[i]->Projection(7); // Dpt gen
		hhreco[i]->SetName(Form("hhreco_%d",i));
		
	
		if (!i){
			hMCpt = (TH1F*)hhMC[0]->Clone("hMCpt");
			hMCpt_reco = (TH1F*)hhreco[0]->Clone("hMCpt_reco");
		}
		else {
			hMCpt->Add(hhMC[i]);
			hMCpt_reco->Add(hhreco[i]);
		}
	}
	hMCpt->Sumw2();
	hMCpt_reco->Sumw2();
	
	hMCpt->SetLineColor(4);
	hMCpt->SetMarkerColor(4);
	hMCpt->SetMarkerStyle(20);

	hMCpt_reco->SetLineColor(2);
	hMCpt_reco->SetMarkerColor(2);
	hMCpt_reco->SetMarkerStyle(20);
	
	
	TH1D * hEff = (TH1D*)hMCpt_reco->Clone("hEff");
	hEff -> Divide(hMCpt_reco,hMCpt,1,1,"b");

	hEff->GetXaxis()->SetRangeUser(ptmin,ptmax);
	
	hMCpt->SetMinimum(10);
	TCanvas *cPt = new TCanvas();
	hMCpt->Draw();
	hMCpt_reco->Draw("same");
	
	if(isPrompt) setHistoDetails(hEff,2,20,1.2);
	else setHistoDetails(hEff,4,25,1.2);

	TCanvas *cEff = new TCanvas();
	hEff->Draw("ep");

	TH1D *hpt_mc_reb = (TH1D*)hMCpt->Rebin(ptbinsDN,"hpt_mc_reb",ptDbins);
	TH1D *hpt_reco_reb = (TH1D*)hMCpt_reco->Rebin(ptbinsDN,"hpt_reco_reb",ptDbins);
	
	TH1D * hEff_reb = (TH1D*)hpt_reco_reb->Clone("hEff_reb");
	hEff_reb -> Divide(hpt_reco_reb,hpt_mc_reb,1,1,"b");
	hEff_reb->GetXaxis()->SetRangeUser(ptmin,ptmax);
	hEff_reb->SetTitle("|#eta_{jet}|<0.5");
	hEff_reb->GetXaxis()->SetTitle("p_{T}^{D0} (GeV/c)");
	hEff_reb->GetYaxis()->SetTitle("acc #times eff");
	
	if(isPrompt) setHistoDetails(hEff_reb,2,20,1.2);
	else setHistoDetails(hEff_reb,4,25,1.2);
	
	TCanvas *cEffReb = new TCanvas();
	hEff_reb->Draw("ep");
	
	TString out = Form("%s/DjetEffRaw_%s_jetpt%d_%d", outDir.Data(), isPrompt ? "prompt" : "nonPrompt", (int)jetptmin,(int)jetptmax );
	SaveCanvas(cEff,out);
	
	out = Form("%s/DjetEff_%s_jetpt%d_%d", outDir.Data(), isPrompt ? "prompt" : "nonPrompt", (int)jetptmin,(int)jetptmax );
	SaveCanvas(cEffReb,out);
	
	TFile *outFile = new TFile(Form("%s/DjetEff_%s_jetpt%d_%d.root", outDir.Data(), isPrompt ? "prompt" : "nonPrompt", (int)jetptmin, (int)jetptmax ),"RECREATE");
	hEff->Write();
	hEff_reb->Write();
	hMCpt->Write();
	hMCpt_reco->Write();
	hpt_mc_reb->Write();
	hpt_reco_reb->Write();
	outFile->Close();	

return; 

}

void setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Size_t size = 0.9, Width_t width=2){
    
    h->SetMarkerStyle(Mstyle);
    h->SetMarkerColor(color);
    h->SetMarkerSize(size);
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetTitle(0);
    h->GetXaxis()->SetTitle("p_{T,D^{0}}(GeV/c)");
    
    return;

}

void SaveCanvas(TCanvas *c, TString name = "tmp"){
    
    c->SaveAs(Form("%s.png",name.Data()));
    c->SaveAs(Form("%s.pdf",name.Data()));
   
}
