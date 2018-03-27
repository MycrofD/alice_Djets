#include <string>
#include <sstream>
#include <iostream>

const int nJetBins = 16;
double ptJetbins[nJetBins] = { 3,4,6,8,10,12,14,16,18,20,22,24,26,30,40,50 };

void JetProperties(int Dptmin = 0, int Dptmax = 100, int Rpar = 4 )
{
    
    char *outDir = "plots";
    char *datafile1 = "/home/basia/Work/alice/analysis/pPb_run2/outData/jetProperties/AnalysisResults_RandomConesUE_Djets_1308.root";
    //char *datafile1 = "/home/basia/Work/alice/analysis/pPb_run2/outData/jetProperties/AnalysisResults_noLeadingJetRejection.root";
  
     gStyle->SetOptStat(0000); //Mean and RMS shown
    
stringstream sst;
sst.clear(); sst.str("");

    double zmin = -2., zmax = 2;
    float jetmin = 0, jetmax = 100;
    double plotmin = 0, plotmax = 50;


    TPaveText *pv2 = new TPaveText(0.15,0.8,0.25,0.9,"brNDC");
    pv2->SetFillStyle(0);
    pv2->SetBorderSize(0);
    pv2->AddText(Form("R=0.%d",Rpar));
    
    char dataFile1[200];
	sprintf(dataFile1,datafile1);
	TFile *File1 = new TFile(dataFile1,"read");
	TDirectoryFile* dir1=(TDirectoryFile*)File1->Get("DmesonsForJetCorrelations");
    TList *histList1 =  (TList*)dir1->Get("histosDStarMBN0");

 //======================= Rho ==================================
     
    TH2F *hRhoMult = (TH2F*)histList1->FindObject("fhRhoMult");
    hRhoMult->SetTitle();
    
    TCanvas *cRhoM = new TCanvas("cRhoM","cRhoM",800,600);
    hRhoMult->Draw("colz");
    
    cRhoM->SaveAs("plotsDjets/RhoVsMult.png");
    
    TProfile *hRhoMultProf = (TProfile*)hRhoMult->ProfileX("hRhoMultProf");
    hRhoMultProf->SetTitle();
    hRhoMultProf->GetYaxis()->SetRangeUser(0,10);
    hRhoMultProf->GetYaxis()->SetTitle("<#rho(GeV/c*rad^{-1})>");

    hRhoMultProf->BuildOptions(0,0,"i");

    TH2F *hRhoLeadPt = (TH2F*)histList1->FindObject("fhRhoLeadPt");
    hRhoLeadPt->SetTitle();
    hRhoLeadPt->GetXaxis()->SetRangeUser(0,80);
    
    TCanvas *cRhoLeadPt = new TCanvas("cRhoLeadPt","cRhoLeadPt",800,600);
    hRhoLeadPt->Draw("colz");
    
    cRhoLeadPt->SaveAs("plotsDjets/RhoVsLeadPt.png");
    
    TProfile *hRhoLeadPtProf = (TProfile*)hRhoLeadPt->ProfileX("hRhoLeadPtProf");
    hRhoLeadPtProf->SetTitle();
    hRhoLeadPtProf->GetYaxis()->SetRangeUser(0,10);
    hRhoLeadPtProf->GetYaxis()->SetTitle("<#rho(GeV/c*rad^{-1})>");
    hRhoLeadPtProf->BuildOptions(0,0,"i");
    
    TH1F *hDLeading = (TH1F*)histList1->FindObject("fhDleadStat");
    hDLeading->SetTitle();
    hDLeading->SetMinimum(0);
    
    TCanvas *cLeadD = new TCanvas("cLeadD","cLeadD",800,600);
    hDLeading->Draw();
    
    cLeadD->SaveAs("plotsDjets/LeadingDJet.png");
    
    TH2F *hDLeading2D = (TH2F*)histList1->FindObject("fhDleadStatJetPt");
    hDLeading2D->SetTitle();
    
    TCanvas *cLeadD2D = new TCanvas("cLeadD2D","cLeadD2D",800,600);
    hDLeading2D->Draw("colz");
    
    cLeadD2D->SaveAs("plotsDjets/LeadingDJetVsPt.png");
  

    TH3F *hRhoLeadPtMult_tmp = (TH3F*)histList1->FindObject("fhRhoLeadPtMult");
    TH3F *hRhoLeadPtMult = (TH3F*)hRhoLeadPtMult_tmp->Clone("hRhoLeadPtMult");
    hRhoLeadPtMult->SetTitle();
    hRhoLeadPtMult->GetXaxis()->SetRangeUser(-20,180);
    hRhoLeadPtMult->GetXaxis()->SetTitle("p_{T,lead}^jet-#rho A (GeV/c)");
   
    TH3F *hRhoLeadPtMult1_tmp = (TH3F*)histList1->FindObject("fhRhoLeadPtMult");
    TH3F *hRhoLeadPtMult1 = (TH3F*)hRhoLeadPtMult1_tmp->Clone("hRhoLeadPtMult1");
    hRhoLeadPtMult1->SetName("hRhoLeadPtMult1");
    hRhoLeadPtMult1->SetTitle();
    hRhoLeadPtMult1->GetXaxis()->SetRangeUser(5,180);
    
    TH2F *hRhoMult_leadpt5 = (TH2F*)hRhoLeadPtMult1->Project3D("zy");
    hRhoMult_leadpt5->SetName("hRhoMult_leadpt5");
    hRhoMult_leadpt5->SetTitle();
    
    TCanvas *cRhoM_pt5 = new TCanvas("cRhoM_pt5","cRhoM_pt5",800,600);
    hRhoMult_leadpt5->Draw("colz");
    
    cRhoM_pt5->SaveAs("plotsDjets/RhoVsMult_pt5.png");
    
    TProfile *hRhoMultProf_leadpt5 = (TProfile*)hRhoMult_leadpt5->ProfileX("hRhoMultProf_leadpt5");
    hRhoMultProf_leadpt5->SetTitle();
    hRhoMultProf_leadpt5->GetYaxis()->SetRangeUser(0,10);
    hRhoMultProf_leadpt5->GetYaxis()->SetTitle("<#rho(GeV/c*rad^{-1})>");
    hRhoMultProf_leadpt5->BuildOptions(0,0,"i");
    
    
    TH3F *hRhoLeadPtMult2_tmp = (TH3F*)histList1->FindObject("fhRhoLeadPtMult");
    TH3F *hRhoLeadPtMult2 = (TH3F*)hRhoLeadPtMult2_tmp->Clone("hRhoLeadPtMult2");
    hRhoLeadPtMult2->SetName("hRhoLeadPtMult2");
    hRhoLeadPtMult2->SetTitle();
    hRhoLeadPtMult2->GetXaxis()->SetRangeUser(10,180);
    
    TH2F *hRhoMult_leadpt10 = (TH2F*)hRhoLeadPtMult2->Project3D("zy");
    hRhoMult_leadpt10->SetName("hRhoMult_leadpt10");
    hRhoMult_leadpt10->SetTitle();
    
    TCanvas *cRhoM_pt10 = new TCanvas("cRhoM_pt10","cRhoM_pt10",800,600);
    hRhoMult_leadpt10->Draw("colz");
    
    cRhoM_pt10->SaveAs("plotsDjets/RhoVsMult_pt10.png");
    
    TProfile *hRhoMultProf_leadpt10 = (TProfile*)hRhoMult_leadpt10->ProfileX("hRhoMultProf_leadpt10");
    hRhoMultProf_leadpt10->SetTitle();
    hRhoMultProf_leadpt10->GetYaxis()->SetRangeUser(0,10);
    hRhoMultProf_leadpt10->GetYaxis()->SetTitle("<#rho(GeV/c*rad^{-1})>");
    hRhoMultProf_leadpt10->BuildOptions(0,0,"i");
    
    //======================= pt bins ==================================
    
    TLegend *leg1 = new TLegend(0.5,0.5,0.88,0.88,"D-jet, p_{T,lead}:");
    leg1->SetBorderSize(0);

    
    const int ptbins = 5;
    int ptleadDown[] = { 0,2,5,10,20 };
    int ptleadUp[] = { 2,5,10,20,80 };
    int colors[] = {1,2,4,6,8,10};
    
    TH3F *hRhoLeadPtMultBinPt[ptbins];
    TH2F *hRhoMultBinsPt[ptbins];
    TProfile *hRhoMultProfBinsPt[ptbins];
    
    for(int i=0; i<ptbins; i++){

    hRhoLeadPtMultBinPt[i] = (TH3F*)hRhoLeadPtMult->Clone(Form("hRhoLeadPtMultBinPt_%d",i));
    hRhoLeadPtMultBinPt[i]->GetXaxis()->SetRangeUser(ptleadDown[i],ptleadUp[i]);
    hRhoMultBinsPt[i] = (TH2F*)hRhoLeadPtMultBinPt[i]->Project3D("zy");
    hRhoMultBinsPt[i]->SetTitle(Form("%d<p_{T}^{lead}<%d",ptleadDown[i],ptleadUp[i]));

    hRhoMultProfBinsPt[i] = (TProfile*)hRhoMultBinsPt[i]->ProfileX(Form("hRhoMultProfBinsPt_%d",i));
    //hRhoMultProf_1->SetTitle(Form("%d<p_{T}^{lead}<%d",ptleadDown[0],ptleadUp[0]));
    hRhoMultProfBinsPt[i]->SetTitle();
    hRhoMultProfBinsPt[i]->SetLineColor(colors[i]);
    hRhoMultProfBinsPt[i]->GetYaxis()->SetRangeUser(0,6);
    hRhoMultProfBinsPt[i]->GetYaxis()->SetTitle("<#rho(GeV/c*rad^{-1})>");
    hRhoMultProfBinsPt[i]->BuildOptions(0,0,"i");
    
    leg1->AddEntry(hRhoMultProfBinsPt[i],Form("%d-%d GeV/c",ptleadDown[i],ptleadUp[i]),"l");
    
    
    }
    
         TH1F* hRho_ptleadbin5 = (TH1F*)hRhoLeadPtMult->ProjectionZ("hRho_ptleadbin5",hRhoLeadPtMult->GetXaxis()->FindBin(5),hRhoLeadPtMult->GetXaxis()->FindBin(180)-1,hRhoLeadPtMult->GetYaxis()->FindBin(0),hRhoLeadPtMult->GetYaxis()->FindBin(100));
    TH1F* hRho_ptleadbin10 = (TH1F*)hRhoLeadPtMult->ProjectionZ("hRho_ptleadbin10",hRhoLeadPtMult->GetXaxis()->FindBin(10),hRhoLeadPtMult->GetXaxis()->FindBin(180)-1,hRhoLeadPtMult->GetYaxis()->FindBin(0),hRhoLeadPtMult->GetYaxis()->FindBin(100));
    
    
     TH1F* hRho_ptleadbin5_cent0_10 = (TH1F*)hRhoLeadPtMult->ProjectionZ("hRho_ptleadbin5_cent0_10",hRhoLeadPtMult->GetXaxis()->FindBin(5),hRhoLeadPtMult->GetXaxis()->FindBin(180)-1,hRhoLeadPtMult->GetYaxis()->FindBin(0),hRhoLeadPtMult->GetYaxis()->FindBin(10)-1);
     TH1F* hRho_ptleadbin5_cent10_20 = (TH1F*)hRhoLeadPtMult->ProjectionZ("hRho_ptleadbin5_cent10_20",hRhoLeadPtMult->GetXaxis()->FindBin(5),hRhoLeadPtMult->GetXaxis()->FindBin(180)-1,hRhoLeadPtMult->GetYaxis()->FindBin(10),hRhoLeadPtMult->GetYaxis()->FindBin(20)-1);
     TH1F* hRho_ptleadbin5_cent20_40 = (TH1F*)hRhoLeadPtMult->ProjectionZ("hRho_ptleadbin5_cent20_40",hRhoLeadPtMult->GetXaxis()->FindBin(5),hRhoLeadPtMult->GetXaxis()->FindBin(180)-1,hRhoLeadPtMult->GetYaxis()->FindBin(20),hRhoLeadPtMult->GetYaxis()->FindBin(40)-1);
     TH1F* hRho_ptleadbin5_cent40_60 = (TH1F*)hRhoLeadPtMult->ProjectionZ("hRho_ptleadbin5_cent40_60",hRhoLeadPtMult->GetXaxis()->FindBin(5),hRhoLeadPtMult->GetXaxis()->FindBin(180)-1,hRhoLeadPtMult->GetYaxis()->FindBin(40),hRhoLeadPtMult->GetYaxis()->FindBin(60)-1);
     TH1F* hRho_ptleadbin5_cent60_100 = (TH1F*)hRhoLeadPtMult->ProjectionZ("hRho_ptleadbin5_cent60_100",hRhoLeadPtMult->GetXaxis()->FindBin(5),hRhoLeadPtMult->GetXaxis()->FindBin(180)-1,hRhoLeadPtMult->GetYaxis()->FindBin(60),hRhoLeadPtMult->GetYaxis()->FindBin(100));
    
    
     //======================= pt lead cut ==================================
    
    TLegend *leg11 = new TLegend(0.5,0.5,0.88,0.88,"D-jet, p_{T,lead}:");
    leg11->SetBorderSize(0);

    TH3F *hRhoLeadPtMultBinPt2[ptbins];
    TH2F *hRhoMultBinsPt2[ptbins];
    TProfile *hRhoMultProfBinsPt2[ptbins];
    
    for(int i=0; i<ptbins; i++){

    hRhoLeadPtMultBinPt2[i] = (TH3F*)hRhoLeadPtMult->Clone(Form("hRhoLeadPtMultBinPt2_%d",i));
    hRhoLeadPtMultBinPt2[i]->GetXaxis()->SetRangeUser(ptleadDown[i],180);
    hRhoMultBinsPt2[i] = (TH2F*)hRhoLeadPtMultBinPt2[i]->Project3D("zy");
    hRhoMultBinsPt2[i]->SetTitle(Form("%d<p_{T}^{lead}<%d",ptleadDown[i],180));

    hRhoMultProfBinsPt2[i] = (TProfile*)hRhoMultBinsPt2[i]->ProfileX(Form("hRhoMultProfBinsPt2_%d",i));
    //hRhoMultProf_1->SetTitle(Form("%d<p_{T}^{lead}<%d",ptleadDown[0],ptleadUp[0]));
    hRhoMultProfBinsPt2[i]->SetTitle();
    hRhoMultProfBinsPt2[i]->SetLineColor(colors[i]);
    hRhoMultProfBinsPt2[i]->GetYaxis()->SetRangeUser(0,6);
    hRhoMultProfBinsPt2[i]->GetYaxis()->SetTitle("<#rho(GeV/c*rad^{-1})>");
    hRhoMultProfBinsPt2[i]->BuildOptions(0,0,"i");
    
    leg11->AddEntry(hRhoMultProfBinsPt2[i],Form(">%d GeV/c",ptleadDown[i]),"l");
    
    
    }
    
    
      
    //======================= cent bins ==================================
    
     TLegend *leg2 = new TLegend(0.6,0.6,0.95,0.88,"D-jet, centrality");
    leg2->SetBorderSize(0);
    
    const int centBins = 5;
    int centDown[] = { 0,10,20,40,60 };
    int centUp[] = { 10,20,40,60,100 };
    int colors[] = {1,2,4,6,8,10};
    
    TH3F *hRhoLeadPtMultBin[centBins];
    TH2F *hRhoMultBins[centBins];
    TProfile *hRhoMultProfBins[centBins];
    
    for(int i=0; i<centBins; i++){

    hRhoLeadPtMultBin[i] = (TH3F*)hRhoLeadPtMult->Clone(Form("hRhoLeadPtMultBin_%d",i));
    hRhoLeadPtMultBin[i]->GetYaxis()->SetRangeUser(centDown[i],centUp[i]);
    hRhoMultBins[i] = (TH2F*)hRhoLeadPtMultBin[i]->Project3D("zx");
    hRhoMultBins[i]->SetTitle(Form("%d<cent<%d",centDown[i],centUp[i]));

    hRhoMultProfBins[i] = (TProfile*)hRhoMultBins[i]->ProfileX(Form("hRhoMultProfBins_%d",i));
    //hRhoMultProf_1->SetTitle(Form("%d<p_{T}^{lead}<%d",ptleadDown[0],ptleadUp[0]));
    hRhoMultProfBins[i]->SetTitle();
    hRhoMultProfBins[i]->SetLineColor(colors[i]);
    hRhoMultProfBins[i]->GetYaxis()->SetRangeUser(0,6);
    hRhoMultProfBins[i]->GetXaxis()->SetRangeUser(0,100);
    hRhoMultProfBins[i]->GetYaxis()->SetTitle("<#rho(GeV/c*rad^{-1})>");
    hRhoMultProfBins[i]->BuildOptions(0,0,"i");
    
     leg2->AddEntry(hRhoMultProfBins[i],Form("%d-%d",centDown[i],centUp[i]),"l");
    
    }
 
 
 //======================= delta pt ==================================
 
 
    // vs centrality
 
    TH2F *hDeltaPtCent = (TH2F*)histList1->FindObject("fDeltaPTCent");
    hDeltaPtCent->SetTitle();
    hDeltaPtCent->GetYaxis()->SetRangeUser(-10,100);
    
    TCanvas *cDetlaPtCent = new TCanvas("cDetlaPtCent","cDetlaPtCent",800,600);
    hDeltaPtCent->Draw("colz");
    
    cDetlaPtCent->SaveAs("plotsDjets/DeltaPtCent.png");
    
    TProfile *hDeltaPtCentProf = (TProfile*)hDeltaPtCent->ProfileX("hDeltaPtCentProf");
    hDeltaPtCentProf->SetTitle();
    //hDeltaPtCentProf->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtCentProf->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtCentProf->BuildOptions(0,0,"i");
 
    TH1F *hDeltaPt_c = (TH1F*)hDeltaPtCent->ProjectionY("hDeltaPt_c");
    hDeltaPt_c->SetTitle();
    hDeltaPt_c->Scale(1./hDeltaPt_c->Integral());
     
    // vs centrality excl. lead. jet 
    
    TH2F *hDeltaPtCent_excleading = (TH2F*)histList1->FindObject("fDeltaPTCent_excl_lead");
    hDeltaPtCent_excleading->SetTitle();
    hDeltaPtCent_excleading->GetYaxis()->SetRangeUser(-10,100);
    
    TCanvas *cDetlaPtCent_excleading = new TCanvas("cDetlaPtCent_excleading","cDetlaPtCent_excleading",800,600);
    hDeltaPtCent_excleading->Draw("colz");
    
    cDetlaPtCent_excleading->SaveAs("plotsDjets/DeltaPtCent_excleading.png");
    
    TProfile *hDeltaPtCentProf_excleading = (TProfile*)hDeltaPtCent_excleading->ProfileX("hDeltaPtCentProf_excleading");
    hDeltaPtCentProf_excleading->SetTitle();
    //hDeltaPtCentProf->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtCentProf_excleading->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtCentProf_excleading->BuildOptions(0,0,"i");
 
    TH1F *hDeltaPt_c_excleading = (TH1F*)hDeltaPtCent_excleading->ProjectionY("hDeltaPt_c_excleading");
    hDeltaPt_c_excleading->SetTitle();
    hDeltaPt_c_excleading->Scale(1./hDeltaPt_c_excleading->Integral());
    

    // vs lead. jet pT

    TH2F *hDeltaPtLeadPt = (TH2F*)histList1->FindObject("fDeltaPTLeadPt");
    hDeltaPtLeadPt->SetTitle();
    hDeltaPtLeadPt->GetYaxis()->SetRangeUser(-10,100);
    
    TCanvas *cDetlaPtLeadPt = new TCanvas("cDetlaPtLeadPt","cDetlaPtLeadPt",800,600);
    hDeltaPtLeadPt->Draw("colz");
    cDetlaPtLeadPt->SaveAs("plotsDjets/DetlaPtLeadPt.png");
    
    TProfile *hDeltaPtLeadPtProf = (TProfile*)hDeltaPtLeadPt->ProfileX("hDeltaPtLeadPtProf");
    hDeltaPtLeadPtProf->SetTitle();
    //hDeltaPtCentProf->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtLeadPtProf->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtLeadPtProf->BuildOptions(0,0,"i");
    
    TH1F *hDeltaPt_ptlead = (TH1F*)hDeltaPtLeadPt->ProjectionY("hDeltaPt_ptlead");
    hDeltaPt_ptlead->SetTitle();
    
    TLegend *legDeltaPt = new TLegend(0.5,0.5,0.88,0.88,"D-jet, p_{T,lead}:");
    legDeltaPt->SetBorderSize(0);
    TH1F *hDeltaPt_ptleadbins[ptbins];
    for(int i=0; i<ptbins; i++){
        
    hDeltaPt_ptleadbins[i] = (TH1F*)hDeltaPtLeadPt->ProjectionY(Form("hDeltaPt_ptleadbins_%d",i),hDeltaPtLeadPt->GetXaxis()->FindBin(ptleadDown[i]),hDeltaPtLeadPt->GetXaxis()->FindBin(ptleadUp[i])-1);
        hDeltaPt_ptleadbins[i]->SetLineColor(colors[i]);
        hDeltaPt_ptleadbins[i]->Scale(1./hDeltaPt_ptleadbins[i]->Integral());
        legDeltaPt->AddEntry(hDeltaPt_ptleadbins[i],Form("%d-%d GeV/c",ptleadDown[i],ptleadUp[i]),"l");
    }
    
   TH1F* hDeltaPt_ptleadbin5 = (TH1F*)hDeltaPtLeadPt->ProjectionY("hDeltaPt_ptleadbin5",hDeltaPtLeadPt->GetXaxis()->FindBin(5),hDeltaPtLeadPt->GetXaxis()->FindBin(100)-1);
    TH1F* hDeltaPt_ptleadbin10 = (TH1F*)hDeltaPtLeadPt->ProjectionY("hDeltaPt_ptleadbin10",hDeltaPtLeadPt->GetXaxis()->FindBin(10),hDeltaPtLeadPt->GetXaxis()->FindBin(100)-1);
    
    // vs lead. jet pT, excl. lead. jet
    
    TH2F *hDeltaPtLeadPt_exlcuding_tmp = (TH2F*)histList1->FindObject("fDeltaPTLeadPt_excl_lead");
    TH2F *hDeltaPtLeadPt_exlcuding = (TH2F*)hDeltaPtLeadPt_exlcuding_tmp->Clone("hDeltaPtLeadPt_exlcuding");
    hDeltaPtLeadPt_exlcuding->SetTitle();
    hDeltaPtLeadPt_exlcuding->GetYaxis()->SetRangeUser(-10,100);
    
    TH2F *hDeltaPtLeadPtFit_exlcuding = (TH2F*)hDeltaPtLeadPt_exlcuding_tmp->Clone("hDeltaPtLeadPtFit_exlcuding");
    hDeltaPtLeadPtFit_exlcuding->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFit_exlcuding->GetYaxis()->SetRangeUser(-2,5);
    hDeltaPtLeadPtFit_exlcuding->RebinX(10);
    //hDeltaPtLeadPtFit_exlcuding->RebinY(2);
    
    TF1 *fit = new TF1("fit","landau",-2,5);
    hDeltaPtLeadPtFit_exlcuding->FitSlicesY(fit,hDeltaPtLeadPtFit_exlcuding->GetXaxis()->FindBin(0),hDeltaPtLeadPtFit_exlcuding->GetXaxis()->FindBin(50),0,"QNR");
    TH1F *hDeltaPtLeadPtFitMean_exlcuding = (TH1F*)gDirectory->Get("hDeltaPtLeadPtFit_exlcuding_1");
    hDeltaPtLeadPtFitMean_exlcuding->SetName("hDeltaPtLeadPtFitMean_exlcuding");
    hDeltaPtLeadPtFitMean_exlcuding->SetTitle();
    hDeltaPtLeadPtFitMean_exlcuding->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFitMean_exlcuding->GetYaxis()->SetRangeUser(-0.6,-0.1);
    TH1F *hDeltaPtLeadPtFitSigma_exlcuding = (TH1F*)gDirectory->Get("hDeltaPtLeadPtFit_exlcuding_2");
    hDeltaPtLeadPtFitSigma_exlcuding->SetName("hDeltaPtLeadPtFitSigma_exlcuding");
    hDeltaPtLeadPtFitSigma_exlcuding->SetTitle();
    hDeltaPtLeadPtFitSigma_exlcuding->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFitSigma_exlcuding->GetYaxis()->SetRangeUser(0.2,0.8);
    
    TH2F *hDeltaPtLeadPt_exlcudingD_tmp = (TH2F*)histList1->FindObject("fDeltaPTLeadPt_excl_Dlead");
    TH2F *hDeltaPtLeadPt_exlcudingD = (TH2F*)hDeltaPtLeadPt_exlcudingD_tmp->Clone("hDeltaPtLeadPt_exlcudingD");
    hDeltaPtLeadPt_exlcudingD->SetTitle();
    hDeltaPtLeadPt_exlcudingD->GetYaxis()->SetRangeUser(-10,100);
    
    TCanvas *cDetlaPtLeadPt_exlcuding = new TCanvas("cDetlaPtLeadPt_exlcuding","cDetlaPtLeadPt_exlcuding",800,600);
    hDeltaPtLeadPt_exlcuding->Draw("colz");
    cDetlaPtLeadPt_exlcuding->SaveAs("plotsDjets/DetlaPtLeadPt_exlcuding.png");
    
    TProfile *hDeltaPtLeadPtProf_excleading = (TProfile*)hDeltaPtLeadPt_exlcuding->ProfileX("hDeltaPtLeadPtProf_excleading");
    hDeltaPtLeadPtProf_excleading->SetTitle();
    //hDeltaPtCentProf->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtLeadPtProf_excleading->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtLeadPtProf_excleading->BuildOptions(0,0,"i");
    
    TH1F *hDeltaPt_ptlead_exlcuding = (TH1F*)hDeltaPtLeadPt_exlcuding->ProjectionY("hDeltaPt_ptlead_exlcuding");
    hDeltaPt_ptlead_exlcuding->SetTitle();
    
    TLegend *legDeltaPt_exlcuding = new TLegend(0.5,0.5,0.88,0.88,"D-jet, p_{T,lead}, exc. lead:");
    legDeltaPt_exlcuding->SetBorderSize(0);
    TH1F *hDeltaPt_ptleadbins_exlcuding[ptbins];
    for(int i=0; i<ptbins; i++){
        
        hDeltaPt_ptleadbins_exlcuding[i] = (TH1F*)hDeltaPtLeadPt_exlcuding->ProjectionY(Form("hDeltaPt_ptleadbins_exlcuding_%d",i),hDeltaPtLeadPt_exlcuding->GetXaxis()->FindBin(ptleadDown[i]),hDeltaPtLeadPt_exlcuding->GetXaxis()->FindBin(ptleadUp[i])-1);
        hDeltaPt_ptleadbins_exlcuding[i]->SetLineColor(colors[i]);
        hDeltaPt_ptleadbins_exlcuding[i]->Scale(1./hDeltaPt_ptleadbins_exlcuding[i]->Integral());
        legDeltaPt_exlcuding->AddEntry(hDeltaPt_ptleadbins_exlcuding[i],Form("%d-%d GeV/c",ptleadDown[i],ptleadUp[i]),"l");
    }
    
       TH1F* hDeltaPt_ptleadbin5_exlcuding = (TH1F*)hDeltaPtLeadPt_exlcuding->ProjectionY("hDeltaPt_ptleadbin5_exlcuding",hDeltaPtLeadPt_exlcuding->GetXaxis()->FindBin(5),hDeltaPtLeadPt_exlcuding->GetXaxis()->FindBin(100)-1);
    TH1F* hDeltaPt_ptleadbin10_exlcuding = (TH1F*)hDeltaPtLeadPt_exlcuding->ProjectionY("hDeltaPt_ptleadbin10_exlcuding",hDeltaPtLeadPt_exlcuding->GetXaxis()->FindBin(10),hDeltaPtLeadPt_exlcuding->GetXaxis()->FindBin(100)-1);
    
    TH1F* hDeltaPt_ptleadbin5_exlcudingD = (TH1F*)hDeltaPtLeadPt_exlcudingD->ProjectionY("hDeltaPt_ptleadbin5_exlcudingD",hDeltaPtLeadPt_exlcudingD->GetXaxis()->FindBin(5),hDeltaPtLeadPt_exlcudingD->GetXaxis()->FindBin(100)-1);
    TH1F* hDeltaPt_ptleadbin10_exlcudingD = (TH1F*)hDeltaPtLeadPt_exlcudingD->ProjectionY("hDeltaPt_ptleadbin10_exlcudingD",hDeltaPtLeadPt_exlcudingD->GetXaxis()->FindBin(10),hDeltaPtLeadPt_exlcudingD->GetXaxis()->FindBin(100)-1);
    
    // === Delta pT no Bkg =====
    
    // vs lead. jet pT

    TH2F *hDeltaPtLeadPtNoBkg = (TH2F*)histList1->FindObject("fDeltaPTLeadPtNoBkg");
    hDeltaPtLeadPtNoBkg->SetTitle();
    hDeltaPtLeadPtNoBkg->GetYaxis()->SetRangeUser(-10,100);
    
    TCanvas *cDetlaPtLeadPtNoBkg = new TCanvas("cDetlaPtLeadPtNoBkg","cDetlaPtLeadPtNoBkg",800,600);
    hDeltaPtLeadPtNoBkg->Draw("colz");
    cDetlaPtLeadPtNoBkg->SaveAs("plotsDjets/DetlaPtLeadPt.png");
    
    TProfile *hDeltaPtLeadPtProfNoBkg = (TProfile*)hDeltaPtLeadPtNoBkg->ProfileX("hDeltaPtLeadPtProfNoBkg");
    hDeltaPtLeadPtProfNoBkg->SetTitle();
    //hDeltaPtCentProfNoBkg->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtLeadPtProfNoBkg->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtLeadPtProfNoBkg->BuildOptions(0,0,"i");
    
    TH1F *hDeltaPt_ptleadNoBkg = (TH1F*)hDeltaPtLeadPtNoBkg->ProjectionY("hDeltaPt_ptleadNoBkg");
    hDeltaPt_ptleadNoBkg->SetTitle();
    
    TLegend *legDeltaPtNoBkg = new TLegend(0.5,0.5,0.89,0.88,"D-jet woBkg, p_{T,lead}:");
    legDeltaPtNoBkg->SetBorderSize(0);
    TH1F *hDeltaPt_ptleadbinsNoBkg[ptbins];
    for(int i=0; i<ptbins; i++){
        
    hDeltaPt_ptleadbinsNoBkg[i] = (TH1F*)hDeltaPtLeadPtNoBkg->ProjectionY(Form("hDeltaPt_ptleadbinsNoBkg_%d",i),hDeltaPtLeadPtNoBkg->GetXaxis()->FindBin(ptleadDown[i]),hDeltaPtLeadPtNoBkg->GetXaxis()->FindBin(ptleadUp[i])-1);
        hDeltaPt_ptleadbinsNoBkg[i]->SetLineColor(colors[i]);
        hDeltaPt_ptleadbinsNoBkg[i]->Scale(1./hDeltaPt_ptleadbinsNoBkg[i]->Integral());
        legDeltaPtNoBkg->AddEntry(hDeltaPt_ptleadbinsNoBkg[i],Form("%d-%d GeV/c",ptleadDown[i],ptleadUp[i]),"l");
    }
    
    // vs lead. jet pT, excl. lead. jet
    
    TH2F *hDeltaPtLeadPtNoBkg_exlcuding = (TH2F*)histList1->FindObject("fDeltaPTLeadPtRawNoBkg_excl_lead");
    hDeltaPtLeadPtNoBkg_exlcuding->SetTitle();
    hDeltaPtLeadPtNoBkg_exlcuding->GetYaxis()->SetRangeUser(-10,100);
    
    TCanvas *cDetlaPtLeadPtNoBkg_exlcuding = new TCanvas("cDetlaPtLeadPtNoBkg_exlcuding","cDetlaPtLeadPtNoBkg_exlcuding",800,600);
    hDeltaPtLeadPtNoBkg_exlcuding->Draw("colz");
    cDetlaPtLeadPtNoBkg_exlcuding->SaveAs("plotsDjets/DetlaPtLeadPtNoBkg_exlcuding.png");
    
    TProfile *hDeltaPtLeadPtProfNoBkg_excleading = (TProfile*)hDeltaPtLeadPtNoBkg_exlcuding->ProfileX("hDeltaPtLeadPtProfNoBkg_excleading");
    hDeltaPtLeadPtProfNoBkg_excleading->SetTitle();
    //hDeltaPtCentProfNoBkg->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtLeadPtProfNoBkg_excleading->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtLeadPtProfNoBkg_excleading->BuildOptions(0,0,"i");
    
    TH1F *hDeltaPt_ptleadNoBkg_exlcuding = (TH1F*)hDeltaPtLeadPtNoBkg_exlcuding->ProjectionY("hDeltaPt_ptleadNoBkg_exlcuding");
    hDeltaPt_ptleadNoBkg_exlcuding->SetTitle();
    
    TLegend *legDeltaPtNoBkg_exlcuding = new TLegend(0.5,0.5,0.88,0.88,"D-jet wo bkg, p_{T,lead}, exc. lead:");
    legDeltaPtNoBkg_exlcuding->SetBorderSize(0);
    TH1F *hDeltaPt_ptleadbinsNoBkg_exlcuding[ptbins];
    for(int i=0; i<ptbins; i++){
        
        hDeltaPt_ptleadbinsNoBkg_exlcuding[i] = (TH1F*)hDeltaPtLeadPtNoBkg_exlcuding->ProjectionY(Form("hDeltaPt_ptleadbinsNoBkg_exlcuding_%d",i),hDeltaPtLeadPtNoBkg_exlcuding->GetXaxis()->FindBin(ptleadDown[i]),hDeltaPtLeadPtNoBkg_exlcuding->GetXaxis()->FindBin(ptleadUp[i])-1);
        hDeltaPt_ptleadbinsNoBkg_exlcuding[i]->SetLineColor(colors[i]);
        hDeltaPt_ptleadbinsNoBkg_exlcuding[i]->Scale(1./hDeltaPt_ptleadbinsNoBkg_exlcuding[i]->Integral());
        legDeltaPtNoBkg_exlcuding->AddEntry(hDeltaPt_ptleadbinsNoBkg_exlcuding[i],Form("%d-%d GeV/c",ptleadDown[i],ptleadUp[i]),"l");
    }
    
       TH1F* hDeltaPt_ptleadbin5NoBkg_exlcuding = (TH1F*)hDeltaPtLeadPtNoBkg_exlcuding->ProjectionY("hDeltaPt_ptleadbin5NoBkg_exlcuding",hDeltaPtLeadPtNoBkg_exlcuding->GetXaxis()->FindBin(5),hDeltaPtLeadPtNoBkg_exlcuding->GetXaxis()->FindBin(100)-1);
       TH1F* hDeltaPt_ptleadbin10NoBkg_exlcuding = (TH1F*)hDeltaPtLeadPtNoBkg_exlcuding->ProjectionY("hDeltaPt_ptleadbin10NoBkg_exlcuding",hDeltaPtLeadPtNoBkg_exlcuding->GetXaxis()->FindBin(10),hDeltaPtLeadPtNoBkg_exlcuding->GetXaxis()->FindBin(100)-1);

    

//======================= nTracks data ==================================

	THnSparseF *sparse1 = (THnSparseF*)histList1->FindObject("hsDphiz");
    //sparse1->GetAxis(2)->SetRangeUser(Dptmin,Dptmax); 
    //sparse1->GetAxis(1)->SetRangeUser(jetmin,jetmax); 
    //sparse1->GetAxis(5)->SetRangeUser(jetmin,jetmax); 
    
    TH2F *hNtracksPtJet = (TH2F*)sparse1->Projection(5,1,"E");
    hNtracksPtJet->SetTitle();
    hNtracksPtJet->Sumw2();
    
    TCanvas *cNtracksPtJet = new TCanvas("cNtracksPtJet","cNtracksPtJet",800,600);
    hNtracksPtJet->Draw("colz");
    hNtracksPtJet->GetXaxis()->SetTitle("p_{T,jet} (GeV/c)");
    hNtracksPtJet->GetYaxis()->SetTitle("#constituents");
    hNtracksPtJet->GetYaxis()->SetRangeUser(0,30);
    cNtracksPtJet->SaveAs("plotsDjets/NtracksVsJetPt.png");
     
    TH1F *hNtracks[16];
    TCanvas *cNtracks = new TCanvas("cNtracks","cNtracks",1600,1600);
    cNtracks->Divide(4,4);
    
    for(int i=0; i<nJetBins-1;i++){
        cNtracks->cd(i+1);
        TPaveText *pv1 = new TPaveText(0.6,0.2,0.9,0.3,"brNDC");
        pv1->SetFillStyle(0);
        pv1->SetBorderSize(0);
        pv1->AddText(Form("%.0f < p_{T,jet} < %.0f GeV/c",ptJetbins[i], ptJetbins[i+1]));
        hNtracks[i] = (TH1F*)hNtracksPtJet->ProjectionY(Form("hNtracks%d",i),hNtracksPtJet->GetXaxis()->FindBin(ptJetbins[i]),hNtracksPtJet->GetXaxis()->FindBin(ptJetbins[i+1]));
        hNtracks[i]->Draw();
        pv1->Draw("same");
       // hNtracks[i]->Write();
    }
    
      cNtracks->SaveAs("plotsDjets/Ntracks.png"); 
      cNtracks->SaveAs("plotsDjets/Ntracks.pdf"); 
   
//======================= nTracks data End ==================================
   

    
    
    //============================== DRAW ========================================
  
  
    TCanvas *cRhoM2 = new TCanvas("cRhoM2","cRhoM2",800,600);
    hRhoMultProf->Draw();
    cRhoM2->SaveAs("plotsDjets/RhoVsMultProf.png");
    
    TCanvas *cRhoLeadPt2 = new TCanvas("cRhoLeadPt2","cRhoLeadPt2",800,600);
    hRhoLeadPtProf->Draw();
    cRhoLeadPt2->SaveAs("plotsDjets/RhoVsLeadPtProf.png");

    TCanvas *cRhoMultProfJetBins = new TCanvas("cRhoMultProfJetBins","cRhoMultProfJetBins",800,600);
    hRhoMultProfBinsPt[0]->Draw();
    for(int j=1;j<ptbins;j++){
        hRhoMultProfBinsPt[j]->Draw("same");
    }
    leg1->Draw("same");

    cRhoMultProfJetBins->SaveAs("plotsDjets/RhoVsMult_leadJetBins.png");
    
     TCanvas *cRhoMultProfJetBins2 = new TCanvas("cRhoMultProfJetBins2","cRhoMultProfJetBins2",800,600);
    hRhoMultProfBinsPt2[0]->Draw();
    for(int j=1;j<ptbins;j++){
        hRhoMultProfBinsPt2[j]->Draw("same");
    }
    leg11->Draw("same");

    cRhoMultProfJetBins2->SaveAs("plotsDjets/RhoVsMult_leadJetCut.png");
    
    
    TLegend *ldelta1 = new TLegend(0.6,0.7,0.85,0.85,"D-jet events");
    ldelta1->SetBorderSize(0);
    ldelta1->AddEntry(hDeltaPt_c,"All jets","l");
    ldelta1->AddEntry(hDeltaPt_c_excleading,"Excl. lead. jet","l");
    
    
    TCanvas *cRhoMultProfCentBins = new TCanvas("cRhoMultProfCentBins","cRhoMultProfCentBins",800,600);
    hRhoMultProfBins[0]->Draw();
    for(int j=1;j<centBins;j++){
        hRhoMultProfBins[j]->Draw("same");
    }
    leg2->Draw("same");

    cRhoMultProfCentBins->SaveAs("plotsDjets/RhoVsLeadJetPt_centBins.png");
  
  // ===== Delta pT ========
  
    // vs cent
  
    // profile
    hDeltaPtCentProf->GetXaxis()->SetRangeUser(-10,60);
    hDeltaPtCentProf->SetLineColor(1);
    hDeltaPtCentProf_excleading->SetLineColor(2);
    TCanvas *cDeltaPtCentProf = new TCanvas("cDeltaPtCentProf","cDeltaPtCentProf",800,600);
    hDeltaPtCentProf->Draw();
    hDeltaPtCentProf_excleading->Draw("same");
    ldelta1->Draw("same");
    cDeltaPtCentProf->SaveAs("plotsDjets/DeltaPtCentProf.png");

   // 1D
    hDeltaPt_c->GetXaxis()->SetRangeUser(-10,60);
    hDeltaPt_c->SetLineColor(1);
    hDeltaPt_c_excleading->SetLineColor(2);
    TCanvas *cDeltaPt_c = new TCanvas("cDeltaPt_c","cDeltaPt_c",800,600);
    cDeltaPt_c->SetLogy();
    hDeltaPt_c->Draw();
    hDeltaPt_c_excleading->Draw("same");
    ldelta1->Draw("same");
    cDeltaPt_c->SaveAs("plotsDjets/DeltaPt_exclComparison.png");

    //==== vs lead jet pT
    
    TCanvas *cDeltaPtLeadPt = new TCanvas("cDeltaPtLeadPt","cDeltaPtLeadPt",800,600);
    hDeltaPt_ptlead->Draw();
    cDeltaPtLeadPt->SaveAs("plotsDjets/DeltaPtLead.png");

    // profile
    hDeltaPtLeadPtProf->Rebin(4);
    hDeltaPtLeadPtProf_excleading->Rebin(4);
    hDeltaPtLeadPtProf->SetLineColor(1);
    hDeltaPtLeadPtProf_excleading->SetLineColor(2);
    hDeltaPtLeadPtProf->GetYaxis()->SetRangeUser(0,5);
    hDeltaPtLeadPtProf->GetXaxis()->SetRangeUser(0,60);
    
    
    TCanvas *cDeltaPtLeadPtProf = new TCanvas("cDeltaPtLeadPtProf","cDeltaPtLeadPtProf",800,600);
    hDeltaPtLeadPtProf->Draw();
    hDeltaPtLeadPtProf_excleading->Draw("same");
    ldelta1->Draw("same");
    cDeltaPtLeadPtProf->SaveAs("plotsDjets/DeltaPtLeadPtProf.png");


    // jet pT bins

    TCanvas *cDeltaPtLeadPtBins = new TCanvas("cDeltaPtLeadPtBins","cDeltaPtLeadPtBins",800,600);
    cDeltaPtLeadPtBins->SetLogy();
    hDeltaPt_ptleadbins[0]->SetMaximum(hDeltaPt_ptleadbins[0]->GetMaximum()*2);
    hDeltaPt_ptleadbins[0]->GetXaxis()->SetRangeUser(-10,60);
    hDeltaPt_ptleadbins[0]->Draw();
    for(int j=1;j<ptbins;j++){
        hDeltaPt_ptleadbins[j]->Draw("same");
    }
    legDeltaPt->Draw("same");
    cDeltaPtLeadPtBins->SaveAs("plotsDjets/DeltaPtLeadPtBins.png");


    TCanvas *cDeltaPtLeadPtBins_exlcuding = new TCanvas("cDeltaPtLeadPtBins_exlcuding","cDeltaPtLeadPtBins_exlcuding",800,600);
    cDeltaPtLeadPtBins_exlcuding->SetLogy();
    hDeltaPt_ptleadbins_exlcuding[0]->SetMaximum(hDeltaPt_ptleadbins_exlcuding[0]->GetMaximum()*2);
    hDeltaPt_ptleadbins_exlcuding[0]->GetXaxis()->SetRangeUser(-10,60);
    hDeltaPt_ptleadbins_exlcuding[0]->Draw();
    for(int j=1;j<centBins;j++){
        hDeltaPt_ptleadbins_exlcuding[j]->Draw("same");
    }
    legDeltaPt_exlcuding->Draw("same");
    cDeltaPtLeadPtBins_exlcuding->SaveAs("plotsDjets/DeltaPtLeadPtBins_exlcuding.png");



    //==== vs lead jet pT, no bkg
    
    TCanvas *cDeltaPtLeadPtNoBkg = new TCanvas("cDeltaPtLeadPtNoBkg","cDeltaPtLeadPtNoBkg",800,600);
    hDeltaPt_ptleadNoBkg->Draw();
    cDeltaPtLeadPtNoBkg->SaveAs("plotsDjets/DeltaPtLeadNoBkg.png");

    // profile
    hDeltaPtLeadPtProfNoBkg->Rebin(4);
    hDeltaPtLeadPtProfNoBkg_excleading->Rebin(4);
    hDeltaPtLeadPtProfNoBkg->SetLineColor(1);
    hDeltaPtLeadPtProfNoBkg_excleading->SetLineColor(2);
    hDeltaPtLeadPtProfNoBkg->GetYaxis()->SetRangeUser(0,5);
    hDeltaPtLeadPtProfNoBkg->GetXaxis()->SetRangeUser(0,60);
    
    TCanvas *cDeltaPtLeadPtProfNoBkg = new TCanvas("cDeltaPtLeadPtProfNoBkg","cDeltaPtLeadPtProfNoBkg",800,600);
    hDeltaPtLeadPtProfNoBkg->Draw();
    hDeltaPtLeadPtProfNoBkg_excleading->Draw("same");
    ldelta1->Draw("same");
    cDeltaPtLeadPtProfNoBkg->SaveAs("plotsDjets/DeltaPtLeadPtProfNoBkg.png");

    // jet pT bins

    TCanvas *cDeltaPtLeadPtBinsNoBkg = new TCanvas("cDeltaPtLeadPtBinsNoBkg","cDeltaPtLeadPtBinsNoBkg",800,600);
    cDeltaPtLeadPtBinsNoBkg->SetLogy();
    hDeltaPt_ptleadbinsNoBkg[0]->SetMaximum(hDeltaPt_ptleadbins[0]->GetMaximum()*2);
    hDeltaPt_ptleadbinsNoBkg[0]->GetXaxis()->SetRangeUser(-2,30);
    hDeltaPt_ptleadbinsNoBkg[0]->Draw();
    for(int j=1;j<ptbins;j++){
        hDeltaPt_ptleadbinsNoBkg[j]->Draw("same");
    }
    legDeltaPt->Draw("same");
    cDeltaPtLeadPtBinsNoBkg->SaveAs("plotsDjets/DeltaPtLeadPtBinsNoBkg.png");


    TCanvas *cDeltaPtLeadPtBinsNoBkg_exlcuding = new TCanvas("cDeltaPtLeadPtBinsNoBkg_exlcuding","cDeltaPtLeadPtBinsNoBkg_exlcuding",800,600);
    cDeltaPtLeadPtBinsNoBkg_exlcuding->SetLogy();
    hDeltaPt_ptleadbinsNoBkg_exlcuding[0]->SetMaximum(hDeltaPt_ptleadbins_exlcuding[0]->GetMaximum()*2);
    hDeltaPt_ptleadbinsNoBkg_exlcuding[0]->GetXaxis()->SetRangeUser(-2,30);
    hDeltaPt_ptleadbinsNoBkg_exlcuding[0]->Draw();
    for(int j=1;j<centBins;j++){
        hDeltaPt_ptleadbinsNoBkg_exlcuding[j]->Draw("same");
    }
    legDeltaPt_exlcuding->Draw("same");
    cDeltaPtLeadPtBinsNoBkg_exlcuding->SaveAs("plotsDjets/DeltaPtLeadPtBinsNoBkg_exlcuding.png");
  
  
    hDeltaPtLeadPtProfNoBkg_excleading->SetLineColor(kGreen);
    hDeltaPtLeadPtProf_excleading->GetXaxis()->SetRangeUser(0,60);
    hDeltaPtLeadPtProf_excleading->GetYaxis()->SetRangeUser(0,5);
    // bkg vs no bkg subtraction
    TCanvas *cDeltaPtLeadPtProf_bkgCOmparison = new TCanvas("cDeltaPtLeadPtProf_bkgCOmparison","cDeltaPtLeadPtProf_bkgCOmparison",800,600);
    hDeltaPtLeadPtProf_excleading->Draw();
    hDeltaPtLeadPtProfNoBkg_excleading->Draw("same");
     TLegend *legBkg = new TLegend(0.6,0.7,0.85,0.82,"lead. jet excl.");
     legBkg->SetBorderSize(0);
    legBkg->AddEntry(hDeltaPtLeadPtProf_excleading,"bkg sub","l");
    legBkg->AddEntry(hDeltaPtLeadPtProfNoBkg_excleading,"bkg not sub","l");
    legBkg->Draw("same");
    cDeltaPtLeadPtProf_bkgCOmparison->SaveAs("plotsDjets/DeltaPtLeadPtProf_bkgComparison.png");
  
  
TFile *fOutDjet = new TFile("fOutDjetProperties.root","RECREATE");
hRhoMult->Write();
hRhoLeadPt->Write();
hDLeading->Write();
hDLeading2D->Write();
hRhoLeadPtMult->Write();
hRho_ptleadbin5->Write();
hRho_ptleadbin10->Write();
hRho_ptleadbin5_cent0_10->Write();
hRho_ptleadbin5_cent10_20->Write();
hRho_ptleadbin5_cent20_40->Write();
hRho_ptleadbin5_cent40_60->Write();
hRho_ptleadbin5_cent60_100->Write();
hRhoMult_leadpt5->Write();
hRhoMultProf_leadpt5->Write();
hRhoMult_leadpt10->Write();
hRhoMultProf_leadpt10->Write();

hDeltaPtCent->Write();
hDeltaPt_c->Write();
hDeltaPtCent_excleading->Write();
hDeltaPt_c_excleading->Write();
hDeltaPtLeadPt->Write();
hDeltaPt_ptlead->Write();
hDeltaPt_ptleadbin5->Write();
hDeltaPt_ptleadbin10->Write();
hDeltaPtLeadPt_exlcuding->Write();
hDeltaPt_ptlead_exlcuding->Write();
hDeltaPt_ptleadbin5_exlcuding->Write();
hDeltaPt_ptleadbin10_exlcuding->Write();
hDeltaPt_ptleadbin5_exlcudingD->Write();
hDeltaPt_ptleadbin10_exlcudingD->Write();

hDeltaPtLeadPtNoBkg->Write();
hDeltaPt_ptleadNoBkg->Write();
hDeltaPtLeadPtNoBkg_exlcuding->Write();
hDeltaPt_ptleadNoBkg_exlcuding->Write();
hDeltaPt_ptleadbin5NoBkg_exlcuding->Write();
hDeltaPt_ptleadbin10NoBkg_exlcuding->Write();
hRhoMultProf->Write();
hRhoLeadPtProf->Write();
hDeltaPtCentProf->Write();
hDeltaPtCentProf_excleading->Write();
hDeltaPtLeadPtProf->Write();
hDeltaPtLeadPtProf_excleading->Write();
hDeltaPtLeadPtProfNoBkg->Write();
hDeltaPtLeadPtProfNoBkg_excleading->Write();

hDeltaPtLeadPtFit_exlcuding->Write();
hDeltaPtLeadPtFitMean_exlcuding->Write();
hDeltaPtLeadPtFitSigma_exlcuding->Write();

    for(int i=0; i<ptbins; i++){
        hDeltaPt_ptleadbins_exlcuding[i]->Write();
    }
    

fOutDjet->Close();
   

}
