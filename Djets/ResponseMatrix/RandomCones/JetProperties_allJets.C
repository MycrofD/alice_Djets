#include <string>
#include <sstream>
#include <iostream>

const int nJetBins = 17;
double ptJetbins[nJetBins] = { 3,4,6,8,10,12,14,16,18,20,22,24,26,30,40,50,70 };

//TString outDir = "plots_DzeroEvt";
TString outDir = "plots_DzeroLeadEvt";
//TString outDir = "plots_AllEvt";

//TString outFile = "DzeroEvt";
TString outFile = "DzeroLeadEvt";
//TString outFile = "AllEvt";

//file: /home/basia/Work/alice/analysis/pPb_run2/D0jet/outRC/AnalysisResults_Dzero.root
//AnalysisResults_allEvents.root
//AnalysisResults_DzeroTrigg.root

void JetProperties_allJets(TString datafile = "/home/basia/Work/alice/analysis/pPb_run2/D0jet/outRC03/AnalysisResults_DzeroTrigg.root", int Dptmin = 0, int Dptmax = 100, int Rpar = 4 )
{

     gStyle->SetOptStat(0000); //Mean and RMS shown
     gStyle->SetPadRightMargin(0.1);
     gSystem->Exec(Form("mkdir %s",outDir.Data()));

    double zmin = -2., zmax = 2;
    float jetmin = 0, jetmax = 100;
    double plotmin = 0, plotmax = 50;

    TPaveText *pv2 = new TPaveText(0.15,0.8,0.25,0.9,"brNDC");
    pv2->SetFillStyle(0);
    pv2->SetBorderSize(0);
    pv2->AddText(Form("R=0.%d",Rpar));

    TFile *File1 = new TFile(datafile,"read");
    TDirectoryFile* dir1=(TDirectoryFile*)File1->Get("DmesonsForJetCorrelations");
    TList *histList1 =  (TList*)dir1->Get("histosD0MBN0");

    TH1F *hPtTrk = (TH1F*)histList1->FindObject("hPtJetTrks");
    hPtTrk->SetTitle();
    hPtTrk->Scale(1./hPtTrk->Integral());

    TH2F *hRhoMult = (TH2F*)histList1->FindObject("fhRhoMult");
    hRhoMult->SetTitle();

    TCanvas *cRhoM = new TCanvas("cRhoM","cRhoM",800,600);
    hRhoMult->Draw("colz");

    SaveCanvas(cRhoM,"RhoVsMult");

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

    SaveCanvas(cRhoLeadPt,"RhoVsLeadPt");

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

    SaveCanvas(cLeadD,"LeadingDJet");

    TH2F *hDLeading2D = (TH2F*)histList1->FindObject("fhDleadStatJetPt");
    hDLeading2D->SetTitle();

    TCanvas *cLeadD2D = new TCanvas("cLeadD2D","cLeadD2D",800,600);
    cLeadD2D->SetLogz();
    hDLeading2D->Draw("colz");

    SaveCanvas(cLeadD2D,"LeadingDJetVsPt");

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

    SaveCanvas(cRhoM_pt5,"RhoVsMult_pt5");

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

    SaveCanvas(cRhoM_pt10,"RhoVsMult_pt10");

    TProfile *hRhoMultProf_leadpt10 = (TProfile*)hRhoMult_leadpt10->ProfileX("hRhoMultProf_leadpt10");
    hRhoMultProf_leadpt10->SetTitle();
    hRhoMultProf_leadpt10->GetYaxis()->SetRangeUser(0,10);
    hRhoMultProf_leadpt10->GetYaxis()->SetTitle("<#rho(GeV/c*rad^{-1})>");
    hRhoMultProf_leadpt10->BuildOptions(0,0,"i");

    //======================= pt bins ==================================

    TLegend *leg1 = new TLegend(0.5,0.5,0.88,0.88,"p_{T,lead}:");
    leg1->SetBorderSize(0);

    const int ptbins = 6;
    int ptleadDown[] = { 0,2,5,10,20,30 };
    int ptleadUp[] = { 2,5,10,20,30,80 };
    const Int_t nptbins = 13;
    Float_t ptlead[nptbins+1] = { 0,3,5,6,8,10,12,16,20,25,30,40,50,80 };
    int colors[] = { 1,2,4,6,8,7 };

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
      hRhoMultProfBinsPt[i]->SetMarkerColor(colors[i]);
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

    TLegend *leg11 = new TLegend(0.5,0.5,0.88,0.88,"p_{T,lead}:");
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
      hRhoMultProfBinsPt2[i]->SetMarkerColor(colors[i]);
      hRhoMultProfBinsPt2[i]->GetYaxis()->SetRangeUser(0,6);
      hRhoMultProfBinsPt2[i]->GetYaxis()->SetTitle("<#rho(GeV/c*rad^{-1})>");
      hRhoMultProfBinsPt2[i]->BuildOptions(0,0,"i");

      leg11->AddEntry(hRhoMultProfBinsPt2[i],Form(">%d GeV/c",ptleadDown[i]),"l");
    }

    //======================= cent bins ==================================

    TLegend *leg2 = new TLegend(0.6,0.6,0.95,0.88,"centrality:");
    leg2->SetBorderSize(0);

    const int centBins = 5;
    int centDown[] = { 0,10,20,40,60 };
    int centUp[] = { 10,20,40,60,100 };

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
      hRhoMultProfBins[i]->SetMarkerColor(colors[i]);
      hRhoMultProfBins[i]->GetYaxis()->SetRangeUser(0,6);
      hRhoMultProfBins[i]->GetXaxis()->SetRangeUser(0,100);
      hRhoMultProfBins[i]->GetYaxis()->SetTitle("<#rho(GeV/c*rad^{-1})>");
      hRhoMultProfBins[i]->BuildOptions(0,0,"i");

      leg2->AddEntry(hRhoMultProfBins[i],Form("%d-%d",centDown[i],centUp[i]),"l");

    }


 //======================= delta pt ==================================


    // vs centrality

  /*  TH2F *hDeltaPtCent = (TH2F*)histList1->FindObject("fDeltaPTCent");
    hDeltaPtCent->SetTitle();
    hDeltaPtCent->GetYaxis()->SetRangeUser(-10,100);

    TCanvas *cDetlaPtCent = new TCanvas("cDetlaPtCent","cDetlaPtCent",800,600);
    hDeltaPtCent->Draw("colz");

    cDetlaPtCent->SaveAs("plotsAlljets/DeltaPtCent.png");

    TProfile *hDeltaPtCentProf = (TProfile*)hDeltaPtCent->ProfileX("hDeltaPtCentProf");
    hDeltaPtCentProf->SetTitle();
    //hDeltaPtCentProf->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtCentProf->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtCentProf->BuildOptions(0,0,"i");

    TH1F *hDeltaPt_c = (TH1F*)hDeltaPtCent->ProjectionY("hDeltaPt_c");
    hDeltaPt_c->SetTitle();
    hDeltaPt_c->Scale(1./hDeltaPt_c->Integral());
*/

    // vs centrality excl. lead. jet

    TH2F *hDeltaPtCent_excleading = (TH2F*)histList1->FindObject("fDeltaPTCent_excl_lead");
    hDeltaPtCent_excleading->SetTitle();
    hDeltaPtCent_excleading->GetYaxis()->SetRangeUser(-10,100);

    TCanvas *cDeltaPtCent_excleading = new TCanvas("cDeltaPtCent_excleading","cDeltaPtCent_excleading",800,600);
    hDeltaPtCent_excleading->Draw("colz");

    SaveCanvas(cDeltaPtCent_excleading,"DeltaPtCent_excleading");

    TProfile *hDeltaPtCentProf_excleading = (TProfile*)hDeltaPtCent_excleading->ProfileX("hDeltaPtCentProf_excleading");
    hDeltaPtCentProf_excleading->SetTitle();
    //hDeltaPtCentProf->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtCentProf_excleading->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtCentProf_excleading->BuildOptions(0,0,"i");

    TH1F *hDeltaPt_c_excleading = (TH1F*)hDeltaPtCent_excleading->ProjectionY("hDeltaPt_c_excleading");
    hDeltaPt_c_excleading->SetTitle();
    hDeltaPt_c_excleading->Scale(1./hDeltaPt_c_excleading->Integral());

    // vs centrality excl. lead. jet -- trans plane

    TH2F *hDeltaPtCent_trans = (TH2F*)histList1->FindObject("fDeltaPTCent_trans");
    hDeltaPtCent_trans->SetTitle();
    hDeltaPtCent_trans->GetYaxis()->SetRangeUser(-10,100);

    TCanvas *cDeltaPtCent_trans = new TCanvas("cDeltaPtCent_trans","cDeltaPtCent_trans",800,600);
    hDeltaPtCent_trans->Draw("colz");

    SaveCanvas(cDeltaPtCent_trans,"DeltaPtCent_trans");

    TProfile *hDeltaPtCentProf_trans = (TProfile*)hDeltaPtCent_trans->ProfileX("hDeltaPtCentProf_trans");
    hDeltaPtCentProf_trans->SetTitle();
    //hDeltaPtCentProf->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtCentProf_trans->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtCentProf_trans->BuildOptions(0,0,"i");

    TH1F *hDeltaPt_c_trans = (TH1F*)hDeltaPtCent_trans->ProjectionY("hDeltaPt_c_trans");
    hDeltaPt_c_trans->SetTitle();
    hDeltaPt_c_trans->Scale(1./hDeltaPt_c_trans->Integral());


/*
    // vs lead. jet pT

    TH2F *hDeltaPtLeadPt = (TH2F*)histList1->FindObject("fDeltaPTLeadPt");
    hDeltaPtLeadPt->SetTitle();
    hDeltaPtLeadPt->GetYaxis()->SetRangeUser(-10,100);

    TCanvas *cDetlaPtLeadPt = new TCanvas("cDetlaPtLeadPt","cDetlaPtLeadPt",800,600);
    hDeltaPtLeadPt->Draw("colz");
    cDetlaPtLeadPt->SaveAs("plotsAlljets/DetlaPtLeadPt.png");

    TProfile *hDeltaPtLeadPtProf = (TProfile*)hDeltaPtLeadPt->ProfileX("hDeltaPtLeadPtProf");
    hDeltaPtLeadPtProf->SetTitle();
    //hDeltaPtCentProf->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtLeadPtProf->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtLeadPtProf->BuildOptions(0,0,"i");

    TH1F *hDeltaPt_ptlead = (TH1F*)hDeltaPtLeadPt->ProjectionY("hDeltaPt_ptlead");
    hDeltaPt_ptlead->SetTitle();

    TLegend *legDeltaPt = new TLegend(0.5,0.5,0.88,0.88,"inc. jet, p_{T,lead}:");
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
*/

    // vs lead. jet pT, excl. lead. jet

    TH2F *hDeltaPtLeadPt_excluding_tmp = (TH2F*)histList1->FindObject("fDeltaPTLeadPt_excl_lead");
    TH2F *hDeltaPtLeadPt_excluding = (TH2F*)hDeltaPtLeadPt_excluding_tmp->Clone("hDeltaPtLeadPt_excluding");
    hDeltaPtLeadPt_excluding->SetTitle();
    hDeltaPtLeadPt_excluding->GetYaxis()->SetRangeUser(-10,100);

  /*  TH2F *hDeltaPtLeadPtFit_excluding = (TH2F*)hDeltaPtLeadPt_excluding_tmp->Clone("hDeltaPtLeadPtFit_excluding");
    hDeltaPtLeadPtFit_excluding->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFit_excluding->GetYaxis()->SetRangeUser(-2,5);
    hDeltaPtLeadPtFit_excluding->RebinX(10);
    //hDeltaPtLeadPtFit_excluding->RebinY(2);

    TF1 *fit = new TF1("fit","landau",-2,5);
    hDeltaPtLeadPtFit_excluding->FitSlicesY(fit,hDeltaPtLeadPtFit_excluding->GetXaxis()->FindBin(0),hDeltaPtLeadPtFit_excluding->GetXaxis()->FindBin(50),0,"QNR");
    TH1F *hDeltaPtLeadPtFitMean_excluding = (TH1F*)gDirectory->Get("hDeltaPtLeadPtFit_excluding_1");
    hDeltaPtLeadPtFitMean_excluding->SetName("hDeltaPtLeadPtFitMean_excluding");
    hDeltaPtLeadPtFitMean_excluding->SetTitle();
    hDeltaPtLeadPtFitMean_excluding->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFitMean_excluding->GetYaxis()->SetRangeUser(-0.6,-0.1);
    TH1F *hDeltaPtLeadPtFitSigma_excluding = (TH1F*)gDirectory->Get("hDeltaPtLeadPtFit_excluding_2");
    hDeltaPtLeadPtFitSigma_excluding->SetName("hDeltaPtLeadPtFitSigma_excluding");
    hDeltaPtLeadPtFitSigma_excluding->SetTitle();
    hDeltaPtLeadPtFitSigma_excluding->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFitSigma_excluding->GetYaxis()->SetRangeUser(0.2,0.8);*/

    TCanvas *cDetlaPtLeadPt_excluding = new TCanvas("cDetlaPtLeadPt_excluding","cDetlaPtLeadPt_excluding",800,600);
    hDeltaPtLeadPt_excluding->Draw("colz");
    SaveCanvas(cDetlaPtLeadPt_excluding,"hDeltaPtLeadPt_excluding");

    TProfile *hDeltaPtLeadPtProf_excleading = (TProfile*)hDeltaPtLeadPt_excluding->ProfileX("hDeltaPtLeadPtProf_excleading");
    hDeltaPtLeadPtProf_excleading->SetTitle();
    //hDeltaPtCentProf->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtLeadPtProf_excleading->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtLeadPtProf_excleading->BuildOptions(0,0,"i");

    TH1F *hDeltaPt_ptlead_excluding = (TH1F*)hDeltaPtLeadPt_excluding->ProjectionY("hDeltaPt_ptlead_excluding");
    hDeltaPt_ptlead_excluding->SetTitle();

    TLegend *legDeltaPt_excluding = new TLegend(0.5,0.5,0.88,0.88,"p_{T,lead}:");
    legDeltaPt_excluding->SetBorderSize(0);
    TH1F *hDeltaPt_ptleadbins_excluding[ptbins];

    for(int i=0; i<ptbins; i++){

        hDeltaPt_ptleadbins_excluding[i] = (TH1F*)hDeltaPtLeadPt_excluding->ProjectionY(Form("hDeltaPt_ptleadbins_excluding_%d",i),hDeltaPtLeadPt_excluding->GetXaxis()->FindBin(ptleadDown[i]),hDeltaPtLeadPt_excluding->GetXaxis()->FindBin(ptleadUp[i])-1);
        hDeltaPt_ptleadbins_excluding[i]->SetLineColor(colors[i]);
        hDeltaPt_ptleadbins_excluding[i]->SetMarkerColor(colors[i]);
        hDeltaPt_ptleadbins_excluding[i]->Scale(1./hDeltaPt_ptleadbins_excluding[i]->Integral());
        legDeltaPt_excluding->AddEntry(hDeltaPt_ptleadbins_excluding[i],Form("%d-%d GeV/c",ptleadDown[i],ptleadUp[i]),"l");
    }

    TH1F *hDeltaPtLeadPtFitMean_excluding = new TH1F("hDeltaPtLeadPtFitMean_excluding","hDeltaPtLeadPtFitMean_excluding",nptbins,ptlead);
    hDeltaPtLeadPtFitMean_excluding->SetName("hDeltaPtLeadPtFitMean_excluding");
    hDeltaPtLeadPtFitMean_excluding->SetTitle();
    hDeltaPtLeadPtFitMean_excluding->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFitMean_excluding->GetYaxis()->SetRangeUser(-0.6,-0.1);
    TH1F *hDeltaPtLeadPtFitSigma_excluding = new TH1F("hDeltaPtLeadPtFitSigma_excluding","hDeltaPtLeadPtFitSigma_excluding",nptbins,ptlead);
    hDeltaPtLeadPtFitSigma_excluding->SetName("hDeltaPtLeadPtFitSigma_excluding");
    hDeltaPtLeadPtFitSigma_excluding->SetTitle();
    hDeltaPtLeadPtFitSigma_excluding->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFitSigma_excluding->GetYaxis()->SetRangeUser(0.2,0.8);
    for(int i=0; i<nptbins; i++){
        TH1F* tmp = (TH1F*)hDeltaPtLeadPt_excluding->ProjectionY(Form("tmp_%d",i),hDeltaPtLeadPt_excluding->GetXaxis()->FindBin(ptlead[i]),hDeltaPtLeadPt_excluding->GetXaxis()->FindBin(ptlead[i+1])-1);
        hDeltaPtLeadPtFitMean_excluding->SetBinContent(i+1,tmp->GetMean());
        hDeltaPtLeadPtFitMean_excluding->SetBinError(i+1,tmp->GetMeanError());
        hDeltaPtLeadPtFitSigma_excluding->SetBinContent(i+1,tmp->GetStdDev());
        hDeltaPtLeadPtFitSigma_excluding->SetBinError(i+1,tmp->GetStdDevError());
    }

    TH1F* hDeltaPt_ptleadbin5_excluding = (TH1F*)hDeltaPtLeadPt_excluding->ProjectionY("hDeltaPt_ptleadbin5_excluding",hDeltaPtLeadPt_excluding->GetXaxis()->FindBin(5),hDeltaPtLeadPt_excluding->GetXaxis()->FindBin(100)-1);
    TH1F* hDeltaPt_ptleadbin10_excluding = (TH1F*)hDeltaPtLeadPt_excluding->ProjectionY("hDeltaPt_ptleadbin10_excluding",hDeltaPtLeadPt_excluding->GetXaxis()->FindBin(10),hDeltaPtLeadPt_excluding->GetXaxis()->FindBin(100)-1);



    // vs lead. jet pT -- trans

    TH2F *hDeltaPtLeadPt_trans_tmp = (TH2F*)histList1->FindObject("fDeltaPTLeadPt_trans");
    TH2F *hDeltaPtLeadPt_trans = (TH2F*)hDeltaPtLeadPt_trans_tmp->Clone("hDeltaPtLeadPt_trans");
    hDeltaPtLeadPt_trans->SetTitle();
    hDeltaPtLeadPt_trans->GetYaxis()->SetRangeUser(-10,100);

    TH2F *hDeltaPtLeadPtFit_trans = (TH2F*)hDeltaPtLeadPt_trans_tmp->Clone("hDeltaPtLeadPtFit_trans");
    hDeltaPtLeadPtFit_trans->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFit_trans->GetYaxis()->SetRangeUser(-2,5);
    hDeltaPtLeadPtFit_trans->RebinX(10);
    //hDeltaPtLeadPtFit_trans->RebinY(2);

  /*  TF1 *fit = new TF1("fit","landau",-2,5);
    hDeltaPtLeadPtFit_trans->FitSlicesY(fit,hDeltaPtLeadPtFit_trans->GetXaxis()->FindBin(0),hDeltaPtLeadPtFit_trans->GetXaxis()->FindBin(50),0,"QNR");
    TH1F *hDeltaPtLeadPtFitMean_trans = (TH1F*)gDirectory->Get("hDeltaPtLeadPtFit_trans_1");
    hDeltaPtLeadPtFitMean_trans->SetName("hDeltaPtLeadPtFitMean_trans");
    hDeltaPtLeadPtFitMean_trans->SetTitle();
    hDeltaPtLeadPtFitMean_trans->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFitMean_trans->GetYaxis()->SetRangeUser(-0.6,-0.1);
    TH1F *hDeltaPtLeadPtFitSigma_trans = (TH1F*)gDirectory->Get("hDeltaPtLeadPtFit_trans_2");
    hDeltaPtLeadPtFitSigma_trans->SetName("hDeltaPtLeadPtFitSigma_trans");
    hDeltaPtLeadPtFitSigma_trans->SetTitle();
    hDeltaPtLeadPtFitSigma_trans->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFitSigma_trans->GetYaxis()->SetRangeUser(0.2,0.8);*/

    TCanvas *cDetlaPtLeadPt_trans = new TCanvas("cDetlaPtLeadPt_trans","cDetlaPtLeadPt_trans",800,600);
    hDeltaPtLeadPt_trans->Draw("colz");
    SaveCanvas(cDetlaPtLeadPt_trans,"hDeltaPtLeadPt_trans");

    TProfile *hDeltaPtLeadPtProf_trans = (TProfile*)hDeltaPtLeadPt_trans->ProfileX("hDeltaPtLeadPtProf_trans");
    hDeltaPtLeadPtProf_trans->SetTitle();
    //hDeltaPtCentProf->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtLeadPtProf_trans->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtLeadPtProf_trans->BuildOptions(0,0,"i");

    TH1F *hDeltaPt_ptlead_trans = (TH1F*)hDeltaPtLeadPt_trans->ProjectionY("hDeltaPt_ptlead_trans");
    hDeltaPt_ptlead_trans->SetTitle();

    TLegend *legDeltaPt_trans = new TLegend(0.5,0.5,0.88,0.88,"perp plane, p_{T,lead}:");
    legDeltaPt_trans->SetBorderSize(0);
    TH1F *hDeltaPt_ptleadbins_trans[ptbins];
    for(int i=0; i<ptbins; i++){

        hDeltaPt_ptleadbins_trans[i] = (TH1F*)hDeltaPtLeadPt_trans->ProjectionY(Form("hDeltaPt_ptleadbins_trans_%d",i),hDeltaPtLeadPt_trans->GetXaxis()->FindBin(ptleadDown[i]),hDeltaPtLeadPt_trans->GetXaxis()->FindBin(ptleadUp[i])-1);
        hDeltaPt_ptleadbins_trans[i]->SetLineColor(colors[i]);
        hDeltaPt_ptleadbins_trans[i]->SetMarkerColor(colors[i]);
        hDeltaPt_ptleadbins_trans[i]->Scale(1./hDeltaPt_ptleadbins_trans[i]->Integral());
        legDeltaPt_trans->AddEntry(hDeltaPt_ptleadbins_trans[i],Form("%d-%d GeV/c",ptleadDown[i],ptleadUp[i]),"l");
    }

    TH1F *hDeltaPtLeadPtFitMean_trans = new TH1F("hDeltaPtLeadPtFitMean_trans","hDeltaPtLeadPtFitMean_trans",nptbins,ptlead);
    hDeltaPtLeadPtFitMean_trans->SetName("hDeltaPtLeadPtFitMean_trans");
    hDeltaPtLeadPtFitMean_trans->SetTitle();
    hDeltaPtLeadPtFitMean_trans->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFitMean_trans->GetYaxis()->SetRangeUser(-0.6,-0.1);
    TH1F *hDeltaPtLeadPtFitSigma_trans = new TH1F("hDeltaPtLeadPtFitSigma_trans","hDeltaPtLeadPtFitSigma_trans",nptbins,ptlead);
    hDeltaPtLeadPtFitSigma_trans->SetName("hDeltaPtLeadPtFitSigma_trans");
    hDeltaPtLeadPtFitSigma_trans->SetTitle();
    hDeltaPtLeadPtFitSigma_trans->GetXaxis()->SetRangeUser(0,50);
    hDeltaPtLeadPtFitSigma_trans->GetYaxis()->SetRangeUser(0.2,0.8);
    for(int i=0; i<nptbins; i++){
        TH1F* tmp = (TH1F*)hDeltaPtLeadPt_trans->ProjectionY(Form("tmp_%d",i),hDeltaPtLeadPt_trans->GetXaxis()->FindBin(ptlead[i]),hDeltaPtLeadPt_trans->GetXaxis()->FindBin(ptlead[i+1])-1);
        hDeltaPtLeadPtFitMean_trans->SetBinContent(i+1,tmp->GetMean());
        hDeltaPtLeadPtFitMean_trans->SetBinError(i+1,tmp->GetMeanError());
        hDeltaPtLeadPtFitSigma_trans->SetBinContent(i+1,tmp->GetStdDev());
        hDeltaPtLeadPtFitSigma_trans->SetBinError(i+1,tmp->GetStdDevError());
    }

    TH1F* hDeltaPt_ptleadbin5_trans = (TH1F*)hDeltaPtLeadPt_trans->ProjectionY("hDeltaPt_ptleadbin5_trans",hDeltaPtLeadPt_trans->GetXaxis()->FindBin(5),hDeltaPtLeadPt_trans->GetXaxis()->FindBin(100)-1);
    TH1F* hDeltaPt_ptleadbin10_trans = (TH1F*)hDeltaPtLeadPt_trans->ProjectionY("hDeltaPt_ptleadbin10_trans",hDeltaPtLeadPt_trans->GetXaxis()->FindBin(10),hDeltaPtLeadPt_trans->GetXaxis()->FindBin(100)-1);


/*
    // === Delta pT no Bkg =====

    // vs lead. jet pT

    TH2F *hDeltaPtLeadPtNoBkg = (TH2F*)histList1->FindObject("fDeltaPTLeadPtNoBkg");
    hDeltaPtLeadPtNoBkg->SetTitle();
    hDeltaPtLeadPtNoBkg->GetYaxis()->SetRangeUser(-10,100);

    TCanvas *cDetlaPtLeadPtNoBkg = new TCanvas("cDetlaPtLeadPtNoBkg","cDetlaPtLeadPtNoBkg",800,600);
    hDeltaPtLeadPtNoBkg->Draw("colz");
    cDetlaPtLeadPtNoBkg->SaveAs("plotsAlljets/DetlaPtLeadPt.png");

    TProfile *hDeltaPtLeadPtProfNoBkg = (TProfile*)hDeltaPtLeadPtNoBkg->ProfileX("hDeltaPtLeadPtProfNoBkg");
    hDeltaPtLeadPtProfNoBkg->SetTitle();
    //hDeltaPtCentProfNoBkg->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtLeadPtProfNoBkg->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtLeadPtProfNoBkg->BuildOptions(0,0,"i");

    TH1F *hDeltaPt_ptleadNoBkg = (TH1F*)hDeltaPtLeadPtNoBkg->ProjectionY("hDeltaPt_ptleadNoBkg");
    hDeltaPt_ptleadNoBkg->SetTitle();

    TLegend *legDeltaPtNoBkg = new TLegend(0.5,0.5,0.89,0.88,"inc. jet woBkg, p_{T,lead}:");
    legDeltaPtNoBkg->SetBorderSize(0);
    TH1F *hDeltaPt_ptleadbinsNoBkg[ptbins];
    for(int i=0; i<ptbins; i++){

    hDeltaPt_ptleadbinsNoBkg[i] = (TH1F*)hDeltaPtLeadPtNoBkg->ProjectionY(Form("hDeltaPt_ptleadbinsNoBkg_%d",i),hDeltaPtLeadPtNoBkg->GetXaxis()->FindBin(ptleadDown[i]),hDeltaPtLeadPtNoBkg->GetXaxis()->FindBin(ptleadUp[i])-1);
        hDeltaPt_ptleadbinsNoBkg[i]->SetLineColor(colors[i]);
        hDeltaPt_ptleadbinsNoBkg[i]->Scale(1./hDeltaPt_ptleadbinsNoBkg[i]->Integral());
        legDeltaPtNoBkg->AddEntry(hDeltaPt_ptleadbinsNoBkg[i],Form("%d-%d GeV/c",ptleadDown[i],ptleadUp[i]),"l");
    }

    // vs lead. jet pT, excl. lead. jet

    TH2F *hDeltaPtLeadPtNoBkg_excluding = (TH2F*)histList1->FindObject("fDeltaPTLeadPtRawNoBkg_excl_lead");
    hDeltaPtLeadPtNoBkg_excluding->SetTitle();
    hDeltaPtLeadPtNoBkg_excluding->GetYaxis()->SetRangeUser(-10,100);

    TCanvas *cDetlaPtLeadPtNoBkg_excluding = new TCanvas("cDetlaPtLeadPtNoBkg_excluding","cDetlaPtLeadPtNoBkg_excluding",800,600);
    hDeltaPtLeadPtNoBkg_excluding->Draw("colz");
    cDetlaPtLeadPtNoBkg_excluding->SaveAs("plotsAlljets/DetlaPtLeadPtNoBkg_excluding.png");

    TProfile *hDeltaPtLeadPtProfNoBkg_excleading = (TProfile*)hDeltaPtLeadPtNoBkg_excluding->ProfileX("hDeltaPtLeadPtProfNoBkg_excleading");
    hDeltaPtLeadPtProfNoBkg_excleading->SetTitle();
    //hDeltaPtCentProfNoBkg->GetYaxis()->SetRangeUser(0,10);
    hDeltaPtLeadPtProfNoBkg_excleading->GetYaxis()->SetTitle("<#Delta p_{T} (GeV/c)>");
    hDeltaPtLeadPtProfNoBkg_excleading->BuildOptions(0,0,"i");

    TH1F *hDeltaPt_ptleadNoBkg_excluding = (TH1F*)hDeltaPtLeadPtNoBkg_excluding->ProjectionY("hDeltaPt_ptleadNoBkg_excluding");
    hDeltaPt_ptleadNoBkg_excluding->SetTitle();

    TLegend *legDeltaPtNoBkg_excluding = new TLegend(0.5,0.5,0.88,0.88,"inc. jet wo bkg, p_{T,lead}, exc. lead:");
    legDeltaPtNoBkg_excluding->SetBorderSize(0);
    TH1F *hDeltaPt_ptleadbinsNoBkg_excluding[ptbins];
    for(int i=0; i<ptbins; i++){

        hDeltaPt_ptleadbinsNoBkg_excluding[i] = (TH1F*)hDeltaPtLeadPtNoBkg_excluding->ProjectionY(Form("hDeltaPt_ptleadbinsNoBkg_excluding_%d",i),hDeltaPtLeadPtNoBkg_excluding->GetXaxis()->FindBin(ptleadDown[i]),hDeltaPtLeadPtNoBkg_excluding->GetXaxis()->FindBin(ptleadUp[i])-1);
        hDeltaPt_ptleadbinsNoBkg_excluding[i]->SetLineColor(colors[i]);
        hDeltaPt_ptleadbinsNoBkg_excluding[i]->Scale(1./hDeltaPt_ptleadbinsNoBkg_excluding[i]->Integral());
        legDeltaPtNoBkg_excluding->AddEntry(hDeltaPt_ptleadbinsNoBkg_excluding[i],Form("%d-%d GeV/c",ptleadDown[i],ptleadUp[i]),"l");
    }

       TH1F* hDeltaPt_ptleadbin5NoBkg_excluding = (TH1F*)hDeltaPtLeadPtNoBkg_excluding->ProjectionY("hDeltaPt_ptleadbin5NoBkg_excluding",hDeltaPtLeadPtNoBkg_excluding->GetXaxis()->FindBin(5),hDeltaPtLeadPtNoBkg_excluding->GetXaxis()->FindBin(100)-1);
       TH1F* hDeltaPt_ptleadbin10NoBkg_excluding = (TH1F*)hDeltaPtLeadPtNoBkg_excluding->ProjectionY("hDeltaPt_ptleadbin10NoBkg_excluding",hDeltaPtLeadPtNoBkg_excluding->GetXaxis()->FindBin(10),hDeltaPtLeadPtNoBkg_excluding->GetXaxis()->FindBin(100)-1);
*/

/*
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
    cNtracksPtJet->SaveAs("plotsAlljets/NtracksVsJetPt.png");

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

      cNtracks->SaveAs("plotsAlljets/Ntracks.png");
      cNtracks->SaveAs("plotsAlljets/Ntracks.pdf");
*/


    //============================== DRAW ========================================


    TCanvas *cRhoM2 = new TCanvas("cRhoM2","cRhoM2",800,600);
    hRhoMultProf->Draw();
    SaveCanvas(cRhoM2,"RhoVsMultProf");

    TCanvas *cRhoLeadPt2 = new TCanvas("cRhoLeadPt2","cRhoLeadPt2",800,600);
    hRhoLeadPtProf->Draw();
    SaveCanvas(cRhoLeadPt2,"RhoVsLeadPtProf");

    TCanvas *cRhoMultProfJetBins = new TCanvas("cRhoMultProfJetBins","cRhoMultProfJetBins",800,600);
    hRhoMultProfBinsPt[0]->Draw();
    for(int j=1;j<ptbins;j++){
        hRhoMultProfBinsPt[j]->Draw("same");
    }
    leg1->Draw("same");
    SaveCanvas(cRhoMultProfJetBins,"RhoVsMult_leadJetBins");

    TCanvas *cRhoMultProfJetBins2 = new TCanvas("cRhoMultProfJetBins2","cRhoMultProfJetBins2",800,600);
    hRhoMultProfBinsPt2[0]->Draw();
    for(int j=1;j<ptbins;j++){
        hRhoMultProfBinsPt2[j]->Draw("same");
    }
    leg11->Draw("same");
    SaveCanvas(cRhoMultProfJetBins2,"RhoVsMult_leadJetCut");

    TCanvas *cRhoMultProfCentBins = new TCanvas("cRhoMultProfCentBins","cRhoMultProfCentBins",800,600);
    hRhoMultProfBins[0]->Draw();
    for(int j=1;j<centBins;j++){
        hRhoMultProfBins[j]->Draw("same");
    }
    leg2->Draw("same");
    SaveCanvas(cRhoMultProfCentBins,"RhoVsLeadJetPt_centBins");

  // ===== Delta pT ========

    TLegend *ldelta1 = new TLegend(0.6,0.7,0.85,0.85);
    ldelta1->SetBorderSize(0);
    //ldelta1->AddEntry(hDeltaPt_c,"All jets","l");
    ldelta1->AddEntry(hDeltaPt_c_excleading,"Excl. lead. jet","l");
    ldelta1->AddEntry(hDeltaPt_c_trans,"perp. plane","l");

    // vs cent

    // profile
  //  hDeltaPtCentProf->GetXaxis()->SetRangeUser(-10,60);
  //  hDeltaPtCentProf->SetLineColor(1);
    hDeltaPtCentProf_excleading->GetXaxis()->SetRangeUser(-10,60);
    hDeltaPtCentProf_excleading->SetLineColor(2);
    hDeltaPtCentProf_trans->SetLineColor(4);
    hDeltaPtCentProf_excleading->SetMarkerColor(2);
    hDeltaPtCentProf_trans->SetMarkerColor(4);
    TCanvas *cDeltaPtCentProf = new TCanvas("cDeltaPtCentProf","cDeltaPtCentProf",800,600);
    //hDeltaPtCentProf->Draw();
    hDeltaPtCentProf_excleading->Draw();
    hDeltaPtCentProf_trans->Draw("same");
    ldelta1->Draw("same");
    SaveCanvas(cDeltaPtCentProf,"DeltaPtCentProf");

   // 1D
    //hDeltaPt_c->GetXaxis()->SetRangeUser(-10,60);
    //hDeltaPt_c->SetLineColor(1);
    hDeltaPt_c_excleading->GetXaxis()->SetRangeUser(-10,60);
    hDeltaPt_c_excleading->SetLineColor(2);
    hDeltaPt_c_trans->SetLineColor(4);
    hDeltaPt_c_excleading->SetMarkerColor(2);
    hDeltaPt_c_trans->SetMarkerColor(4);
    TCanvas *cDeltaPt_c = new TCanvas("cDeltaPt_c","cDeltaPt_c",800,600);
    cDeltaPt_c->SetLogy();
    //hDeltaPt_c->Draw();
    hDeltaPt_c_excleading->Draw();
    hDeltaPt_c_trans->Draw("same");
    ldelta1->Draw("same");
    SaveCanvas(cDeltaPt_c,"DeltaPt_exclComparison");

    hDeltaPt_ptleadbin5_excluding->GetXaxis()->SetRangeUser(-10,60);
    hDeltaPt_ptleadbin5_excluding->SetLineColor(2);
    hDeltaPt_ptleadbin5_trans->SetLineColor(4);
    TCanvas *cDeltaPt_ptleadbin5 = new TCanvas("cDeltaPt_ptleadbin5","cDeltaPt_ptleadbin5",800,600);
    cDeltaPt_ptleadbin5->SetLogy();
    //hDeltaPt_c->Draw();
    hDeltaPt_ptleadbin5_excluding->Draw();
    hDeltaPt_ptleadbin5_trans->Draw("same");
    ldelta1->Draw("same");
    SaveCanvas(cDeltaPt_ptleadbin5,"DeltaPt_exclComparison_ptlead5");

    hDeltaPt_ptleadbin10_excluding->GetXaxis()->SetRangeUser(-10,60);
    hDeltaPt_ptleadbin10_excluding->SetLineColor(2);
    hDeltaPt_ptleadbin10_trans->SetLineColor(4);
    hDeltaPt_ptleadbin10_excluding->SetMarkerColor(2);
    hDeltaPt_ptleadbin10_trans->SetMarkerColor(4);
    TCanvas *cDeltaPt_ptleadbin10 = new TCanvas("cDeltaPt_ptleadbin10","cDeltaPt_ptleadbin10",800,600);
    cDeltaPt_ptleadbin10->SetLogy();
    //hDeltaPt_c->Draw();
    hDeltaPt_ptleadbin10_excluding->Draw();
    hDeltaPt_ptleadbin10_trans->Draw("same");
    ldelta1->Draw("same");
    SaveCanvas(cDeltaPt_ptleadbin10,"DeltaPt_exclComparison_ptlead10");

    //==== vs lead jet pT

//    TCanvas *cDeltaPtLeadPt = new TCanvas("cDeltaPtLeadPt","cDeltaPtLeadPt",800,600);
//    hDeltaPt_ptlead->Draw();
//    cDeltaPtLeadPt->SaveAs("plotsAlljets/DeltaPtLead.png");

    // profile
  /*  hDeltaPtLeadPtProf->Rebin(4);
    hDeltaPtLeadPtProf_excleading->Rebin(4);
    hDeltaPtLeadPtProf->SetLineColor(1);*/

    hDeltaPtLeadPtProf_excleading->SetLineColor(2);
    hDeltaPtLeadPtProf_trans->SetLineColor(4);
    hDeltaPtLeadPtProf_excleading->SetMarkerColor(2);
    hDeltaPtLeadPtProf_trans->SetMarkerColor(4);
    hDeltaPtLeadPtProf_excleading->GetYaxis()->SetRangeUser(0,5);
    hDeltaPtLeadPtProf_excleading->GetXaxis()->SetRangeUser(0,60);


    TCanvas *cDeltaPtLeadPtProf = new TCanvas("cDeltaPtLeadPtProf","cDeltaPtLeadPtProf",800,600);
    //hDeltaPtLeadPtProf->Draw();
    hDeltaPtLeadPtProf_excleading->Draw();
    hDeltaPtLeadPtProf_trans->Draw("same");
    ldelta1->Draw("same");
    SaveCanvas(cDeltaPtLeadPtProf,"DeltaPtLeadPtProf");

    // jet pT bins

  /*  TCanvas *cDeltaPtLeadPtBins = new TCanvas("cDeltaPtLeadPtBins","cDeltaPtLeadPtBins",800,600);
    cDeltaPtLeadPtBins->SetLogy();
    hDeltaPt_ptleadbins[0]->SetMaximum(hDeltaPt_ptleadbins[0]->GetMaximum()*2);
    hDeltaPt_ptleadbins[0]->GetXaxis()->SetRangeUser(-10,60);
    hDeltaPt_ptleadbins[0]->Draw();
    for(int j=1;j<ptbins;j++){
        hDeltaPt_ptleadbins[j]->Draw("same");
    }
    legDeltaPt->Draw("same");
    cDeltaPtLeadPtBins->SaveAs("plotsAlljets/DeltaPtLeadPtBins.png");*/


    TCanvas *cDeltaPtLeadPtBins_excluding = new TCanvas("cDeltaPtLeadPtBins_excluding","cDeltaPtLeadPtBins_excluding",800,600);
    cDeltaPtLeadPtBins_excluding->SetLogy();
    hDeltaPt_ptleadbins_excluding[0]->SetMaximum(hDeltaPt_ptleadbins_excluding[0]->GetMaximum()*2);
    hDeltaPt_ptleadbins_excluding[0]->GetXaxis()->SetRangeUser(-10,60);
    hDeltaPt_ptleadbins_excluding[0]->Draw();
    for(int j=1;j<centBins;j++){
        hDeltaPt_ptleadbins_excluding[j]->Draw("same");
    }
    legDeltaPt_excluding->Draw("same");
    SaveCanvas(cDeltaPtLeadPtBins_excluding,"DeltaPtLeadPtBins_excluding");

    TCanvas *cDeltaPtLeadPtBins_trans = new TCanvas("cDeltaPtLeadPtBins_trans","cDeltaPtLeadPtBins_trans",800,600);
    cDeltaPtLeadPtBins_trans->SetLogy();
    hDeltaPt_ptleadbins_trans[0]->SetMaximum(hDeltaPt_ptleadbins_trans[0]->GetMaximum()*2);
    hDeltaPt_ptleadbins_trans[0]->GetXaxis()->SetRangeUser(-10,60);
    hDeltaPt_ptleadbins_trans[0]->Draw();
    for(int j=1;j<centBins;j++){
        hDeltaPt_ptleadbins_trans[j]->Draw("same");
    }
    legDeltaPt_trans->Draw("same");
    SaveCanvas(cDeltaPtLeadPtBins_trans,"DeltaPtLeadPtBins_trans");

/*
    //==== vs lead jet pT, no bkg

    TCanvas *cDeltaPtLeadPtNoBkg = new TCanvas("cDeltaPtLeadPtNoBkg","cDeltaPtLeadPtNoBkg",800,600);
    hDeltaPt_ptleadNoBkg->Draw();
    cDeltaPtLeadPtNoBkg->SaveAs("plotsAlljets/DeltaPtLeadNoBkg.png");

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
    cDeltaPtLeadPtProfNoBkg->SaveAs("plotsAlljets/DeltaPtLeadPtProfNoBkg.png");

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
    cDeltaPtLeadPtBinsNoBkg->SaveAs("plotsAlljets/DeltaPtLeadPtBinsNoBkg.png");


    TCanvas *cDeltaPtLeadPtBinsNoBkg_excluding = new TCanvas("cDeltaPtLeadPtBinsNoBkg_excluding","cDeltaPtLeadPtBinsNoBkg_excluding",800,600);
    cDeltaPtLeadPtBinsNoBkg_excluding->SetLogy();
    hDeltaPt_ptleadbinsNoBkg_excluding[0]->SetMaximum(hDeltaPt_ptleadbins_excluding[0]->GetMaximum()*2);
    hDeltaPt_ptleadbinsNoBkg_excluding[0]->GetXaxis()->SetRangeUser(-2,30);
    hDeltaPt_ptleadbinsNoBkg_excluding[0]->Draw();
    for(int j=1;j<centBins;j++){
        hDeltaPt_ptleadbinsNoBkg_excluding[j]->Draw("same");
    }
    legDeltaPt_excluding->Draw("same");
    cDeltaPtLeadPtBinsNoBkg_excluding->SaveAs("plotsAlljets/DeltaPtLeadPtBinsNoBkg_excluding.png");


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
    cDeltaPtLeadPtProf_bkgCOmparison->SaveAs("plotsAlljets/DeltaPtLeadPtProf_bkgComparison.png");
*/

TFile *fOut = new TFile(Form("fOut%s.root",outFile.Data()),"RECREATE");
hPtTrk->Write();
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

//hDeltaPtCent->Write();
//hDeltaPt_c->Write();
hDeltaPtCent_excleading->Write();
hDeltaPt_c_excleading->Write();
hDeltaPtCent_trans->Write();
hDeltaPt_c_trans->Write();
//hDeltaPtLeadPt->Write();
//hDeltaPt_ptlead->Write();
//hDeltaPt_ptleadbin5->Write();
//hDeltaPt_ptleadbin10->Write();
hDeltaPtLeadPt_excluding->Write();
hDeltaPt_ptlead_excluding->Write();
hDeltaPt_ptleadbin5_excluding->Write();
hDeltaPt_ptleadbin10_excluding->Write();
hDeltaPtLeadPt_trans->Write();
hDeltaPt_ptlead_trans->Write();
hDeltaPt_ptleadbin5_trans->Write();
hDeltaPt_ptleadbin10_trans->Write();
//hDeltaPtLeadPtNoBkg->Write();
//hDeltaPt_ptleadNoBkg->Write();
//hDeltaPtLeadPtNoBkg_excluding->Write();
//hDeltaPt_ptleadNoBkg_excluding->Write();
//hDeltaPt_ptleadbin5NoBkg_excluding->Write();
//hDeltaPt_ptleadbin10NoBkg_excluding->Write();

hRhoMultProf->Write();
hRhoLeadPtProf->Write();
//hDeltaPtCentProf->Write();
hDeltaPtCentProf_excleading->Write();
hDeltaPtCentProf_trans->Write();
//hDeltaPtLeadPtProf->Write();
hDeltaPtLeadPtProf_excleading->Write();
hDeltaPtLeadPtProf_trans->Write();

//hDeltaPtLeadPtProfNoBkg->Write();
//hDeltaPtLeadPtProfNoBkg_excleading->Write();

//hDeltaPtLeadPtFit_excluding->Write();
hDeltaPtLeadPtFitMean_excluding->Write();
hDeltaPtLeadPtFitSigma_excluding->Write();

//hDeltaPtLeadPtFit_trans->Write();
hDeltaPtLeadPtFitMean_trans->Write();
hDeltaPtLeadPtFitSigma_trans->Write();

    for(int i=0; i<ptbins; i++){
        hDeltaPt_ptleadbins_excluding[i]->Write();
        hDeltaPt_ptleadbins_trans[i]->Write();
    }


fOut->Close();


}

void SaveCanvas(TCanvas *c, TString name = "tmp"){

    c->SaveAs(Form("%s/%s.png",outDir.Data(),name.Data()));
    c->SaveAs(Form("%s/%s.pdf",outDir.Data(),name.Data()));
}
