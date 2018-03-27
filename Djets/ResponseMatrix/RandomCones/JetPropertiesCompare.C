#include <string>
#include <sstream>
#include <iostream>


TString outDir = "plotsComparison";

const int nfiles = 3;
TString datafile[nfiles] = {"fOutDzeroEvt.root","fOutDzeroLeadEvt.root","fOutAllEvt.root"};
TFile *File[nfiles];

const int nhistos = 3;
const int nhistoComp = 4;
TString histosExl[nhistos] = {"hDeltaPt_ptlead_excluding","hDeltaPt_ptleadbin5_excluding","hDeltaPt_ptleadbin10_excluding"};
TString histosTrans[nhistos] = {"hDeltaPt_ptlead_trans","hDeltaPt_ptleadbin5_trans","hDeltaPt_ptleadbin10_trans"};
TString histosComp[nhistoComp] = {"hDeltaPt_ptleadbin5_excluding","hDeltaPt_ptleadbin5_trans","hDeltaPt_ptlead_excluding","hDeltaPt_ptlead_trans"};
TString histosComp2[nhistoComp] = {"hDeltaPt_ptleadbin5_excluding","hDeltaPt_ptleadbin5_trans","hDeltaPt_ptleadbin10_excluding","hDeltaPt_ptleadbin10_trans"};
TString legDsc[] = {"D-jet evt", "D-jet evt, lead. D", "incl.jet evt" };

TString legDsc2[] = {"DjetEvt", "DjetEvt_LeadD", "IncljetEvt" };
TString legDscExl[] = {"All", "p_{T}>5 GeV/#it{c}", "p_{T}>10 GeV/#it{c}" };
TString legDscComp[] = {"exl.lead jet,p_{T}>5 GeV/#it{c}","perp.plane,p_{T}>5 GeV/#it{c}","exl.lead jet","perp.plane"};
TString legDscComp2[] = {"exl.lead jet,p_{T}>5 GeV/#it{c}","perp.plane,p_{T}>5 GeV/#it{c}","exl.lead jet,p_{T}>10 GeV/#it{c}","perp.plane,p_{T}>10 GeV/#it{c}"};

Int_t colors[] = {2,6,4,8};
Float_t deltamax = 40;

void JetPropertiesCompare()
{

  gStyle->SetOptStat(0000); //Mean and RMS shown
  gSystem->Exec(Form("mkdir %s",outDir.Data()));

  if( !setFiles() ) { return ;}


    getHisto("hDeltaPt_ptlead_excluding","#Delta p_{T},exc. lead.",1,1);
    getHisto("hDeltaPt_ptleadbin5_excluding","#Delta p_{T} exc. lead., p_{T}>5 GeV/c",1,1);
    getHisto("hDeltaPt_ptleadbin10_excluding","#Delta p_{T} exc. lead., p_{T}>10 GeV/c",1,1);

    getHisto("hDeltaPt_ptlead_trans","#Delta p_{T},perp.",1,1);
    getHisto("hDeltaPt_ptleadbin5_trans","#Delta p_{T} perp., p_{T}>5 GeV/c",1,1);
    getHisto("hDeltaPt_ptleadbin10_trans","#Delta p_{T} perp., p_{T}>10 GeV/c",1,1);

    for(int i=0; i<nfiles; i++) {
      getHistoFile(i,nhistos,histosExl,"Excl.Lead");
      getHistoFile(i,nhistos,histosTrans,"Perp.Plane");
      getHistoComp(i,nhistoComp,histosComp,"comp");
      getHistoComp2(i,nhistoComp,histosComp2,"com2");
    }


        //  getHisto("hDeltaPt_ptlead","#Delta p_{T}",1,1);
  //  getHisto("hDeltaPt_ptleadbin5","#Delta p_{T}, p_{T}>5 GeV/c",1,1);

  //  getHisto("hDeltaPt_ptleadbin10","#Delta p_{T}, p_{T}>10 GeV/c",1,1);


/*    getHistoDjet("hDeltaPt_ptleadbin5_exlcuding","#Delta p_{T} exc. lead., p_{T}>5 GeV/c",1,1);
    getHistoDjet("hDeltaPt_ptleadbin10_exlcuding","#Delta p_{T} exc. lead., p_{T}>10 GeV/c",1,1);


    getHisto("hDeltaPtLeadPtFitMean_exlcuding","<#Delta p_{T}> exc. lead., mean",0,0);
    getHisto("hDeltaPtLeadPtFitSigma_exlcuding","<#Delta p_{T}> exc. lead., sigma",0,0);


    getHisto("hRho_ptleadbin5","#rho, p_{T,lead}>5 GeV/c",1,1);
    getHisto("hRho_ptleadbin10","#rho, p_{T,lead}>10 GeV/c",1,1);
    getHisto("hRho_ptleadbin5_cent0_10","#rho, p_{T,lead}>5 GeV/c, cent 0-10%",1,1);
    getHisto("hRho_ptleadbin5_cent10_20","#rho, p_{T,lead}>5 GeV/c, cent 10-20%",1,1);
    getHisto("hRho_ptleadbin5_cent20_40","#rho, p_{T,lead}>5 GeV/c, cent 20-40%",1,1);
    getHisto("hRho_ptleadbin5_cent40_60","#rho, p_{T,lead}>5 GeV/c, cent 40-60%",1,1);
    getHisto("hRho_ptleadbin5_cent60_100","#rho, p_{T,lead}>5 GeV/c, cent 60-100%",1,1);

    getProfile("hRhoMultProf","#rho, D-jet/inc.jet ratio",0,0);
    getProfile("hRhoLeadPtProf","#rho, D-jet/inc.jet ratio",0,0,1);
    getProfile("hRhoMultProf_leadpt5","#rho, D-jet/inc.jet ratio,p_{T,lead}>5 GeV/c",0,0,0);
    getProfile("hRhoMultProf_leadpt10","#rho, D-jet/inc.jet ratio,p_{T,lead}>10 GeV/c",0,0,0);

    getProfile("hDeltaPtLeadPtProf_excleading","#delta p_{T}, D-jet/inc.jet ratio",0,0,1);
    getProfile("hDeltaPtLeadPtProfNoBkg_excleading","#delta p_{T}, D-jet/inc.jet ratio",0,0,1);*/



return;

}

int setFiles(){

  for(int i=0; i<nfiles; i++) {
    File[i] = new TFile(datafile[i],"read");
    if(!File[i]) { cout<< "=== File: " << datafile[i] << " does not exist !!!" << endl; return 0; }
  }

  return 1;
}

void getHisto(TString name, TString desc, int scale = 1, int log=0){

  TLegend *leg = new TLegend(0.6,0.63,0.89,0.78);
  leg->SetBorderSize(0);
  TPaveText *pv[nfiles];

  TH1F *hist[nfiles];
  double mean[nfiles], meanErr[nfiles];
  double sigma[nfiles], sigmaErr[nfiles];
  double shift=0;
  for(int i=0; i<nfiles;i++){
    hist[i] = (TH1F*)File[i]->Get(name.Data());
    hist[i]->SetLineColor(colors[i]);
    hist[i]->SetMarkerColor(colors[i]);
    hist[i]->GetXaxis()->SetRangeUser(-10,deltamax);
    leg->AddEntry(hist[i],Form("%s",legDsc[i].Data()),"l");
    if(scale){
        hist[i]->Scale(1./hist[i]->Integral());
    }
    mean[i] = hist[i]->GetMean();
    meanErr[i] = hist[i]->GetMeanError();
    sigma[i] = hist[i]->GetStdDev();
    sigmaErr[i] = hist[i]->GetStdDevError();

    pv[i] = new TPaveText(0.65,0.5-shift,0.89,0.6-shift,"brNDC");
    pv[i]->SetFillStyle(0);
    pv[i]->SetBorderSize(0);
    pv[i]->SetTextColor(colors[i]);
    pv[i]->AddText(Form("mean: %.2f",mean[i]));
    pv[i]->AddText(Form("sigma: %.2f",sigma[i]));
    shift+=0.1;
  }

    TLegend *ll = new TLegend(0.45,0.8,0.89,0.85,desc.Data());
    ll->SetBorderSize(0);

    TCanvas *cc = new TCanvas("cc","cc",800,600);
    cc->SetLogy(log);
    hist[0]->Draw();
    pv[0]->Draw("same");
      for(int i=1; i<nfiles;i++){
        hist[i]->Draw("same");
        pv[i]->Draw("same");
      }
    leg->Draw("same");
    ll->Draw("same");

    cc->SaveAs(Form("plotsComparison/%s.png",name.Data()));

  /*  TH1F *ratio = (TH1F*)hD->Clone("ratio");
    ratio->Divide(hAll);

    ratio->SetLineColor(kGreen+2);
    ratio->SetMarkerColor(kGreen+2);
    ratio->SetMinimum(0.5);
    ratio->SetMaximum(3);

    TCanvas *cc2 = new TCanvas("cc2","cc2",800,400);
    ratio->Draw();
    //leg->Draw("same");
    //ll->Draw("same");
    //if(ispt)line->Draw("same");

    cc2->SaveAs(Form("plotsComparison/%s_ratio.png",name.Data()));*/
}

void getHistoFile(Int_t it, Int_t n, TString *histos, TString desc, int scale = 1, int log=1){

  TLegend *leg = new TLegend(0.6,0.63,0.89,0.78,desc.Data());
  leg->SetBorderSize(0);
  TPaveText *pv[nhistos];

  TH1F *hist[nhistos];
  double mean[nhistos], meanErr[nhistos];
  double sigma[nhistos], sigmaErr[nhistos];
  double shift=0;
  for(int i=0; i<n;i++){
    hist[i] = (TH1F*)File[it]->Get(histos[i].Data());
    hist[i]->SetLineColor(colors[i]);
    hist[i]->SetMarkerColor(colors[i]);
    hist[i]->GetXaxis()->SetRangeUser(-10,deltamax);
    leg->AddEntry(hist[i],Form("%s",legDscExl[i].Data()),"l");
    if(scale){
        hist[i]->Scale(1./hist[i]->Integral());
    }
    mean[i] = hist[i]->GetMean();
    meanErr[i] = hist[i]->GetMeanError();
    sigma[i] = hist[i]->GetStdDev();
    sigmaErr[i] = hist[i]->GetStdDevError();

    pv[i] = new TPaveText(0.65,0.5-shift,0.89,0.6-shift,"brNDC");
    pv[i]->SetFillStyle(0);
    pv[i]->SetBorderSize(0);
    pv[i]->SetTextColor(colors[i]);
    pv[i]->AddText(Form("mean: %.2f",mean[i]));
    pv[i]->AddText(Form("sigma: %.2f",sigma[i]));
    shift+=0.1;
  }

    TLegend *ll = new TLegend(0.55,0.8,0.89,0.85,legDsc[it].Data());
    ll->SetBorderSize(0);

    TCanvas *cc = new TCanvas("cc","cc",800,600);
    cc->SetLogy(log);
    hist[0]->Draw();
    pv[0]->Draw("same");
      for(int i=1; i<n;i++){
        hist[i]->Draw("same");
        pv[i]->Draw("same");
      }
    leg->Draw("same");
    ll->Draw("same");

    cc->SaveAs(Form("plotsComparison/%s_%s.png",legDsc2[it].Data(),desc.Data()));
}

void getHistoComp(Int_t it, Int_t n, TString *histos, TString desc, int scale = 1, int log=1){

  TLegend *leg = new TLegend(0.6,0.63,0.89,0.78);
  leg->SetBorderSize(0);
  TPaveText *pv[nhistoComp];

  TH1F *hist[nhistoComp];
  double mean[nhistoComp], meanErr[nhistoComp];
  double sigma[nhistoComp], sigmaErr[nhistoComp];
  double shift=0;
  for(int i=0; i<n;i++){
    hist[i] = (TH1F*)File[it]->Get(histos[i].Data());
    hist[i]->SetLineColor(colors[i]);
    hist[i]->SetMarkerColor(colors[i]);
    hist[i]->GetXaxis()->SetRangeUser(-10,deltamax);
    leg->AddEntry(hist[i],Form("%s",legDscComp[i].Data()),"l");
    if(scale){
        hist[i]->Scale(1./hist[i]->Integral());
    }
    mean[i] = hist[i]->GetMean();
    meanErr[i] = hist[i]->GetMeanError();
    sigma[i] = hist[i]->GetStdDev();
    sigmaErr[i] = hist[i]->GetStdDevError();

    pv[i] = new TPaveText(0.65,0.5-shift,0.89,0.6-shift,"brNDC");
    pv[i]->SetFillStyle(0);
    pv[i]->SetBorderSize(0);
    pv[i]->SetTextColor(colors[i]);
    pv[i]->AddText(Form("mean: %.2f",mean[i]));
    pv[i]->AddText(Form("sigma: %.2f",sigma[i]));
    shift+=0.1;
  }

    TLegend *ll = new TLegend(0.55,0.8,0.89,0.85,legDsc[it].Data());
    ll->SetBorderSize(0);

    TCanvas *cc = new TCanvas("cc","cc",800,600);
    cc->SetLogy(log);
    hist[0]->Draw();
    pv[0]->Draw("same");
      for(int i=1; i<n;i++){
        hist[i]->Draw("same");
        pv[i]->Draw("same");
      }
    leg->Draw("same");
    ll->Draw("same");

    cc->SaveAs(Form("plotsComparison/%s_%s.png",legDsc2[it].Data(),desc.Data()));
}

void getHistoComp2(Int_t it, Int_t n, TString *histos, TString desc, int scale = 1, int log=1){

  TLegend *leg = new TLegend(0.6,0.63,0.89,0.78);
  leg->SetBorderSize(0);
  TPaveText *pv[nhistoComp];

  TH1F *hist[nhistoComp];
  double mean[nhistoComp], meanErr[nhistoComp];
  double sigma[nhistoComp], sigmaErr[nhistoComp];
  double shift=0;
  for(int i=0; i<n;i++){
    hist[i] = (TH1F*)File[it]->Get(histos[i].Data());
    hist[i]->SetLineColor(colors[i]);
    hist[i]->SetMarkerColor(colors[i]);
    hist[i]->GetXaxis()->SetRangeUser(-10,deltamax);
    leg->AddEntry(hist[i],Form("%s",legDscComp2[i].Data()),"l");
    if(scale){
        hist[i]->Scale(1./hist[i]->Integral());
    }
    mean[i] = hist[i]->GetMean();
    meanErr[i] = hist[i]->GetMeanError();
    sigma[i] = hist[i]->GetStdDev();
    sigmaErr[i] = hist[i]->GetStdDevError();

    pv[i] = new TPaveText(0.65,0.5-shift,0.89,0.6-shift,"brNDC");
    pv[i]->SetFillStyle(0);
    pv[i]->SetBorderSize(0);
    pv[i]->SetTextColor(colors[i]);
    pv[i]->AddText(Form("mean: %.2f",mean[i]));
    pv[i]->AddText(Form("sigma: %.2f",sigma[i]));
    shift+=0.1;
  }

    TLegend *ll = new TLegend(0.55,0.8,0.89,0.85,legDsc[it].Data());
    ll->SetBorderSize(0);

    TCanvas *cc = new TCanvas("cc","cc",800,600);
    cc->SetLogy(log);
    hist[0]->Draw();
    pv[0]->Draw("same");
      for(int i=1; i<n;i++){
        hist[i]->Draw("same");
        pv[i]->Draw("same");
      }
    leg->Draw("same");
    ll->Draw("same");

    cc->SaveAs(Form("plotsComparison/%s_%s.png",legDsc2[it].Data(),desc.Data()));
}


void getProfile(TString name, TString desc, int scale = 1, int log=0, int ispt=0){

    char *datafile1 = "fOutDjetProperties.root";
    char dataFile1[200];
	sprintf(dataFile1,datafile1);
	TFile *File1 = new TFile(dataFile1,"read");

    char *datafile2 = "fOutAlljetProperties.root";
    char dataFile2[200];
	sprintf(dataFile2,datafile2);
	TFile *File2 = new TFile(dataFile2,"read");


    TProfile *hD = (TProfile*)File1->Get(name.Data());
    TProfile *hAll = (TProfile*)File2->Get(name.Data());
    hD->SetLineColor(2);
    hAll->SetLineColor(4);
    if(scale){
        hD->Scale(1./hD->Integral());
        hAll->Scale(1./hAll->Integral());
    }

    TLegend *leg = new TLegend(0.6,0.6,0.89,0.75);
    leg->SetBorderSize(0);
    leg->AddEntry(hD,"D-jet","l");
    leg->AddEntry(hAll,"inc.jet","l");

    TLegend *ll = new TLegend(0.45,0.8,0.89,0.85,desc.Data());
    ll->SetBorderSize(0);

    TCanvas *cc = new TCanvas("cc","cc",800,600);
    cc->SetLogy(log);
    hD->Draw();
    hAll->Draw("same");
    leg->Draw("same");
    ll->Draw("same");

    cc->SaveAs(Form("plotsComparison/%s.png",name.Data()));


    TH1D *h1 = hD->ProjectionX("h1");
    TH1D *h2 = hAll->ProjectionX("h2");

    TH1D *ratio = (TH1D*)h1->Clone("ratio");
    ratio->Divide(h2);

    //TProfile *ratio = (TProfile*)hD->Clone("ratio");
    //ratio->BuildOptions(0,0,"i");
    //ratio->Divide(hAll);

    ratio->SetLineColor(kGreen+2);
    ratio->SetMarkerColor(kGreen+2);
    ratio->SetMaximum(3);
    ratio->SetMinimum(0);
    if(ispt) ratio->GetXaxis()->SetRangeUser(0,50);

    for(int i=0; i<=ratio->GetNbinsX();i++) {

         //ratio->SetBinError(i,0);
    }


    TLine *line = new TLine(5,0.2,5,2.5);
    line->SetLineStyle(2);
    line->SetLineWidth(2);

    if(ispt){
    TF1 *fit = new TF1("fit","[0]",10,80);
    fit->SetLineColor(1);
    fit->SetLineStyle(3);
    ratio->Fit("fit","RM");
    }

    TCanvas *cc2 = new TCanvas("cc2","cc2",800,600);
    ratio->Draw();
    //leg->Draw("same");
    ll->Draw("same");
    if(ispt)line->Draw("same");

    cc2->SaveAs(Form("plotsComparison/%s_ratio.png",name.Data()));


}


void getHistoDjet(TString name, TString desc, int scale = 1, int log=0){

    char *datafile1 = "fOutDjetProperties.root";
    char dataFile1[200];
	sprintf(dataFile1,datafile1);
	TFile *File1 = new TFile(dataFile1,"read");

    char *datafile2 = "fOutAlljetProperties.root";
    char dataFile2[200];
	sprintf(dataFile2,datafile2);
	TFile *File2 = new TFile(dataFile2,"read");


    TH1F *hD = (TH1F*)File1->Get(name.Data());
    TH1F *hD_2 = (TH1F*)File1->Get(Form("%sD", name.Data()));
    TH1F *hAll = (TH1F*)File2->Get(name.Data());
    hD->SetLineColor(2);
    hD_2->SetLineColor(kGreen+1);
    hAll->SetLineColor(4);
    hD->SetMarkerColor(2);
    hD_2->SetMarkerColor(kGreen+1);
    hAll->SetMarkerColor(4);
    if(scale){
        hD->Scale(1./hD->Integral());
        hD_2->Scale(1./hD_2->Integral());
        hAll->Scale(1./hAll->Integral());
    }

    TLegend *leg = new TLegend(0.6,0.6,0.89,0.75);
    leg->SetBorderSize(0);
    leg->AddEntry(hD,"D-jet, excl. lead","l");
    leg->AddEntry(hD_2,"D-jet, excl. lead D","l");
    leg->AddEntry(hAll,"inc.jet, excl. lead","l");

    TLegend *ll = new TLegend(0.45,0.8,0.89,0.85,desc.Data());
    ll->SetBorderSize(0);


    TCanvas *cc = new TCanvas("cc","cc",800,600);
    cc->SetLogy(log);
    hD->Draw();
    hD_2->Draw("same");
    hAll->Draw("same");
    leg->Draw("same");
    ll->Draw("same");

    cc->SaveAs(Form("plotsComparison/%s_exlD.png",name.Data()));


}
