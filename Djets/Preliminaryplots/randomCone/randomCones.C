#include "style.C"
#include <string>
#include <sstream>
#include <iostream>

setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Width_t width);

 double zmin = 0, zmax = 2.;
    double jetmin = 2, jetmax = 40;
    double plotmin = 2, plotmax = 36;


    const int ptbinsJetN = 8;
    float ptJetbins[ptbinsJetN+1] = { 2,4,6,8,10,12,16,24,50 };

    int promptColor = 1;
    int nonpromptColor = kBlue+1;


void randomCones( int Rpar = 4 )
{

    style();
    gStyle->SetOptStat(000);

    gStyle->SetLegendFont(42);
    //gStyle->SetLegendTextSize(0.05);


stringstream sst;
sst.clear(); sst.str("");


   TFile *inFilePrompt = new TFile("/home/basia/Work/alice/analysis/alice_Djets/Djets/ResponseMatrix/BkgRM03/RandCones_BkgM_Djet5Excl.root","read");

   //TH1F *hEffPrompt = (TH1F*)inFilePrompt->Get("hDeltaPt_ptleadbin5_exlcuding");
    TH1F *hEffPrompt = (TH1F*)inFilePrompt->Get("hDeltaPt2");
    //hEffPrompt->Rebin(2);
    //hEffPrompt->Scale(1.,"width");

    hEffPrompt->SetTitle();

            hEffPrompt->SetMarkerColor(promptColor);
            hEffPrompt->SetLineColor(promptColor);
            hEffPrompt->SetMarkerStyle(20);
            hEffPrompt->SetMarkerSize(1.2);

            hEffPrompt->GetXaxis()->SetTitle("#delta #it{p}_{T, ch} (GeV/#it{c})");
            hEffPrompt->GetYaxis()->SetTitle("Probability density");
            hEffPrompt->GetXaxis()->SetLabelSize(0.04);
            hEffPrompt->GetXaxis()->SetTitleSize(0.05);
            hEffPrompt->GetXaxis()->SetTitleOffset(1.);
            hEffPrompt->GetYaxis()->SetLabelSize(0.045);
            hEffPrompt->GetYaxis()->SetTitleSize(0.05);
            hEffPrompt->GetXaxis()->SetRangeUser(-7,25);
            hEffPrompt->SetMaximum(hEffPrompt->GetMaximum()*40);


            TH1F *hEffPromptFill = (TH1F*)hEffPrompt->Clone("hEffPromptFill");
            hEffPromptFill->SetFillStyle(1001);
            hEffPromptFill->SetFillColor(kGray);


/*
    TLegend *leg = new TLegend(0.5,0.25,0.85,0.40);
    leg->SetTextSize(0.045);
    leg->AddEntry(hEffPrompt,"Prompt D^{*+}","p");
    leg->AddEntry(hEffNonPrompt,"Feed-down D^{*+}","p");
   */

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
   // pvEn->AddText("PYTHIA6+HIJING, p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
    pvEn->AddText("p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");

    TPaveText *pvCent= new TPaveText(0.15,0.74,0.8,0.79,"brNDC");
    pvCent->SetFillStyle(0);
    pvCent->SetBorderSize(0);
    pvCent->SetTextFont(42);
    pvCent->SetTextSize(0.045);
    pvCent->SetTextAlign(11);
   // pvEn->AddText("PYTHIA6+HIJING, p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
    pvCent->AddText("Minimum bias");

    double shift = 0.;

     TPaveText *pvEta = new TPaveText(0.6,0.72,0.8,0.82,"brNDC");
    pvEta->SetFillStyle(0);
    pvEta->SetBorderSize(0);
    pvEta->SetTextFont(42);
    pvEta->SetTextSize(0.045);
    pvEta->SetTextAlign(11);
    //pvEta->AddText("|#it{#eta}_{lab}| < 0.9, #it{p}_{T, track} > 0.15 GeV/#it{c}");
    pvEta->AddText("#splitline{|#it{#eta}_{lab}| < 0.9}{#it{p}_{T, track} > 0.15 GeV/#it{c}}");


    TPaveText *pvJet = new TPaveText(0.6,0.6,0.8,0.7,"brNDC");
    pvJet->SetFillStyle(0);
    pvJet->SetBorderSize(0);
    pvJet->SetTextFont(42);
    pvJet->SetTextSize(0.045);
    pvJet->SetTextAlign(11);
    //pvJet->AddText("#it{R} = 0.4");
    pvJet->AddText("#splitline{#it{R} = 0.3}{Random Cones}");


     TPaveText *pvD = new TPaveText(0.6,0.48,0.8,0.59,"brNDC");
    pvD->SetFillStyle(0);
    pvD->SetBorderSize(0);
    pvD->SetTextFont(42);
    pvD->SetTextSize(0.045);
    pvD->SetTextAlign(11);
    pvD->AddText("#splitline{with D^{*+} #rightarrow D^{0}#pi^{+}}{and charge conj.}");



    TCanvas *cEff = new TCanvas("cEff","cEff",1200,900);
    //cEff->SetBatch();
    cEff->SetLogy();
    //cMass->Divide(3,1);
    //cMass->cd(1);
    hEffPromptFill->Draw("hist");
    hEffPrompt->Draw("ep same");



    pvALICE->Draw("same");
    pvEn->Draw("same");
    //pvCent->Draw("same");
    pvJet->Draw("same");
    pvD->Draw("same");

    pvEta->Draw("same");
  //  leg->Draw("same");


    //cEff->SaveAs("DjetEff_Sim.png");
    //cEff->SaveAs("DjetEff_Sim.pdf");
    cEff->Print("deltaPt_prel.pdf");
    cEff->Print("deltaPt_prel.eps");
    cEff->Print("deltaPt_prel.png");




}

void setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Size_t size = 0.9, Width_t width=2, int scale = 0){

    if(scale)h->Scale(1,"width");
    h->SetMarkerStyle(Mstyle);
    h->SetMarkerColor(color);
    h->SetMarkerSize(size);
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetTitle(0);
    h->GetXaxis()->SetTitle("p_{T}^{D^{*+}}(GeV/c)");

    return;

}

void SaveCanvas(TCanvas *c, string name = "tmp"){

    c->SaveAs(Form("%s.png",name.c_str()));
    c->SaveAs(Form("%s.pdf",name.c_str()));

}
