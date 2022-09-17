#include "style.C"
#include <string>
#include <sstream>
#include <iostream>
#include <TPDF.h>

setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Width_t width);

    const int ptbinsDN = 10;
    float ptDbins[ptbinsDN+1] = { 3,4,5,6,7,8,10,12,16,24,36 };

    const int ptbinsJetN = 7;
    double ptJetbins[ptbinsJetN+1] = { 5,6,8,10,14,20,30,50 };

    int massColor = kBlack;
    int signalColor = kRed+1;
    int SBColor = kGreen+3;
    int subColor = kBlue+1;
    double markersize = 2.;
    //double markersize = 1.6;
    //int markerstyle[] = { 20,21,33 };
    int markerstyle[] = { 24,25,27 };

    double ltextsize = 0.053;

void drawSB_jet( int Rpar = 4 )
{

    style();
    gStyle->SetOptStat(000);


    gStyle->SetLegendFont(42);
    //gStyle->SetLegendTextSize(0.05);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.02);
    gStyle->SetPadTopMargin(0.02);


stringstream sst;
sst.clear(); sst.str("");

   //TFile *inFile = new TFile("JetPtSpectra_SB_FASTwoSDD_eff_ptD3_rebin.root","read");
  // TFile *inFile = new TFile("JetPtSpectra_SB_eff.root","read");
    TFile *inFile = new TFile(
//"/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50_ppbinning/signalExtraction/JetPtSpectra_SB_eff.root"
"/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_cutTight/DzeroR03_def_437_old0/Default/signalExtraction/JetPtSpectra_SB_eff.root","read");

    int bin1 = 1, bin2 = 5, bin3 = 7;

    int islog = 1;


    //TH1F* hjetpt[ptbinsDN]=new TH1F("yo","yo",100,4.8,50.2);
    TH1F* hjetpt[ptbinsDN];
    TH1F* hjetpt_s[ptbinsDN];
    TH1F* hjetptsub[ptbinsDN];
    TH1F* hjetptcorr[ptbinsDN];


    for(int i=0; i<ptbinsDN; i++){

            TH1F *hjet = (TH1F*)inFile->Get(Form("hjetpt_%d",i));
            hjetpt[i] = (TH1F*)hjet->Rebin(ptbinsJetN,Form("hjetpt_%d",i),ptJetbins);

            hjetpt[i]->SetTitle();
            hjetpt[i]->SetMarkerColor(signalColor);
            hjetpt[i]->SetLineColor(signalColor);
            hjetpt[i]->SetLineWidth(2);
            //hjetpt[i]->SetMarkerStyle(24);
            hjetpt[i]->SetMarkerStyle(markerstyle[0]);
            hjetpt[i]->SetMarkerSize(markersize);

            hjetpt[i]->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
            hjetpt[i]->GetYaxis()->SetTitle("Entries");
            hjetpt[i]->GetXaxis()->SetLabelSize(0.045);
            hjetpt[i]->GetXaxis()->SetTitleSize(0.055);
            hjetpt[i]->GetYaxis()->SetLabelSize(0.045);
            hjetpt[i]->GetYaxis()->SetTitleSize(0.055);

            hjetpt[i]->GetYaxis()->SetTitleOffset(1.4);

            hjetpt[i]->SetMaximum(hjetpt[i]->GetMaximum()*1.55);
            hjetpt[i]->SetMaximum(1600);
            if(islog) hjetpt[i]->SetMaximum(hjetpt[i]->GetMaximum()*10);
            //hjetpt[i]->SetMaximum(1600);
            hjetpt[i]->SetMinimum(-50);
            if(islog) hjetpt[i]->SetMinimum(1);

            TH1F *hjets = (TH1F*)inFile->Get(Form("hjetpt_sb_%d",i));
            hjetpt_s[i] = (TH1F*)hjets->Rebin(ptbinsJetN,Form("hjetpt_s_%d",i),ptJetbins);
            hjetpt_s[i]->SetTitle();
            hjetpt_s[i]->SetMarkerColor(SBColor);
            hjetpt_s[i]->SetLineColor(SBColor);
            hjetpt_s[i]->SetLineWidth(2);
            //hjetpt_s[i]->SetMarkerStyle(25);
            hjetpt_s[i]->SetMarkerStyle(markerstyle[1]);
            hjetpt_s[i]->SetMarkerSize(markersize);


            TH1F *hjetsub = (TH1F*)inFile->Get(Form("hjetptsub_%d",i));
            hjetptsub[i] = (TH1F*)hjetsub->Rebin(ptbinsJetN,Form("hjetptsub_%d",i),ptJetbins);
            hjetptsub[i]->SetTitle();
            hjetptsub[i]->SetMarkerColor(subColor);
            hjetptsub[i]->SetLineColor(subColor);
            hjetptsub[i]->SetLineWidth(2);
          // hjetptsub[i]->SetMarkerStyle(27);
            hjetptsub[i]->SetMarkerStyle(markerstyle[2]);
            hjetptsub[i]->SetMarkerSize(markersize+0.5);

    }

   // hjetpt[bin1]->SetMaximum(hjetpt[bin1]->GetMaximum()*1.1);
   // hjetpt[bin2]->SetMaximum(1600);
   // hjetpt[bin3]->SetMaximum(1600);
   // hjetpt[bin1]->SetMaximum(1600);

     if(islog) hjetpt[bin3]->SetMaximum(hjetpt[bin3]->GetMaximum()*1000);

//     TH1F *hh = (TH1F*)hjetpt[0]->Clone("hh");
//     hh->SetMarkerSize(0);

    TLegend *legBands1 = new TLegend(0.15,0.76,0.7,0.86);
    legBands1->SetTextSize(ltextsize);
    legBands1->SetFillStyle(0);
    legBands1->SetTextAlign(13);
    //legBands1->AddEntry(hjetpt[0],"#splitline{Signal region}{|#it{M}(K#pi)|<3#sigma}","p");
    legBands1->AddEntry(hjetpt[0],"#splitline{Signal region}{|#it{M}(K#pi)-#it{M}_{D^{0}}|<2#sigma}","p");

    TLegend *legBands2 = new TLegend(0.15,0.65,0.7,0.75);
    legBands2->SetTextSize(ltextsize);
    legBands2->SetFillStyle(0);
    legBands2->SetTextAlign(13);
      legBands2->AddEntry(hjetpt_s[0],"#splitline{Side bands (SB)}{4#sigma<|#it{M}(K#pi)-#it{M}_{D^{0}}|<9#sigma}","p");
  //  legBands2->AddEntry(hjetpt_s[0],"Side Bands (SB)","p");

  //  TLegend *legBands22 = new TLegend(0.28,0.68,0.7,0.72,"4<|#it{M}(K#pi)|<9#sigma");
  //  legBands22->SetTextSize(ltextsize);
  //  legBands22->SetFillStyle(0);
    //legBands22->AddEntry(hh,"|4<(#it{M}(K#pi))<9#sigma|","p");
  //  legBands22->SetTextAlign(13);

    TLegend *legBands23 = new TLegend(0.28,0.63,0.7,0.65,"normalised to signal region");
    legBands23->SetTextSize(ltextsize);
    legBands23->SetFillStyle(0);
    legBands23->SetTextAlign(11);
    //legBands23->AddEntry(hh,"normalised to Signal region","p");

     TLegend *legBands3 = new TLegend(0.15,0.55,0.7,0.6);
    legBands3->SetTextSize(ltextsize);
    legBands3->SetFillStyle(0);
    legBands3->SetTextAlign(13);
    legBands3->AddEntry(hjetptsub[0],"Signal - SB","p");

    //TPaveText *pvALICE = new TPaveText(0.25,0.88,0.6,0.92,"brNDC");
    //TPaveText *pvALICE = new TPaveText(0.247,0.88,0.6,0.92,"brNDC");
    //TPaveText *pvALICE = new TPaveText(0.187,0.88,0.6,0.92,"brNDC");
    TPaveText *pvALICE = new TPaveText(0.187,0.88,0.6,0.92,"brNDC");
    pvALICE->SetFillStyle(0);
    pvALICE->SetBorderSize(0);
    pvALICE->SetTextFont(42);
    pvALICE->SetTextSize(0.055);
    pvALICE->SetTextAlign(11);
    pvALICE->AddText("ALICE Preliminary");

    TPaveText *pvEn = new TPaveText(0.25,0.88,0.8,0.92,"brNDC");
    pvEn->SetFillStyle(0);
    pvEn->SetBorderSize(0);
    pvEn->SetTextFont(42);
    pvEn->SetTextSize(0.053);
    pvEn->SetTextAlign(11);
    //pvEn->AddText("p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, 603M events");
    pvEn->AddText("pp, #sqrt{#it{s}} = 5.02 TeV");

    //TPaveText *pvD = new TPaveText(0.25,0.82,0.55,0.87,"brNDC");
    //TPaveText *pvD = new TPaveText(0.247,0.82,0.55,0.87,"brNDC");
    TPaveText *pvD = new TPaveText(0.187,0.82,0.55,0.87,"brNDC");
    pvD->SetFillStyle(0);
    pvD->SetBorderSize(0);
    pvD->SetTextFont(42);
    pvD->SetTextSize(0.053);
    pvD->SetTextAlign(11);
    pvD->AddText("D^{0} #rightarrow K^{-}#pi^{+} and charge conj.");

    //TPaveText *pvJet = new TPaveText(0.25,0.77,0.55,0.82,"brNDC");
    TPaveText *pvJet = new TPaveText(0.189,0.77,0.55,0.82,"brNDC");
    pvJet->SetFillStyle(0);
    pvJet->SetBorderSize(0);
    pvJet->SetTextFont(42);
    pvJet->SetTextSize(0.053);
    pvJet->SetTextAlign(11);
    pvJet->AddText("in charged jets, anti-#it{k}_{T}, #it{R} = 0.3");

    //TPaveText *pvEta = new TPaveText(0.25,0.72,0.6,0.77,"brNDC");
    //TPaveText *pvEta = new TPaveText(0.247,0.72,0.6,0.77,"brNDC");
    TPaveText *pvEta = new TPaveText(0.187,0.72,0.6,0.77,"brNDC");
    pvEta->SetFillStyle(0);
    pvEta->SetBorderSize(0);
    pvEta->SetTextFont(42);
    pvEta->SetTextSize(0.053);
    pvEta->SetTextAlign(11);
    pvEta->AddText("|#it{#eta}_{lab}^{jet}| < 0.6");
    //pvEta->AddText("|#it{#eta}_{jet}| < 0.6");

    TPaveText *pvpt1 = new TPaveText(0.58,0.66,0.85,0.71,"brNDC");
    pvpt1->SetFillStyle(0);
    pvpt1->SetBorderSize(0);
    pvpt1->SetTextFont(42);
    pvpt1->SetTextSize(0.053);
    pvpt1->AddText(Form("%.0f < #it{p}_{T,D^{0}} < %.0f GeV/#it{c}",ptDbins[bin1],ptDbins[bin1+1]));

    TPaveText *pvpt2 = new TPaveText(0.58,0.66,0.85,0.71,"brNDC");
    //TPaveText *pvpt2 = new TPaveText(0.2,0.88,0.65,0.92,"brNDC");
    pvpt2->SetFillStyle(0);
    pvpt2->SetBorderSize(0);
    pvpt2->SetTextFont(42);
    pvpt2->SetTextSize(0.053);
    pvpt2->AddText(Form("%.0f < #it{p}_{T,D^{0}} < %.0f GeV/#it{c}",ptDbins[bin2],ptDbins[bin2+1]));


    TPaveText *pvpt3 = new TPaveText(0.35,0.88,0.65,0.92,"brNDC");
    TPaveText *pvpt3 = new TPaveText(0.31,0.88,0.65,0.92,"brNDC");
    //TPaveText *pvpt3 = new TPaveText(0.6,0.8,0.83,0.85,"brNDC");
    pvpt3->SetFillStyle(0);
    pvpt3->SetBorderSize(0);
    pvpt3->SetTextFont(42);
    pvpt3->SetTextSize(0.053);
    pvpt3->AddText(Form("%.0f < #it{p}_{T,D^{0}} < %.0f GeV/#it{c}",ptDbins[bin3],ptDbins[bin3+1]));



    //TCanvas *cMass = new TCanvas("cMass","cMass",3000,900);
    //TCanvas *cMass = new TCanvas("cMass","cMass",3000,1800);
    //TCanvas *cMass = new TCanvas("cMass","cMass",3000,1400);
    TCanvas *cMass = new TCanvas("cMass","cMass",1800,840);
//    TCanvas *cMass = new TCanvas("cMass","cMass",2100,980);
    //TCanvas *cMass = new TCanvas("cMass","cMass",1800,600);
    //TCanvas *cMass = new TCanvas("cMass","cMass",900,300);
    //TCanvas *cMass = new TCanvas("cMass","cMass");
    //cMass->SetBatch();
    TCanvas *cMass = new TCanvas("cMass","cMass",2160,1008);



    cMass->Divide(3,1);
    cMass->cd(1);
    gPad->SetLogy(islog);
// the DrawFrame method takes the parameters (x1,y1,x2,y2).
//   TH1F *hr = cMass->cd(1)->DrawFrame(4.8,10,50.2,12);
 //  gr = new TGraph(n,x,y);
    hjetpt[bin1]->Draw();
    hjetpt[bin1]->Draw("same");
    hjetpt[bin1]->Draw("same");
    hjetpt_s[bin1]->Draw("same");
    hjetpt_s[bin1]->Draw("same");
    hjetpt_s[bin1]->Draw("same");
    hjetptsub[bin1]->Draw("same");
    hjetptsub[bin1]->Draw("same");
    hjetptsub[bin1]->Draw("same");


  // change this line and leave out the "A" for axis.
//   gr->Draw("CP");

    pvpt1->Draw("same");
    pvALICE->Draw("same");
    pvD->Draw("same");
    pvJet->Draw("same");
    pvEta->Draw("same");
    //pvsig1->Draw("same");
    //legBands->Draw("same");

    cMass->cd(2);
    gPad->SetLogy(islog);
    hjetpt[bin2]->Draw();
    hjetpt[bin2]->Draw("same");
    hjetpt_s[bin2]->Draw("same");
    hjetpt_s[bin2]->Draw("same");
    hjetptsub[bin2]->Draw("same");
    hjetptsub[bin2]->Draw("same");

    pvpt2->Draw("same");
    pvEn->Draw("same");

    //legBands->Draw("same");


    cMass->cd(3);
    gPad->SetLogy(islog);
    hjetpt[bin3]->Draw();
    hjetpt[bin3]->Draw("same");
    hjetpt_s[bin3]->Draw("same");
    hjetpt_s[bin3]->Draw("same");
    hjetptsub[bin3]->Draw("same");
    hjetptsub[bin3]->Draw("same");

   // pvEn->Draw("same");
    pvpt3->Draw("same");
    //pvEta->Draw("same");

     legBands1->Draw("same");
    legBands2->Draw("same");
    //legBands22->Draw("same");
    legBands23->Draw("same");
    legBands3->Draw("same");


   /* cMass->Print("RawJetPt_Perf_v5.pdf");
    cMass->Print("RawJetPt_Perf_v5.eps");
    cMass->Print("RawJetPt_Perf_v5.root");
    cMass->Print("RawJetPt_Perf_v5.C");
    */

    if(islog) SaveCanvas(cMass,"RawJetPt_Perf_log_2");
    else SaveCanvas(cMass,"RawJetPt_Perf_2");

}

void setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Size_t size = 0.9, Width_t width=2, int scale = 0){

    if(scale)h->Scale(1,"width");
    h->SetMarkerStyle(Mstyle);
    h->SetMarkerColor(color);
    h->SetMarkerSize(size);
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetTitle(0);
    h->GetXaxis()->SetTitle("p_{T}^{#D^0}(GeV/c)");

    return;

}

void SaveCanvas(TCanvas *c, string name = "tmp"){

    c->SaveAs(Form("%s.png",name.c_str()));
    c->SaveAs(Form("%s.pdf",name.c_str()));
    c->SaveAs(Form("%s.eps",name.c_str()));
    c->SaveAs(Form("%s.root",name.c_str()));
    c->SaveAs(Form("%s.C",name.c_str()));

}
