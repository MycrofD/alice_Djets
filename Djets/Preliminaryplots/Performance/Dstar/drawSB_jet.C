#include "style.C"
#include <string>
#include <sstream>
#include <iostream>
#include <TPDF.h>

setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Width_t width);
    
 double zmin = 0, zmax = 2.;
    double jetmin = 0, jetmax = 50;
    double plotmin = 0, plotmax = 50;
    
    const int ptbinsDN = 10;
    float ptDbins[ptbinsDN+1] = { 3,4,5,6,7,8,10,12,16,24,36 };
    
    const int ptbinsJetN = 6;
    double ptJetbins[ptbinsJetN+1] = { 5,6,8,10,14,20,30 };
    
    
    double efficiency[ptbinsDN];// = { 0.0353325, 0.0678059, 0.109101, 0.168871, 0.243708, 0.307365, 0.324496, 0.361858 };
    
    
    int massColor = kBlack;
    int signalColor = kRed+1;
    int SBColor = kGreen+3;
    int subColor = kBlue+1;
    double markersize = 2.2;
    
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
   TFile *inFile = new TFile("JetPtSpectra_SB_FASTwoSDD_eff_ptD3.root","read");
	
    int bin1 = 0, bin2 = 3, bin3 = 7;
    
    int islog = 1;
    
    
    TH1F* hjetpt[ptbinsDN];
    TH1F* hjetpt_s[ptbinsDN];
    TH1F* hjetptsub[ptbinsDN];
    TH1F* hmass_c[ptbinsDN];
   
    TH1F* hjetpt[ptbinsDN];
    TH1F *hjetpt_s[ptbinsDN];
    TH1F *hjetptsub[ptbinsDN];
    TH1F *hjetptcorr[ptbinsDN];
   

    for(int i=0; i<ptbinsDN; i++){
        
            TH1F *hjet = (TH1F*)inFile->Get(Form("hjetpt_%d",i));
            hjetpt[i] = (TH1F*)hjet->Rebin(ptbinsJetN,Form("hjetpt_%d",i),ptJetbins);
            
            hjetpt[i]->SetTitle();
            hjetpt[i]->SetMarkerColor(signalColor);
            hjetpt[i]->SetLineColor(signalColor);
            hjetpt[i]->SetMarkerStyle(24);
            hjetpt[i]->SetMarkerSize(markersize);
            hjetpt[i]->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
            hjetpt[i]->GetYaxis()->SetTitle("Entries");
            hjetpt[i]->GetXaxis()->SetLabelSize(0.04);
            hjetpt[i]->GetXaxis()->SetTitleSize(0.05);
            hjetpt[i]->GetYaxis()->SetLabelSize(0.045);
            hjetpt[i]->GetYaxis()->SetTitleSize(0.05);
            hjetpt[i]->GetYaxis()->SetTitleOffset(1.5);
           
            hjetpt[i]->SetMaximum(hjetpt[i]->GetMaximum()*1.55);
            hjetpt[i]->SetMaximum(1600);
            if(islog) hjetpt[i]->SetMaximum(hjetpt[i]->GetMaximum()*10);
            //hjetpt[i]->SetMaximum(1600);
            hjetpt[i]->SetMinimum(-50);
            if(islog) hjetpt[i]->SetMinimum(1);
            
            TH1F *hjets = (TH1F*)inFile->Get(Form("hjetpt_s_%d",i));
            hjetpt_s[i] = (TH1F*)hjets->Rebin(ptbinsJetN,Form("hjetpt_s_%d",i),ptJetbins);
            hjetpt_s[i]->SetTitle();
            hjetpt_s[i]->SetMarkerColor(SBColor);
            hjetpt_s[i]->SetLineColor(SBColor);
            hjetpt_s[i]->SetMarkerStyle(25);
            hjetpt_s[i]->SetMarkerSize(markersize);
    
            TH1F *hjetsub = (TH1F*)inFile->Get(Form("hjetptsub_%d",i));
            hjetptsub[i] = (TH1F*)hjetsub->Rebin(ptbinsJetN,Form("hjetptsub_%d",i),ptJetbins);
            hjetptsub[i]->SetTitle();
            hjetptsub[i]->SetMarkerColor(subColor);
            hjetptsub[i]->SetLineColor(subColor);
            hjetptsub[i]->SetMarkerStyle(27);
            hjetptsub[i]->SetMarkerSize(markersize+0.5);
           
    }

   // hjetpt[bin1]->SetMaximum(hjetpt[bin1]->GetMaximum()*1.1);
   // hjetpt[bin2]->SetMaximum(1600);
   // hjetpt[bin3]->SetMaximum(1600);
   // hjetpt[bin1]->SetMaximum(1600);
     
     if(islog) hjetpt[bin3]->SetMaximum(hjetpt[bin3]->GetMaximum()*5000);
     
     TH1F *hh = (TH1F*)hjetpt[0]->Clone("hh");
     hh->SetMarkerSize(0);
     
   // TLegend *legBands = new TLegend(0.5,0.55,0.8,0.72);
   /* TLegend *legBands = new TLegend(0.3,0.4,0.88,0.78);
    legBands->SetTextSize(0.045);
    //legBands->AddEntry(hjetpt[0],"Signal region, |#it{M}(K#pi#pi)-#it{M}(K#pi)|<3#sigma","p");
    legBands->AddEntry(hjetpt[0],"#splitline{Signal region}{ |#it{M}(K#pi#pi)-#it{M}(K#pi)|<3#sigma }","p");

    legBands->AddEntry(hjetpt_s[0],"#splitline{Side Bands (SB)}{normalised to Signal region }","p");
    legBands->AddEntry(hh,"#splitline{-8<(#it{M}(K#pi#pi)-#it{M}(K#pi))<-5#sigma}{5<(#it{M}(K#pi#pi)-#it{M}(K#pi))<13#sigma} ","p");
    
   // legBands->AddEntry(hjetpt_s[0],"#splitline{Side Bands (SB)} {|#it{M}(K#pi#pi)-#it{M}(K#pi)|<3#sigma }{normalised to the singal region}","p");
    legBands->AddEntry(hjetptsub[0],"Signal - SB","p");*/
    
     TLegend *legBands1 = new TLegend(0.13,0.8,0.7,0.88);
    legBands1->SetTextSize(0.048);
    legBands1->SetFillStyle(0);
    legBands1->SetTextAlign(13);
    legBands1->AddEntry(hjetpt[0],"#splitline{Signal region}{ |#it{M}(K#pi#pi)-#it{M}(K#pi)|<3#sigma }","p");

    // TLegend *legBands2 = new TLegend(0.15,0.63,0.9,0.73);
     TLegend *legBands2 = new TLegend(0.13,0.71,0.7,0.76);
    legBands2->SetTextSize(0.048);
    legBands2->SetFillStyle(0);
    legBands2->SetTextAlign(13);
     
     //legBands2->AddEntry(hjetpt_s[0],"#splitline{Side Bands (SB)}{normalised to Signal region }","p");
     legBands2->AddEntry(hjetpt_s[0],"Side Bands (SB)","p");
     
     TLegend *legBands22 = new TLegend(0.13,0.6,0.7,0.7);
    legBands22->SetTextSize(0.048);
    legBands22->SetFillStyle(0);
    legBands22->AddEntry(hh,"#splitline{-8<(#it{M}(K#pi#pi)-#it{M}(K#pi))<-5#sigma}{5<(#it{M}(K#pi#pi)-#it{M}(K#pi))<13#sigma} ","p");
    legBands22->SetTextAlign(13);
    
    TLegend *legBands23 = new TLegend(0.13,0.53,0.7,0.58);
    legBands23->SetTextSize(0.048);
    legBands23->SetFillStyle(0);
    legBands23->SetTextAlign(13);
     
    legBands23->AddEntry(hh,"normalised to Signal region","p");


     TLegend *legBands3 = new TLegend(0.13,0.46,0.7,0.51);
    legBands3->SetTextSize(0.048);
    legBands3->SetFillStyle(0);
    legBands3->SetTextAlign(13);
    
    legBands3->AddEntry(hjetptsub[0],"Signal - SB","p");
    
     
    TPaveText *pvALICE = new TPaveText(0.2,0.88,0.6,0.92,"brNDC");
    pvALICE->SetFillStyle(0);
    pvALICE->SetBorderSize(0);
    pvALICE->SetTextFont(42);
    pvALICE->SetTextSize(0.045);
    pvALICE->SetTextAlign(11);
    pvALICE->AddText("ALICE Preliminary");
    
    TPaveText *pvEn = new TPaveText(0.2,0.88,0.8,0.92,"brNDC");
    pvEn->SetFillStyle(0);
    pvEn->SetBorderSize(0);
    pvEn->SetTextFont(42);
    pvEn->SetTextSize(0.045);
    pvEn->SetTextAlign(11);
    //pvEn->AddText("p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, 603M events");
    pvEn->AddText("p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
    
    TPaveText *pvD = new TPaveText(0.2,0.82,0.55,0.87,"brNDC");
    pvD->SetFillStyle(0);
    pvD->SetBorderSize(0);
    pvD->SetTextFont(42);
    pvD->SetTextSize(0.045);
    pvD->SetTextAlign(11);
    pvD->AddText("D^{*+} #rightarrow D^{0}#pi^{+} and charge conj.");
    
    TPaveText *pvJet = new TPaveText(0.2,0.77,0.55,0.82,"brNDC");
    pvJet->SetFillStyle(0);
    pvJet->SetBorderSize(0);
    pvJet->SetTextFont(42);
    pvJet->SetTextSize(0.045);
    pvJet->SetTextAlign(11);
    pvJet->AddText("in Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4");
    
    TPaveText *pvEta = new TPaveText(0.2,0.72,0.4,0.77,"brNDC");
    pvEta->SetFillStyle(0);
    pvEta->SetBorderSize(0);
    pvEta->SetTextFont(42);
    pvEta->SetTextSize(0.045);
    pvEta->SetTextAlign(11);
    pvEta->AddText("|#it{#eta}_{jet}| < 0.5");
    
    TPaveText *pvpt1 = new TPaveText(0.58,0.63,0.85,0.68,"brNDC");
    pvpt1->SetFillStyle(0);
    pvpt1->SetBorderSize(0);
    pvpt1->SetTextFont(42);
    pvpt1->SetTextSize(0.045);
    pvpt1->AddText(Form("%.0f < #it{p}_{T, D^{*+}} < %.0f GeV/#it{c}",ptDbins[bin1],ptDbins[bin1+1]));
    
    TPaveText *pvpt2 = new TPaveText(0.58,0.63,0.85,0.68,"brNDC");
    //TPaveText *pvpt2 = new TPaveText(0.2,0.88,0.65,0.92,"brNDC");
    pvpt2->SetFillStyle(0);
    pvpt2->SetBorderSize(0);
    pvpt2->SetTextFont(42);
    pvpt2->SetTextSize(0.045);
    pvpt2->AddText(Form("%.0f < #it{p}_{T, D^{*+}} < %.0f GeV/#it{c}",ptDbins[bin2],ptDbins[bin2+1]));
 
 
    TPaveText *pvpt3 = new TPaveText(0.65,0.88,0.85,0.92,"brNDC");
    //TPaveText *pvpt3 = new TPaveText(0.6,0.8,0.83,0.85,"brNDC");
    pvpt3->SetFillStyle(0);
    pvpt3->SetBorderSize(0);
    pvpt3->SetTextFont(42);
    pvpt3->SetTextSize(0.045);
    pvpt3->AddText(Form("%.0f < #it{p}_{T, D^{*+}} < %.0f GeV/#it{c}",ptDbins[bin3],ptDbins[bin3+1]));
   
   

    TCanvas *cMass = new TCanvas("cMass","cMass",3000,1000);
    //TCanvas *cMass = new TCanvas("cMass","cMass");
    //cMass->SetBatch();
    cMass->Divide(3,1);
    cMass->cd(1);
    gPad->SetLogy(islog);
    hjetpt[bin1]->Draw();
    hjetpt_s[bin1]->Draw("same");
    hjetptsub[bin1]->Draw("same");

   
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
    hjetpt_s[bin2]->Draw("same");
    hjetptsub[bin2]->Draw("same");
   
    pvpt2->Draw("same");
    pvEn->Draw("same");
   
    //legBands->Draw("same");
   
    
    cMass->cd(3);
    gPad->SetLogy(islog);
    hjetpt[bin3]->Draw();
    hjetpt_s[bin3]->Draw("same");
    hjetptsub[bin3]->Draw("same");
   
   // pvEn->Draw("same");
    pvpt3->Draw("same");
    //pvEta->Draw("same");
    
     legBands1->Draw("same");
    legBands2->Draw("same");
    legBands22->Draw("same");
    legBands23->Draw("same");
    legBands3->Draw("same");

  
   /* cMass->Print("RawJetPt_Perf_v5.pdf");
    cMass->Print("RawJetPt_Perf_v5.eps");
    cMass->Print("RawJetPt_Perf_v5.root");
    cMass->Print("RawJetPt_Perf_v5.C");
    */
    
    if(islog) SaveCanvas(cMass,"RawJetPt_Perf_v5_log");
    else SaveCanvas(cMass,"RawJetPt_Perf_v5");
   
}

void setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Size_t size = 0.9, Width_t width=2, int scale = 0){
    
    if(scale)h->Scale(1,"width");
    h->SetMarkerStyle(Mstyle);
    h->SetMarkerColor(color);
    h->SetMarkerSize(size);
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetTitle(0);
    h->GetXaxis()->SetTitle("p_{T}^{D*}(GeV/c)");
    
    return;

}

void SaveCanvas(TCanvas *c, string name = "tmp"){
    
    c->SaveAs(Form("%s.png",name.c_str()));
    c->SaveAs(Form("%s.pdf",name.c_str()));
    c->SaveAs(Form("%s.eps",name.c_str()));
    c->SaveAs(Form("%s.root",name.c_str()));
    c->SaveAs(Form("%s.C",name.c_str()));
   
}
