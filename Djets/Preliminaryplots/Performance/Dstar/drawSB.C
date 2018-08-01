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
    
    double efficiency[ptbinsDN];// = { 0.0353325, 0.0678059, 0.109101, 0.168871, 0.243708, 0.307365, 0.324496, 0.361858 };
    
    
    int massColor = kBlack;
    int signalColor = kRed+1;
    int SBColor = kGreen+3;
    
void drawSB( int Rpar = 4 )
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
     
     
    TH1F *hmean = (TH1F*)inFile->Get("hmean");
    TH1F *hsigma = (TH1F*)inFile->Get("hsigma");
    TH1F *hsign = (TH1F*)inFile->Get("hsign");
    TH1F *hsb = (TH1F*)inFile->Get("hsb");
    
    TH1F* hmass[ptbinsDN];
    TH1F* hmass_l[ptbinsDN];
    TH1F* hmass_u[ptbinsDN];
    TH1F* hmass_c[ptbinsDN];
    TF1* fullfit[ptbinsDN];
    TH1F* hjetpt[ptbinsDN];
    TH1F *hjetpt_s[ptbinsDN];
    TH1F *hjetptsub[ptbinsDN];
    TH1F *hjetptcorr[ptbinsDN];
   

    for(int i=0; i<ptbinsDN; i++){
            hmass[i] = (TH1F*)inFile->Get(Form("hmass_%d",i));
            hmass[i]->SetTitle();
            hmass[i]->SetMarkerColor(massColor);
            hmass[i]->SetLineColor(massColor);
            hmass[i]->SetMarkerStyle(20);
            hmass[i]->SetMarkerSize(1.2);
            hmass[i]->GetXaxis()->SetTitle("#it{M}(K#pi#pi)-#it{M}(K#pi) (GeV/#it{c^{2}})");
            hmass[i]->GetYaxis()->SetTitle(Form("Entries/%.1f MeV/#it{c^{2}}",hmass[i]->GetBinWidth(1)*1000));
            hmass[i]->GetXaxis()->SetLabelSize(0.04);
            hmass[i]->GetXaxis()->SetTitleSize(0.05);
            hmass[i]->GetYaxis()->SetLabelSize(0.045);
            hmass[i]->GetYaxis()->SetTitleSize(0.05);
           
            hmass[i]->SetMaximum(hmass[i]->GetMaximum()*1.3);
            hmass[i]->SetMinimum(1);
            
            hmass_l[i] = (TH1F*)inFile->Get(Form("hmass_l_%d",i));
            hmass_l[i]->SetTitle();
            hmass_l[i]->SetMarkerColor(massColor);
            hmass_l[i]->SetLineColor(SBColor);
            hmass_l[i]->SetFillColor(SBColor);
            hmass_l[i]->SetFillStyle(3004);
            
            hmass_u[i] = (TH1F*)inFile->Get(Form("hmass_u_%d",i));
            hmass_u[i]->SetTitle();
            hmass_u[i]->SetMarkerColor(massColor);
            hmass_u[i]->SetLineColor(SBColor);
            hmass_u[i]->SetFillColor(SBColor);
            hmass_u[i]->SetFillStyle(3004);
            
            hmass_c[i] = (TH1F*)inFile->Get(Form("hmass_c_%d",i));
            hmass_c[i]->SetTitle();
            hmass_c[i]->SetMarkerColor(massColor);
            hmass_c[i]->SetLineColor(signalColor);
            hmass_c[i]->SetFillColor(signalColor);
            hmass_c[i]->SetFillStyle(3005);
        
            fullfit[i] = (TF1*)inFile->Get(Form("fullfit_%d",i));
            fullfit[i]->SetNpx(150);
            fullfit[i]->SetLineWidth(2);
           
    }

     hmass[0]->SetMaximum(hmass[0]->GetMaximum()*1.1);
     
    TLegend *legBands = new TLegend(0.52,0.55,0.8,0.65);
    legBands->SetTextSize(0.04);
    legBands->AddEntry(hmass_c[0],"Signal region","f");
    legBands->AddEntry(hmass_l[0],"Side Bands","f");
    
     
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
    
    TPaveText *pvpt1 = new TPaveText(0.62,0.68,0.9,0.72,"brNDC");
    pvpt1->SetFillStyle(0);
    pvpt1->SetBorderSize(0);
    pvpt1->SetTextFont(42);
    pvpt1->SetTextSize(0.045);
    pvpt1->AddText(Form("%.0f < #it{p}_{T, D^{*+}} < %.0f GeV/#it{c}",ptDbins[bin1],ptDbins[bin1+1]));
    
    TPaveText *pvpt2 = new TPaveText(0.62,0.8,0.85,0.85,"brNDC");
    pvpt2->SetFillStyle(0);
    pvpt2->SetBorderSize(0);
    pvpt2->SetTextFont(42);
    pvpt2->SetTextSize(0.045);
    pvpt2->AddText(Form("%.0f < #it{p}_{T, D^{*+}} < %.0f GeV/#it{c}",ptDbins[bin2],ptDbins[bin2+1]));
    
    TPaveText *pvpt3 = new TPaveText(0.62,0.88,0.85,0.93,"brNDC");
    pvpt3->SetFillStyle(0);
    pvpt3->SetBorderSize(0);
    pvpt3->SetTextFont(42);
    pvpt3->SetTextSize(0.045);
    pvpt3->AddText(Form("%.0f < #it{p}_{T, D^{*+}} < %.0f GeV/#it{c}",ptDbins[bin3],ptDbins[bin3+1]));
   
    TPaveText *pvsig1 = new TPaveText(0.57,0.60,0.9,0.65,"brNDC");
    pvsig1->SetFillStyle(0);
    pvsig1->SetBorderSize(0);
    pvsig1->SetTextFont(42);
    pvsig1->SetTextSize(0.04);
    pvsig1->AddText(Form("S (3#sigma) = %.1f #pm %.1f",hsign->GetBinContent(hsign->FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2.)),hsign->GetBinError(hsign->FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2.))));
    
    TPaveText *pvsb1 = new TPaveText(0.57,0.55,0.9,0.60,"brNDC");
    pvsb1->SetFillStyle(0);
    pvsb1->SetBorderSize(0);
    pvsb1->SetTextFont(42);
    pvsb1->SetTextSize(0.04);
    pvsb1->AddText(Form("S/B (3#sigma) = %.2f",hsb->GetBinContent(hsb->FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2. ))));
   
    TPaveText *pvmean1 = new TPaveText(0.57,0.50,0.9,0.55,"brNDC");
    pvmean1->SetFillStyle(0);
    pvmean1->SetBorderSize(0);
    pvmean1->SetTextFont(42);
    pvmean1->SetTextSize(0.04);
    pvmean1->AddText(Form("#mu = %.2f #pm %.2f MeV/#it{c^{2}}",hmean->GetBinContent(hmean->FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2. )),hmean->GetBinError(hmean->FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2. ))));
    
    TPaveText *pvsigma1 = new TPaveText(0.57,0.45,0.9,0.50,"brNDC");
    pvsigma1->SetFillStyle(0);
    pvsigma1->SetBorderSize(0);
    pvsigma1->SetTextFont(42);
    pvsigma1->SetTextSize(0.04);
    pvsigma1->AddText(Form("#sigma = %.2f #pm %.2f MeV/#it{c^{2}}",hsigma->GetBinContent(hsigma->FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2. )),hsigma->GetBinError(hsigma->FindBin((ptDbins[bin1]+ptDbins[bin1+1])/2. ))));
    
    TPaveText *pvsig2 = new TPaveText(0.57,0.73,0.9,0.78,"brNDC");
    pvsig2->SetFillStyle(0);
    pvsig2->SetBorderSize(0);
    pvsig2->SetTextFont(42);
    pvsig2->SetTextSize(0.04);
    pvsig2->AddText(Form("S (3#sigma) = %.1f #pm %.1f",hsign->GetBinContent(hsign->FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2.)),hsign->GetBinError(hsign->FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2.))));
    
    TPaveText *pvsb2 = new TPaveText(0.57,0.68,0.9,0.73,"brNDC");
    pvsb2->SetFillStyle(0);
    pvsb2->SetBorderSize(0);
    pvsb2->SetTextFont(42);
    pvsb2->SetTextSize(0.04);
    pvsb2->AddText(Form("S/B (3#sigma) = %.2f",hsb->GetBinContent(hsb->FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2. ))));
    
    
    TPaveText *pvmean2 = new TPaveText(0.57,0.63,0.9,0.68,"brNDC");
    pvmean2->SetFillStyle(0);
    pvmean2->SetBorderSize(0);
    pvmean2->SetTextFont(42);
    pvmean2->SetTextSize(0.04);
    pvmean2->AddText(Form("#mu = %.2f #pm %.2f MeV/#it{c^{2}}",hmean->GetBinContent(hmean->FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2. )),hmean->GetBinError(hmean->FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2. ))));
    
    TPaveText *pvsigma2 = new TPaveText(0.57,0.58,0.9,0.63,"brNDC");
    pvsigma2->SetFillStyle(0);
    pvsigma2->SetBorderSize(0);
    pvsigma2->SetTextFont(42);
    pvsigma2->SetTextSize(0.04);
    pvsigma2->AddText(Form("#sigma = %.2f #pm %.2f MeV/#it{c^{2}}",hsigma->GetBinContent(hsigma->FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2. )),hsigma->GetBinError(hsigma->FindBin((ptDbins[bin2]+ptDbins[bin2+1])/2. ))));
    
    TPaveText *pvsig3 = new TPaveText(0.57,0.82,0.9,0.87,"brNDC");
    pvsig3->SetFillStyle(0);
    pvsig3->SetBorderSize(0);
    pvsig3->SetTextFont(42);
    pvsig3->SetTextSize(0.04);
    pvsig3->AddText(Form("S (3#sigma) = %.1f #pm %.1f",hsign->GetBinContent(hsign->FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2.)),hsign->GetBinError(hsign->FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2.))));
    
    TPaveText *pvsb3 = new TPaveText(0.57,0.77,0.9,0.82,"brNDC");
    pvsb3->SetFillStyle(0);
    pvsb3->SetBorderSize(0);
    pvsb3->SetTextFont(42);
    pvsb3->SetTextSize(0.04);
    pvsb3->AddText(Form("S/B (3#sigma) = %.2f",hsb->GetBinContent(hsb->FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2. ))));
    
    TPaveText *pvmean3 = new TPaveText(0.57,0.72,0.9,0.77,"brNDC");
    pvmean3->SetFillStyle(0);
    pvmean3->SetBorderSize(0);
    pvmean3->SetTextFont(42);
    pvmean3->SetTextSize(0.04);
    pvmean3->AddText(Form("#mu = %.2f #pm %.2f MeV/#it{c^{2}}",hmean->GetBinContent(hmean->FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2. )),hmean->GetBinError(hmean->FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2. ))));
    
    TPaveText *pvsigma3 = new TPaveText(0.57,0.67,0.9,0.72,"brNDC");
    pvsigma3->SetFillStyle(0);
    pvsigma3->SetBorderSize(0);
    pvsigma3->SetTextFont(42);
    pvsigma3->SetTextSize(0.04);
    pvsigma3->AddText(Form("#sigma = %.2f #pm %.2f MeV/#it{c^{2}}",hsigma->GetBinContent(hsigma->FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2. )),hsigma->GetBinError(hsigma->FindBin((ptDbins[bin3]+ptDbins[bin3+1])/2. ))));

    TCanvas *cMass = new TCanvas("cMass","cMass",3000,1000);
    //TCanvas *cMass = new TCanvas("cMass","cMass");
    //cMass->SetBatch();
    cMass->Divide(3,1);
    cMass->cd(1);
    hmass[bin1]->Draw();
    hmass_l[bin1]->Draw("hsame");
    hmass_u[bin1]->Draw("hsame");
    hmass_c[bin1]->Draw("hsame");
    hmass[bin1]->Draw("same");
    fullfit[bin1]->Draw("same");
    pvpt1->Draw("same");
    pvALICE->Draw("same");
    pvD->Draw("same");
    pvJet->Draw("same");
    pvmean1->Draw("same");
    pvsigma1->Draw("same");
    pvsb1->Draw("same");
    pvEta->Draw("same");
    //pvsig1->Draw("same");
    //legBands->Draw("same");
     
    cMass->cd(2);
    hmass[bin2]->Draw();
    hmass_l[bin2]->Draw("hsame");
    hmass_u[bin2]->Draw("hsame");
    hmass_c[bin2]->Draw("hsame");
    hmass[bin2]->Draw("same");
    fullfit[bin2]->Draw("same");
    
    pvpt2->Draw("same");
    pvEn->Draw("same");
    pvmean2->Draw("same");
    pvsigma2->Draw("same");
    pvsb2->Draw("same");
    //pvsig2->Draw("same");
    
    
    cMass->cd(3);
    hmass[bin3]->Draw();
    hmass_l[bin3]->Draw("hsame");
    hmass_u[bin3]->Draw("hsame");
    hmass_c[bin3]->Draw("hsame");
    hmass[bin3]->Draw("same");
    fullfit[bin3]->Draw("same");
    pvpt3->Draw("same");
    //pvEta->Draw("same");
    pvmean3->Draw("same");
    pvsigma3->Draw("same");
    pvsb3->Draw("same");
    //pvsig3->Draw("same");
    legBands->Draw("same");
    
    //cMass->SaveAs("DjetInMass_DPt_Perf.png");
    //cMass->SaveAs("DjetInMass_DPt_Perf.pdf");
    //cMass->SaveAs("DjetInMass_DPt_Perf.eps");
    cMass->Print("DjetInMass_DPt_Perf_v2.pdf");
    cMass->Print("DjetInMass_DPt_Perf_v2.eps");
    cMass->Print("DjetInMass_DPt_Perf_v2.root");
    cMass->Print("DjetInMass_DPt_Perf_v2.C");
   
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
   
}
