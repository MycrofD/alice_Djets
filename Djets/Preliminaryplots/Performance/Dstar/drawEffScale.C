#include "style.C"
#include <string>
#include <sstream>
#include <iostream>

setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Width_t width);
    
 double zmin = 0, zmax = 2.;
    double jetmin = 0, jetmax = 50;
    double plotmin = 0, plotmax = 50;
    
    
    const int ptbinsJetN = 8;
    float ptJetbins[ptbinsJetN+1] = { 3,5,6,8,10,14,20,30,50 };
    
    int massColor = kBlack;
    int signalColor = kRed+1;
    int SBColor = kGreen+3;
    
void drawEffScale( int Rpar = 4 )
{
    
    style();
    gStyle->SetOptStat(000);
        
    gStyle->SetLegendFont(42);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.02);
    gStyle->SetPadTopMargin(0.02);
    //gStyle->SetLegendTextSize(0.05);
  
  
stringstream sst;
sst.clear(); sst.str("");

   //TFile *inFile = new TFile("JetPtSpectra_EffScale_FASTwoSDD_eff_ptD3.root","read");
   TFile *inFile = new TFile("JetPtSpectra_EffScale_FASTwoSDD_eff_ptD3_rebin.root","read");
	
    TH1F *hmean = (TH1F*)inFile->Get("hmean");
    TH1F *hsigma = (TH1F*)inFile->Get("hsigma");
    TH1F *hsign = (TH1F*)inFile->Get("hsign");
    TH1F *hsb = (TH1F*)inFile->Get("hsb");
    TH1F *hchi2 = (TH1F*)inFile->Get("hchi2");
    
    TH1F* hmass[ptbinsJetN];
    TH1F* hmass_l[ptbinsJetN];
    TH1F* hmass_u[ptbinsJetN];
    TH1F* hmass_c[ptbinsJetN];
    TF1* fullfit[ptbinsJetN];
    TH1F* hjetpt[ptbinsJetN];
    TH1F *hjetpt_s[ptbinsJetN];
    TH1F *hjetptsub[ptbinsJetN];
    TH1F *hjetptcorr[ptbinsJetN];
   

    for(int i=0; i<ptbinsJetN; i++){
            hmass[i] = (TH1F*)inFile->Get(Form("hmass_%d",i));
            hmass[i]->SetTitle();
            hmass[i]->SetMarkerColor(massColor);
            hmass[i]->SetLineColor(massColor);
            hmass[i]->SetMarkerStyle(20);
            hmass[i]->SetMarkerSize(1.2);
            hmass[i]->GetXaxis()->SetTitle("#it{M}(K#pi#pi)-#it{M}(K#pi) (GeV/#it{c^{2}})");
            //hmass[i]->GetYaxis()->SetTitle(Form("Entries/%.1f MeV/#it{c^{2}}",hmass[i]->GetBinWidth(1)*1000));
            hmass[i]->GetYaxis()->SetTitle("arb. units");
            hmass[i]->GetXaxis()->SetLabelSize(0.04);
            hmass[i]->GetXaxis()->SetTitleSize(0.045);
            hmass[i]->GetYaxis()->SetTitleOffset(1.55);
            hmass[i]->GetYaxis()->SetLabelSize(0.04);
            hmass[i]->GetYaxis()->SetTitleSize(0.05);
            //hmass[i]->SetMaximum(hmass[i]->GetMaximum()*1.35);
             hmass[i]->SetMinimum(1);
           
            fullfit[i] = (TF1*)inFile->Get(Form("fullfit_%d",i));
            fullfit[i]->SetNpx(150);
            fullfit[i]->SetLineWidth(2);
           
    }
    hmass[1]->GetYaxis()->SetTitleOffset(1.65);
    //hmass[1]->GetYaxis()->SetLabelSize(0.04);
    hmass[1]->SetMinimum(1);
     hmass[1]->SetMaximum(hmass[1]->GetMaximum()*1.5);
     hmass[4]->SetMaximum(hmass[4]->GetMaximum()*1.35);
     hmass[6]->SetMaximum(hmass[6]->GetMaximum()*1.45);
     
    TPaveText *pvALICE = new TPaveText(0.2,0.88,0.6,0.92,"brNDC");
    pvALICE->SetFillStyle(0);
    pvALICE->SetBorderSize(0);
    pvALICE->SetTextFont(42);
    pvALICE->SetTextSize(0.045);
    pvALICE->SetTextAlign(11);
    pvALICE->AddText("ALICE Preliminary");
    
    TPaveText *pvEn= new TPaveText(0.2,0.88,0.8,0.92,"brNDC");
    pvEn->SetFillStyle(0);
    pvEn->SetBorderSize(0);
    pvEn->SetTextFont(42);
    pvEn->SetTextSize(0.045);
    pvEn->SetTextAlign(11);
    pvEn->AddText("p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
    
    TPaveText *pvEff= new TPaveText(0.35,0.17,0.8,0.22,"brNDC");
    pvEff->SetFillStyle(0);
    pvEff->SetBorderSize(0);
    pvEff->SetTextFont(42);
    pvEff->SetTextSize(0.045);
    pvEff->SetTextAlign(11);
    pvEff->AddText("Weighted by the D^{*+} efficiency");
    
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
    
    TPaveText *pvEta = new TPaveText(0.2,0.72,0.5,0.77,"brNDC");
    pvEta->SetFillStyle(0);
    pvEta->SetBorderSize(0);
    pvEta->SetTextFont(42);
    pvEta->SetTextSize(0.045);
    pvEta->SetTextAlign(11);
    pvEta->AddText("|#it{#eta}_{jet}| < 0.5");
    
    TPaveText *pvDpt = new TPaveText(0.2,0.80,0.55,0.86,"brNDC");
    pvDpt->SetFillStyle(0);
    pvDpt->SetBorderSize(0);
    pvDpt->SetTextFont(42);
    pvDpt->SetTextSize(0.045);
    pvDpt->SetTextAlign(11);
    pvDpt->AddText("3 < #it{p}_{T, D^{*+}} < 36 GeV/#it{c}");
    
    TPaveText *pvpt1 = new TPaveText(0.5,0.7,0.8,0.75,"brNDC");
    pvpt1->SetFillStyle(0);
    pvpt1->SetBorderSize(0);
    pvpt1->SetTextFont(42);
    pvpt1->SetTextSize(0.045);
    pvpt1->SetTextAlign(11);
    pvpt1->AddText(Form("%.0f < #it{p}_{T, ch jet} < %.0f GeV/#it{c}",ptJetbins[1],ptJetbins[2]));
    
    TPaveText *pvpt2 = new TPaveText(0.48,0.72,0.8,0.77,"brNDC");
    pvpt2->SetFillStyle(0);
    pvpt2->SetBorderSize(0);
    pvpt2->SetTextFont(42);
    pvpt2->SetTextSize(0.045);
    pvpt2->SetTextAlign(11);
    pvpt2->AddText(Form("%.0f < #it{p}_{T, ch jet} < %.0f GeV/#it{c}",ptJetbins[4],ptJetbins[5]));
    
    TPaveText *pvpt3 = new TPaveText(0.2,0.80,0.5,0.85,"brNDC");
    pvpt3->SetFillStyle(0);
    pvpt3->SetBorderSize(0);
    pvpt3->SetTextFont(42);
    pvpt3->SetTextSize(0.045);
    pvpt3->SetTextAlign(11);
    pvpt3->AddText(Form("%.0f < #it{p}_{T, ch jet} < %.0f GeV/#it{c}",ptJetbins[6],ptJetbins[7]));
   
   
    TPaveText *pvsb1 = new TPaveText(0.57,0.6,0.9,0.65,"brNDC");
    pvsb1->SetFillStyle(0);
    pvsb1->SetBorderSize(0);
    pvsb1->SetTextFont(42);
    pvsb1->SetTextSize(0.04);
    pvsb1->AddText(Form("S/B (3#sigma) = %.2f",hsb->GetBinContent(hsb->FindBin((ptJetbins[1]+ptJetbins[2])/2. ))));
   
    TPaveText *pvmean1 = new TPaveText(0.57,0.55,0.9,0.6,"brNDC");
    pvmean1->SetFillStyle(0);
    pvmean1->SetBorderSize(0);
    pvmean1->SetTextFont(42);
    pvmean1->SetTextSize(0.04);
    pvmean1->AddText(Form("#mu = %.2f #pm %.2f MeV/#it{c^{2}}",hmean->GetBinContent(hmean->FindBin((ptJetbins[1]+ptJetbins[2])/2. )),hmean->GetBinError(hmean->FindBin((ptJetbins[1]+ptJetbins[2])/2. ))));
    
    TPaveText *pvsigma1 = new TPaveText(0.57,0.5,0.9,0.55,"brNDC");
    pvsigma1->SetFillStyle(0);
    pvsigma1->SetBorderSize(0);
    pvsigma1->SetTextFont(42);
    pvsigma1->SetTextSize(0.04);
    pvsigma1->AddText(Form("#sigma = %.2f #pm %.2f MeV/#it{c^{2}}",hsigma->GetBinContent(hsigma->FindBin((ptJetbins[1]+ptJetbins[2])/2. )),hsigma->GetBinError(hsigma->FindBin((ptJetbins[1]+ptJetbins[2])/2. ))));
    
    TPaveText *pvchi1 = new TPaveText(0.57,0.17,0.9,0.22,"brNDC");
    pvchi1->SetFillStyle(0);
    pvchi1->SetBorderSize(0);
    pvchi1->SetTextFont(42);
    pvchi1->SetTextSize(0.04);
    //pvchi1->AddText(Form("#chi^{2}/NDF = %.1f",hchi2->GetBinContent(hchi2->FindBin((ptJetbins[1]+ptJetbins[2])/2. ))));
  
    TPaveText *pvsb2 = new TPaveText(0.57,0.29,0.9,0.34,"brNDC");
    pvsb2->SetFillStyle(0);
    pvsb2->SetBorderSize(0);
    pvsb2->SetTextFont(42);
    pvsb2->SetTextSize(0.04);
    pvsb2->AddText(Form("S/B (3#sigma) = %.2f",hsb->GetBinContent(hsb->FindBin((ptJetbins[4]+ptJetbins[5])/2. ))));
    
    
    TPaveText *pvmean2 = new TPaveText(0.57,0.24,0.9,0.29,"brNDC");
    pvmean2->SetFillStyle(0);
    pvmean2->SetBorderSize(0);
    pvmean2->SetTextFont(42);
    pvmean2->SetTextSize(0.04);
    pvmean2->AddText(Form("#mu = %.2f #pm %.2f MeV/#it{c^{2}}",hmean->GetBinContent(hmean->FindBin((ptJetbins[4]+ptJetbins[5])/2. )),hmean->GetBinError(hmean->FindBin((ptJetbins[4]+ptJetbins[5])/2. ))));
    
    TPaveText *pvsigma2 = new TPaveText(0.57,0.19,0.9,0.24,"brNDC");
    pvsigma2->SetFillStyle(0);
    pvsigma2->SetBorderSize(0);
    pvsigma2->SetTextFont(42);
    pvsigma2->SetTextSize(0.04);
    pvsigma2->AddText(Form("#sigma = %.2f #pm %.2f MeV/#it{c^{2}}",hsigma->GetBinContent(hsigma->FindBin((ptJetbins[4]+ptJetbins[5])/2. )),hsigma->GetBinError(hsigma->FindBin((ptJetbins[4]+ptJetbins[5])/2. ))));
    
    TPaveText *pvchi2 = new TPaveText(0.57,0.17,0.9,0.22,"brNDC");
    pvchi2->SetFillStyle(0);
    pvchi2->SetBorderSize(0);
    pvchi2->SetTextFont(42);
    pvchi2->SetTextSize(0.04);
  //  pvchi2->AddText(Form("#chi^{2}/NDF = %.1f",hchi2->GetBinContent(hchi2->FindBin((ptJetbins[4]+ptJetbins[5])/2. ))));
    
    TPaveText *pvsb3 = new TPaveText(0.57,0.29,0.9,0.34,"brNDC");
    pvsb3->SetFillStyle(0);
    pvsb3->SetBorderSize(0);
    pvsb3->SetTextFont(42);
    pvsb3->SetTextSize(0.04);
    pvsb3->AddText(Form("S/B (3#sigma) = %.2f",hsb->GetBinContent(hsb->FindBin((ptJetbins[6]+ptJetbins[7])/2. ))));
    
    TPaveText *pvmean3 = new TPaveText(0.57,0.24,0.9,0.29,"brNDC");
    pvmean3->SetFillStyle(0);
    pvmean3->SetBorderSize(0);
    pvmean3->SetTextFont(42);
    pvmean3->SetTextSize(0.04);
    pvmean3->AddText(Form("#mu = %.2f #pm %.2f MeV/#it{c^{2}}",hmean->GetBinContent(hmean->FindBin((ptJetbins[6]+ptJetbins[7])/2. )),hmean->GetBinError(hmean->FindBin((ptJetbins[6]+ptJetbins[7])/2. ))));
    
    TPaveText *pvsigma3 = new TPaveText(0.57,0.19,0.9,0.24,"brNDC");
    pvsigma3->SetFillStyle(0);
    pvsigma3->SetBorderSize(0);
    pvsigma3->SetTextFont(42);
    pvsigma3->SetTextSize(0.04);
    pvsigma3->AddText(Form("#sigma = %.2f #pm %.2f MeV/#it{c^{2}}",hsigma->GetBinContent(hsigma->FindBin((ptJetbins[6]+ptJetbins[7])/2. )),hsigma->GetBinError(hsigma->FindBin((ptJetbins[6]+ptJetbins[7])/2. ))));

    TPaveText *pvchi3 = new TPaveText(0.57,0.17,0.9,0.22,"brNDC");
    pvchi3->SetFillStyle(0);
    pvchi3->SetBorderSize(0);
    pvchi3->SetTextFont(42);
    pvchi3->SetTextSize(0.04);
   // pvchi3->AddText(Form("#chi^{2}/NDF = %.1f",hchi2->GetBinContent(hchi2->FindBin((ptJetbins[6]+ptJetbins[7])/2. ))));

    TCanvas *cMass = new TCanvas("cMass","cMass",3000,1000);
    cMass->Divide(3,1);
    cMass->cd(1);
    hmass[1]->Draw();
    fullfit[1]->Draw("same");
    pvpt1->Draw("same");
    pvALICE->Draw("same");
    pvD->Draw("same");
    pvmean1->Draw("same");
    pvsigma1->Draw("same");
    pvsb1->Draw("same");
    pvEta->Draw("same");
    //pvchi1->Draw("same");
    pvJet->Draw("same");
    pvEff->Draw("same");
     
    cMass->cd(2);
    hmass[4]->Draw();
    fullfit[4]->Draw("same");
    
    pvpt2->Draw("same");
    pvEn->Draw("same");
    pvmean2->Draw("same");
    pvsigma2->Draw("same");
    pvsb2->Draw("same");
    pvDpt->Draw("same");
    //pvchi2->Draw("same");
   
    cMass->cd(3);
    hmass[6]->Draw();
    fullfit[6]->Draw("same");
    pvpt3->Draw("same");
    //pvEta->Draw("same");
    pvmean3->Draw("same");
    pvsigma3->Draw("same");
    pvsb3->Draw("same");
    //pvDpt->Draw("same");
    //pvchi3->Draw("same");
    
    
    //cMass->SaveAs("DjetInMass_JetPt_Perf.png");
    //cMass->SaveAs("DjetInMass_JetPt_Perf.pdf");
    cMass->Print("DjetInMass_JetPt_Perf.pdf");
    cMass->Print("DjetInMass_JetPt_Perf.eps");
    cMass->Print("DjetInMass_JetPt_Perf.root");
    cMass->Print("DjetInMass_JetPt_Perf.C");
   
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
