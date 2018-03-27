//
// Macro to extract jet pt spectra from simulation, prompt or non-prompt D
//
// Author: B.Trzeciak (barbara.antonina.trzeciak@cern.ch)
//

#include "config.h"

double plotmin = 0, plotmax = 50;

double BR = 0.0257;
const int APb = 208;
const double simScaling =  APb;
//double pPbscaling = 2.1/208.;

TH1* GetInputHist(TString inFile, string histName,TH1 *hh);

void ScaleHist(TH1 *hh, int full = 0);
void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, int Msize = 1.3, Width_t Lwidth = 2, Style_t Lstyle = 1);
void SaveCanvas(TCanvas *c, TString name = "tmp");
void compareSpectra(TH1D **hFD, const int nFiles, bool isjet, bool quark, TString cName = "canvas");
TH1* GetUpSys(TH1D **hFD, const int nFiles, TH1D *hFD_up);
TH1* GetDownSys(TH1D **hFD, const int nFiles, TH1D *hFD_down);


// quark: 0 --charm, 1 -- beauty
// jet: 1- jet spectrum, 1 - D pT spectrum
// isjetptcut: if D pT spectrum if we apply jet pT range
// isDptcut: if jet pT spectrum if we apply D pt range
// isEff: is jet pT spectrum if the spectrum is scaled by the eff_nonprompt/eff_prompt in D pt bins

void plotSimSpectra(int quark = 1, bool jet = 1, bool isDptcut = 1, bool isEff = 0,
TString outHistDir = "out", TString outPlotDir = "plots",
bool isjetptcut = 0, double jetptmin = 5, double jetptmax = 30 )
{

  /*  gStyle->SetPadBottomMargin(0.1); //margins...
    gStyle->SetPadTopMargin(0.04);
    gStyle->SetPadLeftMargin(0.08);
    gStyle->SetPadRightMargin(0.04);*/
    gStyle->SetOptStat(0000);

    TString outh = outPlotDir;
    outh += "/plots";

    gSystem->Exec(Form("mkdir %s",outh.Data()));

    int simNr = 0; // 0 - central value
    int nFiles = 0;
    if(quark) nFiles = fBsimN;
    else nFiles = fCsimN;
    TH1D **hSpectrum = new TH1D*[nFiles];
    TH1D **hSpectrum_binned = new TH1D*[nFiles];

    for (int nr=simNr; nr<nFiles; nr++){
        TString file = outHistDir;
        if(jet) file += "/JetPt_";
        else file += "/DPt_";
      //  if(quark == 0) file += "charm";
      //  else if(quark == 1) file += "beauty";
      //  file += "_";
        if(quark == 0) file += fRunC[nr];
        else if(quark == 1) file += fRunB[nr];
        if(jet) {
             if(isDptcut) { file += "_Dpt"; file += fptbinsDA[0]; file += "_"; file += fptbinsDA[fptbinsDN];  }
        }
        else {
            if(isjetptcut) { file += "_Jetpt"; file +=  jetptmin; file += "_"; file += jetptmax;  }
        }
        if(jet && isEff) file += "_effScaled";
        if(fDmesonSpecie) file += "_Dstar";
        else file += "_Dzero";
        file += ".root";
        TH1D *htmp;
        htmp = (TH1D*) GetInputHist(file, "hPt", htmp);
        //htmp->Scale(sigma_c[nr]);
        htmp->GetYaxis()->SetTitle("d#sigma/dp_{T} (mb)");
        hSpectrum[nr] = (TH1D*)htmp->Clone(Form("hSpectrum_%d",nr));
        hSpectrum_binned[nr] = (TH1D*)htmp->Rebin(fptbinsJetTrueN,Form("hSpectrum_binned_%d",nr),fptbinsJetTrueA);

    }

    int cent = 0;
    TH1D *htmp = (TH1D*)(hSpectrum[cent]->Clone("htmp"));
    TH1D *hSpectrum_central = (TH1D*)htmp->Clone("hSpectrum_central");
    TH1D *hSpectrum_central_binned = (TH1D*)htmp->Rebin(fptbinsJetTrueN,"hSpectrum_central_binned",fptbinsJetTrueA);

    setHistoDetails(hSpectrum_central,4,24);
    setHistoDetails(hSpectrum_central_binned,4,24);

    hSpectrum_central->Scale(simScaling);
    hSpectrum_central_binned->Scale(simScaling);
    hSpectrum_central_binned->Scale(1,"width");
    if(jet) hSpectrum_central_binned->GetXaxis()->SetTitle("p_{T}^{ch,jet} (GeV/c)");
    else hSpectrum_central_binned->GetXaxis()->SetTitle("p_{T}^{D*} (GeV/c)");

    // ----------------- up/down bands ---------------------
    // get up unc
    TH1D *hSpectrum_up = (TH1D*)hSpectrum_central_binned->Clone("hSpectrum_up");
    hSpectrum_up = (TH1D*)GetUpSys(hSpectrum_binned,nFiles,hSpectrum_up);
    setHistoDetails(hSpectrum_up,4,24,0,2,2);
    hSpectrum_up->Scale(simScaling);
    hSpectrum_up->Scale(1,"width");
    // get down unc
    TH1D *hSpectrum_down = (TH1D*)hSpectrum_central_binned->Clone("hSpectrum_down");
    hSpectrum_down = (TH1D*)GetDownSys(hSpectrum_binned,nFiles,hSpectrum_down);
    setHistoDetails(hSpectrum_down,4,24,0,2,2);
    hSpectrum_down->Scale(simScaling);
    hSpectrum_down->Scale(1,"width");

     // ======================= compare spectra ============================
     if(quark == 0) {
         outh+="/Promptspectra_";
         if(jet) {
             outh += "JetPt";
             if(isDptcut) { outh += "_Dpt"; outh += fptbinsDA[0]; outh += "_"; outh += fptbinsDA[fptbinsDN];  }
             }
         else {
             outh += "DPt";
             if(isjetptcut) { outh += "_Jetpt"; outh +=  jetptmin; outh += "_"; outh += jetptmax;  }
         }

      }
      else if(quark == 1) {
           outh+="/NonPromptspectra_";
           if(jet) {
               outh += "JetPt";
               if(isDptcut) { outh += "_Dpt"; outh += fptbinsDA[0]; outh += "_"; outh += fptbinsDA[fptbinsDN];  }
               if(isEff) outh += "_effScaled";
               }
           else {
               outh += "DPt";
               if(isjetptcut) { outh += "_Jetpt"; outh +=  jetptmin; outh += "_"; outh += jetptmax;  }
           }
      }

      compareSpectra(hSpectrum_binned,nFiles, jet, quark, outh);

     // compare data and central sim with unc
      TCanvas *cSpectra = new TCanvas("cSpectra","cSpectra",1000,800);
      cSpectra->SetLogy();
      hSpectrum_central_binned ->Draw();
      hSpectrum_up->Draw("same");
      hSpectrum_down->Draw("same");

      SaveCanvas(cSpectra,outh+"_un");

    return;
}

TH1* GetUpSys(TH1D **hFD, const int nFiles, TH1D *hFD_up){

        double bin = 0, binerr = 0;
        double max = 0, maxerr = 0;
        for(int j=1; j<=fptbinsJetTrueN; j++ ){
            max = hFD[0]->GetBinContent(j);
            for(int i=1;i<nFiles;i++){
                if(hFD[i]->GetBinContent(j) > max){
                        max = hFD[i]->GetBinContent(j);
                       // maxerr = hFD[j]->GetBinError(j);
                }
            }
            hFD_up->SetBinContent(j,max);
            hFD_up->SetBinError(j,0);
        }
    return hFD_up;
}

TH1* GetDownSys(TH1D **hFD, const int nFiles, TH1D *hFD_down){

        double bin = 0, binerr = 0;
        double max = 0, maxerr = 0;
        for(int j=1; j<=fptbinsJetTrueN; j++ ){
            max = hFD[0]->GetBinContent(j);
            for(int i=1;i<nFiles;i++){
                if(hFD[i]->GetBinContent(j) < max){
                        max = hFD[i]->GetBinContent(j);
                      //  maxerr = hFD[j]->GetBinError(j);
                }
            }
            hFD_down->SetBinContent(j,max);
            hFD_down->SetBinError(j,0);
        }
    return hFD_down;
}

void compareSpectra(TH1D **hFD, const int nFiles, bool isjet, bool quark, TString cName){

    TH1D *hPt[nFiles];

    TCanvas *cB = new TCanvas("cB","cB",1200,800);
    cB->SetLogy();
    TLegend *leg = new TLegend(0.6,0.4,0.85,0.85);
    leg->SetBorderSize(0);
    for(int i=0;i<nFiles;i++){
        hPt[i] = (TH1D*)hFD[i]->Clone(Form("hPt_%d",i));
        //hJetPt_B[i] = (TH1D*)hJetPt_B[i]->Rebin(fptbinsJetTrueN,Form("hJetPt_B_%d",i),fptbinsJetTrueA);
        hPt[i]->Scale(1,"width");

        hPt[i]->SetName(Form("hPt_%d",i));
        hPt[i]->GetYaxis()->SetTitle("d#sigma/dp_{T} (mb)");
        //hJetPt_B[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
        setHistoDetails(hPt[i],fColors[i],fMarkers[i]);
        hPt[i]->SetLineColor(fColors[i]);
        if(isjet) hPt[i]->GetXaxis()->SetTitle("p_{T}^{ch jet} (GeV/c)");
        else hPt[i]->GetXaxis()->SetTitle("p_{T}^{D*} (GeV/c)");

        if(!i)hPt[i]->Draw("h");
        else hPt[i]->Draw("same");

        if(quark == 0) leg->AddEntry(hPt[i],fDescC[i].Data(),"p");
        else if(quark == 1) leg->AddEntry(hPt[i],fDescB[i].Data(),"p");
        leg->Draw("same");
    }

    SaveCanvas(cB,cName);

    TH1D *hCent = (TH1D*)hPt[0]->Clone("hCent");
    hCent->SetName("hCent");
    TH1D *hJetPtRatio[nFiles-1];

    TCanvas *cr = new TCanvas("cr","cr",1200,800);
    if(isjet) { plotmin = fptbinsJetTrueA[0]; plotmax = fptbinsJetTrueA[fptbinsJetTrueN]; }
    else { plotmin = fptbinsDA[0]; plotmax = fptbinsDA[fptbinsDN]; }
    int bin = (plotmax-plotmin);
    TH1D *hh = new TH1D("hh","",bin, plotmin,plotmax);
    hh->GetYaxis()->SetTitle("ratio");
    if(quark == 0) hh->GetYaxis()->SetRangeUser(0.2,3.5);
    else if(quark == 1) hh->GetYaxis()->SetRangeUser(0.5,2);
    if(isjet) hh->GetXaxis()->SetTitle("p_{T}^{ch jet} (GeV/c)");
    else hh->GetXaxis()->SetTitle("p_{T}^{D*} (GeV/c)");
    hh->Draw();

    TLegend *leg2 = new TLegend(0.2,0.85,0.48,0.95);
    leg2->SetBorderSize(0);

    TLegend *leg3 = new TLegend(0.5,0.65,0.85,0.95);
    leg3->SetBorderSize(0);

    for(int i=1;i<nFiles;i++){
        hJetPtRatio[i-1] = (TH1D*)hPt[i]->Clone(Form("hJetPtRatio_%d",i-1));
        hJetPtRatio[i-1]->Divide(hCent);
        hJetPtRatio[i-1]->SetLineStyle(2);
        hJetPtRatio[i-1]->SetLineWidth(3);
        hJetPtRatio[i-1]->Draw("h same");
    }

     for(int i=0;i<2;i++){
        if(quark == 0) leg2->AddEntry(hJetPtRatio[i],fDescC[i+1].Data(),"l");
        else if(quark == 1) leg2->AddEntry(hJetPtRatio[i],fDescB[i+1].Data(),"l");
    }

    for(int i=2;i<nFiles-1;i++) {
        if(quark == 0) leg3->AddEntry(hJetPtRatio[i],fDescC[i+1].Data(),"l");
        else if(quark == 1) leg3->AddEntry(hJetPtRatio[i],fDescB[i+1].Data(),"l");
    }

    TLine *line = new TLine(fptbinsJetTrueA[0],1,fptbinsJetTrueA[fptbinsJetTrueN],1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");

    leg2->Draw("same");
    leg3->Draw("same");

    SaveCanvas(cr,cName+"_ratio");
}

TH1* GetInputHist(TString inFile, string histName,TH1 *hh){

	TFile *jetPtFile = new TFile(inFile,"read");
    hh = (TH1F*)jetPtFile->Get(histName.c_str());

    return hh;

}

void ScaleHist(TH1 *hh, int full){

    if(full){
         hh->Scale(1,"width");
        hh->GetYaxis()->SetTitle("d^{2}#sigma/d#etadp_{T} (mb #it{c}/GeV)");
    }
    else {
        hh->Scale(1,"width");
        hh->GetYaxis()->SetTitle("dN/dp_{T} (#it{c}/GeV)");
    }

}

void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, int Msize, Width_t Lwidth, Style_t Lstyle){

    hh->SetMarkerColor(color);
    hh->SetMarkerStyle(Mstyle);;
    hh->SetLineColor(color);
    hh->SetLineWidth(Lwidth);
    hh->SetMarkerSize(Msize);
    hh->SetLineStyle(Lstyle);

    //hh->SetTitle();
    //hh->GetXaxis()->SetTitle("p_{T}^{ch,jet} (GeV/c)");

}

void SaveCanvas(TCanvas *c, TString name){

    c->SaveAs(Form("%s.png",name.Data()));
    c->SaveAs(Form("%s.pdf",name.Data()));

}
