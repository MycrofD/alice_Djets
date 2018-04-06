//
// Macro to extract jet pt spectra from simulation, prompt or non-prompt D
// 
// Usage:
// .L simSpectra.C
// and then call the needed methods
//
// Author: B.Trzeciak (barbara.antonina.trzeciak@cern.ch)
//



#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TProfile.h"
#include <iostream>
#include "TPaveText.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TH3F.h"
#include "TH1.h"
#include "TH2.h"

#include <iostream>
#include <string>

using namespace std;


int Rpar = 4;
double Dptcut = 3.;
double jetptmin = 5, jetptmax = 50;
double Dptmin = 3, Dptmax = 36;
double jetEta = 0.5;

const int nJetBins = 9;
double ptJetbins[nJetBins+1] = { 3,4,5,6,8,10,14,20,30,50 };

//const int nJetBins = 8;
//double ptJetbins[nJetBins+1] = { 3,5,6,8,10,14,20,30,50 };

const Int_t nDbins = 10;
double ptDbins[nDbins+1] = {3,4,5,6,7,8,10,12,16,24,36};

double plotmin = 0, plotmax = 50;

Int_t colors[] = {1,2,8,4,kOrange-1,6,kGray+1,kCyan+1,kMagenta+2,kViolet+5,kYellow+2};
Int_t markers[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};


double BR = 0.0393; // 0.0257;
const int APb = 208;
const double simScaling =  APb; 

double pPbscaling = 2.1/208.; // not used now
 
TString runC[] = { "1496999831", "1497979153", "1497983041", "1497891573", "1497817179", "1497857391", "1497726041", "1497699194", "1497713588" };
TString desc_C[] = { "central", "m_{c}=1.3", "m_{c}=1.7", "#mu_{R}=2,#mu_{F}=2", "#mu_{R}=1,#mu_{F}=2", "#mu_{R}=2,#mu_{F}=1", "#mu_{R}=0.5,#mu_{F}=0.5", "#mu_{R}=1,#mu_{F}=0.5" ,"#mu_{R}=0.5,#mu_{F}=1" };

TString runB[] = { "1496995621", "1497227464", "1497228262", "1497207414", "1497173624", "1497172782", "1497132057", "1497130041", "1497121734" };
TString desc_B[] = {"central", "m_{b}=4.5", "m_{b}=5", "#mu_{R}=2,#mu_{F}=2", "#mu_{R}=1,#mu_{F}=2", "#mu_{R}=2,#mu_{F}=1", "#mu_{R}=0.5,#mu_{F}=0.5", "#mu_{R}=1,#mu_{F}=0.5" ,"#mu_{R}=0.5,#mu_{F}=1" };


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

void plotSimSpectra(int quark = 1, bool jet = 1, bool isjetptcut = 0, bool isDptcut = 1, bool isEff = 0, TString outPlotDir = "plotsD0", TString outHistName = "outD0/ptSpectrumSim_" ){
    
    
    setStyle();
    
    gStyle->SetPadBottomMargin(0.1); //margins...
    gStyle->SetPadTopMargin(0.04);
    gStyle->SetPadLeftMargin(0.08);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetOptStat(0000);
    

    int simNr = 0; // 0 - central value
    const int nFiles = 9; // 10 - all files
    TH1D *hSpectrum[nFiles];
    TH1D *hSpectrum_binned[nFiles];

    for (int nr=simNr; nr<nFiles; nr++){
        TString file = outHistName;
        if(jet) file += "JetPt_";
        else file += "DPt_";
        if(quark == 0) file += "charm";
        else if(quark == 1) file += "beauty";
        file += "_";
        if(quark == 0) file += runC[nr];
        else if(quark == 1) file += runB[nr];
        if(jet) { 
             if(isDptcut) { file += "_Dpt"; file += Dptmin; file += "_"; file += Dptmax;  }
        }
        else { 
            if(isjetptcut) { file += "_Jetpt"; file +=  jetptmin; file += "_"; file += jetptmax;  }
        }
        if(jet && isEff) file += "_effScaled";
        file += ".root";
        TH1D *htmp;
        htmp = (TH1D*) GetInputHist(file, "hPt", htmp);
        //htmp->Scale(sigma_c[nr]);
        htmp->GetYaxis()->SetTitle("d#sigma/dp_{T} (mb)");
        hSpectrum[nr] = (TH1D*)htmp->Clone(Form("hSpectrum_%d",nr));
        hSpectrum_binned[nr] = (TH1D*)htmp->Rebin(nJetBins,Form("hSpectrum_binned_%d",nr),ptJetbins);
    
    }
    
    int cent = 0;
    TH1D *htmp = (TH1D*)(hSpectrum[cent]->Clone("htmp"));
    TH1D *hSpectrum_central = (TH1D*)htmp->Clone("hSpectrum_central");
    TH1D *hSpectrum_central_binned = (TH1D*)htmp->Rebin(nJetBins,"hSpectrum_central_binned",ptJetbins);
    
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
     TString outh = outPlotDir;
     if(quark == 0) { 
         outh+="/Promptspectra_";  
         if(jet) { 
             outh += "JetPt"; 
             if(isDptcut) { outh += "_Dpt"; outh += Dptmin; outh += "_"; outh += Dptmax;  }
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
             if(isDptcut) { outh += "_Dpt"; outh += Dptmin; outh += "_"; outh += Dptmax;  }
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
        
    
        for(int j=1; j<=nJetBins; j++ ){
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
        
    
        for(int j=1; j<=nJetBins; j++ ){
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
        //hJetPt_B[i] = (TH1D*)hJetPt_B[i]->Rebin(nJetBins,Form("hJetPt_B_%d",i),ptJetbins);
        hPt[i]->Scale(1,"width");
        
        hPt[i]->SetName(Form("hPt_%d",i));
        hPt[i]->GetYaxis()->SetTitle("d#sigma/dp_{T} (mb)");
        //hJetPt_B[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
        setHistoDetails(hPt[i],colors[i],markers[i]);
        hPt[i]->SetLineColor(colors[i]);
        if(isjet) hPt[i]->GetXaxis()->SetTitle("p_{T}^{ch jet} (GeV/c)");
        else hPt[i]->GetXaxis()->SetTitle("p_{T}^{D*} (GeV/c)");
        
        if(!i)hPt[i]->Draw("h");
        else hPt[i]->Draw("same");
        
        if(quark == 0) leg->AddEntry(hPt[i],desc_C[i].Data(),"p");
        else if(quark == 1) leg->AddEntry(hPt[i],desc_B[i].Data(),"p");
        leg->Draw("same");
    }
   
    SaveCanvas(cB,cName);
    
    TH1D *hCent = (TH1D*)hPt[0]->Clone("hCent");
    hCent->SetName("hCent");
    TH1D *hJetPtRatio[nFiles-1];
    
    TCanvas *cr = new TCanvas("cr","cr",1200,800);
    if(isjet) { plotmin = ptJetbins[0]; plotmax = ptJetbins[nJetBins]; }  
    else { plotmin = ptDbins[0]; plotmax = ptDbins[nDbins]; }  
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
        hJetPtRatio[i-1] -> Divide(hCent);
        hJetPtRatio[i-1] -> SetLineStyle(2);
        hJetPtRatio[i-1] -> SetLineWidth(3);
        hJetPtRatio[i-1] -> Draw("h same");
    }
    
     for(int i=0;i<2;i++){
        if(quark == 0) leg2->AddEntry(hJetPtRatio[i],desc_C[i+1].Data(),"l");
        else if(quark == 1) leg2->AddEntry(hJetPtRatio[i],desc_B[i+1].Data(),"l");
    }
    
    for(int i=2;i<nFiles-1;i++) {
        if(quark == 0) leg3->AddEntry(hJetPtRatio[i],desc_C[i+1].Data(),"l");
        else if(quark == 1) leg3->AddEntry(hJetPtRatio[i],desc_B[i+1].Data(),"l");
    }
    
    TLine *line = new TLine(3,1,40,1);
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
        //hh->Scale(pPbscaling);
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
    hh->SetTitle();
   
    //hh->SetTitle();
    //hh->GetXaxis()->SetTitle("p_{T}^{ch,jet} (GeV/c)");
    
}

void SaveCanvas(TCanvas *c, TString name){
    
    c->SaveAs(Form("%s.png",name.Data()));
    c->SaveAs(Form("%s.pdf",name.Data()));
    
}

void setStyle(){

    gStyle->SetOptStat(000);
    gStyle->SetLegendFont(42);
    gStyle->SetTextFont(22) ;
    //gStyle->SetLegendTextSize(0.05);
    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetPadRightMargin(0.02);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetTitleOffset(1.,"x");
    gStyle->SetTitleOffset(0.9.,"y");
    gStyle->SetTitleSize(0.04,"xyz");
    gStyle->SetLabelSize(0.03,"xyz");


}
