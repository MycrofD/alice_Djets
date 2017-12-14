#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <TCanvas.h>

#endif

void FitReflDistr(TString fitType="DoubleGaus"){//DoubleGaus

  TFile *fReflections = TFile::Open("ToyInvariantMassDistributions_pPb_minbias_cent_notopo_StepSmearing0.root");//ToyInvariantMassDistributions_5tev_central_notopo_StepSmearing0.root");
  TH1F *hfitRefl=0x0;
  
  TString fileName = "reflections_fitted_notopo";
  fileName += fitType;
  fileName += ".root";

  TFile *fFitReflection = new TFile(fileName.Data(), "recreate");
 
  TCanvas *cy = new TCanvas("fitCanv", "fitCanv");
  cy->Divide(4,3);

  TCanvas *cy2 = new TCanvas("fitCanv2", "fitCanv2");
  cy2->Divide(4,3);
  
  TCanvas *cyRatio = new TCanvas("fitCanvRatio", "fitCanvRatio");
  cyRatio->Divide(4,3);
  
  for (Int_t iBin=0; iBin<13; iBin++){
    cy->cd(iBin+1);
    hfitRefl= (TH1F*)fReflections->Get(Form("hRflMCScaled_ptBin%d_stepSmearing0", iBin));
    hfitRefl->SetMarkerStyle(1);
    hfitRefl->SetLineStyle(1);
    //  hfitRefl->Rebin(4);
    hfitRefl->Draw();

    //    Double_t minBin = hfitRefl->GetXaxis()->GetLowEdge(1);
    //    Double_t maxBin = hfitRefl->GetXaxis()->GetHighEdge(hfitRefl->GetNbinsX());

    TH1F *hFitReflNewTemp = (TH1F*)hfitRefl->Clone(Form("hRflMCScaledFitted%s_ptBin%d_stepSmearing0", fitType.Data(), iBin));
    TH1F *ratio = (TH1F*)hfitRefl->Clone(Form("ratioRelDistr_%s_bin%d", fitType.Data(), iBin));
    //    TF1 *finput = new TF1("finput", "gaus(0)*gaus(3)");

    TF1 *finput=0x0;
    if (fitType == "DoubleGaus"){//DoubleGaus
      finput= new TF1 ("finput","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]/( TMath::Sqrt(2.*TMath::Pi())*[5])*TMath::Exp(-(x-[4])*(x-[4])/(2.*[5]*[5]))");
      
      finput->SetParameter(0, 1);
      finput->SetParameter(1, 1);
      finput->SetParameter(2, 1);
      finput->SetParameter(3, 1);
      finput->SetParameter(4, 1);
      finput->SetParameter(5, 1);
  
    }
    else if (fitType == "pol3"){
      finput= new TF1 ("finput", "pol3");
      }

    else if (fitType == "pol6"){
      finput= new TF1 ("finput", "pol6");
    }

    else if (fitType == "gaus"){
      finput= new TF1 ("finput", "gaus");
    }
    
    hfitRefl->Fit("finput", "MLF");
    TF1 *fFitRefl = hfitRefl->GetFunction("finput");
    fFitReflection->cd();
//    TCanvas *ciao=new TCanvas(Form("fitPt_%d",iBin), Form("fitPt_%d",iBin));
//    hfitRefl->Draw();
//    fFitRefl->Draw("same");
       
    for(Int_t iBin2=1; iBin2<=hfitRefl->GetNbinsX(); iBin2++){
      cout<< " iBin2 " << iBin2 << endl;
      hFitReflNewTemp->SetBinContent(iBin2, 0.);
      ratio->SetBinContent(iBin2, 0.);

      hFitReflNewTemp->SetBinContent(iBin2, fFitRefl->Eval(hfitRefl->GetBinCenter(iBin2)));
      ratio->SetBinContent(iBin2, (hfitRefl->GetBinContent(iBin2) / fFitRefl->Eval(hfitRefl->GetBinCenter(iBin2))));
 
    }

    cy2->cd(iBin+1);
    hFitReflNewTemp->Draw();
    
    gStyle->SetOptFit(1111);
    cyRatio->cd(iBin+1);
    ratio->GetYaxis()->SetRangeUser(-1.5, 3.);
    ratio->Draw();
    ratio->Fit("pol0", "FM");
   
    //fFitReflection->cd();
    hFitReflNewTemp->Write();
    ratio->Write();

    delete finput;
  }

  //fFitReflection->Close();
}

