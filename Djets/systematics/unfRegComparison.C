
#include <string>
#include <sstream>
#include <iostream>

#include "sys.h"
#include "config.h"


void unfRangesComparison(int reg=3,  TString inDirBase = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts", TString input = "Default_jetMeas3_50_jetTrue3_50", bool isChain = 0, TString int measmin=3, int measmax=50, int truemin=5, int truemax=50)
{

  if(!isChain) {
    if (truemin == 3){
      nJetBins = 9;
      ptJetbins = new double[nJetBins+1];
      for(int i=0;i<nJetBins+1;i++) ptJetbins[i] = ptbins[i];
    }
    else if (truemin == 4){
      nJetBins = 8;
      ptJetbins = new double[nJetBins+1];
      for(int i=0;i<nJetBins+1;i++) ptJetbins[i] = ptbins[i+1];
    }
    else if (truemin == 5){
      nJetBins = 7;
      ptJetbins = new double[nJetBins+1];
      for(int i=0;i<nJetBins+1;i++) ptJetbins[i] = ptbins[i+2];
    }
    else {
      cout << "WRONG true minimum pT !!!!" << endl;
      return;
    }
  }
  else {
    nJetBins = fptbinsJetTrueN;
    ptJetbins = new double[nJetBins+1];
    for(int i=0;i<nJetBins+1;i++) ptJetbins[i] = fptbinsJetTrueA[i];
  }

  TString inDir = inDirBase;
  inDir += "/";
  inDir += input;
  gSystem->Exec(Form("mkdir %s/systematics",inDir.Data()));

    compareReg(inDir);

return;
}

void compareReg(TString inDir = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50")
{

            gStyle->SetOptStat(0000); //Mean and RMS shown
            gSystem->Exec(Form("mkdir %s/systematics",inDir.Data()));

            const int nFiles = 2;
            TString dirName[nFiles];
            int regList[nFiles] = {3,4};

            for (int i=0; i<nFiles; i++){
                dirName[i] = inDir;
                dirName[i] += "/unfolding_Bayes_";
                dirName[i] += regList[i];
            }

            TString desc[nFiles];
            for(int i=0; i<nFiles; i++){
              desc[i] = "reg=";
              desc[i] += regList[i];

            }

            double plotmin = ptJetbins[0], plotmax = ptJetbins[nJetBins];

            TFile *fproj[nFiles];
            for(int i=0; i<nFiles; i++) fproj[i] = new TFile(Form("%s/unfoldedSpectrum_unfoldedJetSpectrum.root",dirName[i].Data()),"READ");

            TCanvas *cspec = new TCanvas("cspec","cspec",800,600);
            cspec->SetLogy();

            TLegend *leg = new TLegend(0.5,0.6,0.85,0.8);
            leg->SetBorderSize(0);

            TH1F *spec[nFiles];
            for(int i=0; i<nFiles; i++) {
                spec[i] = (TH1F*)fproj[i]->Get("unfoldedSpectrum");
                spec[i]->Sumw2();
                spec[i] -> Scale(1,"width");
                spec[i]->SetTitle();
                spec[i]->SetLineColor(colors[i]);
                spec[i]->SetMarkerColor(colors[i]);
                spec[i]->SetMarkerStyle(markers[i]);
                spec[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
                if(!i) spec[i]->Draw();
                else spec[i]->Draw("same");
                leg->AddEntry(spec[i],desc[i].Data());
            }
            leg->Draw("same");

          cspec->SaveAs(Form("%s/systematics/RegComparison_reg.pdf",inDir.Data()));
          cspec->SaveAs(Form("%s/systematics/RegComparison_reg.png",inDir.Data()));
          TLegend *leg2 = new TLegend(0.55,0.75,0.9,0.85);

            leg2->SetBorderSize(0);
            TCanvas *cspec2 = new TCanvas("cspec2","cspec2",800,400);
            TH1F *hratio[nFiles-1];
            for(int i=0; i<nFiles-1; i++){
                hratio[i] = (TH1F*)spec[i+1]->Clone( Form("hratio_%d",i));
                hratio[i]->Divide(spec[0]);
                hratio[i]->SetLineStyle(linestyle[i]);
                hratio[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
                hratio[i]->GetYaxis()->SetRangeUser(0.8,1.15);
                hratio[i]->GetYaxis()->SetTitle(Form("ratio to central (%s)",desc[0].Data()));
                if(!i) hratio[i]->Draw("hist");
                else hratio[i]->Draw("histsame");
                leg2->AddEntry(hratio[i],desc[i+1].Data());
            }
            leg2->Draw("same");

            TLine *line = new TLine(plotmin,1,plotmax,1);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw("same");


            cspec2->SaveAs(Form("%s/systematics/RegComparison_reg_ratio.pdf",inDir.Data()));
            cspec2->SaveAs(Form("%s/systematics/RegComparison_reg_ratio.png",inDir.Data()));


            
            return;

}


void getRMS(const int nFiles, TH1F **hratio, TH1F *hmean, TH1F *hsys)
{

  //TH1D *hsys = new TH1D("hsys","syst. rms; p_{T,ch jet};  sys [%] (rms)",nJetBins,ptJetbins);
  hsys->SetTitle();
  hsys->SetLineColor(1);
  hsys->SetLineWidth(2);
  hsys->SetLineStyle(2);
  hsys->SetMaximum(0.12);
  hsys->GetXaxis()->SetLabelSize(0.05);
  hsys->GetXaxis()->SetTitleSize(0.05);
  hsys->GetYaxis()->SetTitleSize(0.05);
  hsys->GetYaxis()->SetLabelSize(0.05);
  hsys->GetYaxis()->SetTitleOffset(0.8);

  //hmean = (TH1F*)hsys->Clone("hmean");
  hmean->GetYaxis()->SetTitle("mean");
//  hmean->GetYaxis()->SetRangeUser(0.95,1.1);
  hmean->SetMarkerStyle(20);
  hmean->SetLineStyle(1);
  hmean->SetTitle();

  double *rms = new double[nJetBins];
  double *mean = new double[nJetBins];
  for(int i=0; i<nJetBins; i++){
      rms[i] = 0;
       mean[i] = 0;
       for (int j=0; j<nFiles-1; j++){
         mean[i] = mean[i]+ ( hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) );
      //double m = ( hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) );
      //rms[i] = rms[i]+ ( 1-hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) ) * ( 1-hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) ) ;
      }
   mean[i] = mean[i]/(double)(nFiles-1);

   for (int j=0; j<nFiles-1; j++){
      //mean[i] = mean[i]+ ( hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) );
      //double m = ( hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) );
    //  rms[i] = rms[i]+ ( mean[i]-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) ) * ( mean[i]-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) ) ;
      rms[i] = rms[i]+ ( 1-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) ) * ( 1-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) ) ;


  }
      rms[i] = sqrt(rms[i]/(double)(nFiles-1));
      hsys->SetBinContent(i+1,rms[i]*100);

      hmean->SetBinContent(i+1,mean[i]);
      cout << "RMS pT " << (ptJetbins[i]+ptJetbins[i+1])/2. << " GeV/c:\t" << rms[i]*100 << endl;
      cout << "Mean pT " << (ptJetbins[i]+ptJetbins[i+1])/2. << " GeV/c:\t" << mean[i] << endl;
  }

}
