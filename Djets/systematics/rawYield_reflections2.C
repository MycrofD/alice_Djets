
#include <string>
#include <sstream>
#include <iostream>

#include "sys.h"
#include "config.h"

Int_t colors2[] = {1,2,kGreen+3,kMagenta+2,4,6,kCyan+1,8,kOrange-1,kGray+1,kViolet+5,kYellow+2};
Int_t markers2[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};
Int_t linestyle2[] = {1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};


const int nFiles = 3;
TString inDir[nFiles] = {
  "249","249_R05","249_R15"
}
/*
const int nFiles = 3;
TString inDir[nFiles] = {
  "Default_jetMeas3_50_jetTrue3_50_PbPbbinning/unfolding_Bayes_3",
  "Refl05_jetMeas3_50_jetTrue3_50_PbPbbinning/unfolding_Bayes_3",
  "Refl15_jetMeas3_50_jetTrue3_50_PbPbbinning/unfolding_Bayes_3"

}
*/
TString desc[nFiles] = {
  "default",
  "refl, -50%",
  "refl, +50%"
};

double plotmin = 3, plotmax = 50;
const int ptbinsN = 9;
double ptbinsA[ptbinsN+1] = { 3,4,5,6,8,10,14,20,30,50 };

int nJetBins2 = 7;
double ptJetbins2[] = { 5,6,8,10,14,20,30,50 };

void rawYield_reflections2(int reg=4,  TString inDirBase = "/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results/DzeroR03_pPbCuts", TString inName = "systematics/reflections", bool isChain = 0,  int measmin=3, int measmax=50, int truemin=5, int truemax=50)
//void rawYield_reflections(int reg=3,  TString inDirBase = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts", TString inName = "Default_jetMeas3_50_jetTrue3_50_PbPbbinning/systematics", bool isChain = 0,  int measmin=3, int measmax=50, int truemin=5, int truemax=50)
{

  TString inputDir = inDirBase;
//  inDir += "/";
//  inDir += input;
  gSystem->Exec(Form("mkdir %s/%s",inputDir.Data(),inName.Data()));

 yieldSys(inName,inDirBase,measmin,measmax,truemin,truemax);


return;
}

void yieldSys(TString inName, TString inDirBase, int measmin, int measmax, int truemin, int truemax)
{

            TString out = inDirBase;
            out += "/";
            out += inName;
            gStyle->SetOptStat(0000); //Mean and RMS shown
            gSystem->Exec(Form("mkdir %s",out.Data()));

            TString dirName[nFiles];

            for (int i=0; i<nFiles; i++){
                dirName[i] = inDirBase;
                dirName[i] += "/";
                dirName[i] += "Default_";
                dirName[i] += inDir[i];

            }

            TFile *outFile = new TFile(Form("%s/RawYield_reflections.root",out.Data()),"RECREATE");

            TFile *fproj[nFiles];
            for(int i=0; i<nFiles; i++) {
//              fproj[i] = new TFile(Form("%s/unfoldedSpectrum_unfoldedJetSpectrum.root",dirName[i].Data()),"READ");
              fproj[i] = new TFile(Form("%s/unfolding_Bayes_4/unfoldedSpectrum_unfoldedJetSpectrum.root",dirName[i].Data()),"READ");

              if(!fproj[i]) { cout << "NO FILE !!!!!" << endl; return; }
            }

cout<<"================---------------------------------===================="<<endl;
cout<<dirName[0].Data()<<endl;
	    TCanvas *cspec = new TCanvas("cspec","cspec",800,800);
            cspec->SetLogy();

            TLegend *leg = new TLegend(0.5,0.65,0.85,0.85);
            leg->SetBorderSize(0);

            TH1F *spec[nFiles];
            TH1F *specReb[nFiles];

            for(int i=0; i<nFiles; i++) {
                spec[i] = (TH1F*)fproj[i]->Get("unfoldedSpectrum");
                spec[i]->Sumw2();
                spec[i]->Scale(1,"width");
                spec[i]->SetTitle();
                spec[i]->SetLineColor(colors2[i]);
                spec[i]->SetMarkerColor(colors2[i]);
                spec[i]->SetMarkerStyle(markers2[i]);

                specReb[i] = new TH1F(Form("specReb%d",i),"specReb",ptbinsN,ptbinsA);
                for(int j=1;j<specReb[i]->GetNbinsX()+1;j++){
                    double pt = specReb[i]->GetBinCenter(j);
                    int bin = spec[i]->GetXaxis()->FindBin(pt);
                    double value = spec[i]->GetBinContent(bin);
                    double error = spec[i]->GetBinError(bin);
                    specReb[i]->SetBinContent(j,value);
                    specReb[i]->SetBinError(j,error);
                }

                specReb[i]->SetTitle();
                specReb[i]->SetLineColor(colors2[i]);
                specReb[i]->SetMarkerColor(colors2[i]);
                specReb[i]->SetMarkerStyle(markers2[i]);

                spec[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
                if(!i) spec[i]->Draw();
                else spec[i]->Draw("same");
                leg->AddEntry(spec[i],desc[i].Data());
            }
            leg->Draw("same");

          cspec->SaveAs(Form("%s/RawYield_reflections.pdf",out.Data()));
          cspec->SaveAs(Form("%s/RawYield_reflections.png",out.Data()));

            TLegend *leg2 = new TLegend(0.55,0.65,0.85,0.85);
            leg2->SetBorderSize(0);
            TCanvas *cspec2 = new TCanvas("cspec2","cspec2",1200,800);
            TH1F *hratio[nFiles-1];
            for(int i=0; i<nFiles-1; i++){
                hratio[i] = (TH1F*)spec[i+1]->Clone( Form("hratio_%d",i));
                hratio[i]->Divide(spec[0]);
                hratio[i]->SetLineStyle(linestyle2[i]);
                hratio[i]->SetLineWidth(2);
                hratio[i]->GetXaxis()->SetRangeUser(ptbinsA[0],ptbinsA[ptbinsN]);
                hratio[i]->GetYaxis()->SetRangeUser(0.92,1.1);
                hratio[i]->GetYaxis()->SetTitle(Form("ratio to central (%s)",desc[0].Data()));
                if(!i) hratio[i]->Draw("hist");
                else hratio[i]->Draw("samehist");
                leg2->AddEntry(hratio[i],desc[i+1].Data(),"l");
            }
            leg2->Draw("same");

            TLine *line = new TLine(ptbinsA[0],1,ptbinsA[ptbinsN],1);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw("same");

            cspec2->SaveAs(Form("%s/RawYield_reflections_ratio.pdf",out.Data()));
            cspec2->SaveAs(Form("%s/RawYield_reflections_ratio.png",out.Data()));


            //-------------- uncertanties
            TH1D *hsys = (TH1D*)*spec[0]->Clone("hsys");
            hsys->GetYaxis()->SetTitle("Refl. sys. unc [%]");
            hsys->SetLineColor(kViolet+1);
            hsys->SetLineStyle(2);
            hsys->SetLineWidth(2);
            hsys->GetYaxis()->SetRangeUser(0,4);

            for(int j=1; j<spec[0]->GetNbinsX()+1; j++ ){
                    double unc1 = spec[1]->GetBinContent(j) - spec[0]->GetBinContent(j);
                    double unc2 = spec[0]->GetBinContent(j) - spec[2]->GetBinContent(j);
		if (unc1<0)unc1 = -unc1;
		if (unc2<0)unc2 = -unc2;
                    double unc = 0;
                    if(unc1>unc2) unc = unc1;
                    else unc = unc2;
			cout<<"==========="<<unc<<"==="<<unc1<<"==="<<unc2<<endl;
                    unc /= spec[0]->GetBinContent(j);
                    hsys->SetBinContent(j,unc*100);
                    hsys->SetBinError(j,0);
            }
            TCanvas *cspecUnc = new TCanvas("cspecUnc","cspecUnc",800,400);
            hsys->Draw("hist");

            cspecUnc->SaveAs(Form("%s/RawYield_reflections_sys.pdf",out.Data()));
            cspecUnc->SaveAs(Form("%s/RawYield_reflections_sys.png",out.Data()));

            outFile->cd();
            cspec->Write();
            cspec2->Write();
            cspecUnc->Write();
            hsys->Write();

            outFile->Close();

            return;

}

void getRMS(const int nFiles, TH1F **hratio, TH1F *hmean, TH1F *hsys)
{

  //TH1D *hsys = new TH1D("hsys","syst. rms; p_{T,ch jet};  sys [%] (rms)",nJetBins2,ptJetbins2);
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

  double *rms = new double[nJetBins2];
  double *mean = new double[nJetBins2];
  for(int i=0; i<nJetBins2; i++){
      rms[i] = 0;
       mean[i] = 0;
       for (int j=0; j<nFiles-1; j++){
         mean[i] = mean[i]+ ( hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) );
         rms[i] = rms[i]+ ( 1-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) ) * ( 1-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) ) ;

      }
      mean[i] = mean[i]/(double)(nFiles-1);
      rms[i] = sqrt(rms[i]/(double)(nFiles-1));
      hsys->SetBinContent(i+1,rms[i]*100);

      hmean->SetBinContent(i+1,mean[i]);
      cout << "RMS pT " << (ptJetbins2[i]+ptJetbins2[i+1])/2. << " GeV/c:\t" << rms[i]*100 << endl;
      cout << "Mean pT " << (ptJetbins2[i]+ptJetbins2[i+1])/2. << " GeV/c:\t" << mean[i] << endl;
  }

}
