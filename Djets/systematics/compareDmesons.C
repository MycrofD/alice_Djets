
#include <string>
#include <sstream>
#include <iostream>

#include "sys.h"
#include "config.h"

Int_t colors2[] = {1,2,kGreen+3,kMagenta+2,4,6,kCyan+1,8,kOrange-1,kGray+1,kViolet+5,kYellow+2};
Int_t markers2[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};
Int_t linestyle2[] = {1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};


const int nFiles = 2;
TString inDirData[nFiles] = {
  "DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue5_50_ppbinning/unfolding_Bayes_3/finalSpectra",
  "Preliminaryplots"
}

TString dataFile[nFiles] = {
  "JetPtSpectrum_final.root",
  "DstarJet_pPb_Preliminary.root"
}

TString histName[nFiles] = {
  "hData_binned",
  "dataPoints"
}

TString histNameSim[nFiles] = {
  "hPrompt_central_binned",
  "theory"
}

TString desc[nFiles] = {
  "D^{0} p-Pb 5 TeV",
  "D^{*} p-Pb 5 TeV"
};

double plotmin = 5, plotmax = 30;
const int ptbinsN = 6;
double ptbinsA[ptbinsN+1] = { 5,6,8,10,14,20,30 };


void compareDmesons(TString inDirBase = "/home/basia/Work/alice/analysis/pPb_run2/", TString outName = "/home/basia/Work/alice/analysis/pPb_run2/compareDmesons/")
{

  gSystem->Exec(Form("mkdir %s",outName.Data()));

 compareData(outName,inDirBase,"JetSpectraMesonComparison2");

return;

}

void compareData(TString inName, TString inDirBase, TString outHistName)
{

            TString out = inName; //= inDirBase;
          //  out += "/";
          //  out += inName;
            gStyle->SetOptStat(0000); //Mean and RMS shown
            gSystem->Exec(Form("mkdir %s",out.Data()));

            TString dirName[nFiles];
            for (int i=0; i<nFiles; i++){
                dirName[i] = inDirBase;
                dirName[i] += "/";
                dirName[i] += inDirData[i];
            }

            TFile *fproj[nFiles];
            for(int i=0; i<nFiles; i++) {
              fproj[i] = new TFile(Form("%s/%s",dirName[i].Data(),dataFile[i].Data()),"READ");
              if(!fproj[i]) { cout << "NO FILE !!!!!" << endl; return; }
            }

            TCanvas *cspec = new TCanvas("cspec","cspec",900,800);
            cspec->SetLogy();

            TLegend *leg = new TLegend(0.65,0.75,0.85,0.85);
            leg->SetBorderSize(0);

            TH1F *spec[nFiles];
            TH1F *specReb[nFiles];

            TH1F *specSim[nFiles];
            TH1F *specSimReb[nFiles];

            for(int i=0; i<nFiles; i++) {
                spec[i] = (TH1F*)fproj[i]->Get(Form("%s",histName[i].Data()));
                spec[i]->Sumw2();
              //  if(!i)spec[i] -> Scale(1,"width");
                spec[i]->SetTitle();
                spec[i]->SetLineColor(colors2[i]);
                spec[i]->SetMarkerColor(colors2[i]);
                spec[i]->SetMarkerStyle(markers2[i]);
                spec[i]->SetMinimum(0.001);

                specSim[i] = (TH1F*)fproj[i]->Get(Form("%s",histNameSim[i].Data()));
                specSimReb[i] = (TH1F*)specSim[i]->Rebin(ptbinsN,Form("specSimReb_%d",i),ptbinsA);

                specReb[i] = (TH1F*)spec[i]->Rebin(ptbinsN,Form("specReb_%d",i),ptbinsA);
                specReb[i]->SetTitle();
                specReb[i]->SetLineColor(colors2[i]);
                specReb[i]->SetMarkerColor(colors2[i]);
                specReb[i]->SetMarkerStyle(markers2[i]);
                specReb[i]->SetMinimum(0.001);
                //spec[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
              //  if(!i) spec[i]->Draw();
              //else spec[i]->Draw("same");
                leg->AddEntry(spec[i],desc[i].Data());
            }
            spec[1]->Draw();
            spec[0]->Draw("same");
            leg->Draw("same");

          cspec->SaveAs(Form("%s/%s.pdf",out.Data(),outHistName.Data()));
          cspec->SaveAs(Form("%s/%s.png",out.Data(),outHistName.Data()));

            TLegend *leg2 = new TLegend(0.6,0.7,0.85,0.85,"p-Pb @5TeV, D^{*+}/D^{0}");
            leg2->SetBorderSize(0);
            TCanvas *cspec2 = new TCanvas("cspec2","cspec2",1000,600);
            TH1F *hratio[nFiles-1];
            TH1F *hratioSim[nFiles-1];

            for(int i=0; i<nFiles-1; i++){
                hratio[i] = (TH1F*)specReb[i+1]->Clone( Form("hratio_%d",i));
                hratio[i]->Divide(specReb[0]);

                hratioSim[i] = (TH1F*)specSimReb[i+1]->Clone( Form("hratioSim_%d",i));
                hratioSim[i]->Divide(specSimReb[0]);
                hratioSim[i]->SetLineStyle(linestyle2[i]);
                hratioSim[i]->SetLineWidth(2);

                //hratio[i] = (TH1F*)spec[i+1]->Clone(Form("hratio_%d",i));
                //hratio[i]->Divide(spec[0]);
                hratio[i]->SetLineStyle(linestyle2[i]);
                hratio[i]->SetLineWidth(2);
              //  hratio[i]->GetXaxis()->SetRangeUser(ptbinsA[0],ptbinsA[ptbinsN]);
                hratio[i]->GetYaxis()->SetRangeUser(0.,1.4);
                hratio[i]->GetYaxis()->SetTitle("d^{2}#sigma/d#etadp_{T} D^{*+}-jet/D^{0}-jet");
              //  if(!i) hratio[i]->Draw("hist");
              //  else hratio[i]->Draw("samehist");
              if(!i) hratio[i]->Draw();
              else hratio[i]->Draw("same");
              //  leg2->AddEntry(hratio[i],desc[i+1].Data(),"l");
            }
            hratioSim[0]->Draw("same");
            leg2->AddEntry(hratio[0],"data","l");
            leg2->AddEntry(hratioSim[0],"POWHEG+PYTHIA6","l");
            leg2->Draw("same");

            TLine *line = new TLine(ptbinsA[0],0.476,ptbinsA[ptbinsN],0.476);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            //line->Draw("same");


            cspec2->SaveAs(Form("%s/%s_ratio.pdf",out.Data(),outHistName.Data()));
            cspec2->SaveAs(Form("%s/%s_ratio.png",out.Data(),outHistName.Data()));

            return;

}
