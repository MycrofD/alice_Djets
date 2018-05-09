
#include <string>
#include <sstream>
#include <iostream>

#include "sys.h"
#include "config.h"

Int_t colors2[] = {1,2,4,kGreen+3,kMagenta+2,4,6,kCyan+1,8,kOrange-1,kGray+1,kViolet+5,kYellow+2};
Int_t markers2[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};
Int_t linestyle2[] = {1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};


const int nFiles = 3;
/*TString inDirData[nFiles] = {
  "DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50_PbPbbinning/unfolding_Bayes_3",
  "DzeroR03_RefDPt3PythiaEff_Dpt321_BaseCuts/Default_jetMeas3_50_jetTrue3_50_PbPbbinning/unfolding_Bayes_3"
}*/

TString inDirData[nFiles] = {
  "DzeroR03_RefDPt324PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50_ppbinning/unfolding_Bayes_3",
  "DzeroR03_RefDPt424PythiaEff_BaseCuts/Default_jetMeas4_50_jetTrue4_50_ppbinning/unfolding_Bayes_3",
  "DzeroR03_RefDPt224PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50_ppbinning/unfolding_Bayes_3"
}


TString inDirSim[nFiles] = {
  "DzeroR03_RefDPt3PythiaEff_BaseCuts/Simulations/Prompt",
  "DzeroR03_RefDPt3PythiaEff_Dpt320_BaseCuts/Simulations/Prompt",
  "DzeroR03_RefDPt3PythiaEff_Dpt320_BaseCuts/Simulations/Prompt"
}

TString simFile[nFiles] = {
  "JetPt_AnalysisResults_FastSim_powheg+pythia6_charm_1520422975_Dpt3_36_Dzero.root",
  "JetPt_AnalysisResults_FastSim_powheg+pythia6_charm_1520422975_Dpt3_20_Dzero.root",
  "JetPt_AnalysisResults_FastSim_powheg+pythia6_charm_1520422975_Dpt3_20_Dzero.root"
}


TString desc[nFiles] = {
  "3 < D p_{T} < 24 GeV/c",
  "4 < D p_{T} < 24 GeV/c",
  "2 < D p_{T} < 24 GeV/c"
};



double plotmin = 5, plotmax = 30;
//const int ptbinsN = 7;
//double ptbinsA[ptbinsN+1] = { 3,5,10,15,20,25,35,50 };

const int ptbinsN = 6;
double ptbinsA[ptbinsN+1] = { 5,6,8,10,14,20,30 };

void compareResults(TString inDirBase = "/home/basia/Work/alice/analysis/pp5TeV/", TString outName = "/home/basia/Work/alice/analysis/pp5TeV/compareDPtCuts")
{

  gSystem->Exec(Form("mkdir %s",outName.Data()));

 compareData(outName,inDirBase,"JetSpectraComparisonData_Dcut3to4");
 //compareSim(outName,inDirBase,"JetSpectraComparisonSim_Dcut");

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
              fproj[i] = new TFile(Form("%s/unfoldedSpectrum_unfoldedJetSpectrum.root",dirName[i].Data()),"READ");
              if(!fproj[i]) { cout << "NO FILE !!!!!" << endl; return; }
            }

            TCanvas *cspec = new TCanvas("cspec","cspec",800,600);
            cspec->SetLogy();

            TLegend *leg = new TLegend(0.5,0.65,0.85,0.85,"pp data");
            //TLegend *leg = new TLegend(0.5,0.65,0.85,0.85,"p-Pb data");
            leg->SetBorderSize(0);

            TH1F *spec[nFiles];
            TH1F *specReb[nFiles];

            for(int i=0; i<nFiles; i++) {
                spec[i] = (TH1F*)fproj[i]->Get("unfoldedSpectrum");
                spec[i]->Sumw2();
                spec[i] -> Scale(1,"width");
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

          cspec->SaveAs(Form("%s/%s.pdf",out.Data(),outHistName.Data()));
          cspec->SaveAs(Form("%s/%s.png",out.Data(),outHistName.Data()));

          //TLegend *leg2 = new TLegend(0.55,0.75,0.85,0.85,"p-Pb data");
            TLegend *leg2 = new TLegend(0.55,0.75,0.85,0.85,"pp data");
            leg2->SetBorderSize(0);
            TCanvas *cspec2 = new TCanvas("cspec2","cspec2",1200,800);
            TH1F *hratio[nFiles-1];
            for(int i=0; i<nFiles-1; i++){
                hratio[i] = (TH1F*)specReb[i+1]->Clone( Form("hratio_%d",i));
                  hratio[i]->Divide(specReb[0]);
              //  hratio[i] = (TH1F*)spec[i+1]->Clone(Form("hratio_%d",i));
              //  hratio[i]->Divide(spec[0]);
                hratio[i]->SetLineStyle(linestyle2[i]);
                hratio[i]->SetLineWidth(2);
              //  hratio[i]->GetXaxis()->SetRangeUser(ptbinsA[0],ptbinsA[ptbinsN]);
                hratio[i]->GetYaxis()->SetRangeUser(0.6,1.4);
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


            cspec2->SaveAs(Form("%s/%s_ratio.pdf",out.Data(),outHistName.Data()));
            cspec2->SaveAs(Form("%s/%s_ratio.png",out.Data(),outHistName.Data()));


            return;

}


void compareSim(TString inName, TString inDirBase, TString outHistName)
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
                dirName[i] += inDirSim[i];
            }

            TFile *fproj[nFiles];
            for(int i=0; i<nFiles; i++) {
              fproj[i] = new TFile(Form("%s/%s",dirName[i].Data(),simFile[i].Data()),"READ");
              if(!fproj[i]) { cout << "NO FILE !!!!!" << endl; return; }
            }

            TCanvas *cspec = new TCanvas("cspec","cspec",800,600);
            cspec->SetLogy();

            TLegend *leg = new TLegend(0.5,0.65,0.85,0.85,"p-Pb POWHEG+PYTHIA6");
            leg->SetBorderSize(0);

            TH1F *spec[nFiles];
            TH1F *specReb[nFiles];

            for(int i=0; i<nFiles; i++) {
                spec[i] = (TH1F*)fproj[i]->Get("hPt");
                spec[i]->Sumw2();
                spec[i] -> Scale(1,"width");
                spec[i]->SetTitle();
                spec[i]->SetLineColor(colors2[i]);
                spec[i]->SetMarkerColor(colors2[i]);
                spec[i]->SetMarkerStyle(markers2[i]);

                specReb[i] = (TH1F*)spec[i]->Rebin(ptbinsN,Form("specReb_%d",i),ptbinsA);
                specReb[i]->SetTitle();
                specReb[i]->SetLineColor(colors2[i]);
                specReb[i]->SetMarkerColor(colors2[i]);
                specReb[i]->SetMarkerStyle(markers2[i]);

                specReb[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
                if(!i) specReb[i]->Draw();
                else specReb[i]->Draw("same");
                leg->AddEntry(specReb[i],desc[i].Data());
            }
            leg->Draw("same");

          cspec->SaveAs(Form("%s/%s.pdf",out.Data(),outHistName.Data()));
          cspec->SaveAs(Form("%s/%s.png",out.Data(),outHistName.Data()));

            TLegend *leg2 = new TLegend(0.55,0.75,0.85,0.85,"p-Pb POWHEG+PYTHIA6");
            leg2->SetBorderSize(0);
            TCanvas *cspec2 = new TCanvas("cspec2","cspec2",1200,800);
            TH1F *hratio[nFiles-1];
            for(int i=0; i<nFiles-1; i++){
                hratio[i] = (TH1F*)specReb[i+1]->Clone( Form("hratio_%d",i));
                hratio[i]->Divide(specReb[0]);
                hratio[i]->SetLineStyle(linestyle2[i]);
                hratio[i]->SetLineWidth(2);
              //  hratio[i]->GetXaxis()->SetRangeUser(ptbinsA[0],ptbinsA[ptbinsN]);
                hratio[i]->GetYaxis()->SetRangeUser(0.4,1.02);
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


            cspec2->SaveAs(Form("%s/%s_ratio.pdf",out.Data(),outHistName.Data()));
            cspec2->SaveAs(Form("%s/%s_ratio.png",out.Data(),outHistName.Data()));


            return;

}
