
#include <string>
#include <sstream>
#include <iostream>

#include "sys.h"
#include "config.h"

Int_t colors2[] = {2,4,kGreen+2,kMagenta+2,kOrange-1,6,kCyan+1,8,kGray+1,kViolet+5,kYellow+2};
Int_t markers2[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};
Int_t linestyle2[] = {1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};


const int nFiles = 6;
/*TString inDirData[nFiles] = {
  "Default_jetMeas3_50_jetTrue3_50_ppbinning/unfolding_Bayes_3",
  "Default_jetMeas3_50_jetTrue3_50_PbPbbinning/unfolding_Bayes_3",
  "Default_jetMeas3_40_jetTrue5_40_binning2/unfolding_Bayes_3"
}*/

TString inDirData[nFiles] = {
  "Default_jetMeas3_50_jetTrue3_50_ppbinning/unfolding_Bayes_3/finalSpectra",
  "Default_jetMeas3_50_jetTrue3_50_PbPbbinning/unfolding_Bayes_3/finalSpectra",
  "Default_jetMeas3_40_jetTrue3_40_binning2/unfolding_Bayes_3/finalSpectra",
  "Default_jetMeas3_50_jetTrue3_50_binning3/unfolding_Bayes_3/finalSpectra",
  "Default_jetMeas3_50_jetTrue3_50_binning4/unfolding_Bayes_3/finalSpectra",
  "Default_jetMeas3_50_jetTrue3_50_binning5/unfolding_Bayes_3/finalSpectra"
}

TString simFile[nFiles] = {
  "JetPt_AnalysisResults_FastSim_powheg+pythia6_charm_1520422975_Dpt3_36_Dzero.root",
  "JetPt_AnalysisResults_FastSim_powheg+pythia6_charm_1520422975_Dpt3_36_Dzero.root",
  "JetPt_AnalysisResults_FastSim_powheg+pythia6_charm_1520422975_Dpt3_36_Dzero.root"

}


TString desc[nFiles] = {
  "pp binning",
  "Pb-Pb binning",
  "binning 2",
  "binning 3",
  "binning 4",
  "binning 5"
};



double plotmin = 5, plotmax = 50;
const int ptbinsN = 7;
double ptbinsA[ptbinsN+1] = { 3,5,10,15,20,25,35,50 };


void compareSpectra(TString inDirBase = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts")
{

  TString outName = inDirBase;
  outName += "/compareBinning";
  gSystem->Exec(Form("mkdir %s",outName.Data()));

 compareData(outName,inDirBase,"JetSpectra_binningComparison");

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
              //fproj[i] = new TFile(Form("%s/unfoldedSpectrum_unfoldedJetSpectrum.root",dirName[i].Data()),"READ");
              fproj[i] = new TFile(Form("%s/JetPtSpectrum_final.root",dirName[i].Data()),"READ");
              if(!fproj[i]) { cout << "NO FILE !!!!!" << endl; return; }
            }

            TCanvas *cspec = new TCanvas("cspec","cspec",1000,800);
            cspec->SetLogy();
            TLegend *leg = new TLegend(0.5,0.65,0.85,0.85);
            leg->SetBorderSize(0);
            TH1F *spec[nFiles];
            TH1F *specReb[nFiles];
            for(int i=0; i<nFiles; i++) {
              //spec[i] = (TH1F*)fproj[i]->Get("unfoldedSpectrum");
                spec[i] = (TH1F*)fproj[i]->Get("hData_binned");
                spec[i]->Sumw2();
                //spec[i] -> Scale(1,"width");
                spec[i]->SetTitle();
                spec[i]->SetLineColor(colors2[i]);
                spec[i]->SetMarkerColor(colors2[i]);
                spec[i]->SetMarkerStyle(markers2[i]);
                //spec[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
                if(!i) spec[i]->Draw();
                else spec[i]->Draw("same");
                leg->AddEntry(spec[i],desc[i].Data());
            }
            leg->Draw("same");


            int simNr = 0; // 0 - central value
            int cent = 0;
            int nFiles2 = fCsimN;
            TH1D *hPrompt[fCsimN];
            TH1D *hPrompt_binned[fCsimN];
            TString simDir = inDirBase;
            simDir += "/Simulations/Prompt";
            double simScaling = 208./2.;

            for (int nr=simNr; nr<fCsimN; nr++){
                TString file = simDir;
                file += "/JetPt_";
                file += fRunC[nr];
                file += "_Dpt"; file += fptbinsDA[0]; file += "_"; file += fptbinsDA[fptbinsDN];
                if(fDmesonSpecie) file += "_Dstar";
                else file += "_Dzero";
                file += ".root";
                TH1D *htmp;
                htmp = (TH1D*) GetInputHist(file, "hPt", htmp);
                htmp->GetYaxis()->SetTitle("d#sigma/dp_{T} (mb)");
                hPrompt[nr] = (TH1D*)htmp->Clone(Form("hPrompt_%d",nr));
              //  hPrompt_binned[nr] = (TH1D*)htmp->Rebin(fptbinsJetTrueN,Form("hPrompt_binned_%d",nr),fptbinsJetTrueA);

            }

            TH1D *htmp = (TH1D*)(hPrompt[cent]->Clone("htmp"));
            TH1D *hPrompt_central = (TH1D*)htmp->Clone("hPrompt_central");

            // ----------------- prompt syst. (rebinned)---------------------
            // get up unc
              TH1D *hPrompt_up = (TH1D*)hPrompt_central->Clone("hPrompt_up");
              hPrompt_up = (TH1D*)GetUpSys(hPrompt,nFiles2,hPrompt_up);

              // get down unc
              TH1D *hPrompt_down = (TH1D*)hPrompt_central->Clone("hPrompt_down");
              hPrompt_down = (TH1D*)GetDownSys(hPrompt,nFiles2,hPrompt_down);


              hPrompt_central->SetMarkerColor(kGray+1);
              hPrompt_central->SetLineColor(kGray+1);
              hPrompt_central->SetMarkerStyle(24);
              hPrompt_up->SetMarkerColor(kGray+1);
              hPrompt_up->SetLineColor(kGray+1);
              hPrompt_down->SetMarkerColor(kGray+1);
              hPrompt_down->SetLineColor(kGray+1);

              hPrompt_up->Scale(simScaling);
              hPrompt_up->Scale(1,"width");
              hPrompt_up->Scale(1./1.2);
              hPrompt_down->Scale(simScaling);
              hPrompt_down->Scale(1,"width");
              hPrompt_down->Scale(1./1.2);
              hPrompt_central->Scale(simScaling);
              hPrompt_central->Scale(1,"width");
              hPrompt_central->Scale(1./1.2);


              hPrompt_central->Draw("epsame");
              hPrompt_up->Draw("same");
              hPrompt_down->Draw("same");

          cspec->SaveAs(Form("%s/%s.pdf",out.Data(),outHistName.Data()));
          cspec->SaveAs(Form("%s/%s.png",out.Data(),outHistName.Data()));


                      TCanvas *cspecRatio = new TCanvas("cspecRatio","cspecRatio",1200,800);
                      //cspec->SetLogy();
                      TLegend *leg2 = new TLegend(0.15,0.6,0.85,0.85,"ratio to theory");
                      leg2->SetBorderSize(0);
                      TH1F *specR[nFiles];

                      for(int i=0; i<nFiles; i++) {
                        //spec[i] = (TH1F*)fproj[i]->Get("unfoldedSpectrum");
                          specR[i] = (TH1F*)fproj[i]->Get("hData_binned_ratio");
                          specR[i]->Sumw2();
                          //spec[i] -> Scale(1,"width");
                          specR[i]->SetTitle();
                          specR[i]->SetLineColor(colors2[i]);
                          specR[i]->SetMarkerColor(colors2[i]);
                          specR[i]->SetMarkerStyle(markers2[i]);
                          specR[i]->GetYaxis()->SetRangeUser(0.2,4);
                          if(!i) specR[i]->Draw();
                          else specR[i]->Draw("same");
                          leg2->AddEntry(specR[i],desc[i].Data());
                      }
                      leg2->Draw("same");
                      TLine *line = new TLine(5,0.2,5,3.5);
                      line->SetLineStyle(2);
                      line->SetLineWidth(2);
                      line->Draw("same");
                      cspecRatio->SaveAs(Form("%s/%s_theoryRatio.pdf",out.Data(),outHistName.Data()));
                      cspecRatio->SaveAs(Form("%s/%s_theoryRatio.png",out.Data(),outHistName.Data()));

          /*  TLegend *leg2 = new TLegend(0.55,0.75,0.85,0.85,"p-Pb data");
            leg2->SetBorderSize(0);
            TCanvas *cspec2 = new TCanvas("cspec2","cspec2",1200,800);
            TH1F *hratio[nFiles-1];
            for(int i=0; i<nFiles-1; i++){
              //  hratio[i] = (TH1F*)specReb[i+1]->Clone( Form("hratio_%d",i));
                hratio[i] = (TH1F*)spec[i+1]->Clone(Form("hratio_%d",i));
                hratio[i]->Divide(spec[0]);
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
*/

            return;

}



TH1* GetInputHist(TString inFile, TString histName,TH1 *hh){

	TFile *jetPtFile = new TFile(inFile,"read");
    hh = (TH1F*)jetPtFile->Get(histName.Data());

    return hh;

}


TH1* GetUpSys(TH1D **hh, const int nFiles = 11, TH1D *hh_up){


        double bin = 0, binerr = 0;
        double max = 0, maxerr = 0;


      //  for(int j=1; j<fptbinsJetTrueN+1; j++ ){
        for(int j=1; j<hh[0]->GetNbinsX()+1; j++ ){
            max = hh[0]->GetBinContent(j);
            for(int i=1;i<nFiles;i++){
                if(hh[i]->GetBinContent(j) > max){
                        max = hh[i]->GetBinContent(j);
                        maxerr = hh[i]->GetBinError(j);
                }

            }
            hh_up->SetBinContent(j,max);
            hh_up->SetBinError(j,0);

        }

    return hh_up;
}

TH1* GetDownSys(TH1D **hh, const int nFiles = 11, TH1D *hh_down){


        double bin = 0, binerr = 0;
        double max = 0, maxerr = 0;

        //for(int j=1; j<fptbinsJetTrueN+1; j++ ){
        for(int j=1; j<hh[0]->GetNbinsX()+1; j++ ){
            max = hh[0]->GetBinContent(j);
            for(int i=1;i<nFiles;i++){
                if(hh[i]->GetBinContent(j) < max){
                        max = hh[i]->GetBinContent(j);
                        maxerr = hh[i]->GetBinError(j);
                }

            }
            hh_down->SetBinContent(j,max);
            hh_down->SetBinError(j,0);

        }

    return hh_down;
}
