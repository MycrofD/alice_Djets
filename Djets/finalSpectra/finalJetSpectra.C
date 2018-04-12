//
//
// Author: B.Trzeciak (barbara.antonina.trzeciak@cern.ch)
//

#include "config.h"

double jetEta = 0.9 - fRpar;
double dy = 2*jetEta;
//const int fptbinsJetTrueN = 7;
//const double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 5,6,8,10,14,20,30,50 };

double *systuncP;
//double systuncP[] = { 18,11,10,10,12,20,23 };
const double sysG = 3.8;

const int Naxis = fptbinsJetTrueN;
Double_t *xAxis;

Int_t colors[] = {1,2,8,4,6,kOrange-1,kGray+2,kCyan+1,kMagenta+2,kViolet+5,,kYellow+2};
Int_t markers[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};

TGraphAsymmErrors *grsystheory;
TGraphAsymmErrors *grsys;
TGraphAsymmErrors *grsysRatio;
TGraphAsymmErrors *grsystheoryratio;
TH1D *hData_binned;
TH1D *hData_binned_ratio;
TH1D *hPrompt_central_binned;
TH1D *hPrompt_up;
TH1D *hPrompt_down;

void ScaleHist(TH1 *hh, int full = 0);
void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, int Msize = 1.1, Width_t Lwidth = 2, Style_t Lstyle = 1);
void SaveCanvas(TCanvas *c, TString name = "tmp");

TH1D *CentralPointsStatisticalUncertainty__4;
TH1D *TH1D *GeneratorLevel_JetPtSpectrum__3 ;

bool isSimSys = 1;

void finalJetSpectra(TString dataFile = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50_ppbinning/unfolding_Bayes_3/unfoldedSpectrum_unfoldedJetSpectrum.root",
TString dataAnalysisFile = "/home/basia/Work/alice/analysis/pPb_run2/D0jet/outData/AnalysisResults_LHC16R03.root",
TString simDir = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Simulations/Prompt",
TString outSpectraDir = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50_ppbinning/finalSpectra",
TString histBase = "unfoldedSpectrum",
bool simsys = 1 )
{

  isSimSys = simsys;
  gStyle->SetOptStat(0000);

  TString outPlotDir = outSpectraDir;
  outPlotDir+="/plots";
  gSystem->Exec(Form("mkdir %s",outSpectraDir.Data()));
  gSystem->Exec(Form("mkdir %s",outPlotDir.Data()));

  xAxis = new double[Naxis+1];
  for(int k=0; k<Naxis+1; k++) xAxis[k] = fptbinsJetTrueA[k];
  systuncP = new double[fptbinsJetTrueN+1];
  for(int k=0; k<fptbinsJetTrueN+1; k++) systuncP[k] = 15;

  TFile *File = new TFile(dataAnalysisFile,"read");
  TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
  TList *histList;
  if(fDmesonSpecie) histList =  (TList*)dir->Get("histosDStarMBN0");
  else histList =  (TList*)dir->Get("histosD0MBN0");
  TH1F* hEvents = (TH1F*)histList->FindObject("hstat");
  double nEvSel = hEvents->GetBinContent(2);

  double dataLum = 1.0444*nEvSel/(sigma_in*1000) ;//Luminosity in mbar
  double simScaling = APb/2.;
  double dataScaling;
  if(fDmesonSpecie) dataScaling = 1. /(BRDstar * dataLum)/2.;
  else dataScaling = 1. /(BRDzero * dataLum)/2.;

    bool dopp = 0;

    // ----------------- prompt simulation ---------------------

    int cent = 0;
    bool jet = 1;

    int simNr = 0; // 0 - central value
    int nFiles = fCsimN;
    TH1D *hPrompt[fCsimN];
    TH1D *hPrompt_binned[fCsimN];

    for (int nr=simNr; nr<fCsimN; nr++){
        TString file = simDir;
        file += "/JetPt_";
        file += fRunC[nr];
        file += "_Dpt"; file += fptbinsDA[0]; file += "_"; file += fptbinsDA[fptbinsDN];
        if(fDmesonSpecie) file += "_Dstar";
        else file += "_Dzero";
        file += ".root";
        if(!file) { cout << "sim file doesn't exist !!! " << file << endl; return; }
      
        TH1D *htmp;
        htmp = (TH1D*) GetInputHist(file, "hPt", htmp);
        //htmp->Scale(sigma_c[nr]);
        htmp->GetYaxis()->SetTitle("d#sigma/dp_{T} (mb)");
        hPrompt[nr] = (TH1D*)htmp->Clone(Form("hPrompt_%d",nr));
        hPrompt_binned[nr] = (TH1D*)htmp->Rebin(fptbinsJetTrueN,Form("hPrompt_binned_%d",nr),fptbinsJetTrueA);

    }

    TH1D *htmp = (TH1D*)(hPrompt[cent]->Clone("htmp"));
    TH1D *hPrompt_central = (TH1D*)htmp->Clone("hPrompt_central");
    hPrompt_central_binned = (TH1D*)htmp->Rebin(fptbinsJetTrueN,"hPrompt_central_binned",fptbinsJetTrueA);

    setHistoDetails(hPrompt_central,4,24);
    setHistoDetails(hPrompt_central_binned,4,24);

    hPrompt_central->Scale(simScaling);
    hPrompt_central_binned->Scale(simScaling);
    hPrompt_central_binned->Scale(1,"width");
    hPrompt_central_binned->Scale(1./dy);


    if(isSimSys){
    // ----------------- prompt syst. (rebinned)---------------------
    // get up unc
      hPrompt_up = (TH1D*)hPrompt_central_binned->Clone("hPrompt_up");
      hPrompt_up = (TH1D*)GetUpSys(hPrompt_binned,nFiles,hPrompt_up);
      setHistoDetails(hPrompt_up,4,24,0,2,2);
      hPrompt_up->Scale(simScaling);
      hPrompt_up->Scale(1,"width");
      hPrompt_up->Scale(1./dy);
      // get down unc
      hPrompt_down = (TH1D*)hPrompt_central_binned->Clone("hPrompt_down");
      hPrompt_down = (TH1D*)GetDownSys(hPrompt_binned,nFiles,hPrompt_down);
      setHistoDetails(hPrompt_down,4,24,0,2,2);
      hPrompt_down->Scale(simScaling);
      hPrompt_down->Scale(1,"width");
      hPrompt_down->Scale(1./dy);
    }


    // ----------------- data ---------------------
    //TH1D *hData_binned;
    //hData_binned = (TH1D*)GetInputHist(dataFile, histBase, hData_binned);
    TH1D *hData_binned2;
    hData_binned2 = (TH1D*)GetInputHist(dataFile, histBase, hData_binned2);
    //TH1D *hData_binned = (TH1D*)hData_binned2->Rebin(fptbinsJetTrueN,"hData_binned", fptbinsJetTrueA);
    hData_binned = (TH1D*)hData_binned2->Clone("hData_binned");
    hData_binned->Scale(1,"width");
    hData_binned->Scale(dataScaling);
    hData_binned->Scale(1./dy);
    hData_binned->SetTitle();
    //hData_binned->SetMinimum(1);
    hData_binned->SetMaximum(hData_binned->GetMaximum()*2);
    hData_binned->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}d#eta (mb)");



    //_______________________________________________ drawing and saving
     // ------------------ draw sim and data output ---------------
    // ------------------ before B spectra folded ---------------

    if(dopp) {
         hPrompt_central_binned->SetMarkerColor(kRed+1);
         hPrompt_central_binned->SetLineColor(kRed+1);
    }

    // compare data and central sim with unc
  /*  TCanvas *cSpectra = new TCanvas("cSpectra","cSpectra",900,1000);
    //cSpectra->Divide(1,2);
    //cSpectra->cd(1);
     TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        pad1->SetBottomMargin(0);
        pad1->SetLogy();
         pad1->Draw();
         pad1->cd();
    hData_binned->SetLineColor(kRed+1);
    hData_binned->SetMarkerColor(kRed+1);
    hData_binned->GetXaxis()->SetRangeUser(plotmin,plotmax);

    hData_binned->Draw();
    hPrompt_central_binned->Draw("same p e0 x0");

    hData_binned->GetYaxis()->SetLabelSize(0.);
    TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
     axis->SetLabelFont(43);
      axis->SetLabelSize(15);
      axis->Draw();*/

    Double_t sysunc[fptbinsJetTrueN];
    Double_t sysuncAbs[fptbinsJetTrueN];
    Double_t statunc[fptbinsJetTrueN];
    Double_t value[fptbinsJetTrueN];
    Double_t ptval[fptbinsJetTrueN];
    Double_t ptvalunc[fptbinsJetTrueN];
    cout << endl;
    for(int j=0; j<fptbinsJetTrueN; j++){
            ptval[j] = (fptbinsJetTrueA[j]+fptbinsJetTrueA[j+1]) / 2.;
            ptvalunc[j] = (fptbinsJetTrueA[j+1]-fptbinsJetTrueA[j]) / 2.;
            value[j] = hData_binned->GetBinContent(hData_binned->GetXaxis()->FindBin(ptval[j]));
            error = hData_binned->GetBinError(hData_binned->GetXaxis()->FindBin(ptval[j]));
            sysunc[j] =  systuncP[j] * 0.01;
            sysuncAbs[j] = value[j] * systuncP[j] * 0.01;
            statunc[j] = error/ value[j] *100;
            cout << j << " sys: " << sysunc[j] << "\t sys2: " << sysunc[j]*100/value[j] << "\t\t stat: " << statunc[j] << endl;
    }

    grsys = new TGraphAsymmErrors(fptbinsJetTrueN,ptval,value,ptvalunc,ptvalunc,sysuncAbs,sysuncAbs);

//   grsys->SetMarkerColor(4);
//   grsys->SetMarkerStyle(21);
//   grsys->SetFillColor(kGray+3);
//   grsys->SetFillStyle(3003);
  // if(!dopp) grsys->Draw("2");

    Double_t sysuncTheory[fptbinsJetTrueN];
    Double_t ptvaltheory[fptbinsJetTrueN];
    Double_t ptvalunctheory[fptbinsJetTrueN];
    Double_t valuetheory[fptbinsJetTrueN];
    Double_t valuetheoryerrup[fptbinsJetTrueN];
    Double_t valuetheoryerrdown[fptbinsJetTrueN];
      for(int j=0; j<fptbinsJetTrueN; j++){
            ptvaltheory[j] = (fptbinsJetTrueA[j]+fptbinsJetTrueA[j+1]) / 2.;
            ptvalunctheory[j] = (fptbinsJetTrueA[j+1]-fptbinsJetTrueA[j]) / 2.;
            valuetheory[j] = hPrompt_central_binned->GetBinContent(hPrompt_central_binned->GetXaxis()->FindBin(ptvaltheory[j]));
            if(isSimSys){
              valuetheoryerrup[j] = hPrompt_up->GetBinContent(hPrompt_up->GetXaxis()->FindBin(ptvaltheory[j])) - valuetheory[j];
              valuetheoryerrdown[j] = valuetheory[j] - hPrompt_down->GetBinContent(hPrompt_up->GetXaxis()->FindBin(ptvaltheory[j]));
            }
    }

if(isSimSys){
  grsystheory = new TGraphAsymmErrors(fptbinsJetTrueN,ptvaltheory,valuetheory,ptvalunctheory,ptvalunctheory,valuetheoryerrdown,valuetheoryerrup);
}

//   grsystheory->SetMarkerColor(kBlue+2);
//   grsystheory->SetLineColor(kBlue+2);
   //grsystheory->SetMarkerStyle(21);
  // grsystheory->SetFillColor(0);
  // grsystheory->SetFillStyle(3003);
   //if(!dopp) grsystheory->Draw("2");
//   if(!dopp) grsystheory->Draw("Psame");


 //  if(isSimSys && !dopp) hPrompt_up->Draw("same");
 //   if(isSimSys && !dopp) hPrompt_down->Draw("same");

   //======= Ratios to powheg ======


   hData_binned_ratio = (TH1D*)hData_binned->Clone("hData_binned_ratio");
   TH1D *hPrompt_central_binned_ratio = (TH1D*)hPrompt_central_binned->Clone("hPrompt_central_binned_ratio");
   hPrompt_central_binned_ratio->Divide(hPrompt_central_binned);
   TH1D *hPrompt_down_ratio;
   TH1D *hPrompt_up_ratio;
   if(isSimSys){
     hPrompt_up_ratio = (TH1D*)hPrompt_up->Clone("hPrompt_up_ratio");
     hPrompt_down_ratio = (TH1D*)hPrompt_down->Clone("hPrompt_down_ratio");
     hPrompt_up_ratio->Divide(hPrompt_central_binned);
     hPrompt_down_ratio->Divide(hPrompt_central_binned);

   }

     Double_t sysuncRatio[fptbinsJetTrueN];
     Double_t valRatio[fptbinsJetTrueN];

    for(int j=0; j<fptbinsJetTrueN; j++){

            double pt = (fptbinsJetTrueA[j]+fptbinsJetTrueA[j+1]) / 2.;
            double val = hData_binned->GetBinContent(hData_binned->GetXaxis()->FindBin(pt));
            double valPred = hPrompt_central_binned->GetBinContent(hPrompt_central_binned->GetXaxis()->FindBin(pt));
            valRatio[j] = val / valPred;
            double err = hData_binned->GetBinError(hData_binned->GetXaxis()->FindBin(pt)) / val;
            err = err * valRatio[j];
            sysuncRatio[j] = sysunc[j]*valRatio[j];

            hData_binned_ratio->SetBinContent(hData_binned_ratio->GetXaxis()->FindBin(pt),valRatio[j]);
            hData_binned_ratio->SetBinError(hData_binned_ratio->GetXaxis()->FindBin(pt),err);

    }

     grsysRatio = new TGraphAsymmErrors(fptbinsJetTrueN,ptval,valRatio,ptvalunc,ptvalunc,sysuncRatio,sysuncRatio);


    Double_t sysuncTheoryratio[fptbinsJetTrueN];
    Double_t ptvaltheoryratio[fptbinsJetTrueN];
    Double_t ptvalunctheoryratio[fptbinsJetTrueN];
    Double_t valuetheoryratio[fptbinsJetTrueN];
    Double_t valuetheoryerrupratio[fptbinsJetTrueN];
    Double_t valuetheoryerrdownratio[fptbinsJetTrueN];
    if(isSimSys){
      for(int j=0; j<fptbinsJetTrueN; j++){
            ptvaltheoryratio[j] = (fptbinsJetTrueA[j]+fptbinsJetTrueA[j+1]) / 2.;
            ptvalunctheoryratio[j] = (fptbinsJetTrueA[j+1]-fptbinsJetTrueA[j]) / 2.;
            valuetheoryratio[j] = hPrompt_central_binned_ratio->GetBinContent(hPrompt_central_binned_ratio->GetXaxis()->FindBin(ptvaltheory[j]));
            valuetheoryerrupratio[j] = hPrompt_up_ratio->GetBinContent(hPrompt_up_ratio->GetXaxis()->FindBin(ptvaltheory[j])) - valuetheoryratio[j];
            valuetheoryerrdownratio[j] = valuetheoryratio[j] - hPrompt_down_ratio->GetBinContent(hPrompt_down_ratio->GetXaxis()->FindBin(ptvaltheory[j]));

          }
    }
    if(isSimSys) grsystheoryratio = new TGraphAsymmErrors(fptbinsJetTrueN,ptvaltheoryratio,valuetheoryratio,ptvalunctheoryratio,ptvalunctheoryratio,valuetheoryerrdownratio,valuetheoryerrupratio);

//   grsystheoryratio->SetMarkerColor(kBlue+2);
//   grsystheoryratio->SetLineColor(kBlue+2);
//


    //hData_binned_up->Draw("same");
    //hData_binned_down->Draw("same");
    //if(!isSimSys) hPrompt_central_binned ->Draw("same");
   // if(is&& !dopp) hPrompt_central_binned->Draw("psame");

drawFinal(outPlotDir);

TFile *ofile = new TFile(Form("%s/JetPtSpectrum_final.root",outSpectraDir.Data()),"RECREATE");
hData_binned->Write();
hData_binned_ratio->Write();
hPrompt_central_binned->Write();
if(isSimSys){
hPrompt_up->Write();
hPrompt_down->Write();

}
if(isSimSys) grsystheory->Write();
grsys->Write();
grsysRatio->Write();
if(isSimSys) grsystheoryratio->Write();

ofile->Close();

return;

    TH1D *hdatapp;
    TH1D *hR;
    if(dopp) {

        hdatapp = drawPP();

        hR = (TH1D*)CentralPointsStatisticalUncertainty__4->Clone("hR");
        for(int i=1; i<fptbinsJetTrueN; i++){
         double rr =  hData_binned->GetBinContent(i+1) / CentralPointsStatisticalUncertainty__4->GetBinContent(i);
       // cout << "==== RpPb: " << rr << endl;

        }

        TH1D *hratioData = (TH1D*)hData_binned->Clone("hratioData");
        hratioData->GetYaxis()->SetTitle("ratio");
        for(int j= 1; j<hData_binned->GetNbinsX();j++){
                double val = hData_binned->GetBinContent(j);
                double err = hData_binned->GetBinError(j);
                err = err / val;
                double valpp = CentralPointsStatisticalUncertainty__4->GetBinContent(CentralPointsStatisticalUncertainty__4->GetXaxis()->FindBin(hData_binned->GetBinCenter(j)));
                double errpp = CentralPointsStatisticalUncertainty__4->GetBinError(CentralPointsStatisticalUncertainty__4->GetXaxis()->FindBin(hData_binned->GetBinCenter(j)));
                if(valpp) errpp = errpp/valpp;
                else errpp = 0;
                //cout << "pp val: " << valpp << '\t';
               // cout << "bin cen: " << hData_binned->GetBinCenter(j) << endl;


                double ratio = 0;
                if(valpp) ratio = val / valpp;
                double errratio = ratio * (sqrt(err*err+errpp*errpp));

                hratioData->SetBinContent(j,ratio);
                hratioData->SetBinError(j,errratio);
        }

        TH1D *hratioSim = (TH1D*)hPrompt_central_binned->Clone("hratioSim");
         hratioSim->SetMarkerColor(kBlue+1);
         hratioSim->SetMarkerStyle(20);
        hratioSim->SetLineColor(kBlue+1);
        for(int j= 1; j< hPrompt_central_binned->GetNbinsX()+1;j++){
                double val = hPrompt_central_binned->GetBinContent(j);
                double err = 0;
                err = err / val;
                double valpp = GeneratorLevel_JetPtSpectrum__3->GetBinContent(GeneratorLevel_JetPtSpectrum__3->GetXaxis()->FindBin(hPrompt_central_binned->GetBinCenter(j)));
                double errpp = 0;
                if(valpp) errpp = errpp/valpp;
                else errpp = 0;
                cout << "pp val: " << valpp << '\t';
                cout << "bin cen: " << hPrompt_central_binned->GetBinCenter(j) << endl;

                double ratio = 0;
                if(valpp) ratio = val / valpp;
                double errratio = ratio * (sqrt(err*err+errpp*errpp));

                hratioSim->SetBinContent(j,ratio);
                hratioSim->SetBinError(j,errratio);
        }


    }
    TLegend *leg1 = new TLegend(0.6,0.7,0.9,0.9);
    leg1->SetBorderSize(0);


         leg1->AddEntry(hData_binned,"prompt D^{*+}-jet","p");
        //leg1->AddEntry(hData_binned_sub,"data-FD","p");
        if(isSimSys) leg1->AddEntry(hPrompt_up,"POWHEG+PYTHIA","l");
          leg1->Draw("same");




   // TCanvas *cratioPred = new TCanvas("cratioPred","cratioPred",600,300);
   // cSpectra->cd(2);
   cSpectra->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    if(isSimSys) {
      hPrompt_up_ratio->GetYaxis()->SetTitle("data/theory");
      hPrompt_up_ratio->GetYaxis()->SetRangeUser(0,2.6);
    }
     hData_binned_ratio->GetYaxis()->SetTitle("data/theory");
    hData_binned_ratio->GetYaxis()->SetRangeUser(0,2.6);
    //hPrompt_central_binned_ratio->Draw();
   // hPrompt_up_ratio->Draw("h");
  //  hPrompt_down_ratio->Draw("hsame");
    hData_binned_ratio->Draw();

    if(isSimSys){
      hPrompt_up_ratio->GetYaxis()->SetTitle("data/theory");
     hPrompt_up_ratio->GetYaxis()->SetNdivisions(505);
     hPrompt_up_ratio->GetYaxis()->SetTitleSize(20);
     hPrompt_up_ratio->GetYaxis()->SetTitleFont(43);
     hPrompt_up_ratio->GetYaxis()->SetTitleOffset(1.2);
     hPrompt_up_ratio->GetYaxis()->SetLabelFont(43);
      hPrompt_up_ratio->GetYaxis()->SetLabelSize(15);

      hPrompt_up_ratio->GetXaxis()->SetTitleSize(20);
     hPrompt_up_ratio->GetXaxis()->SetTitleFont(43);
     hPrompt_up_ratio->GetXaxis()->SetTitleOffset(4);
     hPrompt_up_ratio->GetXaxis()->SetLabelFont(43);
      hPrompt_up_ratio->GetXaxis()->SetLabelSize(15);
  }
    grsysRatio->SetMarkerColor(4);
   grsysRatio->SetMarkerStyle(21);
   grsysRatio->SetFillColor(kGray+3);
   grsysRatio->SetFillStyle(3003);
   grsysRatio->Draw("2");

  if(isSimSys)  grsystheoryratio->Draw("same");

    TLine *ll = new TLine(5,1,30,1);
    ll->SetLineStyle(3);
    ll->SetLineWidth(2);
    ll->Draw("same");

    if(dopp) {
         cSpectra->SaveAs(Form("%s/JetPtSpectra_final_wpp.png",outPlotDir.Data()));
         cratio->SaveAs(Form("%s/JetPtSpectra_ppratio.png",outPlotDir.Data()));
    }
    else cSpectra->SaveAs(Form("%s/JetPtSpectra_final.png",outPlotDir.Data()));

     if(dopp) {
         cSpectra->SaveAs(Form("%s/JetPtSpectra_final_wpp.pdf",outPlotDir.Data()));
         cratio->SaveAs(Form("%s/JetPtSpectra_ppratio.pdf",outPlotDir.Data()));
    }
    else cSpectra->SaveAs(Form("%s/JetPtSpectra_final.pdf",outPlotDir.Data()));



    return;
}

TH1* GetUpSys(TH1D **hh, const int nFiles = 11, TH1D *hh_up){


        double bin = 0, binerr = 0;
        double max = 0, maxerr = 0;


        for(int j=1; j<fptbinsJetTrueN+1; j++ ){
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

        for(int j=1; j<fptbinsJetTrueN+1; j++ ){
      //for(int j=1; j<hh[0]->GetNbinsX()+1; j++ ){
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



TH1* GetInputHist(TString inFile, TString histName,TH1 *hh){

	TFile *jetPtFile = new TFile(inFile,"read");
    hh = (TH1F*)jetPtFile->Get(histName.Data());

    return hh;

}



void ScaleHist(TH1 *hh, int full){


    if(full){
        //hh->Scale(1,"width");
        //hh->Scale(pPbscaling);
        //hh->Scale(scalingF);
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
   // hh->SetName(name.c_str());

    hh->SetTitle();
    hh->GetXaxis()->SetTitle("p_{T}^{ch,jet} (GeV/c)");

}

void SaveCanvas(TCanvas *c, TString name){

    c->SaveAs(Form("%s.png",name.Data()));
    c->SaveAs(Form("%s.pdf",name.Data()));
    //c->SaveAs(Form("%s.png",name.c_str()));
    //c->SaveAs(Form("%s.pdf",name.c_str()));

}

void drawFinal(TString outPlotDir){



  //=========Macro generated from canvas: FinalSpectrum/FinalSpectrum
//=========  (Mon Jun 12 17:00:30 2017) by ROOT version5.34/30
   TCanvas *FinalSpectrum = new TCanvas("FinalSpectrum", "FinalSpectrum",0,45,700,700);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   FinalSpectrum->SetHighLightColor(2);
   FinalSpectrum->Range(0,0,1,1);
   FinalSpectrum->SetFillColor(0);
   FinalSpectrum->SetBorderMode(0);
   FinalSpectrum->SetBorderSize(2);
   FinalSpectrum->SetFrameBorderMode(0);

// ------------>Primitives in pad: FinalSpectrum_1
   TPad *FinalSpectrum_1 = new TPad("FinalSpectrum_1", "FinalSpectrum_1",0,0.35,1,1);
   FinalSpectrum_1->Draw();
   FinalSpectrum_1->cd();
   FinalSpectrum_1->Range(-1.986821e-07,-4.69897,33.33333,0.3499945);
   FinalSpectrum_1->SetFillColor(0);
   FinalSpectrum_1->SetBorderMode(0);
   FinalSpectrum_1->SetBorderSize(2);
   FinalSpectrum_1->SetLogy();
   FinalSpectrum_1->SetTickx(1);
   FinalSpectrum_1->SetTicky(1);
   FinalSpectrum_1->SetLeftMargin(0.15);
   FinalSpectrum_1->SetBottomMargin(0);
   FinalSpectrum_1->SetFrameBorderMode(0);
   FinalSpectrum_1->SetFrameBorderMode(0);



   TH1D *CentralPointsStatisticalUncertainty__1 = new TH1D("CentralPointsStatisticalUncertainty__1","Central Values",Naxis, xAxis);
   if(fSystem){
     CentralPointsStatisticalUncertainty__1->SetMinimum(2.e-04);
     CentralPointsStatisticalUncertainty__1->SetMaximum(500);
   }
   else{
     CentralPointsStatisticalUncertainty__1->SetMinimum(2.e-05);
     CentralPointsStatisticalUncertainty__1->SetMaximum(0.5);
  }
   CentralPointsStatisticalUncertainty__1->SetEntries(8);
   CentralPointsStatisticalUncertainty__1->SetDirectory(0);
   CentralPointsStatisticalUncertainty__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__1->SetLineColor(ci);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__1->SetMarkerColor(ci);
   CentralPointsStatisticalUncertainty__1->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__1->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__1->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__1->GetXaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__1->GetXaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__1->GetXaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__1->GetXaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__1->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]");
   CentralPointsStatisticalUncertainty__1->GetYaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__1->GetYaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__1->GetYaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__1->GetYaxis()->SetTitleOffset(1.6);
   CentralPointsStatisticalUncertainty__1->GetYaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__1->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__1->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__1->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__1->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__1->Draw("axis");

   // dat syst. unc.
   TGraphAsymmErrors *grae = (TGraphAsymmErrors*) grsys->Clone("grae"); // = new TGraphAsymmErrors(6);
   grae->SetName("CentralPointsSystematicUncertainty_copy");
   grae->SetTitle("Bayes, iter=4, prior=ResponseTruth Systematics");

   ci = TColor::GetColor("#cccccc");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#cccccc");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#cccccc");
/*   grae->SetMarkerColor(ci);
   grae->SetPoint(0,5.5,0.9776732);
   grae->SetPointError(0,0.5,0.5,0.146651,0.146651);
   grae->SetPoint(1,7,0.4808065);
   grae->SetPointError(1,1,1,0.06250484,0.06250484);
   grae->SetPoint(2,9,0.2106707);
   grae->SetPointError(2,1,1,0.0294939,0.0294939);
   grae->SetPoint(3,12,0.09992121);
   grae->SetPointError(3,2,2,0.01199054,0.01199054);
   grae->SetPoint(4,17,0.02028462);
   grae->SetPointError(4,3,3,0.002839846,0.002839846);
   grae->SetPoint(5,25,0.004155956);
   grae->SetPointError(5,5,5,0.000664953,0.000664953);*/

   //=== data uncertantity from grae
   TH1F *Graph_central_syst_unc1 = new TH1F("Graph_central_syst_unc1","Bayes, iter=4, prior=ResponseTruth Systematics",100,2.5,52.5);
   Graph_central_syst_unc1->SetMinimum(4.779682e-05);
   Graph_central_syst_unc1->SetMaximum(0.02142993);
   Graph_central_syst_unc1->SetDirectory(0);
   Graph_central_syst_unc1->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_central_syst_unc1->SetLineColor(ci);
   Graph_central_syst_unc1->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   Graph_central_syst_unc1->GetXaxis()->SetLabelFont(42);
   Graph_central_syst_unc1->GetXaxis()->SetLabelSize(0.035);
   Graph_central_syst_unc1->GetXaxis()->SetTitleSize(0.035);
   Graph_central_syst_unc1->GetXaxis()->SetTitleFont(42);
   Graph_central_syst_unc1->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]");
   Graph_central_syst_unc1->GetYaxis()->SetLabelFont(42);
   Graph_central_syst_unc1->GetYaxis()->SetLabelSize(0.035);
   Graph_central_syst_unc1->GetYaxis()->SetTitleSize(0.035);
   Graph_central_syst_unc1->GetYaxis()->SetTitleFont(42);
   Graph_central_syst_unc1->GetZaxis()->SetLabelFont(42);
   Graph_central_syst_unc1->GetZaxis()->SetLabelSize(0.035);
   Graph_central_syst_unc1->GetZaxis()->SetTitleSize(0.035);
   Graph_central_syst_unc1->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_central_syst_unc1);

   grae->Draw("2");


   // Central data
   TH1D *CentralPointsStatisticalUncertainty__2 = (TH1D*) hData_binned->Clone("CentralPointsStatisticalUncertainty__2"); // = new TH1D("CentralPointsStatisticalUncertainty__2","Central Values",6, xAxis2);

   /*CentralPointsStatisticalUncertainty__2->SetBinContent(0,1.834125);
   CentralPointsStatisticalUncertainty__2->SetBinContent(1,0.9776732);
   CentralPointsStatisticalUncertainty__2->SetBinContent(2,0.4808065);
   CentralPointsStatisticalUncertainty__2->SetBinContent(3,0.2106707);
   CentralPointsStatisticalUncertainty__2->SetBinContent(4,0.09992121);
   CentralPointsStatisticalUncertainty__2->SetBinContent(5,0.02028462);
   CentralPointsStatisticalUncertainty__2->SetBinContent(6,0.004155956);
   CentralPointsStatisticalUncertainty__2->SetBinContent(7,0.0002553179);
   CentralPointsStatisticalUncertainty__2->SetBinError(0,0.06249125);
   CentralPointsStatisticalUncertainty__2->SetBinError(1,0.05513667);
   CentralPointsStatisticalUncertainty__2->SetBinError(2,0.03286703);
   CentralPointsStatisticalUncertainty__2->SetBinError(3,0.02035879);
   CentralPointsStatisticalUncertainty__2->SetBinError(4,0.01021899);
   CentralPointsStatisticalUncertainty__2->SetBinError(5,0.00376896);
   CentralPointsStatisticalUncertainty__2->SetBinError(6,0.001104688);
   CentralPointsStatisticalUncertainty__2->SetBinError(7,0.0002005154);*/
   CentralPointsStatisticalUncertainty__2->SetMinimum(2e-05);
   CentralPointsStatisticalUncertainty__2->SetMaximum(4);
   CentralPointsStatisticalUncertainty__2->SetEntries(8);
   CentralPointsStatisticalUncertainty__2->SetDirectory(0);
   CentralPointsStatisticalUncertainty__2->SetStats(0);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__2->SetLineColor(ci);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__2->SetMarkerColor(ci);
   CentralPointsStatisticalUncertainty__2->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__2->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__2->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__2->GetXaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__2->GetXaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__2->GetXaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__2->GetXaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__2->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]");
   CentralPointsStatisticalUncertainty__2->GetYaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__2->GetYaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__2->GetYaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__2->GetYaxis()->SetTitleOffset(1.4);
   CentralPointsStatisticalUncertainty__2->GetYaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__2->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__2->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__2->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__2->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__2->Draw("same p e0 x0");

   // central theory
   TH1D *GeneratorLevel_JetPtSpectrum__3 = (TH1D*)hPrompt_central_binned->Clone("GeneratorLevel_JetPtSpectrum__3"); //  new TH1D("GeneratorLevel_JetPtSpectrum__3","D0_MCTruth_Charged_R040_JetPtDPtSpectrum",6, xAxis3);
  /* GeneratorLevel_JetPtSpectrum__3->SetBinContent(0,1.586436);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(1,0.5048265);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(2,0.2791936);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(3,0.1260778);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(4,0.04664339);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(5,0.01253488);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(6,0.002788148);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(7,0.00081422);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(0,0.004003268);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(1,0.002258262);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(2,0.00118752);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(3,0.0007980088);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(4,0.0003432167);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(5,0.000145274);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(6,5.307149e-05);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(7,2.867968e-05);*/
   GeneratorLevel_JetPtSpectrum__3->SetEntries(316731);
   GeneratorLevel_JetPtSpectrum__3->SetDirectory(0);
   GeneratorLevel_JetPtSpectrum__3->SetStats(0);

   ci = TColor::GetColor("#000099");
   GeneratorLevel_JetPtSpectrum__3->SetLineColor(ci);

   ci = TColor::GetColor("#000099");
   GeneratorLevel_JetPtSpectrum__3->SetMarkerColor(ci);
   GeneratorLevel_JetPtSpectrum__3->SetMarkerStyle(24);
   GeneratorLevel_JetPtSpectrum__3->SetMarkerSize(1.2);
   GeneratorLevel_JetPtSpectrum__3->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   GeneratorLevel_JetPtSpectrum__3->GetXaxis()->SetLabelFont(42);
   GeneratorLevel_JetPtSpectrum__3->GetXaxis()->SetLabelSize(0.035);
   GeneratorLevel_JetPtSpectrum__3->GetXaxis()->SetTitleSize(0.035);
   GeneratorLevel_JetPtSpectrum__3->GetXaxis()->SetTitleFont(42);
   GeneratorLevel_JetPtSpectrum__3->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}} #times #Delta#it{p}_{T} (mb)");
   GeneratorLevel_JetPtSpectrum__3->GetYaxis()->SetLabelFont(42);
   GeneratorLevel_JetPtSpectrum__3->GetYaxis()->SetLabelSize(0.035);
   GeneratorLevel_JetPtSpectrum__3->GetYaxis()->SetTitleSize(0.035);
   GeneratorLevel_JetPtSpectrum__3->GetYaxis()->SetTitleFont(42);
   GeneratorLevel_JetPtSpectrum__3->GetZaxis()->SetLabelFont(42);
   GeneratorLevel_JetPtSpectrum__3->GetZaxis()->SetLabelSize(0.035);
   GeneratorLevel_JetPtSpectrum__3->GetZaxis()->SetTitleSize(0.035);
   GeneratorLevel_JetPtSpectrum__3->GetZaxis()->SetTitleFont(42);
   GeneratorLevel_JetPtSpectrum__3->Draw("same p e0 x0");



   // theory syst unc
   if(isSimSys){
     grae = (TGraphAsymmErrors*) grsystheory->Clone("grae"); //  new TGraphAsymmErrors(6);
     grae->SetName("theorySyst_copy");
     grae->SetTitle("Graph");
     grae->SetFillColor(1);
     grae->SetFillStyle(0);

     ci = TColor::GetColor("#000099");
     grae->SetLineColor(ci);
     grae->SetLineWidth(2);
  }
  /* grae->SetPoint(0,5.5,0.5048265);
   grae->SetPointError(0,0.5,0.5,0.3307398,0.6140599);
   grae->SetPoint(1,7,0.2791936);
   grae->SetPointError(1,1,1,0.1709185,0.3091237);
   grae->SetPoint(2,9,0.1260778);
   grae->SetPointError(2,1,1,0.07066835,0.1294046);
   grae->SetPoint(3,12,0.04664339);
   grae->SetPointError(3,2,2,0.02322777,0.04466855);
   grae->SetPoint(4,17,0.01253488);
   grae->SetPointError(4,3,3,0.005750655,0.01229934);
   grae->SetPoint(5,25,0.002788148);
   grae->SetPointError(5,5,5,0.001190245,0.002366597);*/


   TH1F *Graph_theorySyst_copy2 = new TH1F("Graph_theorySyst_copy2","Graph",100,2.5,52.5);
   Graph_theorySyst_copy2->SetMinimum(2.864129e-05);
   Graph_theorySyst_copy2->SetMaximum(0.01825038);
   Graph_theorySyst_copy2->SetDirectory(0);
   Graph_theorySyst_copy2->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_theorySyst_copy2->SetLineColor(ci);
   Graph_theorySyst_copy2->GetXaxis()->SetLabelFont(42);
   Graph_theorySyst_copy2->GetXaxis()->SetLabelSize(0.035);
   Graph_theorySyst_copy2->GetXaxis()->SetTitleSize(0.035);
   Graph_theorySyst_copy2->GetXaxis()->SetTitleFont(42);
   Graph_theorySyst_copy2->GetYaxis()->SetLabelFont(42);
   Graph_theorySyst_copy2->GetYaxis()->SetLabelSize(0.035);
   Graph_theorySyst_copy2->GetYaxis()->SetTitleSize(0.035);
   Graph_theorySyst_copy2->GetYaxis()->SetTitleFont(42);
   Graph_theorySyst_copy2->GetZaxis()->SetLabelFont(42);
   Graph_theorySyst_copy2->GetZaxis()->SetLabelSize(0.035);
   Graph_theorySyst_copy2->GetZaxis()->SetTitleSize(0.035);
   Graph_theorySyst_copy2->GetZaxis()->SetTitleFont(42);
   if(isSimSys){
     grae->SetHistogram(Graph_theorySyst_copy2);
     grae->Draw("2");
  }

   TLegend *leg = new TLegend(0.5,0.35,0.8,0.61,NULL,"NB NDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(43);
   leg->SetTextSize(23);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("CentralPointsStatisticalUncertainty","Data","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#990000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(0.9);
   entry->SetTextFont(43);
   entry=leg->AddEntry("CentralPointsSystematicUncertainty_copy","Syst. Unc. (data)","f");

   ci = TColor::GetColor("#cccccc");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#cccccc");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(43);
   entry=leg->AddEntry("GeneratorLevel_JetPtSpectrum","POWHEG+PYTHIA6","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#000099");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(43);
   entry=leg->AddEntry("theorySyst_copy","Syst. Unc. (theory)","f");
   entry->SetFillColor(1);

   ci = TColor::GetColor("#000099");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(43);
   leg->Draw();

   TPaveText *pt = new TPaveText(0.16,0.64,0.55,0.9,"NB NDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(13);
   pt->SetTextFont(43);
   pt->SetTextSize(22);
   TText *text = pt->AddText("ALICE Preliminary");
   if(fSystem) text = pt->AddText("p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
   else text = pt->AddText("pp, #sqrt{#it{s}} = 5.02 TeV");
   text = pt->AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.3, |#eta_{jet}| < 0.6");
   text = pt->AddText("with D^{0}, #it{p}_{T,D} > 3 GeV/#it{c}");
   pt->Draw();

   // does nothing
   TH1D *CentralPointsStatisticalUncertainty__4 = new TH1D("CentralPointsStatisticalUncertainty__4","Central Values",Naxis, xAxis);
   CentralPointsStatisticalUncertainty__4->SetMinimum(2e-05);
   CentralPointsStatisticalUncertainty__4->SetMaximum(0.7);
   CentralPointsStatisticalUncertainty__4->SetEntries(8);
   CentralPointsStatisticalUncertainty__4->SetDirectory(0);
   CentralPointsStatisticalUncertainty__4->SetStats(0);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__4->SetLineColor(ci);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__4->SetMarkerColor(ci);
   CentralPointsStatisticalUncertainty__4->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__4->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__4->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__4->GetXaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__4->GetXaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__4->GetXaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__4->GetXaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__4->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]");
   CentralPointsStatisticalUncertainty__4->GetYaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__4->GetYaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__4->GetYaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__4->GetYaxis()->SetTitleOffset(1.6);
   CentralPointsStatisticalUncertainty__4->GetYaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__4->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__4->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__4->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__4->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__4->Draw("sameaxis");
   FinalSpectrum_1->Modified();
   FinalSpectrum->cd();


// ------------>Primitives in pad: FinalSpectrum_2
   FinalSpectrum_2 = new TPad("FinalSpectrum_2", "FinalSpectrum_2",0,0,1,0.35);
   FinalSpectrum_2->Draw();
   FinalSpectrum_2->cd();
   FinalSpectrum_2->Range(-1.986821e-07,-0.9209589,33.33333,2.49);
   FinalSpectrum_2->SetFillColor(0);
   FinalSpectrum_2->SetBorderMode(0);
   FinalSpectrum_2->SetBorderSize(2);
   FinalSpectrum_2->SetGridy();
   FinalSpectrum_2->SetTickx(1);
   FinalSpectrum_2->SetTicky(1);
   FinalSpectrum_2->SetLeftMargin(0.15);
   FinalSpectrum_2->SetTopMargin(0);
   FinalSpectrum_2->SetBottomMargin(0.27);
   FinalSpectrum_2->SetFrameBorderMode(0);
   FinalSpectrum_2->SetFrameBorderMode(0);

   // central points data (values don't really needed)
   TH1D *CentralPointsStatisticalUncertainty__5 = new TH1D("CentralPointsStatisticalUncertainty__5","Central Values",Naxis, xAxis);
   CentralPointsStatisticalUncertainty__5->SetMinimum(0);
   CentralPointsStatisticalUncertainty__5->SetMaximum(2.49);
   CentralPointsStatisticalUncertainty__5->SetEntries(8);
   CentralPointsStatisticalUncertainty__5->SetDirectory(0);
   CentralPointsStatisticalUncertainty__5->SetStats(0);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__5->SetLineColor(ci);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__5->SetMarkerColor(ci);
   CentralPointsStatisticalUncertainty__5->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__5->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__5->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__5->GetXaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__5->GetXaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__5->GetXaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__5->GetXaxis()->SetTitleOffset(2.9);
   CentralPointsStatisticalUncertainty__5->GetXaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__5->GetYaxis()->SetTitle("data / theory");
   CentralPointsStatisticalUncertainty__5->GetYaxis()->SetNdivisions(509);
   CentralPointsStatisticalUncertainty__5->GetYaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__5->GetYaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__5->GetYaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__5->GetYaxis()->SetTitleOffset(1.4);
   CentralPointsStatisticalUncertainty__5->GetYaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__5->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__5->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__5->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__5->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__5->Draw("axis");



   // data syst unc ratio
   grae =  (TGraphAsymmErrors*) grsysRatio->Clone("grae"); //   new TGraphAsymmErrors(6);
   grae->SetName("ratioSyst");
   grae->SetTitle("Bayes, iter=4, prior=ResponseTruth Systematics");

   ci = TColor::GetColor("#cccccc");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#cccccc");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#cccccc");
   grae->SetMarkerColor(ci);
   /*grae->SetPoint(0,5.5,1.936652);
   grae->SetPointError(0,0.5,0.5,0.2904978,0.2904978);
   grae->SetPoint(1,7,1.722126);
   grae->SetPointError(1,1,1,0.2238763,0.2238763);
   grae->SetPoint(2,9,1.670958);
   grae->SetPointError(2,1,1,0.2339341,0.2339341);
   grae->SetPoint(3,12,2.142237);
   grae->SetPointError(3,2,2,0.2570685,0.2570685);
   grae->SetPoint(4,17,1.618254);
   grae->SetPointError(4,3,3,0.2265555,0.2265555);
   grae->SetPoint(5,25,1.49058);
   grae->SetPointError(5,5,5,0.2384927,0.2384927);*/

   TH1F *Graph_ratioSyst3 = new TH1F("Graph_ratioSyst3","Bayes, iter=4, prior=ResponseTruth Systematics",100,2.5,52.5);
   Graph_ratioSyst3->SetMinimum(0.5839274);
   Graph_ratioSyst3->SetMaximum(2.28012);
   Graph_ratioSyst3->SetDirectory(0);
   Graph_ratioSyst3->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_ratioSyst3->SetLineColor(ci);
   Graph_ratioSyst3->GetXaxis()->SetLabelFont(42);
   Graph_ratioSyst3->GetXaxis()->SetLabelSize(0.035);
   Graph_ratioSyst3->GetXaxis()->SetTitleSize(0.035);
   Graph_ratioSyst3->GetXaxis()->SetTitleFont(42);
   Graph_ratioSyst3->GetYaxis()->SetLabelFont(42);
   Graph_ratioSyst3->GetYaxis()->SetLabelSize(0.035);
   Graph_ratioSyst3->GetYaxis()->SetTitleSize(0.035);
   Graph_ratioSyst3->GetYaxis()->SetTitleFont(42);
   Graph_ratioSyst3->GetZaxis()->SetLabelFont(42);
   Graph_ratioSyst3->GetZaxis()->SetLabelSize(0.035);
   Graph_ratioSyst3->GetZaxis()->SetTitleSize(0.035);
   Graph_ratioSyst3->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_ratioSyst3);

   grae->Draw("2");

   // data central ratio
   TH1D *CentralPointsStatisticalUncertainty__6 =  (TH1D*) hData_binned_ratio->Clone("CentralPointsStatisticalUncertainty__6"); // new TH1D("CentralPointsStatisticalUncertainty__6","Central Values",6, xAxis6);
   /*CentralPointsStatisticalUncertainty__6->SetBinContent(0,1.83412);
   CentralPointsStatisticalUncertainty__6->SetBinContent(1,1.936652);
   CentralPointsStatisticalUncertainty__6->SetBinContent(2,1.722126);
   CentralPointsStatisticalUncertainty__6->SetBinContent(3,1.670958);
   CentralPointsStatisticalUncertainty__6->SetBinContent(4,2.142237);
   CentralPointsStatisticalUncertainty__6->SetBinContent(5,1.618254);
   CentralPointsStatisticalUncertainty__6->SetBinContent(6,1.49058);
   CentralPointsStatisticalUncertainty__6->SetBinContent(7,0.0002553179);
   CentralPointsStatisticalUncertainty__6->SetBinError(0,0.06249125);
   CentralPointsStatisticalUncertainty__6->SetBinError(1,0.109219);
   CentralPointsStatisticalUncertainty__6->SetBinError(2,0.1177213);
   CentralPointsStatisticalUncertainty__6->SetBinError(3,0.161478);
   CentralPointsStatisticalUncertainty__6->SetBinError(4,0.2190877);
   CentralPointsStatisticalUncertainty__6->SetBinError(5,0.3006778);
   CentralPointsStatisticalUncertainty__6->SetBinError(6,0.3962086);
   CentralPointsStatisticalUncertainty__6->SetBinError(7,0.0002005154);*/
   CentralPointsStatisticalUncertainty__6->SetMinimum(1.122659e-05);
   CentralPointsStatisticalUncertainty__6->SetMaximum(0.06069752);
   CentralPointsStatisticalUncertainty__6->SetEntries(14);
   CentralPointsStatisticalUncertainty__6->SetDirectory(0);
   CentralPointsStatisticalUncertainty__6->SetStats(0);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__6->SetLineColor(ci);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__6->SetMarkerColor(ci);
   CentralPointsStatisticalUncertainty__6->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__6->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__6->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__6->GetXaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__6->GetXaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__6->GetXaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__6->GetXaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__6->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]");
   CentralPointsStatisticalUncertainty__6->GetYaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__6->GetYaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__6->GetYaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__6->GetYaxis()->SetTitleOffset(1.4);
   CentralPointsStatisticalUncertainty__6->GetYaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__6->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__6->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__6->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__6->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__6->Draw("same p e0 x0");

   // theory syst ratio
   if(isSimSys) {
     grae = (TGraphAsymmErrors*)  grsystheoryratio->Clone("grae"); // new TGraphAsymmErrors(6);
     grae->SetName("ratioTheorySyst");
     grae->SetTitle("Graph");
     grae->SetFillColor(1);
     grae->SetFillStyle(0);

     ci = TColor::GetColor("#000099");
     grae->SetLineColor(ci);
     grae->SetLineWidth(2);
  }
/*     grae->SetPoint(0,5.5,1);
   grae->SetPointError(0,0.5,0.5,0.6551554,1.216378);
   grae->SetPoint(1,7,1);
   grae->SetPointError(1,1,1,0.6121863,1.107202);
   grae->SetPoint(2,9,1);
   grae->SetPointError(2,1,1,0.5605137,1.026387);
   grae->SetPoint(3,12,1);
   grae->SetPointError(3,2,2,0.4979863,0.957661);
   grae->SetPoint(4,17,1);
   grae->SetPointError(4,3,3,0.4587722,0.9812092);
   grae->SetPoint(5,25,1);
   grae->SetPointError(5,5,5,0.4268945,0.848806);*/

   TH1F *Graph_ratioTheorySyst4 = new TH1F("Graph_ratioTheorySyst4","Graph",100,2.5,32.5);
   Graph_ratioTheorySyst4->SetMinimum(0.2058133);
   Graph_ratioTheorySyst4->SetMaximum(1.968171);
   Graph_ratioTheorySyst4->SetDirectory(0);
   Graph_ratioTheorySyst4->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_ratioTheorySyst4->SetLineColor(ci);
   Graph_ratioTheorySyst4->GetXaxis()->SetLabelFont(42);
   Graph_ratioTheorySyst4->GetXaxis()->SetLabelSize(0.035);
   Graph_ratioTheorySyst4->GetXaxis()->SetTitleSize(0.035);
   Graph_ratioTheorySyst4->GetXaxis()->SetTitleFont(42);
   Graph_ratioTheorySyst4->GetYaxis()->SetLabelFont(42);
   Graph_ratioTheorySyst4->GetYaxis()->SetLabelSize(0.035);
   Graph_ratioTheorySyst4->GetYaxis()->SetTitleSize(0.035);
   Graph_ratioTheorySyst4->GetYaxis()->SetTitleFont(42);
   Graph_ratioTheorySyst4->GetZaxis()->SetLabelFont(42);
   Graph_ratioTheorySyst4->GetZaxis()->SetLabelSize(0.035);
   Graph_ratioTheorySyst4->GetZaxis()->SetTitleSize(0.035);
   Graph_ratioTheorySyst4->GetZaxis()->SetTitleFont(42);
   if(isSimSys) {
     grae->SetHistogram(Graph_ratioTheorySyst4);
     grae->Draw("2");
  }
   // just for graphics, data not really needed
   TH1D *CentralPointsStatisticalUncertainty__7 = new TH1D("CentralPointsStatisticalUncertainty__7","Central Values",Naxis, xAxis);
   CentralPointsStatisticalUncertainty__7->SetMinimum(0);
   CentralPointsStatisticalUncertainty__7->SetMaximum(2.49);
   CentralPointsStatisticalUncertainty__7->SetEntries(8);
   CentralPointsStatisticalUncertainty__7->SetDirectory(0);
   CentralPointsStatisticalUncertainty__7->SetStats(0);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__7->SetLineColor(ci);

   ci = TColor::GetColor("#990000");
   CentralPointsStatisticalUncertainty__7->SetMarkerColor(ci);
   CentralPointsStatisticalUncertainty__7->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__7->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__7->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__7->GetXaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__7->GetXaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__7->GetXaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__7->GetXaxis()->SetTitleOffset(2.9);
   CentralPointsStatisticalUncertainty__7->GetXaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__7->GetYaxis()->SetTitle("data / theory");
   CentralPointsStatisticalUncertainty__7->GetYaxis()->SetNdivisions(509);
   CentralPointsStatisticalUncertainty__7->GetYaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__7->GetYaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__7->GetYaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__7->GetYaxis()->SetTitleOffset(1.4);
   CentralPointsStatisticalUncertainty__7->GetYaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__7->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__7->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__7->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__7->GetZaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__7->Draw("sameaxig");

   FinalSpectrum_2->Modified();
   FinalSpectrum->cd();
   FinalSpectrum->Modified();
   FinalSpectrum->cd();
   FinalSpectrum->SetSelected(FinalSpectrum);

   FinalSpectrum->SaveAs(Form("%s/JetPtSpectra_final.png",outPlotDir.Data()));
   FinalSpectrum->SaveAs(Form("%s/JetPtSpectra_final.pdf",outPlotDir.Data()));
   FinalSpectrum->SaveAs(Form("%s/JetPtSpectra_final.eps",outPlotDir.Data()));
}


TH1* drawPP(){




 Double_t xAxis4[7] = {5, 6, 8, 10, 14, 20, 30};

   CentralPointsStatisticalUncertainty__4 = new TH1D("CentralPointsStatisticalUncertainty__4","Central Values",6, xAxis4);
   CentralPointsStatisticalUncertainty__4->SetBinContent(0,0.04368426);
   CentralPointsStatisticalUncertainty__4->SetBinContent(1,0.0176198);
   CentralPointsStatisticalUncertainty__4->SetBinContent(2,0.009234982);
   CentralPointsStatisticalUncertainty__4->SetBinContent(3,0.004266254);
   CentralPointsStatisticalUncertainty__4->SetBinContent(4,0.0007877781);
   CentralPointsStatisticalUncertainty__4->SetBinContent(5,0.0002866535);
   CentralPointsStatisticalUncertainty__4->SetBinContent(6,6.840787e-05);
   CentralPointsStatisticalUncertainty__4->SetBinContent(7,1.640951e-05);
   CentralPointsStatisticalUncertainty__4->SetBinError(0,0.003617416);
   CentralPointsStatisticalUncertainty__4->SetBinError(1,0.001885998);
   CentralPointsStatisticalUncertainty__4->SetBinError(2,0.0009230143);
   CentralPointsStatisticalUncertainty__4->SetBinError(3,0.0006312643);
   CentralPointsStatisticalUncertainty__4->SetBinError(4,0.0002221694);
   CentralPointsStatisticalUncertainty__4->SetBinError(5,8.947761e-05);
   CentralPointsStatisticalUncertainty__4->SetBinError(6,3.312407e-05);
   CentralPointsStatisticalUncertainty__4->SetBinError(7,8.561932e-06);
   CentralPointsStatisticalUncertainty__4->SetMinimum(2e-05);
   CentralPointsStatisticalUncertainty__4->SetMaximum(0.7);
   CentralPointsStatisticalUncertainty__4->SetEntries(8);
   CentralPointsStatisticalUncertainty__4->SetDirectory(0);
   CentralPointsStatisticalUncertainty__4->SetStats(0);

   //ci = TColor::GetColor("#990000");
   //CentralPointsStatisticalUncertainty__4->SetLineColor(ci);
   CentralPointsStatisticalUncertainty__4->SetLineColor(kBlue+1);
   CentralPointsStatisticalUncertainty__4->SetMarkerColor(kBlue+1);
   //CentralPointsStatisticalUncertainty__4->SetMarkerColor(ci);

   CentralPointsStatisticalUncertainty__4->SetMarkerStyle(20);
   CentralPointsStatisticalUncertainty__4->SetMarkerSize(0.9);
   CentralPointsStatisticalUncertainty__4->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   CentralPointsStatisticalUncertainty__4->GetXaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__4->GetXaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__4->GetXaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__4->GetXaxis()->SetTitleFont(42);
   CentralPointsStatisticalUncertainty__4->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]");
   CentralPointsStatisticalUncertainty__4->GetYaxis()->SetLabelFont(43);
   CentralPointsStatisticalUncertainty__4->GetYaxis()->SetLabelSize(22);
   CentralPointsStatisticalUncertainty__4->GetYaxis()->SetTitleSize(26);
   CentralPointsStatisticalUncertainty__4->GetYaxis()->SetTitleOffset(1.6);
   CentralPointsStatisticalUncertainty__4->GetYaxis()->SetTitleFont(43);
   CentralPointsStatisticalUncertainty__4->GetZaxis()->SetLabelFont(42);
   CentralPointsStatisticalUncertainty__4->GetZaxis()->SetLabelSize(0.035);
   CentralPointsStatisticalUncertainty__4->GetZaxis()->SetTitleSize(0.035);
   CentralPointsStatisticalUncertainty__4->GetZaxis()->SetTitleFont(42);

  // CentralPointsStatisticalUncertainty__4->Scale(208/2.1);

   CentralPointsStatisticalUncertainty__4->Draw("same");

      Double_t xAxis3[7] = {5, 6, 8, 10, 14, 20, 30};

   GeneratorLevel_JetPtSpectrum__3 = new TH1D("GeneratorLevel_JetPtSpectrum__3","D0_MCTruth_Charged_R040_JetPtDPtSpectrum",6, xAxis3);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(0,0.02642345);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(1,0.009111114);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(2,0.005175957);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(3,0.002394825);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(4,0.0009202893);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(5,0.0002493489);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(6,5.509451e-05);
   GeneratorLevel_JetPtSpectrum__3->SetBinContent(7,1.998239e-05);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(0,4.240729e-05);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(1,2.490184e-05);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(2,1.327169e-05);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(3,9.027502e-06);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(4,3.957108e-06);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(5,1.681798e-06);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(6,6.123502e-07);
   GeneratorLevel_JetPtSpectrum__3->SetBinError(7,3.687818e-07);
   GeneratorLevel_JetPtSpectrum__3->SetEntries(440507);
   GeneratorLevel_JetPtSpectrum__3->SetDirectory(0);
   GeneratorLevel_JetPtSpectrum__3->SetStats(0);

   int ci = TColor::GetColor("#000099");
   GeneratorLevel_JetPtSpectrum__3->SetLineColor(ci);

   ci = TColor::GetColor("#000099");
   GeneratorLevel_JetPtSpectrum__3->SetMarkerColor(ci);
   GeneratorLevel_JetPtSpectrum__3->SetMarkerStyle(24);
   GeneratorLevel_JetPtSpectrum__3->SetMarkerSize(1.2);
   GeneratorLevel_JetPtSpectrum__3->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   GeneratorLevel_JetPtSpectrum__3->GetXaxis()->SetLabelFont(42);
   GeneratorLevel_JetPtSpectrum__3->GetXaxis()->SetLabelSize(0.035);
   GeneratorLevel_JetPtSpectrum__3->GetXaxis()->SetTitleSize(0.035);
   GeneratorLevel_JetPtSpectrum__3->GetXaxis()->SetTitleFont(42);
   GeneratorLevel_JetPtSpectrum__3->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}} #times #Delta#it{p}_{T} (mb)");
   GeneratorLevel_JetPtSpectrum__3->GetYaxis()->SetLabelFont(42);
   GeneratorLevel_JetPtSpectrum__3->GetYaxis()->SetLabelSize(0.035);
   GeneratorLevel_JetPtSpectrum__3->GetYaxis()->SetTitleSize(0.035);
   GeneratorLevel_JetPtSpectrum__3->GetYaxis()->SetTitleFont(42);
   GeneratorLevel_JetPtSpectrum__3->GetZaxis()->SetLabelFont(42);
   GeneratorLevel_JetPtSpectrum__3->GetZaxis()->SetLabelSize(0.035);
   GeneratorLevel_JetPtSpectrum__3->GetZaxis()->SetTitleSize(0.035);
   GeneratorLevel_JetPtSpectrum__3->GetZaxis()->SetTitleFont(42);
   GeneratorLevel_JetPtSpectrum__3->Draw("same p e0 x0");

   grae = new TGraphAsymmErrors(6);
   grae->SetName("theorySyst_copy");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetFillStyle(0);

   ci = TColor::GetColor("#000099");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);
   grae->SetLineStyle(2);
   grae->SetPoint(0,5.5,0.009111114);
   grae->SetPointError(0,0.5,0.5,0.005897839,0.007483031);
   grae->SetPoint(1,7,0.005175957);
   grae->SetPointError(1,1,1,0.003157424,0.003836946);
   grae->SetPoint(2,9,0.002394825);
   grae->SetPointError(2,1,1,0.001348567,0.001666138);
   grae->SetPoint(3,12,0.0009202893);
   grae->SetPointError(3,2,2,0.0004774547,0.0005949356);
   grae->SetPoint(4,17,0.0002493489);
   grae->SetPointError(4,3,3,0.0001182088,0.000151889);
   grae->SetPoint(5,25,5.509451e-05);
   grae->SetPointError(5,5,5,2.327085e-05,3.399796e-05);

   TH1F *Graph_theorySyst_copy2 = new TH1F("Graph_theorySyst_copy2","Graph",100,2.5,32.5);
   Graph_theorySyst_copy2->SetMinimum(2.864129e-05);
   Graph_theorySyst_copy2->SetMaximum(0.01825038);
   Graph_theorySyst_copy2->SetDirectory(0);
   Graph_theorySyst_copy2->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_theorySyst_copy2->SetLineColor(ci);
   Graph_theorySyst_copy2->SetLineStyle(2);
   Graph_theorySyst_copy2->GetXaxis()->SetLabelFont(42);
   Graph_theorySyst_copy2->GetXaxis()->SetLabelSize(0.035);
   Graph_theorySyst_copy2->GetXaxis()->SetTitleSize(0.035);
   Graph_theorySyst_copy2->GetXaxis()->SetTitleFont(42);
   Graph_theorySyst_copy2->GetYaxis()->SetLabelFont(42);
   Graph_theorySyst_copy2->GetYaxis()->SetLabelSize(0.035);
   Graph_theorySyst_copy2->GetYaxis()->SetTitleSize(0.035);
   Graph_theorySyst_copy2->GetYaxis()->SetTitleFont(42);
   Graph_theorySyst_copy2->GetZaxis()->SetLabelFont(42);
   Graph_theorySyst_copy2->GetZaxis()->SetLabelSize(0.035);
   Graph_theorySyst_copy2->GetZaxis()->SetTitleSize(0.035);
   Graph_theorySyst_copy2->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_theorySyst_copy2);

   grae->Draw("2");

   return;


   }
