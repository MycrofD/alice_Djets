//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include "TObject.h"
#include "TMath.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TDatabasePDG.h"
#include "AliNormalizationCounter.h"
#include "TStyle.h"
#include "AliHFInvMassFitter.h"
#include "TProfile.h"
#include "AliAnalysisTaskDmesonJets.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TLegendEntry.h"
#include "TPad.h"
#include "../Preliminaryplots/FD/style.C"

static Double_t dy;
static Bool_t pdf;
static Int_t *jetpTbins;
static Int_t *DpTbins[2];
static Int_t zBin;
static Double_t plotRanges[4]{0,2,10e-9,1};


void ScaleHist(TH1 *hh, int full = 0);
void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, Size_t Msize = 1.1f, Width_t Lwidth = 2, Style_t Lstyle = 1);
void SaveCanvas(TCanvas *c, TString name = "tmp");
void getSystematics(TString inDir, TString outPlotDir);
TH1D* GetUpSys(TH1D *central, TH1D **hh, Int_t nFiles = 0);
TH1D* GetDownSys(TH1D *central, TH1D **hh, Int_t nFiles = 0);
TH1* GetInputHist(TString inFile, TString histName);
TGraphAsymmErrors* GetInputGraph(TString inFile, TString histName);
void ScaleHist(TH1 *hh, int full);
void drawFinal(TString outPlotDir);
//TH1D* GetRatio(TH1D *h1, TH1D *h2);
//TH1D* GetRatio(TH1D *h1, TH1D *h2, Double_t shift);


static TString fFiles[12] = {
    "/media/jackbauer/data/z_out/R_02_finaltry/unfolding/final/Jetpt5_7JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root"
    ,"/media/jackbauer/data/z_out/R_02_finaltry/unfolding/final/Jetpt7_10JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root"
    ,"/media/jackbauer/data/z_out/R_02_finaltry/unfolding/final/Jetpt10_15JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root"
    ,"/media/jackbauer/data/z_out/R_02_finaltry/unfolding/final/Jetpt15_50JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root"
    ,"/media/jackbauer/data/z_out/R_04_finaltry/unfolding/final/Jetpt5_7JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root"
    ,"/media/jackbauer/data/z_out/R_04_finaltry/unfolding/final/Jetpt7_10JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root"
    ,"/media/jackbauer/data/z_out/R_04_finaltry/unfolding/final/Jetpt10_15JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root"
    ,"/media/jackbauer/data/z_out/R_04_finaltry/unfolding/final/Jetpt15_50JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root"
    ,"/media/jackbauer/data/z_out/R_06_finaltry/unfolding/final/Jetpt5_7JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root"
    ,"/media/jackbauer/data/z_out/R_06_finaltry/unfolding/final/Jetpt7_10JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root"
    ,"/media/jackbauer/data/z_out/R_06_finaltry/unfolding/final/Jetpt10_15JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root"
    ,"/media/jackbauer/data/z_out/R_06_finaltry/unfolding/final/Jetpt15_50JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root"
    };

/*
static TString fFiles[12] = {
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_z02F/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW2/JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_z02F/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW3/JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_z02F/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW4/JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_z02F/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW5/JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_z04F/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW2/JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_z04F/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW3/JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_z04F/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW4/JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_z04F/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW5/JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_z06F/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW2/JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_z06F/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW3/JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_z06F/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW4/JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_z06F/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW5/JetPtSpectrum_final_PDF_noneGlobal_addedCUTandJES.root",
    };
*/

static TString fFiles2[9] = {
    "/home/kvapil/Desktop/auro/5TeV_forJakub_15Dec_v0/JetPtSpectrum_final_fullGlobal_addedCUTandJES_5TeVR02.root",
    "/home/kvapil/Desktop/auro/5TeV_forJakub_15Dec_v0/JetPtSpectrum_final_fullGlobal_addedCUTandJES_5TeVR04.root",
    "/home/kvapil/Desktop/auro/5TeV_forJakub_15Dec_v0/JetPtSpectrum_final_fullGlobal_addedCUTandJES_5TeVR06.root",
    };

static TString fFiles13TeVforR[9] = {
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_X02/Default_AnalysisResults_Run2w18b.root/unfolding_Bayes_5/finalSpectraNEW/JetPtSpectrum_final_noneGlobal_separateCUTnoneJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_X04/Default_AnalysisResults_Run2w18b.root/unfolding_Bayes_5/finalSpectraNEW/JetPtSpectrum_final_noneGlobal_separateCUTnoneJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_X06/Default_AnalysisResults_Run2w18b.root/unfolding_Bayes_5/finalSpectraNEW/JetPtSpectrum_final_noneGlobal_separateCUTnoneJES.root",
    };

static TString fFiles5TeVforR[9] = {
    "/home/kvapil/Desktop/auro/5TeV_forJakub_15Dec_sysGlobal_1/JetPtSpectrum_final_noneGlobal_separateCUTnoneJES_sysGlo1_R02_5TeV.root",
    "/home/kvapil/Desktop/auro/5TeV_forJakub_15Dec_sysGlobal_1/JetPtSpectrum_final_noneGlobal_separateCUTnoneJES_sysGlo1_R04_5TeV.root",
    "/home/kvapil/Desktop/auro/5TeV_forJakub_15Dec_sysGlobal_1/JetPtSpectrum_final_noneGlobal_separateCUTnoneJES_sysGlo1_R06_5TeV.root"
};

static TString fFiles13TeVforE[9] = {
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_X02/Default_AnalysisResults_Run2w18b.root/unfolding_Bayes_5/finalSpectraNEW/JetPtSpectrum_final_noBRUnc_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_X04/Default_AnalysisResults_Run2w18b.root/unfolding_Bayes_5/finalSpectraNEW/JetPtSpectrum_final_noBRUnc_addedCUTandJES.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_X06/Default_AnalysisResults_Run2w18b.root/unfolding_Bayes_5/finalSpectraNEW/JetPtSpectrum_final_noBRUnc_addedCUTandJES.root",
    };

static TString fFiles5TeVforE[9] = {
    "/home/kvapil/Desktop/auro/5_TeV_forJakub_15Dec_sysGlobal_2/JetPtSpectrum_final_noBRUnc_addedCUTandJES_R02.root",
    "/home/kvapil/Desktop/auro/5_TeV_forJakub_15Dec_sysGlobal_2/JetPtSpectrum_final_noBRUnc_addedCUTandJES_R04.root",
    "/home/kvapil/Desktop/auro/5_TeV_forJakub_15Dec_sysGlobal_2/JetPtSpectrum_final_noBRUnc_addedCUTandJES_R06.root"
};



/*
static TString fFiles[9] = {
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_Z02/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW25/JetPtSpectrum_final.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_04zF/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW25/JetPtSpectrum_final.root",
    "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_Z06/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW25/JetPtSpectrum_final.root",
    };
*/

std::tuple<TGraphAsymmErrors*, TGraphAsymmErrors*, TH1D**> GetSim(TString hash, TString dataFile,TString simname, Int_t nVar);
std::tuple<TGraphAsymmErrors*, TH1D*, TH1D*, TH1D*> GetDataSimRatio(TString simname, TH1D *data_cent, TH1D *sim_cent, TH1D *sim_up, TH1D *sim_down, UInt_t nBins, Double_t *xBins);
std::tuple<TGraphAsymmErrors*, TH1D*,TGraphAsymmErrors*, TH1D*> GetData(TString dataFile, TString histBase, Double_t dataScaling, UInt_t nBins, Double_t *xBins, Double_t *systUncD_down, Double_t *systUncD_up);
std::tuple<TCanvas*, TPad*, TPad*,  TPad*, TH1D*, TH1D*, TH1D*> PrepareCanvas(UInt_t xAxisBins, Double_t *xAxis);
void PlaceOnPadData(TPad* pad,TGraphAsymmErrors* histo1, TH1D* histo2, Size_t markersize,Style_t markerstyle,Color_t histo2col, Double_t scale, TString mode ="*");
void PlaceOnPad(TPad* pad, TGraphAsymmErrors* histo, Color_t ci, Style_t style, Style_t linestyle, TString opt,Double_t scale, TString mode ="*");
void PlaceOnPad(TPad* pad, TH1D* histo, Color_t ci, Style_t style, Style_t linestyle, TString opt,Double_t scale, TString mode ="*");
void TerminateCanvas(TPad* pad1,TPad* pad2,TPad* pad3,TH1D* histo1,TH1D* histo2,TH1D* histo3);
std::tuple<TCanvas*, TPad**, TH1D**> PrepareCanvas8Pad(UInt_t xAxisBins, Double_t *xAxisl, Double_t *xAxisr);
std::tuple<TCanvas*, TPad**, TH1D**> PrepareCanvas6Pad(UInt_t xAxisBins, Double_t *xAxisl, Double_t *xAxisr);
std::tuple<TCanvas*, TPad**, TH1D**> PrepareCanvas4Pad(UInt_t xAxisBins, Double_t *xAxis);
std::tuple<TCanvas*, TPad**, TH1D**> PrepareCanvas3x8Pad(UInt_t xAxisBins, Double_t *xAxisl, Double_t *xAxisr);
std::tuple<TCanvas*, TPad**, TH1D**> PrepareCanvas4x6Pad(UInt_t xAxisBins, Double_t *xAxisl, Double_t *xAxisr);
void TerminateCanvas3x2Pad(TPad **pad, TH1D** histo);

//std::tuple<TH1D*, TH1D*, TH1D*> GetRatio(TH1D *sim_cent, TH1D *sim_up, TH1D *sim_down,TH1D *sim2_cent, TH1D *sim2_up, TH1D *sim2_down);
std::tuple<TH1D*, TGraphAsymmErrors*, TH1D*,  TGraphAsymmErrors*> GetData(TString dataFile, TString histBase1, TString histBase2, TString histBase3, TString histBase4);
std::tuple<TH1D*, TGraphAsymmErrors*, TGraph*, TH1D*,  TGraphAsymmErrors*> GetData(TString dataFile, TString histBase1, TString histBase2, TString histBase3, TString histBase4, TString histBase5);
TH1D** RebinVar(UInt_t xAxisBins_rebin, Double_t *xAxis_rebin, TH1D **var1, UInt_t vars);

std::tuple<TH1D*, TGraphAsymmErrors*> CalculateDataRRatio(TH1D *stat1,TGraphAsymmErrors *sys1,TGraph *sysCut1, TH1D *stat2,TGraphAsymmErrors *sys2, TGraph *sysCut2, Double_t *sysJES, Double_t shift);
TGraphAsymmErrors* CalculateSimRRatio(TGraphAsymmErrors *hist1,TH1D **var1,TGraphAsymmErrors *hist2,TH1D **var2, UInt_t vars);
TGraphAsymmErrors* Divide(TGraphAsymmErrors *hist1,TH1D *hist2);
void TerminateCanvas4x6Pad(TPad **pad, TH1D** histo);

std::tuple<TH1D*, TGraphAsymmErrors*> Rebin(UInt_t xAxisBins_rebin, Double_t *xAxis_rebin, TH1D* histo1, TGraphAsymmErrors* histo2);

void finalJetSpectraInvUpdatedRatioZ(
)

{

    // ----------------------------------------------------------------
    // ------------------- ENABLE SIMULATIONS HERE --------------------
    // ----------------------------------------------------------------
    Bool_t ePowhegPythia8 = true;
    Bool_t ePythia8 = true;
    Bool_t ePythia8SoftMode2 = true;


    // ----------------------------------------------------------------

    UInt_t xAxisBins = 5;
  //  UInt_t xAxisBins_rebin = 7;
 //   Double_t *xAxis = new Double_t[xAxisBins+1]{5,6,8,10,12,14,20,30,50};
 //   Double_t *xAxis_rebin = new Double_t[xAxisBins+1]{5,6,8,10,14,20,30,50};
    Double_t *xAxisl = new Double_t[xAxisBins+1]{0.4,0.6,0.7,0.8,0.9,0.9999};
    Double_t *xAxisr = new Double_t[xAxisBins+1]{0.4,0.6,0.7,0.8,0.9,1.0};
  //  Double_t *xAxis_rebinl = new Double_t[xAxisBins+1]{5,6,8,10,14,20,30,49.99};
  //  Double_t *xAxis_rebinr = new Double_t[xAxisBins+1]{5,6,8,10,14,20,30,50};
  //  UInt_t xAxisBins = 5;
   // Double_t *xAxis = new Double_t[xAxisBins+1]{0.4,0.6,0.7,0.8,0.9,1.0};
   // plotRanges[0]=0; plotRanges[1]=2.1; plotRanges[2]=0.1; plotRanges[3]=3.5;
     plotRanges[0]=0.4; plotRanges[1]=1.6; plotRanges[2]=0.0; plotRanges[3]=5;


    // ----------------- Initiate Canvas ---------------------

    TCanvas *canvas;
    TPad **Pad;
    TH1D **placeholder;
    //std::tie(canvas,Pad,placeholder) = PrepareCanvas3x8Pad(xAxisBins,xAxisl, xAxisr);
    std::tie(canvas,Pad,placeholder) =PrepareCanvas4x6Pad(xAxisBins,xAxisl, xAxisr);
    //Int_t markers[3]{20,21,47};
    Int_t markers[3]{20,22,47};
    Double_t markersize[3]{1,1.1,1.3};

    // ----------------- data ---------------------
    TH1D *s13_hData_binned[12],*s13_hData_binned_ratio[12];
    TGraphAsymmErrors *s13_hData_binned_syst[12],*s13_hData_binned_ratio_syst[12];
    for(Int_t iradius = 0; iradius < 12;iradius++){
        std::tie(s13_hData_binned[iradius], s13_hData_binned_syst[iradius], s13_hData_binned_ratio[iradius], s13_hData_binned_ratio_syst[iradius])
                =GetData(fFiles[iradius],"hData_binned","haeData_binned_syst","hData_binned_ratio","haeData_binned_syst_ratio");
        //PlaceOnPadData(Pad[2*iradius-iradius%3],s13_hData_binned_syst[iradius],s13_hData_binned[iradius],0.9,20,-999);
        //PlaceOnPadData(Pad[2*iradius-iradius%3+3],s13_hData_binned_ratio_syst[iradius],s13_hData_binned_ratio[iradius],0,20,-999);
       // s13_hData_binned_syst[iradius/4]->SetMarkerStyle(markers)
        PlaceOnPadData(Pad[2*iradius-iradius%4],s13_hData_binned_syst[iradius],s13_hData_binned[iradius],markersize[iradius/4], markers[iradius/4],kRed+2,-999);
        PlaceOnPadData(Pad[2*iradius-iradius%4+4],s13_hData_binned_ratio_syst[iradius],s13_hData_binned_ratio[iradius],0, markers[iradius/4],kRed+2,-999);
    }


    TGraphAsymmErrors *s13_PowhegPythia8[12],*s13_PowhegPythia8_ratio[12];
    TH1D **s13_PowhegPythia8_var[12];

    TString dummy = "";
    if(ePowhegPythia8){
       for(Int_t iradius = 0; iradius < 12;iradius++){
           std::tie(s13_PowhegPythia8[iradius], s13_PowhegPythia8_ratio[iradius], s13_PowhegPythia8_var[iradius])
                   =GetSim(dummy+"13TeV"+Form("%d",iradius), fFiles[iradius],"PowhegPythia8",8);
           PlaceOnPad(Pad[2*iradius-iradius%4],s13_PowhegPythia8[iradius],-1,-1,-1,"same2p",-999);
           PlaceOnPad(Pad[2*iradius-iradius%4+4],s13_PowhegPythia8_ratio[iradius],-1,-1,-1,"same2p",-999);
       }
    }


    TGraphAsymmErrors *s13_Pythia8[12],*s13_Pythia8_ratio[12];
    TH1D **s13_Pythia8_var[12];
    if(ePythia8){
        for(Int_t iradius = 0; iradius < 12;iradius++){
           std::tie(s13_Pythia8[iradius], s13_Pythia8_ratio[iradius], s13_Pythia8_var[iradius])
                   =GetSim(dummy+"13TeV"+Form("%d",iradius), fFiles[iradius],"Pythia8",0);
           PlaceOnPad(Pad[2*iradius-iradius%4],s13_Pythia8[iradius],kViolet+2,-1,6,"same1",-999);
           PlaceOnPad(Pad[2*iradius-iradius%4+4],s13_Pythia8_ratio[iradius],kViolet+2,-1,6,"same1",-999);
        }
    }


    TGraphAsymmErrors *s13_Pythia8SoftMode2[12],*s13_Pythia8SoftMode2_ratio[12];
    TH1D **s13_Pythia8SoftMode2_var[12];
    if(ePythia8SoftMode2){
        for(Int_t iradius = 0; iradius < 12;iradius++){
           std::tie(s13_Pythia8SoftMode2[iradius], s13_Pythia8SoftMode2_ratio[iradius], s13_Pythia8SoftMode2_var[iradius])
                   =GetSim(dummy+"13TeV"+Form("%d",iradius), fFiles[iradius],"Pythia8Soft2",0);
           PlaceOnPad(Pad[2*iradius-iradius%4],s13_Pythia8SoftMode2[iradius],kGreen+2,-1,7,"same1",-999);
           PlaceOnPad(Pad[2*iradius-iradius%4+4],s13_Pythia8SoftMode2_ratio[iradius],kGreen+2,-1,7,"same1",-999);
        }
    }

    TPaveText *ptAlice = new TPaveText(0.25,0.7,0.85,0.9,"NB NDC");
    ptAlice->SetBorderSize(0);
    ptAlice->SetFillStyle(0);
    ptAlice->SetTextAlign(13);
    ptAlice->SetTextFont(43);
    ptAlice->SetTextSize(18);
    ptAlice->AddText("ALICE");

    TPaveText *ptInfo = new TPaveText(0.25,0.64,0.85,0.75,"NB NDC");
    ptInfo->SetBorderSize(0);
    ptInfo->SetFillStyle(0);
    ptInfo->SetTextAlign(13);
    ptInfo->SetTextFont(43);
    ptInfo->SetTextSize(18);
    ptInfo->AddText("Charged jets, Anti-#it{k}_{T}");

    TPaveText *ptRadius[3];
    ptRadius[0] = new TPaveText(0.2,0.4,0.95,0.8,"NB NDC");
    ptRadius[1] = new TPaveText(0.2,0.7,0.95,0.8,"NB NDC");
    ptRadius[2] = new TPaveText(0.2,0.7,0.95,0.8,"NB NDC");
    for(Int_t s = 0; s<3;s++){
        ptRadius[s]->SetBorderSize(0);
        ptRadius[s]->SetFillStyle(0);
        ptRadius[s]->SetTextAlign(13);
        ptRadius[s]->SetTextFont(43);
        ptRadius[s]->SetTextSize(19);
    }
    ptRadius[0]->AddText("R=0.2");
    ptRadius[1]->AddText("R=0.4");
    ptRadius[2]->AddText("R=0.6");

    TPaveText *ptRadius2[12];
    ptRadius2[0] = new TPaveText(0.25,0.4,0.95,0.8,"NB NDC");
    ptRadius2[1] = new TPaveText(0.09,0.75,0.95,0.8,"NB NDC");
    ptRadius2[2] = new TPaveText(0.09,0.75,0.95,0.8,"NB NDC");
    ptRadius2[3] = new TPaveText(0.09,0.75,0.95,0.8,"NB NDC");
    ptRadius2[4] = new TPaveText(0.25,0.85,0.95,0.95,"NB NDC");
    ptRadius2[5] = new TPaveText(0.09,0.85,0.95,0.95,"NB NDC");
    ptRadius2[6] = new TPaveText(0.09,0.85,0.95,0.95,"NB NDC");
    ptRadius2[7] = new TPaveText(0.09,0.85,0.95,0.95,"NB NDC");
    ptRadius2[8] = new TPaveText(0.25,0.85,0.95,0.95,"NB NDC");
    ptRadius2[9] = new TPaveText(0.09,0.85,0.95,0.95,"NB NDC");
    ptRadius2[10]= new TPaveText(0.09,0.85,0.95,0.95,"NB NDC");
    ptRadius2[11]= new TPaveText(0.09,0.85,0.95,0.95,"NB NDC");
    for(Int_t s = 0; s<12;s++){
        ptRadius2[s]->SetBorderSize(0);
        ptRadius2[s]->SetFillStyle(0);
        ptRadius2[s]->SetTextAlign(13);
        ptRadius2[s]->SetTextFont(43);
        ptRadius2[s]->SetTextSize(15);
    }

    ptRadius2[0]->AddText("5 < #it{p}_{T,jet} < 7 GeV/#it{c}, #it{p}_{T,D} > 2 GeV/#it{c}");
    ptRadius2[1]->AddText("7 < #it{p}_{T,jet} < 10 GeV/#it{c}, #it{p}_{T,D} > 4 GeV/#it{c}");
    ptRadius2[2]->AddText("10 < #it{p}_{T,jet} < 15 GeV/#it{c}, #it{p}_{T,D} > 5 GeV/#it{c}");
    ptRadius2[3]->AddText("15 < #it{p}_{T,jet} < 50 GeV/#it{c}, #it{p}_{T,D} > 10 GeV/#it{c}");

    ptRadius2[4]->AddText("5 < #it{p}_{T,jet} < 7 GeV/#it{c}, #it{p}_{T,D} > 2 GeV/#it{c}");
    ptRadius2[5]->AddText("7 < #it{p}_{T,jet} < 10 GeV/#it{c}, #it{p}_{T,D} > 3 GeV/#it{c}");
    ptRadius2[6]->AddText("10 < #it{p}_{T,jet} < 15 GeV/#it{c}, #it{p}_{T,D} > 5 GeV/#it{c}");
    ptRadius2[7]->AddText("15 < #it{p}_{T,jet} < 50 GeV/#it{c}, #it{p}_{T,D} > 5 GeV/#it{c}");

    ptRadius2[8]->AddText("5 < #it{p}_{T,jet} < 7 GeV/#it{c}, #it{p}_{T,D} > 2 GeV/#it{c}");
    ptRadius2[9]->AddText("7 < #it{p}_{T,jet} < 10 GeV/#it{c}, #it{p}_{T,D} > 3 GeV/#it{c}");
    ptRadius2[10]->AddText("10 < #it{p}_{T,jet} < 15 GeV/#it{c}, #it{p}_{T,D} > 5 GeV/#it{c}");
    ptRadius2[11]->AddText("15 < #it{p}_{T,jet} < 50 GeV/#it{c}, #it{p}_{T,D} > 5 GeV/#it{c}");

    TPaveText *ptEnergy[3];
    ptEnergy[0] = new TPaveText(0.45,0.72,0.95,0.9,"NB NDC");
    ptEnergy[1] = new TPaveText(0.55,0.85,0.95,0.9,"NB NDC");
    ptEnergy[2] = new TPaveText(0.2,0.85,0.85,0.9,"NB NDC");
    for(Int_t s = 0; s<3;s++){
        ptEnergy[s]->SetBorderSize(0);
        ptEnergy[s]->SetFillStyle(0);
        ptEnergy[s]->SetTextAlign(13);
        ptEnergy[s]->SetTextFont(43);
        ptEnergy[s]->SetTextSize(18);
    }
    //ptEnergy[0]->AddText("pp, #sqrt{#it{s}} = 13 TeV");
    ptEnergy[0]->AddText("pp, #sqrt{#it{s}} = 5.02 TeV");

    TLegend *leg =nullptr;
    TLegend *leg2 =nullptr;
    TLegend *leg3 =nullptr;
    Double_t shift = 0.06*(2+3);
    leg = new TLegend(0.08,0.72-shift,0.72,0.675,nullptr,"NB NDC");
    leg2= new TLegend(0.08,0.9-shift/2.,0.72,0.855,nullptr,"NB NDC");
    leg3= new TLegend(0.08,0.9-shift/2.,0.72,0.855,nullptr,"NB NDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(43); leg2->SetTextFont(43); leg3->SetTextFont(43);
    leg->SetTextSize(14); leg2->SetTextSize(14); leg3->SetTextSize(14);
    leg->SetLineColor(1); leg2->SetLineColor(1); leg3->SetLineColor(1);
    leg->SetLineStyle(1); leg2->SetLineStyle(1); leg3->SetLineStyle(1);
    leg->SetLineWidth(1); leg2->SetLineWidth(1); leg3->SetLineWidth(1);
    leg->SetFillColor(0); leg2->SetFillColor(0); leg3->SetFillColor(0);
    leg->SetFillStyle(0); leg2->SetFillStyle(0); leg3->SetFillStyle(0);
    leg->AddEntry(s13_hData_binned_syst[0],"R=0.2","fp");
    leg2->AddEntry(s13_hData_binned_syst[4],"R=0.4","fp");
    leg3->AddEntry(s13_hData_binned_syst[8],"R=0.6","fp");
    leg->AddEntry(s13_PowhegPythia8[0],"POWHEG hvq + PYTHIA 8","pf");
    leg->AddEntry(s13_Pythia8[0],"PYTHIA 8 Monash 2013","l");
    leg->AddEntry(s13_Pythia8SoftMode2[0],"PYTHIA 8 Monash 2013 Soft mode 2","l");


    Pad[0]->cd();  ptRadius2[0]->Draw();
    Pad[1]->cd();  ptRadius2[1]->Draw();
    Pad[2]->cd();  ptRadius2[2]->Draw();
    Pad[3]->cd();  ptRadius2[3]->Draw();
    Pad[8]->cd();  ptRadius2[4]->Draw();
    Pad[9]->cd();  ptRadius2[5]->Draw();
    Pad[10]->cd(); ptRadius2[6]->Draw();
    Pad[11]->cd(); ptRadius2[7]->Draw();
    Pad[16]->cd(); ptRadius2[8]->Draw();
    Pad[17]->cd(); ptRadius2[9]->Draw();
    Pad[18]->cd(); ptRadius2[10]->Draw();
    Pad[19]->cd(); ptRadius2[11]->Draw();

    Pad[0]->cd();
    ptAlice->Draw();
    ptInfo->Draw();
    ptEnergy[0]->Draw();
    Pad[3]->cd();
    leg->Draw();
    Pad[11]->cd();
    leg2->Draw();
    Pad[19]->cd();
    leg3->Draw();

    TerminateCanvas4x6Pad(Pad,placeholder);
    canvas->SaveAs("PDFF.png");
    canvas->SaveAs("PDFF.pdf");


    return;
}

std::tuple<TH1D*, TGraphAsymmErrors*> CalculateDataRRatio(TH1D *stat1,TGraphAsymmErrors *sys1,TGraph *sysCut1, TH1D *stat2,TGraphAsymmErrors *sys2, TGraph *sysCut2, Double_t *sysJES, Double_t shift){
    TString name1 = stat1->GetName();
    TString name2 = stat2->GetName();
    TH1D *tmp = dynamic_cast<TH1D*>(stat1->Clone(name1+"_"+name2+"_ratio"));
    tmp->Divide(stat2);

  /*  if(shift > -900){
        for(Int_t i=1;i<=tmp->GetNbinsX();i++){
            tmp->SetBinContent(i,tmp->GetBinContent(i)+shift);
        }
    }*/

    UInt_t nBins = static_cast<UInt_t>(sys1->GetN());
    Double_t *sysuncAbs_down = new Double_t[nBins];
    Double_t *sysuncAbs_up = new Double_t[nBins];
    Double_t *statunc = new Double_t[nBins];
    Double_t *value = new Double_t[nBins];
    Double_t *xval = new Double_t[nBins];
    Double_t *xvalwidth_down = new Double_t[nBins];
    Double_t *xvalwidth_up = new Double_t[nBins];

   // TGraphAsymmErrors *grsysl = new TGraphAsymmErrors(static_cast<Int_t>(nBins),xval,value,xvalwidth,xvalwidth,sysuncAbs_down,sysuncAbs_up);

    for(Int_t bin = 0; bin < sys1->GetN();bin++){
        Double_t x1 = 0,y1=0;
        Double_t x2 = 0,y2=0;
        Double_t x1c = 0,y1c=0;
        Double_t x2c = 0,y2c=0;
        Double_t JES = 0;
        sys1->GetPoint(bin,x1,y1);
        sys2->GetPoint(bin,x2,y2);
        if(sysCut1)sysCut1->GetPoint(bin,x1c,y1c);
        if(sysCut2)sysCut2->GetPoint(bin,x2c,y2c);
        if(sysJES) JES = sysJES[bin];


        if(bin ==0){
        std::cout<<"hist1:"<<std::endl;
        std::cout<<x1<<" "<<y1<<std::endl;//<<" "<<hData_binned4->GetBinError(bin)/hData_binned4->GetBinContent(bin)<<" "<<hDataR24->GetBinError(bin)/hDataR24->GetBinContent(bin)<<std::endl;
        std::cout<<sys1->GetErrorXlow(bin)<<" "<<sys1->GetErrorXhigh(bin)<<" "<<sys1->GetErrorYlow(bin)/y1<<" "<<sys1->GetErrorYhigh(bin)/y1<<"cut: "<<x1c<<" "<<y1c<<" jes "<<JES<<std::endl;

        std::cout<<"hist2:"<<std::endl;
        std::cout<<sys2->GetErrorXlow(bin)<<" "<<sys2->GetErrorXhigh(bin)<<" "<<sys2->GetErrorYlow(bin)/y2<<" "<<sys2->GetErrorYhigh(bin)/y2<<"cut: "<<x2c<<" "<<y2c<<" jes "<<JES<<std::endl;
        std::cout<<"merge:"<<std::endl;
        std::cout<<"lq: "<<sys1->GetErrorYlow(bin)/y1<<" "<<sys2->GetErrorYlow(bin)/y2<<std::endl;
        std::cout<<"lq: "<<sys1->GetErrorYlow(bin)/y1*sys1->GetErrorYlow(bin)/y1<<" "<<sys2->GetErrorYlow(bin)/y2*sys2->GetErrorYlow(bin)/y2<<std::endl;
        std::cout<<"lq: "<<sys1->GetErrorYlow(bin)/y1*sys1->GetErrorYlow(bin)/y1+sys2->GetErrorYlow(bin)/y2*sys2->GetErrorYlow(bin)/y2<<std::endl;
       std::cout<<"lq: "<<TMath::Sqrt(sys1->GetErrorYlow(bin)/y1*sys1->GetErrorYlow(bin)/y1  +sys2->GetErrorYlow(bin)/y2*sys2->GetErrorYlow(bin)/y2)<<std::endl;
       std::cout<<"uq: "<<TMath::Sqrt(sys1->GetErrorYhigh(bin)/y1*sys1->GetErrorYhigh(bin)/y1 + sys2->GetErrorYhigh(bin)/y2*sys2->GetErrorYhigh(bin)/y2)<<std::endl;
}
       Double_t sysUP = TMath::Sqrt(sys1->GetErrorYhigh(bin)/y1*sys1->GetErrorYhigh(bin)/y1 + sys2->GetErrorYhigh(bin)/y2*sys2->GetErrorYhigh(bin)/y2 + (y1c+y1c)*(y1c+y1c)/4 + JES*JES);
       Double_t sysDown = TMath::Sqrt(sys1->GetErrorYlow(bin)/y1*sys1->GetErrorYlow(bin)/y1 + sys2->GetErrorYlow(bin)/y2*sys2->GetErrorYlow(bin)/y2 + (y1c+y1c)*(y1c+y1c)/4 + JES*JES) ;

       xval[bin] = x1;
       xvalwidth_down[bin] = sys1->GetErrorXlow(bin);
       xvalwidth_up[bin] = sys1->GetErrorXhigh(bin);
       value[bin] = y1/y2;
      // std::cout<<j<<" "<<xval[j]<<" "<<xvalwidth[j]<<" "<<value[j]<<" "<<error<<std::endl;
       sysuncAbs_down[bin] = value[bin] * sysDown;
       sysuncAbs_up[bin] = value[bin] * sysUP;

     /*  if(shift > -900){
          value[bin]+=shift;
          //sysuncAbs_down[bin]+=shift;
          //sysuncAbs_up[bin]+=shift;
       }*/
       //if(value[j]>1E-9)statunc[j] = error/ value[j] *100;
std::cout<<"sum: bin: "<<bin<<" xval: "<<x1<<" xwidth low: "<<xvalwidth_down[bin]<<" xwidth high: "<<xvalwidth_down[bin]<<" value: "<<value[bin]<<" sysdown: "<<sysuncAbs_down[bin]<<" sysup: "<<sysuncAbs_up[bin]<<std::endl;
    }
    TGraphAsymmErrors *grsysl = new TGraphAsymmErrors(static_cast<Int_t>(nBins),xval,value,xvalwidth_down,xvalwidth_up,sysuncAbs_down,sysuncAbs_up);
    grsysl->SetName(name1+"_"+name2+"_ratio_sys");
    return std::make_tuple(tmp,grsysl);
}

TGraphAsymmErrors* CalculateSimRRatio(TGraphAsymmErrors *hist1,TH1D **var1,TGraphAsymmErrors *hist2,TH1D **var2, UInt_t vars){
    std::cout<<"GETTING SIMS"<<std::endl;
    TString name1 = hist1->GetName();
    TString name2 = hist2->GetName();
    std::cout<<name1+name2<<std::endl;
    TH1D **var1var2 = new TH1D*[static_cast<ULong_t>(vars)];
    TH1D * hPrompt_up = nullptr;
    TH1D * hPrompt_down = nullptr;
    TH1D * central1 = nullptr;
    TH1D * central2 = nullptr;
std::cout<<"a"<<std::endl;
    if(var1!=nullptr){
        for(UInt_t ivar = 0; ivar < vars;ivar++){
            std::cout<<ivar<<std::endl;
            if(!var1[ivar])std::cout<<"ivar not exist: "<<ivar<<std::endl;
            var1var2[ivar] = dynamic_cast<TH1D*>(var1[ivar]->Clone(Form(name1+"_"+name2+"_%d",ivar)));
            var1var2[ivar]->Divide(var2[ivar]);
            std::cout<<"VAR: " <<ivar<<" v1: "<<var1[ivar]->GetBinContent(1)<<" v2: "<<var2[ivar]->GetBinContent(1)<<" v1v2: "<<var1var2[ivar]->GetBinContent(1)<<std::endl;

        }
        std::cout<<"a"<<std::endl;
        central1 = dynamic_cast<TH1D*>(var1[0]->Clone(name1+"_"+name2+"_central1"));
        central2 = dynamic_cast<TH1D*>(var1[0]->Clone(name1+"_"+name2+"_central2"));
std::cout<<"a"<<std::endl;
        for(UInt_t bin=0; bin<hist1->GetN(); bin++){
            Double_t x1 = 0,y1=0;
            Double_t x2 = 0,y2=0;
            hist1->GetPoint(bin,x1,y1);
            hist2->GetPoint(bin,x2,y2);
            central1->SetBinContent(bin+1,y1);
            central2->SetBinContent(bin+1,y2);
            central1->SetBinError(bin+1,0);
            central2->SetBinError(bin+1,0);
            std::cout<<"central1: "<<central1->GetBinCenter(bin+1)<<" "<<y1<<" central2: "<<central2->GetBinCenter(bin+1)<<" "<<y2<<" div "<<y1/y2<<std::endl;
        }
        std::cout<<"a"<<std::endl;
        central1->Divide(central2);
        for(UInt_t bin=1; bin<=central1->GetNbinsX(); bin++)std::cout<<"central val: "<<central1->GetBinContent(bin)<<std::endl;

        hPrompt_up = GetUpSys(central1,var1var2,vars);
        hPrompt_up->SetName(name1+"_"+name2+"_up");
        // get down unc
        hPrompt_down = GetDownSys(central1,var1var2,vars);
        hPrompt_down->SetName(name1+"_"+name2+"_down");

        std::cout<<" UP: "<<hPrompt_up->GetBinContent(1)<<" DOWN: "<<hPrompt_down->GetBinContent(1)<<std::endl;
    }
    std::cout<<"b"<<std::endl;

    UInt_t nBin = hist1->GetN();
    Double_t *xval = new Double_t[nBin];
    Double_t *xvalwidth_down = new Double_t[nBin];
    Double_t *xvalwidth_up = new Double_t[nBin];
    Double_t *value = new Double_t[nBin];
    Double_t *valuetheoryerrup = new Double_t[nBin];
    Double_t *valuetheoryerrdown = new Double_t[nBin];
    for(UInt_t bin=0; bin<nBin; bin++){
        Double_t x1 = 0,y1=0;
        Double_t x2 = 0,y2=0;
        hist1->GetPoint(bin,x1,y1);
        hist2->GetPoint(bin,x2,y2);
        std::cout<<"SIM: "<<y1<<" "<<y2<< " "<<y1/y2<<std::endl;
        xval[bin] = x1;
        xvalwidth_down[bin] = hist1->GetErrorXlow(bin);
        xvalwidth_up[bin] = hist1->GetErrorXhigh(bin);
        value[bin] = y1/y2;
        if(var1!=nullptr){
            valuetheoryerrup[bin] = hPrompt_up->GetBinContent(hPrompt_up->GetXaxis()->FindBin(xval[bin])) - value[bin];
            valuetheoryerrdown[bin] = value[bin] - hPrompt_down->GetBinContent(hPrompt_up->GetXaxis()->FindBin(xval[bin]));
            std::cout<<"UP: "<<valuetheoryerrup[bin]<<" DOWN: "<<valuetheoryerrdown[bin]<<std::endl;
        }
        else{
            valuetheoryerrup[bin] = 0;
            valuetheoryerrdown[bin] = 0;
        }

    }
    TGraphAsymmErrors *grsystheory = new TGraphAsymmErrors(static_cast<Int_t>(nBin),xval,value,xvalwidth_down,xvalwidth_up,valuetheoryerrdown,valuetheoryerrup);
    grsystheory->SetName(name1+"_"+name2+"_ratio_sys");
    return grsystheory;

}

TH1D** RebinVar(UInt_t xAxisBins_rebin, Double_t *xAxis_rebin, TH1D **var1, UInt_t vars){
    TString nametmp = "";
    std::cout<<"REBIN"<<std::endl;
    TH1D **var1_rebin = new TH1D*[static_cast<ULong_t>(vars)];
    TString name1 = var1[0]->GetName();
    for(UInt_t ivar = 0; ivar < vars;ivar++){
        std::cout<<ivar<<std::endl;
        if(!var1[ivar])std::cout<<"ivar not exist: "<<ivar<<std::endl;
        var1_rebin[ivar] = new TH1D(nametmp+var1[0]->GetName()+"_rebinned"+Form("%d",ivar),nametmp+var1[0]->GetName()+"_rebinned"+Form("%d",ivar),xAxisBins_rebin,xAxis_rebin);
        for(Int_t bin = 1; bin <= var1_rebin[ivar]->GetNbinsX(); bin++){
            if(bin <= 3){
                var1_rebin[ivar]->SetBinContent(bin,var1[ivar]->GetBinContent(bin));
                var1_rebin[ivar]->SetBinError(bin,var1[ivar]->GetBinError(bin));
            }
            else if(bin == 4){
                Double_t w1 = 1./var1[ivar]->GetBinContent(bin);
                Double_t w2 = 1./var1[ivar]->GetBinContent(bin+1);
                Double_t x1 = var1[ivar]->GetBinContent(bin);
                Double_t x2 = var1[ivar]->GetBinContent(bin+1);
                Double_t e1 = var1[ivar]->GetBinError(bin);
                Double_t e2 = var1[ivar]->GetBinError(bin+1);
                Double_t val = (w1*x1+w2*x2)/(w1+w2);
                Double_t err = (w1*e1+w2*e2)/(w1+w2);
                var1_rebin[ivar]->SetBinContent(bin,val);
                var1_rebin[ivar]->SetBinError(bin,err);
            }
            else if(bin >= 5){
                var1_rebin[ivar]->SetBinContent(bin,var1[ivar]->GetBinContent(bin+1));
                var1_rebin[ivar]->SetBinError(bin,var1[ivar]->GetBinError(bin+1));
            }
        }
    }
return var1_rebin;
}

void TerminateCanvas(TPad* pad1, TPad* pad2, TPad *pad3, TH1D* histo1, TH1D* histo2, TH1D *histo3){
    pad1->cd();
    histo1->Draw("sameaxis");
    pad2->cd();
    histo2->Draw("sameaxis");
    histo2->Draw("sameaxig");
    pad3->cd();
    histo3->Draw("sameaxis");
    histo3->Draw("sameaxig");
}


void TerminateCanvas3x2Pad(TPad **pad, TH1D** histo){
    std::cout<<"ipad0"<<std::endl;
    pad[0]->cd();
    histo[0]->Draw("sameaxis");
    std::cout<<"ipad1"<<std::endl;
    pad[1]->cd();
    histo[1]->Draw("sameaxis");
    pad[2]->cd();
    histo[2]->Draw("sameaxis");
    for(Int_t ipad = 3; ipad<6;ipad++){
        pad[ipad]->cd();
        histo[ipad]->Draw("sameaxis");
        histo[ipad]->Draw("sameaxig");
    }
}

void TerminateCanvas4x6Pad(TPad **pad, TH1D** histo){
    std::cout<<"ipad0"<<std::endl;
    TLine *line = new TLine(0.4,1,1,1);
    line->SetLineStyle(2);
    for(Int_t ipad = 0; ipad<24;ipad++){
        if(ipad%8==4 || ipad%8==5 || ipad%8==6 || ipad%8==7){
            pad[ipad]->cd();
            histo[ipad]->Draw("sameaxis");
            histo[ipad]->Draw("sameaxig");
            line->Draw("same");
        }
        else{
            pad[ipad]->cd();
            histo[ipad]->Draw("sameaxis");
        }
    }
}





void PlaceOnPad(TPad* pad,TGraphAsymmErrors* histo, Color_t ci, Style_t style, Style_t linestyle, TString opt,Double_t scale, TString mode){
    pad->cd();
    histo->SetFillColor(1);
    histo->SetFillStyle(0);
    if(ci!=-1)histo->SetLineColor(ci);
    if(ci!=-1)histo->SetMarkerColor(ci);
   if(style!=-1) histo->SetMarkerStyle(style);
    if(linestyle!=-1)histo->SetLineStyle(linestyle);
   // histo->SetFillStyle(3005);
   // histo->SetMarkerSize(markersize); //add up
    histo->SetLineWidth(1);
    if(scale>-1){
        TString tmp = "";
        TGraphAsymmErrors* histo_tmp = dynamic_cast<TGraphAsymmErrors*>(histo->Clone(tmp+histo->GetName()+"_plotCp"));
        if(mode == "+"){
            for (int i=0;i<histo_tmp->GetN();i++) histo_tmp->GetY()[i] += scale;
        }else if(mode =="*"){
            for (int i=0;i<histo_tmp->GetN();i++) histo_tmp->GetY()[i] *= scale;
            for (int i=0;i<histo_tmp->GetN();i++) histo_tmp->GetEYhigh()[i] *= scale;
            for (int i=0;i<histo_tmp->GetN();i++) histo_tmp->GetEYlow()[i] *= scale;
        }
        if(opt !="") histo_tmp->Draw(opt);
        else histo_tmp->Draw("psame");
    }
    else if(opt !="") histo->Draw(opt);
    else histo->Draw("psame");
}

void PlaceOnPad(TPad* pad,TH1D* histo, Color_t ci, Style_t style, Style_t linestyle, TString opt,Double_t scale, TString mode){
    pad->cd();
    histo->SetFillColor(1);
    histo->SetFillStyle(0);
    histo->SetLineColor(ci);
    histo->SetMarkerColor(ci);
    histo->SetMarkerStyle(style);
    histo->SetLineStyle(linestyle);
   // histo->SetFillStyle(3005);
   // histo->SetMarkerSize(markersize); //add up
    histo->SetLineWidth(2);
    if(scale>-1){
        TString tmp = "";
        TH1D* histo_tmp = dynamic_cast<TH1D*>(histo->Clone(tmp+histo->GetName()+"_plotCp"));
        if(mode == "+"){
            for (int i=1;i<=histo_tmp->GetNbinsX();i++) histo_tmp->SetBinContent(i, histo_tmp->GetBinContent(i)+scale);
        }else if(mode =="*"){
            histo_tmp->Scale(scale);
        }
        if(opt !="") histo_tmp->Draw(opt);
        else histo_tmp->Draw("psame");
    }
    else if(opt !="") histo->Draw(opt);
    else histo->Draw("psame");
}

void PlaceOnPadData(TPad* pad,TGraphAsymmErrors* histo1, TH1D* histo2, Size_t markersize, Style_t markerstyle, Color_t histo2col, Double_t scale,  TString mode ){
    pad->cd();
    Color_t ci;// = static_cast<Color_t>(TColor::GetColor("#990000"));
    //histo1->SetLineColor(ci);
    //histo1->SetMarkerColor(ci);
    histo1->SetMarkerStyle(markerstyle);
    histo1->SetMarkerSize(markersize);
    ci = static_cast<Color_t>(TColor::GetColor("#cccccc"));
    histo1->SetFillColor(ci);
    histo1->SetLineColor(ci);
    histo1->SetLineWidth(0);
    if(scale>-999){
        TString tmp = "";
        TGraphAsymmErrors* histo_tmp = dynamic_cast<TGraphAsymmErrors*>(histo1->Clone(tmp+histo1->GetName()+"_plotCp"));
        if(mode == "+"){
            for (int i=0;i<histo_tmp->GetN();i++) histo_tmp->GetY()[i] += scale;
        }else if(mode =="*"){
            for (int i=0;i<histo_tmp->GetN();i++) histo_tmp->GetY()[i] *= scale;
            for (int i=0;i<histo_tmp->GetN();i++) histo_tmp->GetEYhigh()[i] *= scale;
            for (int i=0;i<histo_tmp->GetN();i++) histo_tmp->GetEYlow()[i] *= scale;
        }
        histo_tmp->Draw("2p same");
    }    else histo1->Draw("2p same");
    //data central w stat. unc.
    histo2->SetMarkerStyle(markerstyle);
    histo2->SetMarkerSize(markersize);
    histo2->SetMarkerColor(histo2col);
    histo2->SetLineColor(histo2col);
    if(scale>-999){
        TString tmp = "";
        TH1D* histo_tmp = dynamic_cast<TH1D*>(histo2->Clone(tmp+histo2->GetName()+"_plotCp"));
        if(mode == "+"){
            for (int i=1;i<=histo_tmp->GetNbinsX();i++) histo_tmp->SetBinContent(i, histo_tmp->GetBinContent(i)+scale);
        }else if(mode =="*"){
            histo_tmp->Scale(scale);
        }
        //for(Int_t i = 1 ; i<= histo_tmp->GetNbinsX();i++){
        //    histo_tmp->SetBinContent(i,histo_tmp->GetBinContent(i)*scale);
        //}
        histo_tmp->Draw("same p  e0 x0");
    }    else histo2->Draw("same p  e0 x0");
};


std::tuple<TCanvas*, TPad**, TH1D**> PrepareCanvas3x8Pad(UInt_t xAxisBins, Double_t *xAxisl, Double_t *xAxisr){
    style();
    //prepare main canvas
    TCanvas *FinalSpectrum = new TCanvas("FinalSpectrum4", "FinalSpectrum4",0,45,1200,1650);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    FinalSpectrum->SetHighLightColor(2);
    FinalSpectrum->Range(0,0,1,1);
    FinalSpectrum->SetFillColor(0);
    FinalSpectrum->SetBorderMode(0);
    FinalSpectrum->SetBorderSize(2);
    FinalSpectrum->SetFrameBorderMode(0);
    TPad **pad = new TPad*[static_cast<ULong_t>(24)];
    //0  1  2 norm
    //3  4  5 ratio
    //6 7 8   norm
    //9 10 11 ratio

    Double_t marginLeft=0.07;
    Double_t marginRight=0.02;
    Double_t innerPadWidth=(1-marginLeft-marginRight)/3.; //0.91/2 = 0.455
    Double_t marginLeftForXAxis=0.02;

    Double_t marginTop=0.04;
    Double_t marginBottom=0.06;
    Double_t innerPadHeight=(1-marginTop-marginBottom)/16.;
    Double_t marginBottomForYAxis=0.0;

    auto xbaseL = [innerPadWidth,marginLeft,marginLeftForXAxis](int m) ->Double_t {return (m*innerPadWidth+marginLeft-marginLeftForXAxis);};
    auto xbaseR = [innerPadWidth,marginLeft](int m) ->Double_t {return (m*innerPadWidth+marginLeft);};
    auto ybaseL = [innerPadHeight,marginBottom,marginBottomForYAxis](int m) ->Double_t {return (m*innerPadHeight+marginBottom-marginBottomForYAxis);};
    auto ybaseR = [innerPadHeight,marginBottom](int m) ->Double_t {return (m*innerPadHeight+marginBottom);};

    //pads for the plots
    for(Int_t ipad = 0;ipad < 24;ipad++){
        std::cout<<"padID: "<<ipad<<std::endl;
        if(ipad ==0) pad[ipad] = new TPad("pad0", "pad0",0,        ybaseL(13)  ,xbaseR(1)    ,1);
        if(ipad ==1) pad[ipad] = new TPad("pad1", "pad1",xbaseL(1) ,  ybaseL(13)  ,xbaseR(2)  ,1);
        if(ipad ==2) pad[ipad] = new TPad("pad2", "pad2",xbaseL(2), ybaseL(13)  ,1         ,1);
        if(ipad ==3) pad[ipad] = new TPad("pad3", "pad3",0,        ybaseL(12)  ,xbaseR(1)    ,ybaseR(13));
        if(ipad ==4) pad[ipad] = new TPad("pad4", "pad4",xbaseL(1) ,  ybaseL(12)  ,xbaseR(2)  ,ybaseR(13));
        if(ipad ==5) pad[ipad] = new TPad("pad5", "pad5",xbaseL(2), ybaseL(12)  ,1         ,ybaseR(13));

        if(ipad ==6) pad[ipad] = new TPad("pad6", "pad6",0,        ybaseL(9)  ,xbaseR(1)    ,ybaseR(12));
        if(ipad ==7) pad[ipad] = new TPad("pad7", "pad7",xbaseL(1) ,  ybaseL(9)  ,xbaseR(2)  ,ybaseR(12));
        if(ipad ==8) pad[ipad] = new TPad("pad8", "pad8",xbaseL(2), ybaseL(9)  ,1         ,ybaseR(12));
        if(ipad ==9) pad[ipad] = new TPad("pad9", "pad9",0,        ybaseL(8)  ,xbaseR(1)    ,ybaseR(9));
        if(ipad ==10) pad[ipad] = new TPad("pad10", "pad10",xbaseL(1) ,  ybaseL(8)  ,xbaseR(2)  ,ybaseR(9));
        if(ipad ==11) pad[ipad] = new TPad("pad11", "pad11",xbaseL(2), ybaseL(8)  ,1         ,ybaseR(9));

        if(ipad ==12) pad[ipad] = new TPad("pad12", "pad12",0,        ybaseL(5)  ,xbaseR(1)    ,ybaseR(8));
        if(ipad ==13) pad[ipad] = new TPad("pad13", "pad13",xbaseL(1) ,  ybaseL(5)  ,xbaseR(2)  ,ybaseR(8));
        if(ipad ==14) pad[ipad] = new TPad("pad14", "pad14",xbaseL(2), ybaseL(5)  ,1         ,ybaseR(8));
        if(ipad ==15) pad[ipad] = new TPad("pad15", "pad15",0,        ybaseL(4)  ,xbaseR(1)    ,ybaseR(5));
        if(ipad ==16) pad[ipad] = new TPad("pad16", "pad16",xbaseL(1) ,  ybaseL(4)  ,xbaseR(2)  ,ybaseR(5));
        if(ipad ==17) pad[ipad] = new TPad("pad17", "pad17",xbaseL(2), ybaseL(4)  ,1         ,ybaseR(5));

        if(ipad ==18) pad[ipad] = new TPad("pad18", "pad18",0,        ybaseL(1)  ,xbaseR(1)    ,ybaseR(4));
        if(ipad ==19) pad[ipad] = new TPad("pad19", "pad19",xbaseL(1) ,  ybaseL(1)  ,xbaseR(2)  ,ybaseR(4));
        if(ipad ==20) pad[ipad] = new TPad("pad20", "pad20",xbaseL(2), ybaseL(1)  ,1         ,ybaseR(4));
        if(ipad ==21) pad[ipad] = new TPad("pad21", "pad21",0,        0  ,xbaseR(1)    ,ybaseR(1));
        if(ipad ==22) pad[ipad] = new TPad("pad22", "pad22",xbaseL(1) ,  0  ,xbaseR(2)  ,ybaseR(1));
        if(ipad ==23) pad[ipad] = new TPad("pad23", "pad23",xbaseL(2), 0  ,1         ,ybaseR(1));


        pad[ipad]->Draw();
        pad[ipad]->cd();

        //left column
        if(ipad%3 ==0){
            pad[ipad]->SetLeftMargin(static_cast<Float_t>(marginLeft/(marginLeft+innerPadWidth))); //0.07/(0.07+0.455) = 0.13333
            pad[ipad]->SetRightMargin(0);
            pad[ipad]->Modified();
            pad[ipad]->SetFillStyle(0);
        }
        //right column
        else if (ipad%3 ==1){
            pad[ipad]->SetLeftMargin(static_cast<Float_t>(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis)));//0.02/(0.02+0.455) = 0.0421
            pad[ipad]->SetRightMargin(0);
            pad[ipad]->Modified();
            pad[ipad]->SetFillStyle(0);
        }
        else{
            pad[ipad]->SetLeftMargin(static_cast<Float_t>(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight)));//0.02/(0.02+0.455) = 0.0421
            pad[ipad]->SetRightMargin(static_cast<Float_t>(marginRight/(innerPadWidth+marginRight)));
            pad[ipad]->Modified();
            pad[ipad]->SetFillStyle(0);
        }
        //up 2 rows
        if(ipad < 3){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(3*innerPadHeight+marginBottomForYAxis+marginTop));
            pad[ipad]->SetTopMargin(marginTop/(3*innerPadHeight+marginTop));
        }
        else if(ipad < 6){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
            pad[ipad]->SetTopMargin(0.);
        }
        else if(ipad < 9){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(3*innerPadHeight+marginBottomForYAxis));
            pad[ipad]->SetTopMargin(0.);
        }
        else if(ipad < 12){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
            pad[ipad]->SetTopMargin(0.);
        }
        else if(ipad < 15){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(3*innerPadHeight+marginBottomForYAxis));
            pad[ipad]->SetTopMargin(0.);
        }
        else if(ipad < 18){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
            pad[ipad]->SetTopMargin(0.);
        }
        else if(ipad < 21){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(3*innerPadHeight+marginBottomForYAxis));
            pad[ipad]->SetTopMargin(0.);
        }
        else{
            pad[ipad]->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
            pad[ipad]->SetTopMargin(0.);
        }

        pad[ipad]->SetFillStyle(0);
        pad[ipad]->Modified();
        pad[ipad]->Update();
        FinalSpectrum->cd();
    }

    //ratio Y-axis title placeholder
    TH1D **placeholder = new TH1D*[static_cast<ULong_t>(24)];
    //other axis placeholders
    for(Int_t ipad = 0;ipad <24;ipad++){
            pad[ipad]->cd();
            if(ipad%3 ==0 || ipad%3 ==1)placeholder[ipad] = new TH1D(Form("placeholder2%d",ipad),"Central Values",static_cast<Int_t>(xAxisBins), xAxisl);
            else placeholder[ipad] = new TH1D(Form("placeholder2%d",ipad),"Central Values",static_cast<Int_t>(xAxisBins), xAxisr);
            //set ranges
            //left
            if(ipad/3 == 1|| ipad/3 == 3 || ipad/3 == 5 || ipad/3 == 7){
               placeholder[ipad]->SetMinimum(0.2);
               placeholder[ipad]->SetMaximum(1.8);
            }
            //right
            else{
               placeholder[ipad]->SetMinimum(0.1);
               placeholder[ipad]->SetMaximum(7.8);
            }

            if(ipad > 20)placeholder[ipad]->GetXaxis()->SetTitle("z_{||}");//parallel
            if(ipad%6==0)placeholder[ipad]->GetYaxis()->SetTitle("1/N_{jets}dN/dz_{||}");
            if(ipad%6==3){placeholder[ipad]->GetYaxis()->SetTitle("MC/Data");
              placeholder[ipad]->GetYaxis()->CenterTitle();
            }

            placeholder[ipad]->GetYaxis()->SetTitleFont(43);
            placeholder[ipad]->GetXaxis()->SetTitleFont(43);
            placeholder[ipad]->GetYaxis()->SetLabelFont(43);
            placeholder[ipad]->GetXaxis()->SetLabelFont(43);
            if(ipad/3 == 1|| ipad/3 == 3 || ipad/3 == 5 || ipad/3 == 7)placeholder[ipad]->GetYaxis()->SetTitleSize(12);
            else placeholder[ipad]->GetYaxis()->SetTitleSize(19);
            placeholder[ipad]->GetXaxis()->SetTitleSize(19);
            placeholder[ipad]->GetYaxis()->SetLabelSize(19);
            placeholder[ipad]->GetXaxis()->SetLabelSize(19);

            placeholder[ipad]->GetYaxis()->SetTitleOffset(4.f);//0.95*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
            placeholder[ipad]->GetXaxis()->SetTitleOffset(75*innerPadHeight/0.4);// not number before, no number was optimized before font 43 was introduced
            if(ipad%3!=0){
            placeholder[ipad]->GetYaxis()->SetTitleSize(0);
            placeholder[ipad]->GetYaxis()->SetLabelSize(0);
            }
            else placeholder[ipad]->GetYaxis()->SetDecimals();
            if(ipad/3 == 1|| ipad/3 == 3 || ipad/3 == 5 || ipad/3 == 7)placeholder[ipad]->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);
            //because root is so advanced it decided to put different ticks lenght for 2 particular plots - lets just force it to look same
           // if(ipad==4 || ipad==5) placeholder[ipad]->GetYaxis()->SetTickLength(0.04f);

            placeholder[ipad]->Draw();
        }

    return std::make_tuple(FinalSpectrum,pad, placeholder);

}

std::tuple<TCanvas*, TPad**, TH1D**> PrepareCanvas4x6Pad(UInt_t xAxisBins, Double_t *xAxisl, Double_t *xAxisr){
    style();
    //prepare main canvas
    TCanvas *FinalSpectrum = new TCanvas("FinalSpectrum4", "FinalSpectrum4",0,45,1250,1650);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    FinalSpectrum->SetHighLightColor(2);
    FinalSpectrum->Range(0,0,1,1);
    FinalSpectrum->SetFillColor(0);
    FinalSpectrum->SetBorderMode(0);
    FinalSpectrum->SetBorderSize(2);
    FinalSpectrum->SetFrameBorderMode(0);
    TPad **pad = new TPad*[static_cast<ULong_t>(24)];
    //0  1  2  3 norm
    //4  5  6  7 ratio
    //8  9  10 11     norm
    //12 13 14 15 ratio

    Double_t marginLeft=0.06;
    Double_t marginRight=0.02;
    Double_t innerPadWidth=(1-marginLeft-marginRight)/4.; //0.91/2 = 0.455
    Double_t marginLeftForXAxis=0.015;

    Double_t marginTop=0.04;
    Double_t marginBottom=0.06;
    Double_t innerPadHeight=(1-marginTop-marginBottom)/12.;
    Double_t marginBottomForYAxis=0.0;

    auto xbaseL = [innerPadWidth,marginLeft,marginLeftForXAxis](int m) ->Double_t {return (m*innerPadWidth+marginLeft-marginLeftForXAxis);};
    auto xbaseR = [innerPadWidth,marginLeft](int m) ->Double_t {return (m*innerPadWidth+marginLeft);};
    auto ybaseL = [innerPadHeight,marginBottom,marginBottomForYAxis](int m) ->Double_t {return (m*innerPadHeight+marginBottom-marginBottomForYAxis);};
    auto ybaseR = [innerPadHeight,marginBottom](int m) ->Double_t {return (m*innerPadHeight+marginBottom);};

    //pads for the plots
    for(Int_t ipad = 0;ipad < 24;ipad++){
        std::cout<<"padID: "<<ipad<<std::endl;
        if(ipad ==0) pad[ipad] = new TPad("pad0", "pad0",0,         ybaseL(9)  ,xbaseR(1) ,1);
        if(ipad ==1) pad[ipad] = new TPad("pad1", "pad1",xbaseL(1), ybaseL(9)  ,xbaseR(2) ,1);
        if(ipad ==2) pad[ipad] = new TPad("pad2", "pad2",xbaseL(2), ybaseL(9)  ,xbaseR(3)         ,1);
        if(ipad ==3) pad[ipad] = new TPad("pad3", "pad3",xbaseL(3), ybaseL(9)  ,1         ,1);
        if(ipad ==4) pad[ipad] = new TPad("pad4", "pad4",0 ,        ybaseL(8)  ,xbaseR(1)  ,ybaseR(9));
        if(ipad ==5) pad[ipad] = new TPad("pad5", "pad5",xbaseL(1), ybaseL(8)  ,xbaseR(2)  ,ybaseR(9));
        if(ipad ==6) pad[ipad] = new TPad("pad6", "pad6",xbaseL(2), ybaseL(8)   ,xbaseR(3) ,ybaseR(9));
        if(ipad ==7) pad[ipad] = new TPad("pad7", "pad7",xbaseL(3), ybaseL(8)   ,1         ,ybaseR(9));

        if(ipad ==8) pad[ipad] = new TPad("pad8", "pad8",0,            ybaseL(5)  ,xbaseR(1) ,ybaseR(8));
        if(ipad ==9) pad[ipad] = new TPad("pad9", "pad9",xbaseL(1),    ybaseL(5)  ,xbaseR(2) ,ybaseR(8));
        if(ipad ==10) pad[ipad] = new TPad("pad10", "pad10",xbaseL(2), ybaseL(5)  ,xbaseR(3) ,ybaseR(8));
        if(ipad ==11) pad[ipad] = new TPad("pad11", "pad11",xbaseL(3), ybaseL(5)  ,1         ,ybaseR(8));
        if(ipad ==12) pad[ipad] = new TPad("pad12", "pad12",0,         ybaseL(4)  ,xbaseR(1) ,ybaseR(5));
        if(ipad ==13) pad[ipad] = new TPad("pad13", "pad13",xbaseL(1), ybaseL(4)  ,xbaseR(2) ,ybaseR(5));
        if(ipad ==14) pad[ipad] = new TPad("pad14", "pad14",xbaseL(2), ybaseL(4)  ,xbaseR(3) ,ybaseR(5));
        if(ipad ==15) pad[ipad] = new TPad("pad15", "pad15",xbaseL(3), ybaseL(4)  ,1         ,ybaseR(5));

        if(ipad ==16) pad[ipad] = new TPad("pad16", "pad16",0 ,        ybaseL(1)  ,xbaseR(1)  ,ybaseR(4));
        if(ipad ==17) pad[ipad] = new TPad("pad17", "pad17",xbaseL(1), ybaseL(1)  ,xbaseR(2)  ,ybaseR(4));
        if(ipad ==18) pad[ipad] = new TPad("pad18", "pad18",xbaseL(2), ybaseL(1)  ,xbaseR(3)  ,ybaseR(4));
        if(ipad ==19) pad[ipad] = new TPad("pad19", "pad19",xbaseL(3), ybaseL(1)  ,1          ,ybaseR(4));
        if(ipad ==20) pad[ipad] = new TPad("pad20", "pad20",0,         0          ,xbaseR(1)  ,ybaseR(1));
        if(ipad ==21) pad[ipad] = new TPad("pad21", "pad21",xbaseL(1), 0          ,xbaseR(2)  ,ybaseR(1));
        if(ipad ==22) pad[ipad] = new TPad("pad22", "pad22",xbaseL(2) ,0          ,xbaseR(3)  ,ybaseR(1));
        if(ipad ==23) pad[ipad] = new TPad("pad23", "pad23",xbaseL(3), 0          ,1          ,ybaseR(1));


        pad[ipad]->Draw();
        pad[ipad]->cd();

        //left column
        if(ipad%4 ==0){
            pad[ipad]->SetLeftMargin(static_cast<Float_t>(marginLeft/(marginLeft+innerPadWidth))); //0.07/(0.07+0.455) = 0.13333
            pad[ipad]->SetRightMargin(0);
            pad[ipad]->Modified();
            pad[ipad]->SetFillStyle(0);
        }
        //right column
        else if (ipad%4 ==1 || ipad%4 ==2){
            pad[ipad]->SetLeftMargin(static_cast<Float_t>(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis)));//0.02/(0.02+0.455) = 0.0421
            pad[ipad]->SetRightMargin(0);
            pad[ipad]->Modified();
            pad[ipad]->SetFillStyle(0);
        }
        else if (ipad%4 ==3){
            pad[ipad]->SetLeftMargin(static_cast<Float_t>(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight)));//0.02/(0.02+0.455) = 0.0421
            pad[ipad]->SetRightMargin(static_cast<Float_t>(marginRight/(innerPadWidth+marginRight)));
            pad[ipad]->Modified();
            pad[ipad]->SetFillStyle(0);
        }
        //up 2 rows
        if(ipad < 4){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(3*innerPadHeight+marginBottomForYAxis+marginTop));
            pad[ipad]->SetTopMargin(marginTop/(3*innerPadHeight+marginTop));
        }
        else if(ipad < 8){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
            pad[ipad]->SetTopMargin(0.);
        }
        else if(ipad < 12){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(3*innerPadHeight+marginBottomForYAxis));
            pad[ipad]->SetTopMargin(0.);
        }
        else if(ipad < 16){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
            pad[ipad]->SetTopMargin(0.);
        }
        else if(ipad < 20){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(3*innerPadHeight+marginBottomForYAxis));
            pad[ipad]->SetTopMargin(0.);
        }
        else{
            pad[ipad]->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
            pad[ipad]->SetTopMargin(0.);
        }

        pad[ipad]->SetFillStyle(0);
        pad[ipad]->Modified();
        pad[ipad]->Update();
        FinalSpectrum->cd();
    }

    //ratio Y-axis title placeholder
    TH1D **placeholder = new TH1D*[static_cast<ULong_t>(24)];
    //other axis placeholders
    for(Int_t ipad = 0;ipad <24;ipad++){
            pad[ipad]->cd();
            if(ipad%4 ==0 || ipad%4 ==1|| ipad%4 ==2)placeholder[ipad] = new TH1D(Form("placeholder2%d",ipad),"Central Values",static_cast<Int_t>(xAxisBins), xAxisl);
            else placeholder[ipad] = new TH1D(Form("placeholder2%d",ipad),"Central Values",static_cast<Int_t>(xAxisBins), xAxisr);
            //set ranges
            //left
            if(ipad/4 == 1){
               placeholder[ipad]->SetMinimum(0.2);
               placeholder[ipad]->SetMaximum(1.8);
            }
            else if(ipad/4 == 3){
               placeholder[ipad]->SetMinimum(0.2);
               placeholder[ipad]->SetMaximum(1.9);
            }
            else if(ipad/4 == 5){
               placeholder[ipad]->SetMinimum(0.2);
               placeholder[ipad]->SetMaximum(2.2);
            }
            //right
            else if(ipad == 0 ||ipad == 1 ||ipad == 2 ||ipad == 3){
               placeholder[ipad]->SetMinimum(0.1);
               placeholder[ipad]->SetMaximum(7.8);
            }
            else if(ipad == 8 ||ipad == 9 ||ipad == 10 ||ipad == 11){
               placeholder[ipad]->SetMinimum(0.1);
               placeholder[ipad]->SetMaximum(5.3);
            }
            else if(ipad == 16 ||ipad == 17 ||ipad == 18 ||ipad == 19){
               placeholder[ipad]->SetMinimum(0.1);
               placeholder[ipad]->SetMaximum(5.3);
            }

            if(ipad > 19)placeholder[ipad]->GetXaxis()->SetTitle("z_{||}");
            if(ipad%8==0)placeholder[ipad]->GetYaxis()->SetTitle("1/N_{jets}dN/dz_{||}");//parallel
            if(ipad%8==4){placeholder[ipad]->GetYaxis()->SetTitle("MC/Data");
              placeholder[ipad]->GetYaxis()->CenterTitle();
            }

            placeholder[ipad]->GetYaxis()->SetTitleFont(43);
            placeholder[ipad]->GetXaxis()->SetTitleFont(43);
            placeholder[ipad]->GetYaxis()->SetLabelFont(43);
            placeholder[ipad]->GetXaxis()->SetLabelFont(43);
            if(ipad/4 == 1|| ipad/4 == 3 || ipad/4 == 5)placeholder[ipad]->GetYaxis()->SetTitleSize(12);
            else placeholder[ipad]->GetYaxis()->SetTitleSize(19);
            placeholder[ipad]->GetXaxis()->SetTitleSize(19);
            placeholder[ipad]->GetYaxis()->SetLabelSize(19);
            placeholder[ipad]->GetXaxis()->SetLabelSize(19);

            placeholder[ipad]->GetYaxis()->SetTitleOffset(4.f);//0.95*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
            if(ipad%8==4)placeholder[ipad]->GetYaxis()->SetTitleOffset(6.f);
            placeholder[ipad]->GetXaxis()->SetTitleOffset(45*innerPadHeight/0.4);// not number before, no number was optimized before font 43 was introduced
            if(ipad%4!=0){
            placeholder[ipad]->GetYaxis()->SetTitleSize(0);
            placeholder[ipad]->GetYaxis()->SetLabelSize(0);
            }
            else placeholder[ipad]->GetYaxis()->SetDecimals();
            if(ipad/4 == 1|| ipad/4 == 3 || ipad/4 == 5)placeholder[ipad]->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);
            //because root is so advanced it decided to put different ticks lenght for 2 particular plots - lets just force it to look same
            if(ipad/4 == 5) placeholder[ipad]->GetYaxis()->SetTickLength(0.055f);

            placeholder[ipad]->Draw();
        }

    return std::make_tuple(FinalSpectrum,pad, placeholder);

}


std::tuple<TCanvas*, TPad**, TH1D**> PrepareCanvas4Pad(UInt_t xAxisBins, Double_t *xAxis){
    style();
    //prepare main canvas
    TCanvas *FinalSpectrum = new TCanvas("FinalSpectrum2", "FinalSpectrum2",0,45,700,700);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    FinalSpectrum->SetHighLightColor(2);
    FinalSpectrum->Range(0,0,1,1);
    FinalSpectrum->SetFillColor(0);
    FinalSpectrum->SetBorderMode(0);
    FinalSpectrum->SetBorderSize(2);
    FinalSpectrum->SetFrameBorderMode(0);
    TPad **pad = new TPad*[static_cast<ULong_t>(10)];
    //0  1
    //2  3
    //4  5

    Double_t marginLeft=0.175;
    Double_t marginRight=0.02;
    Double_t innerPadWidth=(1-marginLeft-marginRight); //0.91/2 = 0.455
    Double_t marginLeftForXAxis=0.02;

    Double_t marginTop=0.04;
    Double_t marginBottom=0.08;
    Double_t innerPadHeight=(1-marginTop-marginBottom)/6.;
    Double_t marginBottomForYAxis=0.0;

    for(Int_t ipad = 0;ipad <4;ipad++){
        Double_t shift = 0.002;
        std::cout<<"padID: "<<ipad<<std::endl;
        if(ipad ==0) pad[ipad] = new TPad("pad0", "pad0",0,     3*innerPadHeight+marginBottom-marginBottomForYAxis  ,1 ,1);
        if(ipad ==1) pad[ipad] = new TPad("pad1", "pad1",0,     2*innerPadHeight+marginBottom-marginBottomForYAxis ,1 ,3*innerPadHeight+marginBottom);
        if(ipad ==2) pad[ipad] = new TPad("pad2", "pad2",0,     innerPadHeight+marginBottom-marginBottomForYAxis  ,1 ,2*innerPadHeight+marginBottom);
        if(ipad ==3) pad[ipad] = new TPad("pad3", "pad3",0,     0.0    ,1 , innerPadHeight+marginBottom);
        pad[ipad]->Draw();
        pad[ipad]->cd();

        //left column
        pad[ipad]->SetLeftMargin(static_cast<Float_t>(marginLeft/(marginLeft+innerPadWidth))); //0.07/(0.07+0.455) = 0.13333
        pad[ipad]->SetRightMargin(static_cast<Float_t>(marginRight/(innerPadWidth+marginRight)));
        if(ipad == 0) {
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(3*innerPadHeight+marginBottomForYAxis+marginTop));
            pad[ipad]->SetTopMargin(marginTop/(3*innerPadHeight+marginTop));
        }
        else if (ipad == 1 || ipad ==2){
            pad[ipad]->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
            pad[ipad]->SetTopMargin(0.);
        }
        else{
            pad[ipad]->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
            pad[ipad]->SetTopMargin(0.);
        }

        FinalSpectrum->cd();
    }

    //TH1D *placeholder[8];
    TH1D **placeholder = new TH1D*[static_cast<ULong_t>(10)];

    for(Int_t ipad = 0;ipad <4;ipad++){
       // if(ipad==6 || ipad==7) continue;
        pad[ipad]->cd();
        //if(ipad==0)pad[ipad]->cd()->SetLogy();
        placeholder[ipad] = new TH1D(Form("placeholder2%d",ipad),"Central Values",static_cast<Int_t>(xAxisBins), xAxis);
        if(ipad<1){
            placeholder[ipad]->SetMinimum(1);
            placeholder[ipad]->SetMaximum(10);
        }
        else{
            placeholder[ipad]->SetMinimum(0.1);
            placeholder[ipad]->SetMaximum(1.99);
        }
        if(ipad==3)placeholder[ipad]->GetXaxis()->SetTitle("#it{p}_{T,ch. jet} (GeV/#it{c})");
        if(ipad==0)placeholder[ipad]->GetYaxis()->SetTitle("#frac{d^{2}#it{#sigma^{13TeV}}}{d#it{p}_{T,ch. jet}d#it{#eta_{T,ch. jet}}}/#frac{d^{2}#it{#sigma^{5TeV}}}{d#it{p}_{T,ch. jet}d#it{#eta_{T,ch. jet}}}");
        if(ipad==2)placeholder[ipad]->GetYaxis()->SetTitle("MC/Data");
        placeholder[ipad]->GetYaxis()->SetTitleFont(43);
        placeholder[ipad]->GetXaxis()->SetTitleFont(43);
        placeholder[ipad]->GetYaxis()->SetLabelFont(43);
        placeholder[ipad]->GetXaxis()->SetLabelFont(43);
        placeholder[ipad]->GetYaxis()->SetTitleSize(19);
        placeholder[ipad]->GetXaxis()->SetTitleSize(19);
        placeholder[ipad]->GetYaxis()->SetLabelSize(19);
        placeholder[ipad]->GetXaxis()->SetLabelSize(19);

        placeholder[ipad]->GetYaxis()->SetTitleOffset(2.2f);//0.95*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
        placeholder[ipad]->GetXaxis()->SetTitleOffset(7*innerPadHeight/0.4);// not number before, no number was optimized before font 43 was introduced

        if(ipad>0)placeholder[ipad]->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);
        //because root is so advanced it decided to put different ticks lenght for 2 particular plots - lets just force it to look same
        if(ipad==1 || ipad ==2) placeholder[ipad]->GetYaxis()->SetTickLength(0.0275f);
        if(ipad==3) placeholder[ipad]->GetYaxis()->SetTickLength(0.043f);

        placeholder[ipad]->Draw();
    }

    return std::make_tuple(FinalSpectrum,pad, placeholder);

}


std::tuple<TCanvas*, TPad*, TPad*,  TPad*, TH1D*, TH1D*, TH1D*> PrepareCanvas(UInt_t xAxisBins, Double_t *xAxis){
    style();
    //prepare main canvas
    TCanvas *FinalSpectrum = new TCanvas("FinalSpectrum2", "FinalSpectrum2",0,45,700,700);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    FinalSpectrum->SetHighLightColor(2);
    FinalSpectrum->Range(0,0,1,1);
    FinalSpectrum->SetFillColor(0);
    FinalSpectrum->SetBorderMode(0);
    FinalSpectrum->SetBorderSize(2);
    FinalSpectrum->SetFrameBorderMode(0);

    //Set primitives in upper pad
    TPad *pad_top = new TPad("pad_top", "pad_top",0,0.5,1,1);
    pad_top->Draw();
    pad_top->cd();
    pad_top->Range(-1.986821e-07,-4.69897,33.33333,0.3499945);
    pad_top->SetFillColor(0);
    pad_top->SetBorderMode(0);
    pad_top->SetBorderSize(2);
//    if(zBin==0)pad_top->SetLogy();
    pad_top->SetTickx(1);
    pad_top->SetTicky(1);
    pad_top->SetLeftMargin(0.18f);
    pad_top->SetBottomMargin(0);
    pad_top->SetFrameBorderMode(0);
    pad_top->SetFrameBorderMode(0);

    FinalSpectrum->cd();



    //Set primitives in bottom pad
    TPad *pad_bottom = new TPad("pad_bottom", "pad_bottom",0,0.275,1,0.5);
    pad_bottom->Draw();
    pad_bottom->cd();
    pad_bottom->Range(-1.986821e-07,-0.9209589,33.33333,2.49);
    pad_bottom->SetFillColor(0);
    pad_bottom->SetBorderMode(0);
    pad_bottom->SetBorderSize(2);
    pad_bottom->SetGridy();
    pad_bottom->SetTickx(1);
    pad_bottom->SetTicky(1);
    pad_bottom->SetLeftMargin(0.18f);
    pad_bottom->SetTopMargin(0);
    pad_bottom->SetBottomMargin(0);
    pad_bottom->SetFrameBorderMode(0);
    pad_bottom->SetFrameBorderMode(0);

    FinalSpectrum->cd();

    //Set primitives in bottom pad
    TPad *pad_bottom2 = new TPad("pad_bottom2", "pad_bottom2",0,0,1,0.275);
    pad_bottom2->Draw();
    pad_bottom2->cd();
    pad_bottom2->Range(-1.986821e-07,-0.9209589,33.33333,2.49);
    pad_bottom2->SetFillColor(0);
    pad_bottom2->SetBorderMode(0);
    pad_bottom2->SetBorderSize(2);
    pad_bottom2->SetGridy();
    pad_bottom2->SetTickx(1);
    pad_bottom2->SetTicky(1);
    pad_bottom2->SetLeftMargin(0.18f);
    pad_bottom2->SetTopMargin(0);
    pad_bottom2->SetBottomMargin(0.27f);
    pad_bottom2->SetFrameBorderMode(0);
    pad_bottom2->SetFrameBorderMode(0);

    pad_top->cd();
    TH1D *hEmpty_up = new TH1D("hEmpty_up","Central Values",static_cast<Int_t>(xAxisBins), xAxis);
    hEmpty_up->SetMinimum(plotRanges[2]);
    hEmpty_up->SetMaximum(plotRanges[3]);
    hEmpty_up->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
    hEmpty_up->GetXaxis()->SetLabelFont(43);
    hEmpty_up->GetXaxis()->SetLabelSize(0.035f);
    hEmpty_up->GetXaxis()->SetTitleSize(0.035f);
    hEmpty_up->GetXaxis()->SetTitleFont(42);
    hEmpty_up->GetYaxis()->SetRangeUser(plotRanges[2],plotRanges[3]);
  //  hEmpty_up->GetYaxis()->SetTitle("Ratio of cross-sections");
    hEmpty_up->GetYaxis()->SetTitle("Ratio of probability densities");
 //   if(zBin ==2 && !pdf)hEmpty_up->GetYaxis()->SetTitle("#frac{d^{2}#it{#sigma}}{d#it{p}_{T}d#it{#eta}} (mb#upointGeV^{#minus1}#upoint#it{c})");
 //   else if(!pdf)hEmpty_up->GetYaxis()->SetTitle("#frac{d^{2}#it{#sigma}}{d#it{z}_{#parallel}d#it{#eta}} (mb)");
 //   else hEmpty_up->GetYaxis()->SetTitle("Probability density");
    hEmpty_up->GetYaxis()->SetLabelFont(43);
    hEmpty_up->GetYaxis()->SetLabelSize(19);
    hEmpty_up->GetYaxis()->SetTitleSize(26);
    hEmpty_up->GetYaxis()->SetLabelOffset(0.015f);
    hEmpty_up->GetYaxis()->SetTitleOffset(2.f);
    hEmpty_up->GetYaxis()->SetTitleFont(43);
    hEmpty_up->GetYaxis()->SetDecimals();
    hEmpty_up->Draw("axis");
    pad_bottom->cd();
    TH1D *hEmpty_bottom = new TH1D("hEmpty_bottom","Central Values",static_cast<Int_t>(xAxisBins), xAxis);
    hEmpty_bottom->SetMinimum(0.4);
    hEmpty_bottom->SetMaximum(plotRanges[1]);
   // if(zBin ==0)hEmpty_bottom->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
   // else hEmpty_bottom->GetXaxis()->SetTitle("z_{#parallel}");
    hEmpty_bottom->GetXaxis()->SetLabelFont(43);
    hEmpty_bottom->GetXaxis()->SetLabelSize(19);
    hEmpty_bottom->GetXaxis()->SetTitleSize(26);
    hEmpty_bottom->GetXaxis()->SetTitleOffset(3.f);
    hEmpty_bottom->GetXaxis()->SetTitleFont(43);
    hEmpty_bottom->GetYaxis()->SetTitle("MC/Data");
    hEmpty_bottom->GetYaxis()->SetNdivisions(509);
    hEmpty_bottom->GetYaxis()->CenterTitle();
    hEmpty_bottom->GetYaxis()->SetDecimals();
    hEmpty_bottom->GetYaxis()->SetLabelOffset(0.015f);
    hEmpty_bottom->GetXaxis()->SetLabelOffset(0.02f);
    hEmpty_bottom->GetYaxis()->SetLabelFont(43);
    hEmpty_bottom->GetYaxis()->SetLabelSize(19);
    hEmpty_bottom->GetYaxis()->SetTitleSize(26);
    hEmpty_bottom->GetYaxis()->SetTitleOffset(2.f);
    hEmpty_bottom->GetYaxis()->SetTitleFont(43);
    hEmpty_bottom->Draw("axis");

    pad_bottom2->cd();
    TH1D *hEmpty_bottom2 = new TH1D("hEmpty_bottom2","Central Values",static_cast<Int_t>(xAxisBins), xAxis);
    hEmpty_bottom2->SetMinimum(plotRanges[0]);
    hEmpty_bottom2->SetMaximum(1.6);
    hEmpty_bottom2->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
    //else hEmpty_bottom->GetXaxis()->SetTitle("z_{#parallel}");
    hEmpty_bottom2->GetXaxis()->SetLabelFont(43);
    hEmpty_bottom2->GetXaxis()->SetLabelSize(19);
    hEmpty_bottom2->GetXaxis()->SetTitleSize(26);
    hEmpty_bottom2->GetXaxis()->SetTitleOffset(3.f);
    hEmpty_bottom2->GetXaxis()->SetTitleFont(43);
    hEmpty_bottom2->GetYaxis()->SetTitle("MC/Data");
    hEmpty_bottom2->GetYaxis()->SetNdivisions(509);
    hEmpty_bottom2->GetYaxis()->CenterTitle();
    hEmpty_bottom2->GetYaxis()->SetDecimals();
    hEmpty_bottom2->GetYaxis()->SetLabelOffset(0.015f);
    hEmpty_bottom2->GetXaxis()->SetLabelOffset(0.02f);
    hEmpty_bottom2->GetYaxis()->SetLabelFont(43);
    hEmpty_bottom2->GetYaxis()->SetLabelSize(19);
    hEmpty_bottom2->GetYaxis()->SetTitleSize(26);
    hEmpty_bottom2->GetYaxis()->SetTitleOffset(2.f);
    hEmpty_bottom2->GetYaxis()->SetTitleFont(43);
    hEmpty_bottom2->Draw("axis");

    return std::make_tuple(FinalSpectrum,pad_top, pad_bottom, pad_bottom2, hEmpty_up,hEmpty_bottom,hEmpty_bottom2);

}


std::tuple<TGraphAsymmErrors*, TH1D*,TGraphAsymmErrors*, TH1D*> GetData(TString dataFile, TString histBase, Double_t dataScaling, UInt_t nBins, Double_t *xBins, Double_t *systUncD_down, Double_t *systUncD_up){
    std::cout<<"get DATA"<<std::endl;
    TH1D *htmp = dynamic_cast<TH1D*>(GetInputHist(dataFile, histBase));
    std::cout<<"getting "<<dataFile<<" "<<histBase<<std::endl;
    if(!htmp)std::cout<<"not histo"<<std::endl;
    TH1D *hData_binned = dynamic_cast<TH1D*>(htmp->Rebin(static_cast<Int_t>(nBins),"hData_binned", xBins));
    //Double_t data_int = hData_binned->Integral();

    if(pdf){
      hData_binned->Scale(1./hData_binned->Integral());
      hData_binned->Scale(1,"width");
    }else{
        hData_binned->Scale(1,"width");
        hData_binned->Scale(dataScaling);
        hData_binned->Scale(1./dy);
    }

    hData_binned->SetTitle("");
    hData_binned->SetMaximum(hData_binned->GetMaximum()*2);
    //hData_binned->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}d#it{#eta} (mb)");
    if(!pdf)hData_binned->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#it{#eta}} mb (GeV/#it{c})^{-1}");
    else hData_binned->GetYaxis()->SetTitle("Probability density");

    Double_t *sysuncAbs_down = new Double_t[nBins];
    Double_t *sysuncAbs_up = new Double_t[nBins];
    Double_t *statunc = new Double_t[nBins];
    Double_t *value = new Double_t[nBins];
    Double_t *xval = new Double_t[nBins];
    Double_t *xvalwidth = new Double_t[nBins];
    for(UInt_t j=0; j<nBins; j++){
          xval[j] = (xBins[j]+xBins[j+1]) / 2.;
          xvalwidth[j] = (xBins[j+1]-xBins[j]) / 2.;
          value[j] = hData_binned->GetBinContent(hData_binned->GetXaxis()->FindBin(xval[j]));
          Double_t error = hData_binned->GetBinError(hData_binned->GetXaxis()->FindBin(xval[j]));
          std::cout<<j<<" "<<xval[j]<<" "<<xvalwidth[j]<<" "<<value[j]<<" "<<error<<std::endl;
          sysuncAbs_down[j] = value[j] * systUncD_down[j];
          sysuncAbs_up[j] = value[j] * systUncD_up[j];
          if(value[j]>1E-9)statunc[j] = error/ value[j] *100;
    }
    TGraphAsymmErrors *grsysl = new TGraphAsymmErrors(static_cast<Int_t>(nBins),xval,value,xvalwidth,xvalwidth,sysuncAbs_down,sysuncAbs_up);
    grsysl->SetName("haeData_binned_syst");

    TH1D *hData_binned_ratio = dynamic_cast<TH1D*>(hData_binned->Clone("hData_binned_ratio"));
      double *sysuncRatio_up = new double[nBins];
      double *sysuncRatio_down = new double[nBins];
      double *valRatio = new double[nBins];
      for(UInt_t j=0; j<nBins; j++){
          xval[j] = (xBins[j]+xBins[j+1]) / 2.;
          xvalwidth[j] = (xBins[j+1]-xBins[j]) / 2.;
              double valPred = hData_binned->GetBinContent(hData_binned->GetXaxis()->FindBin(xval[j]));
              if(valPred<1E-9)valPred=1000;
              if(valPred>1E-9)valRatio[j] = 1.0;
              double err = 0;
              if(valPred>1E-9) err=hData_binned->GetBinError(hData_binned->GetXaxis()->FindBin(xval[j])) / valPred;

              sysuncRatio_down[j] = valRatio[j] * systUncD_down[j];
              sysuncRatio_up[j] = valRatio[j] * systUncD_up[j];

              hData_binned_ratio->SetBinContent(hData_binned_ratio->GetXaxis()->FindBin(xval[j]),valRatio[j]);
              hData_binned_ratio->SetBinError(hData_binned_ratio->GetXaxis()->FindBin(xval[j]),err);
      }
      TGraphAsymmErrors *grsysRatio = new TGraphAsymmErrors(static_cast<Int_t>(nBins),xval,valRatio,xvalwidth,xvalwidth,sysuncRatio_down,sysuncRatio_up);
      grsysRatio->SetName("haeData_binned_syst_ratio");
    return std::make_tuple(grsysl, hData_binned, grsysRatio,hData_binned_ratio);
}

TGraphAsymmErrors* Divide(TGraphAsymmErrors *hist1,TH1D *hist2){
    TString name1 = hist1->GetName();
    TString name2 = hist2->GetName();

    UInt_t nBin = hist1->GetN();
    Double_t *xval = new Double_t[nBin];
    Double_t *xvalwidth_down = new Double_t[nBin];
    Double_t *xvalwidth_up = new Double_t[nBin];
    Double_t *value = new Double_t[nBin];
    Double_t *valuetheoryerrup = new Double_t[nBin];
    Double_t *valuetheoryerrdown = new Double_t[nBin];
    for(UInt_t bin=0; bin<nBin; bin++){
        Double_t x1 = 0,y1=0;
        hist1->GetPoint(bin,x1,y1);
        xval[bin] = x1;
        xvalwidth_down[bin] = hist1->GetErrorXlow(bin);
        xvalwidth_up[bin] = hist1->GetErrorXhigh(bin);
        value[bin] = y1/hist2->GetBinContent(bin+1);
        std::cout<<"DIVIDE: "<<x1<<" "<<y1<<" "<<hist2->GetBinCenter(bin+1)<<" "<<hist2->GetBinContent(bin+1)<<" "<<value[bin]<<std::endl;
        valuetheoryerrup[bin] = hist1->GetErrorYhigh(bin)/hist2->GetBinContent(bin+1);
        valuetheoryerrdown[bin] = hist1->GetErrorYlow(bin)/hist2->GetBinContent(bin+1);
    }
    TGraphAsymmErrors *grsystheory = new TGraphAsymmErrors(static_cast<Int_t>(nBin),xval,value,xvalwidth_down,xvalwidth_up,valuetheoryerrdown,valuetheoryerrup);
    grsystheory->SetName(name1+"_"+name2+"_ratio_sys");
    return grsystheory;
}


std::tuple<TGraphAsymmErrors*, TGraphAsymmErrors*, TH1D**> GetSim(TString hash, TString dataFile,TString simname, Int_t nVar){
    TH1D **hPrompt = nullptr;
    if(nVar>0)hPrompt= new TH1D*[static_cast<ULong_t>(nVar)];
    TFile *jetPtFile = new TFile(dataFile,"read");
    TString name = "";
     TGraphAsymmErrors *h1 = dynamic_cast<TGraphAsymmErrors*>(jetPtFile->Get(name+"haesim"+simname));
     TGraphAsymmErrors *h2 = dynamic_cast<TGraphAsymmErrors*>(jetPtFile->Get(name+"haesim"+simname+"_ratio"));
     h1->SetName(name+"haesim"+simname+hash);
     h2->SetName(name+"haesim"+simname+hash+"_ratio");
     for (int nr=0; nr<nVar; nr++){
         hPrompt[nr]= dynamic_cast<TH1D*>(jetPtFile->Get(name+"hsim"+simname+Form("_var%d",nr+1)));
         std::cout<<"aaaaaaaaa"<<name+hash+"hsim"+simname+Form("_var%d",nr+1)<<std::endl;
         hPrompt[nr]->SetName(name+hash+"hsim"+simname+Form("_var%d",nr+1));
     }
     return std::make_tuple(h1, h2,hPrompt);
}

std::tuple<TGraphAsymmErrors*, TH1D*, TH1D*, TH1D*> GetDataSimRatio(TString simname, TH1D *data_cent, TH1D *sim_cent, TH1D *sim_up, TH1D *sim_down, UInt_t nBins, Double_t *xBins){
    TString name = sim_cent->GetName();
    TString name_up = sim_up->GetName();
    TString name_down = sim_down->GetName();
    TH1D *hPrompt_central_binned_ratio = dynamic_cast<TH1D*>(sim_cent->Clone(name+"_ratio"));
    hPrompt_central_binned_ratio->Divide(data_cent);
    TH1D *hPrompt_down_ratio = dynamic_cast<TH1D*>(sim_down->Clone(name_down+"_ratio"));
    TH1D *hPrompt_up_ratio = dynamic_cast<TH1D*>(sim_up->Clone(name_up+"_ratio"));
    hPrompt_up_ratio->Divide(data_cent);
    hPrompt_down_ratio->Divide(data_cent);

    Double_t *ptvaltheoryratio = new Double_t[nBins];
    Double_t *ptvalunctheoryratio = new Double_t[nBins];
    Double_t *valuetheoryratio = new Double_t[nBins];
    Double_t *valuetheoryerrupratio = new Double_t[nBins];
    Double_t *valuetheoryerrdownratio = new  Double_t[nBins];
    for(UInt_t j=0; j<nBins; j++){
          ptvaltheoryratio[j] = (xBins[j]+xBins[j+1]) / 2.;
          ptvalunctheoryratio[j] = (xBins[j+1]-xBins[j]) / 2.;
          valuetheoryratio[j] = hPrompt_central_binned_ratio->GetBinContent(hPrompt_central_binned_ratio->GetXaxis()->FindBin(ptvaltheoryratio[j]));
          valuetheoryerrupratio[j] = hPrompt_up_ratio->GetBinContent(hPrompt_up_ratio->GetXaxis()->FindBin(ptvaltheoryratio[j])) - valuetheoryratio[j];
          valuetheoryerrdownratio[j] = valuetheoryratio[j] - hPrompt_down_ratio->GetBinContent(hPrompt_down_ratio->GetXaxis()->FindBin(ptvaltheoryratio[j]));
    }
    TGraphAsymmErrors *grsystheoryratio = new TGraphAsymmErrors(static_cast<Int_t>(nBins),ptvaltheoryratio,valuetheoryratio,ptvalunctheoryratio,ptvalunctheoryratio,valuetheoryerrdownratio,valuetheoryerrupratio);
    TString tmp = "hae";
    tmp+=simname;
    grsystheoryratio->SetName(tmp+"_ratio");
    return std::make_tuple(grsystheoryratio, hPrompt_central_binned_ratio, hPrompt_up_ratio,hPrompt_down_ratio);
}

std::tuple<TH1D*, TGraphAsymmErrors*> Rebin(UInt_t xAxisBins_rebin, Double_t *xAxis_rebin, TH1D* histo1, TGraphAsymmErrors* histo2){
    TString nametmp = "";
    std::cout<<"REBIN"<<std::endl;
    TH1D *rebin_tmp = nullptr;
    if(histo1){
        rebin_tmp = new TH1D(nametmp+histo1->GetName()+"_rebinned",nametmp+histo1->GetName()+"_rebinned",xAxisBins_rebin,xAxis_rebin);
        for(Int_t bin = 1; bin <= rebin_tmp->GetNbinsX(); bin++){
            if(bin <= 3){
                rebin_tmp->SetBinContent(bin,histo1->GetBinContent(bin));
                rebin_tmp->SetBinError(bin,histo1->GetBinError(bin));
            }
            else if(bin == 4){
                Double_t w1 = 1./histo1->GetBinContent(bin);
                Double_t w2 = 1./histo1->GetBinContent(bin+1);
                Double_t x1 = histo1->GetBinContent(bin);
                Double_t x2 = histo1->GetBinContent(bin+1);
                Double_t e1 = histo1->GetBinError(bin);
                Double_t e2 = histo1->GetBinError(bin+1);
                Double_t val = (w1*x1+w2*x2)/(w1+w2);
                Double_t err = (w1*e1+w2*e2)/(w1+w2);
                rebin_tmp->SetBinContent(bin,val);
                rebin_tmp->SetBinError(bin,err);
            }
            else if(bin >= 5){
                rebin_tmp->SetBinContent(bin,histo1->GetBinContent(bin+1));
                rebin_tmp->SetBinError(bin,histo1->GetBinError(bin+1));
            }
        }
    }

    Double_t *ptvaltheory = new Double_t[xAxisBins_rebin];
    Double_t *ptvalunctheory = new Double_t[xAxisBins_rebin];
    Double_t *valuetheory = new Double_t[xAxisBins_rebin];
    Double_t *valuetheoryerrup = new Double_t[xAxisBins_rebin];
    Double_t *valuetheoryerrdown = new  Double_t[xAxisBins_rebin];
    Double_t x1,y1;
    Double_t x2,y2,e1l,e1h,e2l,e2h,w1,w2;
    for(UInt_t j=0; j<xAxisBins_rebin; j++){
        if(j < 3){
            histo2->GetPoint(j,x1,y1);
            ptvaltheory[j] = (xAxis_rebin[j]+xAxis_rebin[j+1]) / 2.;
            ptvalunctheory[j] = (xAxis_rebin[j+1]-xAxis_rebin[j]) / 2.;
            valuetheory[j] = y1;
            valuetheoryerrup[j] = histo2->GetErrorYhigh(j);
            valuetheoryerrdown[j] = histo2->GetErrorYlow(j);
        }
        else if(j == 3){
            histo2->GetPoint(j,x1,y1);
            histo2->GetPoint(j+1,x2,y2);
            w1 = 1./y1;
            w2 = 1./y2;
            e1l = histo2->GetErrorYlow(j);
            e2l = histo2->GetErrorYlow(j+1);
            e1h = histo2->GetErrorYhigh(j);
            e2h = histo2->GetErrorYhigh(j+1);
            ptvaltheory[j] = (xAxis_rebin[j]+xAxis_rebin[j+1]) / 2.;
            ptvalunctheory[j] = (xAxis_rebin[j+1]-xAxis_rebin[j]) / 2.;
            valuetheory[j] = (w1*y1+w2*y2)/(w1+w2);;
            valuetheoryerrup[j] = (w1*e1h+w2*e2h)/(w1+w2);
            valuetheoryerrdown[j] = (w1*e1l+w2*e2l)/(w1+w2);

        }
        else if(j >= 4){
            histo2->GetPoint(j+1,x1,y1);
            ptvaltheory[j] = (xAxis_rebin[j]+xAxis_rebin[j+1]) / 2.;
            ptvalunctheory[j] = (xAxis_rebin[j+1]-xAxis_rebin[j]) / 2.;
            valuetheory[j] = y1;
            valuetheoryerrup[j] = histo2->GetErrorYhigh(j+1);
            valuetheoryerrdown[j] = histo2->GetErrorYlow(j+1);
        }

    }
    TGraphAsymmErrors *grsystheory = new TGraphAsymmErrors(static_cast<Int_t>(xAxisBins_rebin),ptvaltheory,valuetheory,ptvalunctheory,ptvalunctheory,valuetheoryerrdown,valuetheoryerrup);
    TString tmp = "hae";
    tmp+=histo2->GetName();
    grsystheory->SetName(tmp+"_rebin");
        //std::cout<<bin<<" "<<rebin_tmp->GetBinCenter(bin)<<" "<<rebin_tmp->GetBinContent(bin)<<" "<<histo1->GetBinCenter(bin)<<" "<<histo1->GetBinContent(bin)<<std::endl;
    return std::make_tuple(rebin_tmp, grsystheory);
}
/*
TH1D* GetRatio(TH1D *h1, TH1D *h2){
    TString name1 = h1->GetName();
    TString name2 = h1->GetName();
    TH1D *hRatio = dynamic_cast<TH1D*>(h1->Clone(name1+"_"+name2+"_ratio"));
    hRatio->Divide(h2);
    return hRatio;
}

TH1D* GetRatio(TH1D *h1, TH1D *h2, Double_t shift){
    TString name1 = h1->GetName();
    TString name2 = h1->GetName();
    TH1D *hRatio = dynamic_cast<TH1D*>(h1->Clone(name1+"_"+name2+"_ratio"));
    hRatio->Divide(h2);
    for(Int_t bin = 1; bin <= hRatio->GetNbinsX();bin++){
        hRatio->SetBinContent(bin,hRatio->GetBinContent(bin)+shift);
    }

    return hRatio;
}
*/
/*
std::tuple<TH1D*, TH1D*, TH1D*> GetRatio(TH1D *sim_cent, TH1D *sim_up, TH1D *sim_down,TH1D *sim2_cent, TH1D *sim2_up, TH1D *sim2_down){
    std::cout<<"A"<<std::endl;
    TString name = sim_cent->GetName();
    TString name_up = sim_up->GetName();
    TString name_down = sim_down->GetName();
    TString name2 = sim2_cent->GetName();
    TString name2_up = sim2_up->GetName();
    TString name2_down = sim2_down->GetName();
    std::cout<<"A"<<std::endl;
    TH1D *hPrompt_central_binned_ratio = dynamic_cast<TH1D*>(sim_cent->Clone(name+name2+"_ratio"));
    TH1D *hPrompt_down_ratio = dynamic_cast<TH1D*>(sim_down->Clone(name_down+name2_down+"_ratio"));
    TH1D *hPrompt_up_ratio = dynamic_cast<TH1D*>(sim_up->Clone(name_up+name2_up+"_ratio"));
    std::cout<<"A"<<std::endl;
    hPrompt_central_binned_ratio->Divide(sim2_cent);
    hPrompt_down_ratio->Divide(sim2_down);
    hPrompt_up_ratio->Divide(sim2_up);
    std::cout<<"A"<<std::endl;
    return std::make_tuple(hPrompt_central_binned_ratio,hPrompt_up_ratio,hPrompt_down_ratio);

}
*/
TH1D* GetUpSys(TH1D *central, TH1D **hh, Int_t nFiles){
    double max = 0, maxerr = 0;
    TString name = central->GetName();
    TH1D *hh_up = dynamic_cast<TH1D*>(central->Clone(name +"_up"));

    for(int iBin=1; iBin<central->GetNbinsX()+1; iBin++ ){
      //  std::cout<<"iBin: "<<iBin<<std::endl;
        max = central->GetBinContent(iBin);
        for(Int_t iFile=1; iFile < nFiles; iFile++){
           // std::cout<<"iFile: "<<iFile<<" nfiles "<<nFiles<<" binc "<<hh[iFile]->GetBinContent(iBin)<<std::endl;
            if(hh[iFile]->GetBinContent(iBin) > max){
                    max = hh[iFile]->GetBinContent(iBin);
                    maxerr = hh[iFile]->GetBinError(iBin);
            }
        }
        hh_up->SetBinContent(iBin,max);
        hh_up->SetBinError(iBin,0);
    }

    return hh_up;
}

TH1D* GetDownSys(TH1D *central, TH1D **hh, Int_t nFiles){
    double max = 0, maxerr = 0;
    TString name = central->GetName();
    TH1D *hh_down = dynamic_cast<TH1D*>(central->Clone(name +"_down"));

    for(int iBin=1; iBin<central->GetNbinsX()+1; iBin++ ){
        max = central->GetBinContent(iBin);
        for(Int_t iFile=1; iFile < nFiles; iFile++){
            if(hh[iFile]->GetBinContent(iBin) < max){
                    max = hh[iFile]->GetBinContent(iBin);
                    maxerr = hh[iFile]->GetBinError(iBin);
            }
        }
        hh_down->SetBinContent(iBin,max);
        hh_down->SetBinError(iBin,0);
    }

    return hh_down;
}

TH1* GetInputHist(TString inFile, TString histName){
    TFile *jetPtFile = new TFile(inFile,"read");
    TH1 *hh = dynamic_cast<TH1*>(jetPtFile->Get(histName.Data()));
    return hh;
}

TGraphAsymmErrors* GetInputGraph(TString inFile, TString histName){
    TFile *jetPtFile = new TFile(inFile,"read");
    TGraphAsymmErrors *hh = dynamic_cast<TGraphAsymmErrors*>(jetPtFile->Get(histName.Data()));
    return hh;
}

std::tuple<TH1D*, TGraphAsymmErrors*, TH1D*,  TGraphAsymmErrors*> GetData(TString dataFile, TString histBase1, TString histBase2, TString histBase3, TString histBase4){
    TFile *jetPtFile = new TFile(dataFile,"read");
    TH1D *h1 = dynamic_cast<TH1D*>(jetPtFile->Get(histBase1.Data()));
     TGraphAsymmErrors *h2 = dynamic_cast<TGraphAsymmErrors*>(jetPtFile->Get(histBase2.Data()));
    TH1D *h3 = dynamic_cast<TH1D*>(jetPtFile->Get(histBase3.Data()));
     TGraphAsymmErrors *h4 = dynamic_cast<TGraphAsymmErrors*>(jetPtFile->Get(histBase4.Data()));
     return std::make_tuple(h1, h2, h3,h4);
}

std::tuple<TH1D*, TGraphAsymmErrors*, TGraph*, TH1D*,  TGraphAsymmErrors*> GetData(TString dataFile, TString histBase1, TString histBase2, TString histBase3, TString histBase4, TString histBase5){
    TFile *jetPtFile = new TFile(dataFile,"read");
    TH1D *h1 = dynamic_cast<TH1D*>(jetPtFile->Get(histBase1.Data()));
     TGraphAsymmErrors *h2 = dynamic_cast<TGraphAsymmErrors*>(jetPtFile->Get(histBase2.Data()));
     TGraph *h3 = dynamic_cast<TGraph*>(jetPtFile->Get(histBase3.Data()));
    TH1D *h4 = dynamic_cast<TH1D*>(jetPtFile->Get(histBase4.Data()));
     TGraphAsymmErrors *h5 = dynamic_cast<TGraphAsymmErrors*>(jetPtFile->Get(histBase5.Data()));
     return std::make_tuple(h1, h2, h3,h4,h5);
}



void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, Size_t Msize, Width_t Lwidth, Style_t Lstyle){
    hh->SetMarkerColor(color);
    hh->SetMarkerStyle(Mstyle);;
    hh->SetLineColor(color);
    hh->SetLineWidth(Lwidth);
    hh->SetMarkerSize(Msize);
    hh->SetLineStyle(Lstyle);
    hh->SetTitle("");
    hh->GetXaxis()->SetTitle("p_{T}^{ch,jet} (GeV/c)");
}

void SaveCanvas(TCanvas *c, TString name){
    c->SaveAs(Form("%s.png",name.Data()));
    c->SaveAs(Form("%s.pdf",name.Data()));
}
