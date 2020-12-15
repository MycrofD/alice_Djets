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

//Double_t fRpar = 0.4;
//double jetEta = 0.9 - fRpar;
//double dy = 2*jetEta;
//Bool_t pdf = false;

static Double_t dy;
static Bool_t pdf;
static Int_t *jetpTbins;
static Int_t *DpTbins[2];
static Int_t zBin;
static Double_t plotRanges[4]{0,2,10e-9,1};


//Color_t colors[] = {1,8,4,2,6,kOrange-1,kGray+2,kCyan+1,kMagenta+2,kViolet+5,kYellow+2};
//Int_t markers[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};
//Style_t linestyle[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};

void ScaleHist(TH1 *hh, int full = 0);
void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, Size_t Msize = 1.1f, Width_t Lwidth = 2, Style_t Lstyle = 1);
void SaveCanvas(TCanvas *c, TString name = "tmp");
void getSystematics(TString inDir, TString outPlotDir);
TH1D* GetUpSys(TH1D **hh, Int_t nFiles = 0);
TH1D* GetDownSys(TH1D **hh, Int_t nFiles = 0);
TH1* GetInputHist(TString inFile, TString histName);
void ScaleHist(TH1 *hh, int full);
void drawFinal(TString outPlotDir);

static TString fPowhegPythia6[9] = {
     "AnalysisResults_FastSim_powheg+pythia6_charm_central"
    ,"AnalysisResults_FastSim_powheg+pythia6_charm_m13_1536595965"
    ,"AnalysisResults_FastSim_powheg+pythia6_charm_m17_1536655729"
    ,"AnalysisResults_FastSim_powheg+pythia6_charm_F1R05_1536598175"
    ,"AnalysisResults_FastSim_powheg+pythia6_charm_F05R1_1536604800"
    ,"AnalysisResults_FastSim_powheg+pythia6_charm_F2R1_1535916012"
    ,"AnalysisResults_FastSim_powheg+pythia6_charm_F1R2_1536594271"
    ,"AnalysisResults_FastSim_powheg+pythia6_charm_F05R05_1535894261"
    ,"AnalysisResults_FastSim_powheg+pythia6_charm_F2R2_1535895146"
    //"AnalysisResults_FastSim_powheg+pythia6_charm_mc13",
    //"AnalysisResults_FastSim_powheg+pythia6_charm_mc17",
    //"AnalysisResults_FastSim_powheg+pythia6_charm_uF1uR05",
    //"AnalysisResults_FastSim_powheg+pythia6_charm_uF05uR1",
    //"AnalysisResults_FastSim_powheg+pythia6_charm_uF2uR1",
    //"AnalysisResults_FastSim_powheg+pythia6_charm_uF1uR2",
    //"AnalysisResults_FastSim_powheg+pythia6_charm_uF05uR05",
    //"AnalysisResults_FastSim_powheg+pythia6_charm_uF2uR2"
    };
static TString fPowhegPythia8[9] = {
     "AnalysisResults_FastSim_powheg+pythia8_charm_central_1593017599"
    ,"AnalysisResults_FastSim_powheg+pythia8_charm_m13_1593091148"
    ,"AnalysisResults_FastSim_powheg+pythia8_charm_m17_1593090378"
    ,"AnalysisResults_FastSim_powheg+pythia8_charm_F1R05_1593338244"
    ,"AnalysisResults_FastSim_powheg+pythia8_charm_F05R1_1593337904"
    ,"AnalysisResults_FastSim_powheg+pythia8_charm_F2R1_1593341637"
    ,"AnalysisResults_FastSim_powheg+pythia8_charm_F1R2_1593339352"
    ,"AnalysisResults_FastSim_powheg+pythia8_charm_F05R05_1593334134"
    ,"AnalysisResults_FastSim_powheg+pythia8_charm_F2R2_1593354950"
    //"AnalysisResults_FastSim_powheg+pythia8_charm_central",
    //"AnalysisResults_FastSim_powheg+pythia8_charm_mc13",
    //"AnalysisResults_FastSim_powheg+pythia8_charm_mc17",
    //"AnalysisResults_FastSim_powheg+pythia8_charm_uF1uR05",
    //"AnalysisResults_FastSim_powheg+pythia8_charm_uF05uR1",
    //"AnalysisResults_FastSim_powheg+pythia8_charm_uF2uR1",
    //"AnalysisResults_FastSim_powheg+pythia8_charm_uF1uR2",
    //"AnalysisResults_FastSim_powheg+pythia8_charm_uF05uR05",
    //"AnalysisResults_FastSim_powheg+pythia8_charm_uF2uR2"
    };


static TString fPythia6[9] = {
   //"AnalysisResults_FastSim_pythia6_charm"
   "AnalysisResults_FastSim_pythia6_charm_1594235514"
};

static TString fPythia8[9] = {
   //"AnalysisResults_FastSim_pythia8_charm"
   "AnalysisResults_FastSim_pythia8_charm_1594137862"
};

static TString fPythia8SoftMode2[9] = {
   //"AnalysisResults_FastSim_pythia8_charm_Color2Soft50"
   "AnalysisResults_FastSim_pythia8_charm_soft2color_1598639645"
};
static TString fPowhegPythia6dijet[9] = {
   "AnalysisResults_FastSim_powheg+pythia8_dijetHadi"
};
static TString fPowhegPythia8dijet[9] = {
   "AnalysisResults_FastSim_powheg+pythia8_dijetSalvatore"
};

std::tuple<TGraphAsymmErrors*, TH1D*, TH1D*, TH1D*, TH1D**> GetSim(TString simname, Int_t type, Int_t nFiles, TString *names, TString simDir, Double_t simScaling, UInt_t nBins, Double_t *xBins);
std::tuple<TGraphAsymmErrors*, TH1D*, TH1D*, TH1D*> GetDataSimRatio(TString simname, TH1D *data_cent, TH1D *sim_cent, TH1D *sim_up, TH1D *sim_down, UInt_t nBins, Double_t *xBins);
std::tuple<TGraphAsymmErrors*, TH1D*,TGraphAsymmErrors*, TH1D*> GetData(TString dataFile, TString histBase, Double_t dataScaling, UInt_t nBins, Double_t *xBins, Double_t *systUncD_down, Double_t *systUncD_up);
std::tuple<TCanvas*, TPad*, TPad*, TH1D*, TH1D*> PrepareCanvas(UInt_t xAxisBins, Double_t *xAxis);
void PlaceOnPadData(TPad* pad,TGraphAsymmErrors* histo1, TH1D* histo2, Size_t markersize);
void PlaceOnPadSim(TPad* pad, TGraphAsymmErrors* histo, Color_t ci, Style_t style, Style_t linestyle);
void TerminateCanvas(TPad* pad1,TPad* pad2,TH1D* histo1,TH1D* histo2);


void finalJetSpectraInvUpdated(
Int_t type = 0, //0 - pt x-section, 1 - z x-section, 2 - z PDF
Int_t radius = 4, // 2, 4 or 6
Int_t z = 0, //which z jet pT bin 1,2,3,4,5
Int_t sysGlobal = 0, //0 - (defualt )add LumiUnc, BRUnc and Tracking Unc. + add cut var and JES;  for x-section
                     //1 - (R comparison) no global Unc; cut var separate and no JES - available only for pT x-section (type 0)
                     //2 - (energy comparison) add LumiUnc and Tracking Unc cut var and JES
TString dataFile = "/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR04_paperCuts/Default/unfolding_Bayes_4/unfoldedSpectrum_unfoldedJetSpectrum.root",
TString dataAnalysisFile = "/mnt/hgfs/vmware/data_R04_050219/data/AnalysisResults_Run2w18b.root",
TString simDir = "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts/Simulations/Prompt/AnalysisResults_Run2w18b.root",
TString outSpectraDir = "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW",
TString histBase = "unfoldedSpectrumKineEff"
)
/*
void finalJetSpectraInvUpdated(
Int_t type = 0,
Int_t radius = 4,
Int_t z = 0,
TString dataFile = "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_ztest/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/unfoldedSpectrum_unfoldedJetSpectrum.root",
TString dataAnalysisFile = "/mnt/hgfs/vmware/data_R04_050219/data/AnalysisResults_Run2w18b.root",
TString simDir = "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_ztest/Simulations/Prompt/AnalysisResults_Run2w18b.root",
TString outSpectraDir = "/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts_ztest/Default_AnalysisResults_Run2w18b.root/unfolding_2D_5/finalSpectraNEW",
TString histBase = "unfoldedSpectrumKineEff"
)*/
{
bool fivetev = 1;
radius = 4;
dataFile = Form("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0%d_paperCuts/Default/unfolding_Bayes_5/unfoldedSpectrum_unfoldedJetSpectrum.root",(int)radius);
dataAnalysisFile="/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_";//745_R04_pp_5cuts.root";
if (radius == 2){dataAnalysisFile+="744_R02_pp_5cuts.root";}
else if (radius == 4){dataAnalysisFile+="745_R04_pp_5cuts.root";}
else if (radius == 6){dataAnalysisFile+="751_R06_pp_5cuts.root";}
simDir=Form("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0%d_paperCuts/Simulations/Prompt/",(int)radius);
outSpectraDir = Form("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR0%d_paperCuts/Default/unfolding_Bayes_5/final/",(int)radius);
histBase="unfoldedSpectrum";
    TString outPlotDir = outSpectraDir;
    outPlotDir+="/plots";
 //   gSystem->Exec(Form("mkdir %s",outSpectraDir.Data()));
 //   gSystem->Exec(Form("mkdir %s",outPlotDir.Data()));

    Bool_t addBRDzeroUnc = false;
    Bool_t addDtrackingUnc = false;
    Bool_t addLumiUnc = false;
    Bool_t separateCUTUnc = false;
    Bool_t addJESUnc = false;



    UInt_t xAxisBins = 10;
    Double_t *xAxis =nullptr;// = new Double_t[xAxisBins];
    Double_t *xAxisC =nullptr;

    Double_t *systUncD_up = nullptr;
    Double_t *systUncD_down = nullptr;
    Double_t *systUncD_JES = nullptr;
    Double_t *systUncD_CUTS = nullptr;

    zBin=z;
    if(type==0) zBin = 0;

    // ----------------------------------------------------------------
    // ------------------- ENABLE SIMULATIONS HERE --------------------
    // ----------------------------------------------------------------
    Bool_t ePowhegPythia6 = false;
    Bool_t ePowhegPythia8 = true;
    Bool_t ePythia6 = false;
    Bool_t ePythia8 = true;
    Bool_t ePythia8SoftMode2 = true;
    Bool_t ePowhegPythia6dijet = false;
    Bool_t ePowhegPythia8dijet = false;

    // ----------------------------------------------------------------
    // --------------- SET PARAMETERS + BINNING HERE ------------------
    // ----------------------------------------------------------------
    Double_t sigma_in = 0.0578;
    const Double_t  BRDzero = 0.0389;
    Double_t BRDzeroUnc = 0.0004;
    Double_t DtrackingUnc = 0.05;
    Double_t LumiUnc = 0.05;
    if(fivetev){
        sigma_in = 0.05077;
        DtrackingUnc = 0.03;
        LumiUnc=0.021;
    }
    if(type ==0){//x-section
        xAxisBins = 8;
        xAxis = new Double_t[xAxisBins+1]{5,6,8,10,12,14,20,30,50};
        xAxisC = new Double_t[xAxisBins]{5.5,7,9,11,13,17,25,40};
        if(fivetev){
            xAxisBins = 7;
            xAxis = new Double_t[xAxisBins+1]{5,6,8,10,14,20,30,50};
            xAxisC = new Double_t[xAxisBins]{5.5,7,9,12,17,25,40};
        }
    }
    else if(type ==1 || type ==2){ //Z x-section or Z PDF
        xAxisBins = 5;
        xAxis = new Double_t[xAxisBins+1]{0.4,0.6,0.7,0.8,0.9,1.0};
        jetpTbins = new Int_t[6]{2,5,7,10,15,50};
        if(radius ==4 || radius ==6 || radius ==3)DpTbins[0] = new Int_t[5]{2,2,3,5,5};
        if(radius ==4 || radius ==6 || radius ==3)DpTbins[1] = new Int_t[5]{5,7,10,15,36};
        if(radius ==2)DpTbins[0] = new Int_t[5]{2,2,4,5,10};
        if(radius ==2)DpTbins[1] = new Int_t[5]{5,7,10,15,36};
        simDir+=Form("%d",zBin);
        histBase+=Form("%d",zBin);
        outSpectraDir+=Form("%d",zBin);
        outPlotDir = outSpectraDir+"/plots";
    }
    if(type ==0 || type ==1){//x-section for [jetpt or z]
        if(sysGlobal ==0){
            addBRDzeroUnc = true;
            addDtrackingUnc = true;
            addLumiUnc = true;
            separateCUTUnc = false;
            addJESUnc = true;
        }
        else if(sysGlobal ==1){//R comparison, no global unc; cut var separate and no JES
            addBRDzeroUnc = false;
            addDtrackingUnc = false;
            addLumiUnc = false;
            separateCUTUnc = true;
            addJESUnc = false;
        }
        else if(sysGlobal ==2){//energy comparison
            addBRDzeroUnc = false;
            addDtrackingUnc = true;
            addLumiUnc = true;
            separateCUTUnc = false;
            addJESUnc = true;
        }
    }
    else if(type ==2){
        addBRDzeroUnc = false;
        addDtrackingUnc = false;
        addLumiUnc = false;
        separateCUTUnc = false;
        addJESUnc = true; //z-already have it in total
    }

    if(type ==2){plotRanges[2]=0; plotRanges[3]=10;}

   // std::cout<<DpTbins[0][2]<<" "<<DpTbins[1][4]<<std::endl;
    gSystem->Exec(Form("mkdir %s",outSpectraDir.Data()));
    gSystem->Exec(Form("mkdir %s",outPlotDir.Data()));

   // gSystem->Exec(Form("mkdir %s",outSpectraDir.Data()));

    // ----------------------------------------------------------------
    // ------------ SET SYSTEMATICS AND PLOT RANGES HERE --------------
    // ----------------------------------------------------------------
    if(type ==0){//x-section
        pdf = false;
        if(radius == 2){
            dy = 2*(0.9 - 0.2);
            plotRanges[0]=0; plotRanges[1]=2.1; plotRanges[2]=0.000004; plotRanges[3]=2.5;
            if(fivetev){
                //total
                //systUncD_up = new Double_t[xAxisBins]{0.0891, 0.0912, 0.0721, 0.0989, 0.078, 0.0827, 0.1145, 0.2137};
                //systUncD_down = new Double_t[xAxisBins]{0.1009, 0.1034,	0.0907, 0.1161,	0.098, 0.1144,	0.1386, 0.2486};
                //JES and CUTS not in systUncD_up and systUncD_down
                systUncD_up = new Double_t[xAxisBins]{0.065, 0.067, 0.062, 0.072, 0.109, 0.171, 0.225};
                systUncD_down = new Double_t[xAxisBins]{0.083, 0.084, 0.082, 0.093,	0.135, 0.193, 0.283};
                systUncD_JES = new Double_t[xAxisBins]{0.011, 0.016, 0.022, 0.031, 0.047, 0.071, 0.117};
                systUncD_CUTS = new Double_t[xAxisBins]{0.014, 0.017, 0.021, 0.026, 0.036, 0.052, 0.081};
            }
            else{
                //total
                //systUncD_up = new Double_t[xAxisBins]{0.0891, 0.0912, 0.0721, 0.0989, 0.078, 0.0827, 0.1145, 0.2137};
                //systUncD_down = new Double_t[xAxisBins]{0.1009, 0.1034,	0.0907, 0.1161,	0.098, 0.1144,	0.1386, 0.2486};
                //JES and CUTS not in systUncD_up and systUncD_down
                systUncD_up = new Double_t[xAxisBins]{0.0883, 0.0891, 0.0690, 0.0893, 0.0699, 0.0727, 0.0921, 0.1790};
                systUncD_down = new Double_t[xAxisBins]{0.1002, 0.1015,	0.0883, 0.1080,	0.0917, 0.1074,	0.1208, 0.2195};
                systUncD_JES = new Double_t[xAxisBins]{0.005, 0.009, 0.013, 0.014, 0.017, 0.022, 0.03, 0.066};
                systUncD_CUTS = new Double_t[xAxisBins]{0.011, 0.018, 0.017, 0.04, 0.03, 0.032,	0.061, 0.096};
            }
        }
        else if(radius == 3){
            dy = 2*(0.9 - 0.3);
            plotRanges[0]=0; plotRanges[1]=2.1; plotRanges[2]=0.000004; plotRanges[3]=2.5;
            if(fivetev){
                //total
                //systUncD_up = new Double_t[xAxisBins]{0.0891, 0.0912, 0.0721, 0.0989, 0.078, 0.0827, 0.1145, 0.2137};
                //systUncD_down = new Double_t[xAxisBins]{0.1009, 0.1034,	0.0907, 0.1161,	0.098, 0.1144,	0.1386, 0.2486};
                //JES and CUTS not in systUncD_up and systUncD_down
                systUncD_up = new Double_t[xAxisBins]{0.078, 0.067, 0.08, 0.093, 0.115, 0.174, 0.254};
                systUncD_down = new Double_t[xAxisBins]{0.094, 0.088, 0.105, 0.123,	0.153, 0.215, 0.274};
                systUncD_JES = new Double_t[xAxisBins]{0.016, 0.02, 0.024, 0.031, 0.043, 0.061, 0.096};
                systUncD_CUTS = new Double_t[xAxisBins]{0.034, 0.043, 0.056, 0.075, 0.107, 0.158, 0.253};
            }
            else{
                //total
                return;
            }
        }
        else if(radius == 4){
            dy = 2*(0.9 - 0.4);
            plotRanges[0]=0; plotRanges[1]=2.1; plotRanges[2]=0.000004; plotRanges[3]=2.5;
            if(fivetev){
                //total
                //systUncD_up = new Double_t[xAxisBins]{0.0891, 0.0912, 0.0721, 0.0989, 0.078, 0.0827, 0.1145, 0.2137};
                //systUncD_down = new Double_t[xAxisBins]{0.1009, 0.1034,	0.0907, 0.1161,	0.098, 0.1144,	0.1386, 0.2486};
                //JES and CUTS not in systUncD_up and systUncD_down
                systUncD_up = new Double_t[xAxisBins]{0.075, 0.062, 0.076, 0.099, 0.155, 0.169, 0.307};
                systUncD_down = new Double_t[xAxisBins]{0.091, 0.082, 0.104, 0.133,	0.203, 0.232, 0.362};
                systUncD_JES = new Double_t[xAxisBins]{0.016, 0.02, 0.024, 0.031, 0.043, 0.061, 0.096};
                systUncD_CUTS = new Double_t[xAxisBins]{0.034, 0.043, 0.056, 0.075, 0.107, 0.158, 0.253};
            }
            else{
                //total
                //systUncD_up = new Double_t[xAxisBins]{0.0762, 0.0704, 0.0781, 0.0972, 0.1156, 0.1195, 0.1531, 0.2121};
                //systUncD_down = new Double_t[xAxisBins]{0.0883, 0.0853,	0.0989, 0.1224,	0.1425, 0.1582,	0.2135, 0.2780};
                //JES and CUTS not in systUncD_up and systUncD_down
                systUncD_up = new Double_t[xAxisBins]{0.0668, 0.0669, 0.0695, 0.0881, 0.1049, 0.1002, 0.1279, 0.1664};
                systUncD_down = new Double_t[xAxisBins]{0.0803, 0.0824,	0.0922, 0.1153,	0.1339, 0.1442,	0.1962, 0.2450};
                systUncD_JES = new Double_t[xAxisBins]{0.008, 0.019, 0.021, 0.035, 0.039, 0.042, 0.064, 0.097};
                systUncD_CUTS = new Double_t[xAxisBins]{0.036, 0.012, 0.029, 0.021,	0.030, 0.050, 0.054, 0.088};
            }
        }
        else if(radius == 6){
            dy = 2*(0.9 - 0.6);
            plotRanges[0]=0; plotRanges[1]=2.7; plotRanges[2]=0.000004; plotRanges[3]=2.5;
            if(fivetev){
                //total
                //systUncD_up = new Double_t[xAxisBins]{0.0891, 0.0912, 0.0721, 0.0989, 0.078, 0.0827, 0.1145, 0.2137};
                //systUncD_down = new Double_t[xAxisBins]{0.1009, 0.1034,	0.0907, 0.1161,	0.098, 0.1144,	0.1386, 0.2486};
                //JES and CUTS not in systUncD_up and systUncD_down
                systUncD_up = new Double_t[xAxisBins]{0.091, 0.085, 0.067, 0.116, 0.156, 0.233, 0.35};
                systUncD_down = new Double_t[xAxisBins]{0.101, 0.10, 0.087, 0.141,	0.195, 0.342, 0.436};
                systUncD_JES = new Double_t[xAxisBins]{0.003, 0.009, 0.017, 0.029, 0.05, 0.082, 0.143};
                systUncD_CUTS = new Double_t[xAxisBins]{0.048, 0.055, 0.064, 0.078, 0.101, 0.138, 0.208};
            }
            else{
                //total
                //systUncD_up = new Double_t[xAxisBins]{0.097, 0.073, 0.09, 0.103, 0.119, 0.136, 0.189, 0.358};
                //systUncD_down = new Double_t[xAxisBins]{0.103, 0.085,	0.104, 0.119,	0.142, 0.165,	0.221, 0.43};
                //JES and CUTS not in systUncD_up and systUncD_down
                systUncD_up = new Double_t[xAxisBins]{0.074, 0.068, 0.071, 0.087, 0.100, 0.108, 0.126, 0.202};
                systUncD_down = new Double_t[xAxisBins]{0.082, 0.081, 0.088, 0.106,	0.126, 0.142, 0.170, 0.313};
                systUncD_JES = new Double_t[xAxisBins]{0.001, 0.006, 0.029, 0.038, 0.052, 0.064, 0.075, 0.151};
                systUncD_CUTS = new Double_t[xAxisBins]{0.063, 0.026, 0.048, 0.040,	0.040, 0.053, 0.119, 0.254};
            }
        }
    }
    else if(type ==1 || type ==2){//Z x-section
        pdf = false;
        if(radius == 2){
            dy = 2*(0.9 - 0.2);
            if(zBin ==2){
                systUncD_up = new Double_t[xAxisBins]{0,0,0,0,0};
                systUncD_down = new Double_t[xAxisBins]{0,0,0,0,0};
            }
            else if(zBin ==3){
                systUncD_up = new Double_t[xAxisBins]{0,0,0,0,0};
                systUncD_down = new Double_t[xAxisBins]{0,0,0,0,0};
            }
            else if(zBin ==4){
                systUncD_up = new Double_t[xAxisBins]{0,0,0,0,0};
                systUncD_down = new Double_t[xAxisBins]{0,0,0,0,0};
            }
            else if(zBin ==5){
                systUncD_up = new Double_t[xAxisBins]{0,0,0,0,0};
                systUncD_down = new Double_t[xAxisBins]{0,0,0,0,0};
            }
        }
        else if(radius == 4){
            dy = 2*(0.9 - 0.4);
            if(zBin ==2){
                plotRanges[0]=0; plotRanges[1]=2.1; plotRanges[2]=0.001; plotRanges[3]=0.4;
                if(type == 2)plotRanges[0]=0; plotRanges[1]=2.1; plotRanges[2]=0.001; plotRanges[3]=10;
                systUncD_up = new Double_t[xAxisBins]{0.14, 0.08, 0.08,	0.09, 0.09};
                systUncD_down = new Double_t[xAxisBins]{0.15, 0.09, 0.09, 0.10,	0.09};
            }
            else if(zBin ==3){
                systUncD_up = new Double_t[xAxisBins]{0.14, 0.08, 0.08,	0.09, 0.09};
                systUncD_down = new Double_t[xAxisBins]{0.15, 0.09, 0.09, 0.10,	0.09};
            }
            else if(zBin ==4){
                systUncD_up = new Double_t[xAxisBins]{0.14, 0.08, 0.08,	0.09, 0.09};
                systUncD_down = new Double_t[xAxisBins]{0.15, 0.09, 0.09, 0.10,	0.09};
            }
            else if(zBin ==5){
                systUncD_up = new Double_t[xAxisBins]{0.14, 0.08, 0.08,	0.09, 0.09};
                systUncD_down = new Double_t[xAxisBins]{0.15, 0.09, 0.09, 0.10,	0.09};
            };
        }
        else if(radius == 6){
            dy = 2*(0.9 - 0.6);
            if(zBin ==2){
                systUncD_up = new Double_t[xAxisBins]{0,0,0,0,0};
                systUncD_down = new Double_t[xAxisBins]{0,0,0,0,0};
            }
            else if(zBin ==3){
                systUncD_up = new Double_t[xAxisBins]{0,0,0,0,0};
                systUncD_down = new Double_t[xAxisBins]{0,0,0,0,0};
            }
            else if(zBin ==4){
                systUncD_up = new Double_t[xAxisBins]{0,0,0,0,0};
                systUncD_down = new Double_t[xAxisBins]{0,0,0,0,0};
            }
            else if(zBin ==5){
                systUncD_up = new Double_t[xAxisBins]{0,0,0,0,0};
                systUncD_down = new Double_t[xAxisBins]{0,0,0,0,0};
            }
        }
    }

    // ----------------------------------------------------------------
    // add uncertainties
    Int_t bit = -1;
    std::cout<<"sys bitmap: ";
    std::cout<<(bit=addBRDzeroUnc?1:0);
    std::cout<<(bit=addDtrackingUnc?1:0);
    std::cout<<(bit=addLumiUnc?1:0);
    std::cout<<(bit=(addJESUnc && systUncD_JES)?1:0);
    std::cout<<(bit=(!separateCUTUnc && systUncD_CUTS)?1:0);
    std::cout<<std::endl;
    for (UInt_t bin = 0; bin < xAxisBins; bin++){
        Double_t totalUnc = 0;
        if(addBRDzeroUnc) totalUnc += BRDzeroUnc*BRDzeroUnc;
        if(addDtrackingUnc) totalUnc += DtrackingUnc*DtrackingUnc;
        if(addLumiUnc) totalUnc += LumiUnc*LumiUnc;
        if(addJESUnc && systUncD_JES) totalUnc += systUncD_JES[bin]*systUncD_JES[bin];
        if(!separateCUTUnc && systUncD_CUTS) totalUnc += systUncD_CUTS[bin]*systUncD_CUTS[bin];
        systUncD_up[bin] = TMath::Sqrt(systUncD_up[bin]*systUncD_up[bin]+totalUnc);
        systUncD_down[bin] = TMath::Sqrt(systUncD_down[bin]*systUncD_down[bin]+totalUnc);
        std::cout<<systUncD_up[bin]<<" ";
    }
    std::cout<<std::endl;
    // in case of separate Cut sys create separate plot
    // JES is not added only for ratios plots wjich must be calculated outside manually anyway, thus is not added
    TGraph *gDataCutSys = nullptr;
    if(separateCUTUnc && systUncD_CUTS) gDataCutSys = new TGraph(xAxisBins, xAxisC, systUncD_CUTS);

    //get Data n Events
    std::cout<<"get Data n Events"<<std::endl;
    TFile *File = new TFile(dataAnalysisFile,"read");
    TDirectoryFile* dir = dynamic_cast<TDirectoryFile*>(File->Get("PWG3_D2H_DmesonsForJetCorrelationsMBN0"));
    AliNormalizationCounter *c = dynamic_cast<AliNormalizationCounter*>(dir->Get("NormalizationCounter"));
    double nEv = c->GetNEventsForNorm();
    double dataLum = nEv/(sigma_in*1000) ;//Luminosity in mbar
    double simScaling = 0.5; //pp 0.5
    double dataScaling = 1. /(BRDzero * dataLum)/2.;
//    std::cout<<"scalings "<<jetEta<<" "<<dy<<" "<<dataScaling<<" "<<BRDzero<<" "<<dataLum<<" "<<nEv<<" "<<sigma_in<<std::endl;
    std::cout<<"scalings: "<<dy<<" "<<dataScaling<<" "<<BRDzero<<" "<<dataLum<<" "<<nEv<<" "<<sigma_in<<std::endl;



    // ----------------- Initiate Canvas ---------------------
    TCanvas *canvas;
    TPad *upPad, *dowmPad;
    TH1D *placeholder_up, *placeholder_down;
    std::tie(canvas,upPad,dowmPad,placeholder_up,placeholder_down) = PrepareCanvas(xAxisBins,xAxis);

    // ----------------- data ---------------------
    TGraphAsymmErrors *hDataSys, *hDataSysRatio;
    TH1D *hData_binned, *hData_binned_ratio;
    std::tie(hDataSys, hData_binned, hDataSysRatio, hData_binned_ratio) = GetData(dataFile, histBase, dataScaling, xAxisBins, xAxis, systUncD_down, systUncD_up);

    PlaceOnPadData(upPad,hDataSys,hData_binned,0.9);
    PlaceOnPadData(dowmPad,hDataSysRatio,hData_binned_ratio,0);


    // ----------------- prompt simulation ---------------------
    TH1D **hBlackHole = nullptr; //just dump

    TH1D **simPowhegPythia6var = nullptr;
    TGraphAsymmErrors *simPowhegPythia6 = nullptr;
    TH1D* simPowhegPythia6_cent = nullptr, *simPowhegPythia6_down = nullptr, *simPowhegPythia6_up = nullptr;
    TGraphAsymmErrors *simPowhegPythia6_R = nullptr;
    TH1D* simPowhegPythia6_cent_R = nullptr, *simPowhegPythia6_down_R = nullptr, *simPowhegPythia6_up_R = nullptr;
    if(ePowhegPythia6){
        std::cout<<"get POWHWG+PYTHIA6"<<std::endl;
        std::tie(simPowhegPythia6, simPowhegPythia6_cent, simPowhegPythia6_up, simPowhegPythia6_down,simPowhegPythia6var) = GetSim("simPowhegPythia6",type, 9, fPowhegPythia6, simDir, simScaling, xAxisBins, xAxis);
        std::tie(simPowhegPythia6_R, simPowhegPythia6_cent_R, simPowhegPythia6_up_R, simPowhegPythia6_down_R) = GetDataSimRatio("simPowhegPythia6",hData_binned,simPowhegPythia6_cent, simPowhegPythia6_up, simPowhegPythia6_down, xAxisBins, xAxis);
        PlaceOnPadSim(upPad,simPowhegPythia6,static_cast<Color_t>(TColor::GetColor("#000099")),24,1);
        PlaceOnPadSim(dowmPad,simPowhegPythia6_R,static_cast<Color_t>(TColor::GetColor("#000099")),24,1);
    }

    TH1D **simPowhegPythia8var = nullptr;
    TGraphAsymmErrors *simPowhegPythia8 = nullptr;
    TH1D* simPowhegPythia8_cent = nullptr, *simPowhegPythia8_down = nullptr, *simPowhegPythia8_up = nullptr;
    TGraphAsymmErrors *simPowhegPythia8_R = nullptr;
    TH1D* simPowhegPythia8_cent_R = nullptr, *simPowhegPythia8_down_R = nullptr, *simPowhegPythia8_up_R = nullptr;
    if(ePowhegPythia8){
        std::cout<<"get POWHWG+PYTHIA8"<<std::endl;
        std::tie(simPowhegPythia8, simPowhegPythia8_cent, simPowhegPythia8_up, simPowhegPythia8_down,simPowhegPythia8var) = GetSim("simPowhegPythia8",type, 9, fPowhegPythia8, simDir, simScaling, xAxisBins, xAxis);
        std::tie(simPowhegPythia8_R, simPowhegPythia8_cent_R, simPowhegPythia8_up_R, simPowhegPythia8_down_R) = GetDataSimRatio("simPowhegPythia8",hData_binned,simPowhegPythia8_cent, simPowhegPythia8_up, simPowhegPythia8_down, xAxisBins, xAxis);
        PlaceOnPadSim(upPad,simPowhegPythia8,static_cast<Color_t>(TColor::GetColor("#000099")),24,1);
        PlaceOnPadSim(dowmPad,simPowhegPythia8_R,static_cast<Color_t>(TColor::GetColor("#000099")),24,1);
    }

    TGraphAsymmErrors *simPythia6 = nullptr;
    TH1D* simPythia6_cent = nullptr, *simPythia6_down = nullptr, *simPythia6_up = nullptr;
    TGraphAsymmErrors *simPythia6_R = nullptr;
    TH1D* simPythia6_cent_R = nullptr, *simPythia6_down_R = nullptr, *simPythia6_up_R = nullptr;
    if(ePythia6){
        std::cout<<"get PYTHIA6"<<std::endl;
        std::tie(simPythia6, simPythia6_cent, simPythia6_up, simPythia6_down,hBlackHole) = GetSim("simPythia6",type, 1,fPythia6, simDir, simScaling, xAxisBins, xAxis);
        std::tie(simPythia6_R, simPythia6_cent_R, simPythia6_up_R, simPythia6_down_R) = GetDataSimRatio("simPythia6",hData_binned,simPythia6_cent, simPythia6_up, simPythia6_down, xAxisBins, xAxis);
        PlaceOnPadSim(upPad,simPythia6,static_cast<Color_t>(TColor::GetColor("#009933")),25,2);
        PlaceOnPadSim(dowmPad,simPythia6_R,static_cast<Color_t>(TColor::GetColor("#009933")),25,2);
    }

    TGraphAsymmErrors *simPythia8 = nullptr;
    TH1D* simPythia8_cent = nullptr, *simPythia8_down = nullptr, *simPythia8_up = nullptr;
    TGraphAsymmErrors *simPythia8_R = nullptr;
    TH1D* simPythia8_cent_R = nullptr, *simPythia8_down_R = nullptr, *simPythia8_up_R = nullptr;
    if(ePythia8){
        std::cout<<"get PYTHIA8"<<std::endl;
        std::tie(simPythia8, simPythia8_cent, simPythia8_up, simPythia8_down,hBlackHole) = GetSim("simPythia8",type, 1,fPythia8, simDir, simScaling, xAxisBins, xAxis);
        std::tie(simPythia8_R, simPythia8_cent_R, simPythia8_up_R, simPythia8_down_R) = GetDataSimRatio("simPythia8",hData_binned,simPythia8_cent, simPythia8_up, simPythia8_down, xAxisBins, xAxis);
        PlaceOnPadSim(upPad,simPythia8,kViolet+2,27,3);
        PlaceOnPadSim(dowmPad,simPythia8_R,kViolet+2,27,3);
    }

    TGraphAsymmErrors *simPythia8Soft2 = nullptr;
    TH1D* simPythia8Soft2_cent = nullptr, *simPythia8Soft2_down = nullptr, *simPythia8Soft2_up = nullptr;
    TGraphAsymmErrors *simPythia8Soft2_R = nullptr;
    TH1D* simPythia8Soft2_cent_R = nullptr, *simPythia8Soft2_down_R = nullptr, *simPythia8Soft2_up_R = nullptr;
    if(ePythia8SoftMode2){
        std::cout<<"get PYTHIA8 soft mode2"<<std::endl;
        std::tie(simPythia8Soft2, simPythia8Soft2_cent, simPythia8Soft2_up, simPythia8Soft2_down,hBlackHole) = GetSim("simPythia8Soft2",type, 1,fPythia8SoftMode2, simDir, simScaling, xAxisBins, xAxis);
        std::tie(simPythia8Soft2_R, simPythia8Soft2_cent_R, simPythia8Soft2_up_R, simPythia8Soft2_down_R) = GetDataSimRatio("simPythia8Soft2",hData_binned,simPythia8Soft2_cent, simPythia8Soft2_up, simPythia8Soft2_down, xAxisBins, xAxis);
        PlaceOnPadSim(upPad,simPythia8Soft2,kOrange+2,28,4);
        PlaceOnPadSim(dowmPad,simPythia8Soft2_R,kOrange+2,28,4);
    }

    TH1D **simPowhegPythia6dijetvar = nullptr;
    TGraphAsymmErrors *simPowhegPythia6dijet = nullptr;
    TGraphAsymmErrors *simPowhegPythia6dijet_R = nullptr;
    TH1D* simPowhegPythia6dijet_cent = nullptr, *simPowhegPythia6dijet_down = nullptr, *simPowhegPythia6dijet_up = nullptr;
    TH1D* simPowhegPythia6dijet_cent_R = nullptr, *simPowhegPythia6dijet_down_R = nullptr, *simPowhegPythia6dijet_up_R = nullptr;
    if(ePowhegPythia6dijet){
        std::cout<<"get POWHWG+PYTHIA6 dijet"<<std::endl;
        std::tie(simPowhegPythia6dijet, simPowhegPythia6dijet_cent, simPowhegPythia6dijet_up, simPowhegPythia6dijet_down,simPowhegPythia6dijetvar) = GetSim("simPowhegPythia6dijet",type, 1, fPowhegPythia6dijet, simDir, simScaling, xAxisBins, xAxis);
        std::tie(simPowhegPythia6dijet_R, simPowhegPythia6dijet_cent_R, simPowhegPythia6dijet_up_R, simPowhegPythia6dijet_down_R) = GetDataSimRatio("simPowhegPythia6dijet",hData_binned,simPowhegPythia6dijet_cent, simPowhegPythia6dijet_up, simPowhegPythia6dijet_down, xAxisBins, xAxis);
        PlaceOnPadSim(upPad,simPowhegPythia6dijet,static_cast<Color_t>(TColor::GetColor("#000099")),28,1);
        PlaceOnPadSim(dowmPad,simPowhegPythia6dijet_R,static_cast<Color_t>(TColor::GetColor("#000099")),28,1);
    }

    TH1D **simPowhegPythia8dijetvar = nullptr;
    TGraphAsymmErrors *simPowhegPythia8dijet = nullptr;
    TGraphAsymmErrors *simPowhegPythia8dijet_R = nullptr;
    TH1D* simPowhegPythia8dijet_cent = nullptr, *simPowhegPythia8dijet_down = nullptr, *simPowhegPythia8dijet_up = nullptr;
    TH1D* simPowhegPythia8dijet_cent_R = nullptr, *simPowhegPythia8dijet_down_R = nullptr, *simPowhegPythia8dijet_up_R = nullptr;
    if(ePowhegPythia8dijet){
        std::cout<<"get POWHWG+PYTHIA8 dijet"<<std::endl;
        std::tie(simPowhegPythia8dijet, simPowhegPythia8dijet_cent, simPowhegPythia8dijet_up, simPowhegPythia8dijet_down,simPowhegPythia8dijetvar) = GetSim("simPowhegPythia8dijet",type, 1, fPowhegPythia8dijet, simDir, simScaling, xAxisBins, xAxis);
        std::tie(simPowhegPythia8dijet_R, simPowhegPythia8dijet_cent_R, simPowhegPythia8dijet_up_R, simPowhegPythia8dijet_down_R) = GetDataSimRatio("simPowhegPythia8dijet",hData_binned,simPowhegPythia8dijet_cent, simPowhegPythia8dijet_up, simPowhegPythia8dijet_down, xAxisBins, xAxis);
        PlaceOnPadSim(upPad,simPowhegPythia8dijet,static_cast<Color_t>(TColor::GetColor("#009999")),28,1);
        PlaceOnPadSim(dowmPad,simPowhegPythia8dijet_R,static_cast<Color_t>(TColor::GetColor("#009999")),28,1);
    }

    // ----------------- Legend and text PAD---------------------
    TLegend *leg =nullptr;
    Double_t shift = 0.06*(ePowhegPythia6+ePowhegPythia8+ePythia6+ePythia8+ePythia8SoftMode2+ePowhegPythia6dijet);
    if(type==0)leg = new TLegend(0.35,0.45,0.65,0.7,nullptr,"NB NDC");
    if(type==1 || type ==2)leg = new TLegend(0.22,0.7-shift,0.5,0.7,nullptr,"NB NDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(43);
    leg->SetTextSize(21);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hDataSys,"Data","fp");

    if(ePowhegPythia6)leg->AddEntry(simPowhegPythia6,"POWHEG hvq + PYTHIA 6","pf");
    if(ePowhegPythia8)leg->AddEntry(simPowhegPythia8,"POWHEG hvq + PYTHIA 8","pf");
    if(ePythia6)leg->AddEntry(simPythia6,"PYTHIA 6 Perugia 2011","l");
    if(ePythia8)leg->AddEntry(simPythia8,"PYTHIA 8 Monash 2013","l");
    if(ePythia8SoftMode2)leg->AddEntry(simPythia8Soft2,"PYTHIA 8 Monash 2013 mode 2","l");
    if(ePowhegPythia6dijet)leg->AddEntry(simPowhegPythia6dijet,"POWHEG dijet + PYTHIA 8 hadi","pf");
    if(ePowhegPythia8dijet)leg->AddEntry(simPowhegPythia8dijet,"POWHEG dijet + PYTHIA 8","l");
    TPaveText *pt[2];
    pt[0] = new TPaveText(0.2,0.9,0.85,0.95,"NB NDC");
    if(type==0)pt[1] = new TPaveText(0.2,0.75,0.85,0.9,"NB NDC");
    if(type==1 ||type==2)pt[1] = new TPaveText(0.2,0.75,0.85,0.9,"NB NDC");
    for(Int_t s = 0; s<2;s++){
        pt[s]->SetBorderSize(0);
        pt[s]->SetFillStyle(0);
        pt[s]->SetTextAlign(13);
        pt[s]->SetTextFont(43);
        pt[s]->SetTextSize(22);
    }
    //pt[0]->AddText("ALICE Preliminary"); //uncomment
    pt[0]->AddText(Form("ALICE, pp, #sqrt{#it{s}} = %s TeV",fivetev?"5.02":"13"));
    pt[1]->AddText(Form("Charged Jets, anti-#it{k}_{T}, #it{R} = 0.%d, |#it{#eta}_{lab}^{jet}| < 0.%d",radius,9-radius));
    if(type==0)pt[1]->AddText(Form ("with D^{0}, %d < #it{p}_{T,D^{0}} < %d GeV/#it{c}",static_cast<int>(2),static_cast<int>(36)));
   // if(type==1 ||type==2)pt[1]->AddText(Form ("%d < #it{p}_{T,jet} < %d GeV/#it{c}",jetpTbins[zBin-1],jetpTbins[zBin]));
   // if(type==1 ||type==2)pt[1]->AddText(Form ("with D^{0}, %d < #it{p}_{T,D^{0}} < %d GeV/#it{c}",DpTbins[0][zBin-1],DpTbins[1][zBin-1]));
    if(type==1 ||type==2)pt[1]->AddText(Form ("%d < #it{p}_{T,jet} < %d GeV/#it{c} with D^{0}, %d < #it{p}_{T,D^{0}} < %d GeV/#it{c}",jetpTbins[zBin-1],jetpTbins[zBin],DpTbins[0][zBin-1],DpTbins[1][zBin-1]));
    //if(type==1 ||type==2)pt[1]->AddText(Form ("",DpTbins[0][zBin-1],DpTbins[1][zBin-1]));

    upPad->cd();
    leg->Draw();
    pt[0]->Draw();
    pt[1]->Draw();
std::cout<<"a"<<std::endl;
    // ----------------- Terminate Canvas ---------------------
    TerminateCanvas(upPad,dowmPad,placeholder_up,placeholder_down);
    //canvas->SaveAs(outPlotDir+"/finalSpectra.png");
    TString sysmode = "";
    if((type ==0 || type ==1) && sysGlobal==0) sysmode = "_fullGlobal_addedCUTandJES";     //default
    if((type ==0 || type ==1) && sysGlobal==1) sysmode = "_noneGlobal_separateCUTnoneJES"; //for R ratio
    if((type ==0 || type ==1) && sysGlobal==2) sysmode = "_noBRUnc_addedCUTandJES";        //for E ratio
    if(type ==2) sysmode = "_PDF_noneGlobal_addedCUTandJES";

    canvas->SaveAs(outPlotDir+Form("/finalSpectra_%s.png",sysmode.Data()));

    TFile *ofile = new TFile(Form("%s/JetPtSpectrum_final%s.root",outSpectraDir.Data(),sysmode.Data()),"RECREATE");
    hData_binned->Write();
    hDataSys->Write();
    if(separateCUTUnc && gDataCutSys) gDataCutSys->Write("gDataCutSys");
    hData_binned_ratio->Write();
    hDataSysRatio->Write();
    if(ePowhegPythia6){
        simPowhegPythia6_cent->Write();
        simPowhegPythia6_up->Write();
        simPowhegPythia6_down->Write();
        simPowhegPythia6->Write();
        simPowhegPythia6_cent_R->Write();
        simPowhegPythia6_up_R->Write();
        simPowhegPythia6_down_R->Write();
        simPowhegPythia6_R->Write();
        for(Int_t ivar = 1; ivar < 9; ivar++){
            simPowhegPythia6var[ivar]->Write();
        }
    }

    if(ePowhegPythia8){
        simPowhegPythia8_cent->Write();
        simPowhegPythia8_up->Write();
        simPowhegPythia8_down->Write();
        simPowhegPythia8->Write();
        simPowhegPythia8_cent_R->Write();
        simPowhegPythia8_up_R->Write();
        simPowhegPythia8_down_R->Write();
        simPowhegPythia8_R->Write();
        for(Int_t ivar = 1; ivar < 9; ivar++){
            simPowhegPythia8var[ivar]->Write();
        }
    }

    if(ePythia6){
        simPythia6_cent->Write();
        simPythia6_up->Write();
        simPythia6_down->Write();
        simPythia6->Write();
        simPythia6_cent_R->Write();
        simPythia6_up_R->Write();
        simPythia6_down_R->Write();
        simPythia6_R->Write();
    }

    if(ePythia8){
        simPythia8_cent->Write();
        simPythia8_up->Write();
        simPythia8_down->Write();
        simPythia8->Write();
        simPythia8_cent_R->Write();
        simPythia8_up_R->Write();
        simPythia8_down_R->Write();
        simPythia8_R->Write();
    }

    if(ePythia8SoftMode2){
        simPythia8Soft2_cent->Write();
        simPythia8Soft2_up->Write();
        simPythia8Soft2_down->Write();
        simPythia8Soft2->Write();
        simPythia8Soft2_cent_R->Write();
        simPythia8Soft2_up_R->Write();
        simPythia8Soft2_down_R->Write();
        simPythia8Soft2_R->Write();
    }

    if(ePowhegPythia6dijet){
        simPowhegPythia6dijet_cent->Write();
        simPowhegPythia6dijet_up->Write();
        simPowhegPythia6dijet_down->Write();
        simPowhegPythia6dijet->Write();
        simPowhegPythia6dijet_cent_R->Write();
        simPowhegPythia6dijet_up_R->Write();
        simPowhegPythia6dijet_down_R->Write();
        simPowhegPythia6dijet_R->Write();
        for(Int_t ivar = 1; ivar < 9; ivar++){
            simPowhegPythia6dijetvar[ivar]->Write();
        }
    }
std::cout<<"a"<<std::endl;
    if(ePowhegPythia8dijet){
        simPowhegPythia8dijet_cent->Write();
      //  simPowhegPythia8dijet_up->Write();
      //  simPowhegPythia8dijet_down->Write();
        simPowhegPythia8dijet->Write();
        simPowhegPythia8dijet_cent_R->Write();
     //   simPowhegPythia8dijet_up_R->Write();
     //   simPowhegPythia8dijet_down_R->Write();
        simPowhegPythia8dijet_R->Write();
     //   for(Int_t ivar = 1; ivar < 9; ivar++){
    //        simPowhegPythia8dijetvar[ivar]->Write();
     //   }
    }
std::cout<<"a"<<std::endl;
    ofile->Close();
    return;
}

void TerminateCanvas(TPad* pad1,TPad* pad2,TH1D* histo1,TH1D* histo2){
    pad1->cd();
    histo1->Draw("sameaxis");
    pad2->cd();
    histo2->Draw("sameaxis");
    histo2->Draw("sameaxig");
}

void PlaceOnPadSim(TPad* pad,TGraphAsymmErrors* histo, Color_t ci, Style_t style, Style_t linestyle){
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
    TString odraaw = "";
    if(linestyle ==1) odraaw = "2p";
    else{
        histo->SetLineWidth(3);
        histo->SetMarkerSize(0);
        odraaw = "E";
    }
    histo->Draw(odraaw);
}

void PlaceOnPadData(TPad* pad,TGraphAsymmErrors* histo1, TH1D* histo2, Size_t markersize){
    pad->cd();
    Color_t ci = static_cast<Color_t>(TColor::GetColor("#990000"));
    ci = kBlack;
    histo1->SetLineColor(ci);
    histo1->SetMarkerColor(ci);
    histo1->SetMarkerStyle(20);
    histo1->SetMarkerSize(markersize);
    ci = static_cast<Color_t>(TColor::GetColor("#cccccc"));
    histo1->SetFillColor(ci);
    histo1->SetLineColor(ci);
    histo1->Draw("2p");
    //data central w stat. unc.
    //ci = static_cast<Color_t>(TColor::GetColor("#990000"));
    ci = kBlack;
    histo2->SetLineColor(ci);
    histo2->SetMarkerColor(ci);
    histo2->SetMarkerStyle(20);
    histo2->SetMarkerSize(markersize);
    histo2->Draw("same p  e0 x0");
};

std::tuple<TCanvas*, TPad*, TPad*, TH1D*, TH1D*> PrepareCanvas(UInt_t xAxisBins, Double_t *xAxis){
    style();
    //prepare main canvas
    TCanvas *FinalSpectrum = new TCanvas("FinalSpectrum", "FinalSpectrum",0,45,700,700);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    FinalSpectrum->SetHighLightColor(2);
    FinalSpectrum->Range(0,0,1,1);
    FinalSpectrum->SetFillColor(0);
    FinalSpectrum->SetBorderMode(0);
    FinalSpectrum->SetBorderSize(2);
    FinalSpectrum->SetFrameBorderMode(0);

    //Set primitives in upper pad
    TPad *pad_top = new TPad("pad_top", "pad_top",0,0.35,1,1);
    pad_top->Draw();
    pad_top->cd();
    pad_top->Range(-1.986821e-07,-4.69897,33.33333,0.3499945);
    pad_top->SetFillColor(0);
    pad_top->SetBorderMode(0);
    pad_top->SetBorderSize(2);
    if(zBin==0 || !pdf)pad_top->SetLogy();
    pad_top->SetTickx(1);
    pad_top->SetTicky(1);
    pad_top->SetLeftMargin(0.18f);
    pad_top->SetBottomMargin(0);
    pad_top->SetFrameBorderMode(0);
    pad_top->SetFrameBorderMode(0);

    FinalSpectrum->cd();

    //Set primitives in bottom pad
    TPad *pad_bottom = new TPad("pad_bottom", "pad_bottom",0,0,1,0.35);
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
    pad_bottom->SetBottomMargin(0.27f);
    pad_bottom->SetFrameBorderMode(0);
    pad_bottom->SetFrameBorderMode(0);

    pad_top->cd();
    TH1D *hEmpty_up = new TH1D("hEmpty_up","Central Values",static_cast<Int_t>(xAxisBins), xAxis);
    hEmpty_up->SetMinimum(plotRanges[2]);
    hEmpty_up->SetMaximum(plotRanges[3]);
    hEmpty_up->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
    hEmpty_up->GetXaxis()->SetLabelFont(43);
    hEmpty_up->GetXaxis()->SetLabelSize(0.035f);
    hEmpty_up->GetXaxis()->SetTitleSize(0.035f);
    hEmpty_up->GetXaxis()->SetTitleFont(42);
    if(zBin==0)hEmpty_up->GetYaxis()->SetTitle("#frac{d^{2}#it{#sigma}}{d#it{p}_{T}d#it{#eta}} [mb (GeV/it{c})^{#minus1}]");
    else if(!pdf)hEmpty_up->GetYaxis()->SetTitle("#frac{d^{2}#it{#sigma}}{d#it{z}_{#parallel}d#it{#eta}} (mb)");
    else hEmpty_up->GetYaxis()->SetTitle("Probability density");
    hEmpty_up->GetYaxis()->SetLabelFont(43);
    hEmpty_up->GetYaxis()->SetLabelSize(22);
    hEmpty_up->GetYaxis()->SetTitleSize(26);
    hEmpty_up->GetYaxis()->SetLabelOffset(0.015f);
    hEmpty_up->GetYaxis()->SetTitleOffset(2.f);
    hEmpty_up->GetYaxis()->SetTitleFont(43);
    hEmpty_up->GetYaxis()->SetDecimals();
    hEmpty_up->Draw("axis");
    pad_bottom->cd();
    TH1D *hEmpty_bottom = new TH1D("hEmpty_bottom","Central Values",static_cast<Int_t>(xAxisBins), xAxis);
    hEmpty_bottom->SetMinimum(plotRanges[0]);
    hEmpty_bottom->SetMaximum(plotRanges[1]);
    if(zBin ==0)hEmpty_bottom->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
    else hEmpty_bottom->GetXaxis()->SetTitle("z_{#parallel}");
    hEmpty_bottom->GetXaxis()->SetLabelFont(43);
    hEmpty_bottom->GetXaxis()->SetLabelSize(22);
    hEmpty_bottom->GetXaxis()->SetTitleSize(26);
    hEmpty_bottom->GetXaxis()->SetTitleOffset(3.f);
    hEmpty_bottom->GetXaxis()->SetTitleFont(43);
    hEmpty_bottom->GetYaxis()->SetTitle("theory / data");
    hEmpty_bottom->GetYaxis()->SetNdivisions(509);
    hEmpty_bottom->GetYaxis()->CenterTitle();
    hEmpty_bottom->GetYaxis()->SetDecimals();
    hEmpty_bottom->GetYaxis()->SetLabelOffset(0.015f);
    hEmpty_bottom->GetXaxis()->SetLabelOffset(0.02f);
    hEmpty_bottom->GetYaxis()->SetLabelFont(43);
    hEmpty_bottom->GetYaxis()->SetLabelSize(22);
    hEmpty_bottom->GetYaxis()->SetTitleSize(26);
    hEmpty_bottom->GetYaxis()->SetTitleOffset(2.f);
    hEmpty_bottom->GetYaxis()->SetTitleFont(43);
    hEmpty_bottom->Draw("axis");

    return std::make_tuple(FinalSpectrum,pad_top, pad_bottom, hEmpty_up,hEmpty_bottom);

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
         // std::cout<<j<<" "<<xval[j]<<" "<<xvalwidth[j]<<" "<<value[j]<<" "<<error<<std::endl;
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


std::tuple<TGraphAsymmErrors*, TH1D*, TH1D*, TH1D*, TH1D**> GetSim(TString simname, Int_t type, Int_t nFiles, TString *names, TString simDir, Double_t simScaling, UInt_t nBins, Double_t *xBins){
    TH1D **hPrompt = new TH1D*[static_cast<ULong_t>(nFiles)];
    TH1D **hPrompt_binned = new TH1D*[static_cast<ULong_t>(nFiles)];
    //TH1D *hPrompt_binned[static_cast<ULong_t>(nFiles)];
    //std::cout<<"type "<<type<<std::endl;
    TString tmp = "h";
    tmp+=simname;

    for (int nr=0; nr<nFiles; nr++){
        TString file = simDir;
        file += "/JetPt_";
        file += names[nr];
        if(type ==0){
            file += "_Dpt"; file += 2; file += "_"; file += 36;
        }
        if(type ==1 || type ==2){
            file += "_Dpt"; file += DpTbins[0][zBin-1]; file += "_"; file += DpTbins[1][zBin-1];
            file += "_Jetpt"; file += jetpTbins[zBin-1]; file += "_"; file += jetpTbins[zBin];
        }
        file += "_Dzero";
        file += ".root";
        TH1D *htmp = nullptr;
        if(type ==0){
            htmp = dynamic_cast<TH1D*>(GetInputHist(file, "hPt"));
            if(!htmp) std::cout<<"no histogram in file: "<<file<<std::endl;
            htmp->GetYaxis()->SetTitle("d#sigma/dp_{T} (mb)");
        }
        if(type ==1 || type ==2){
            htmp = dynamic_cast<TH1D*>(GetInputHist(file, "hz"));
            if(!htmp) std::cout<<"no histogram in file: "<<file<<std::endl;
            htmp->GetYaxis()->SetTitle("d#sigma/dz (mb)");
        }
        hPrompt[nr] = dynamic_cast<TH1D*>(htmp->Clone(Form("hPrompt_%d",nr)));
        hPrompt_binned[nr] = dynamic_cast<TH1D*>(htmp->Rebin(static_cast<Int_t>(nBins),Form("hPrompt_binned_%d",nr),xBins));
        hPrompt_binned[nr]->SetName(tmp+"_var"+Form("%d",nr));
        if(pdf){
            hPrompt_binned[nr]->Scale(1./hPrompt_binned[nr]->Integral());
            hPrompt_binned[nr]->Scale(1,"width");
        }
    }


    TH1D *hPrompt_central_binned = dynamic_cast<TH1D*>(hPrompt_binned[0]->Clone(tmp+"_central"));
    setHistoDetails(hPrompt_central_binned,4,24);

    // get up unc
    TH1D * hPrompt_up = GetUpSys(hPrompt_binned,nFiles);
    hPrompt_up->SetName(tmp+"_up");
    setHistoDetails(hPrompt_up,4,24,0,2,2);
    // get down unc
    TH1D *hPrompt_down = GetDownSys(hPrompt_binned,nFiles);
    hPrompt_down->SetName(tmp+"_down");
    setHistoDetails(hPrompt_down,4,24,0,2,2);

    if(!pdf){
        hPrompt_central_binned->Scale(simScaling);
        hPrompt_central_binned->Scale(1,"width");
        hPrompt_central_binned->Scale(1./dy);//2*jetEta;
        hPrompt_up->Scale(simScaling);
        hPrompt_up->Scale(1,"width");
        hPrompt_up->Scale(1./dy);
        hPrompt_down->Scale(simScaling);
        hPrompt_down->Scale(1,"width");
        hPrompt_down->Scale(1./dy);
        for (int nr=0; nr<nFiles; nr++){
            hPrompt_binned[nr]->Scale(simScaling);
            hPrompt_binned[nr]->Scale(1,"width");
            hPrompt_binned[nr]->Scale(1./dy);
        }
    }

    Double_t *xval = new Double_t[nBins];
    Double_t *xvalwidth = new Double_t[nBins];
    Double_t *valuetheory = new Double_t[nBins];
    Double_t *valuetheoryerrup = new Double_t[nBins];
    Double_t *valuetheoryerrdown = new Double_t[nBins];
    for(UInt_t j=0; j<nBins; j++){
          xval[j] = (xBins[j]+xBins[j+1]) / 2.;
          xvalwidth[j] = (xBins[j+1]-xBins[j]) / 2.;
          valuetheory[j] = hPrompt_central_binned->GetBinContent(hPrompt_central_binned->GetXaxis()->FindBin(xval[j]));
          valuetheoryerrup[j] = hPrompt_up->GetBinContent(hPrompt_up->GetXaxis()->FindBin(xval[j])) - valuetheory[j];
          valuetheoryerrdown[j] = valuetheory[j] - hPrompt_down->GetBinContent(hPrompt_up->GetXaxis()->FindBin(xval[j]));

    }
    TGraphAsymmErrors *grsystheory = new TGraphAsymmErrors(static_cast<Int_t>(nBins),xval,valuetheory,xvalwidth,xvalwidth,valuetheoryerrdown,valuetheoryerrup);
    tmp = "hae";
    tmp+=simname;
    grsystheory->SetName(tmp);
    return std::make_tuple(grsystheory, hPrompt_central_binned, hPrompt_up,hPrompt_down,hPrompt_binned);
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

TH1D* GetUpSys(TH1D **hh, Int_t nFiles){
    double max = 0, maxerr = 0;
    TString name = hh[0]->GetName();
    TH1D *hh_up = dynamic_cast<TH1D*>(hh[0]->Clone(name +"_up"));

    for(int iBin=1; iBin<hh[0]->GetNbinsX()+1; iBin++ ){
      //  std::cout<<"iBin: "<<iBin<<std::endl;
        max = hh[0]->GetBinContent(iBin);
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

TH1D* GetDownSys(TH1D **hh, Int_t nFiles){
    double max = 0, maxerr = 0;
    TString name = hh[0]->GetName();
    TH1D *hh_down = dynamic_cast<TH1D*>(hh[0]->Clone(name +"_down"));

    for(int iBin=1; iBin<hh[0]->GetNbinsX()+1; iBin++ ){
        max = hh[0]->GetBinContent(iBin);
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
