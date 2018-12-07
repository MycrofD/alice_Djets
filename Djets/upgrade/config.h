//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------
#include <iostream>
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
//#include "AliHFInvMassFitter.h"
/*
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBinByBin.h"*/


    TString OUTDIRECTORY="/home/basia/Work/alice/analysis/upgradeProjections/results_Bsys/";
    TString OUTDIRECTORYPLOTS="/home/basia/Work/alice/analysis/upgradeProjections/results_Bsys/plots/";

    TString fFileSignal_prompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_charm";
    TString fFileSignal_nonprompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_beauty";
    TString fFileSignalBkg_prompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_charmBkg";
    TString fFileBkg = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_Hijing.root";

    TString fFileEff_prompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_charm";
    TString fFileEff_nprompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_beauty";

    TString fFilePowhegPrompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/POWHEG_charm_central";
    TString fFilePowhegPrompt_sys = "/home/basia/Work/alice/analysis/upgradeProjections/out/POWHEG_beautyEvGen_muR05muF1";
    TString fFilePowhegNonPrompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/POWHEG_beautyEvGen_central";

	TFile 			*fOutFile;
    TString         fSMBBin = "_1.root";
    Bool_t          fIsHBEffWeighting = 1;
    Double_t        fHBrejection = 2;
    TString         fOutNamePostfix = "";
    Int_t           findex = 2;
    
    Bool_t 			fBRaaSys = 0;


    //TString fFilePowhegPrompt = "/media/basia/Disk2/Work/Djets/POWHEGSimulations/fastSim_pp5TeV/RAW_CHARM_POWHEG.root";
    //TString fFilePowhegNonPrompt = "/media/basia/Disk2/Work/Djets/POWHEGSimulations/fastSim_pp5TeV/AnalysisResults_FastSim_powheg+pythia6_beauty_150593961473.root";

    TH1D *fHEff_prompt, *fHEff_nprompt;
    TH1D *fHPowhegDPrompt, *fHPowhegDNonPrompt;
    TH1D *fHPowhegPrompt, *fHPowhegPromptEffScaled, *fHPowhegPromptEffScaledRaa, *fHPowhegPromptRaa;
    TH1D *fHPowhegNonPrompt, *fHPowhegNonPromptEffScaled, *fHPowhegNonPromptEffScaledRaa, *fHPowhegNonPromptEffRatioScaled, *fHPowhegNonPromptEffRatioScaledRaa, *fHPowhegNonPromptRaa;
    TH1D *fHPowhegNonPromptEffScaledRaa_powhegSys, *fHPowhegNonPromptEffScaledRaa_RaaSys;
    TH1D *fHJetPtSpectrum, *fHJetPtSpectrumUnc, *fHJetPtSpectrumPrompt;
    TH1D *fHJetPtSpectrumUnf, *fHJetPtSpectrumUncUnf, *fHJetPtSpectrumPromptUnf;
    TH1D *fHSB;
    TF1  *fSBfun;
    TH2D *fHDettMatrix;

    const double    fDataLum = 10e9;      // Luminosity 10nb-1
    const double    fSigmaIn = 0.7750;       // Sigma inelastic in bars (7.775b), for 10% centrality
    double          fDataEv = fDataLum*fSigmaIn;         // Number of events
    const double    BRDstar = 0.0257;
    const double    BRDzero = 0.0389;
    const int       APb = 208;
    const double    Taa = 23.2; // in mb-1, for 0-10%
    double          fSimScalingC = BRDzero*Taa;
    double          fSimScalingB = BRDzero*Taa;
    double          *fRaaC, *fRaaB, *fRaaBSys;
    double          *fSB;

    const int       fHardBinsN = 6;
    double          fHBWeightsC[5], fHBWeightsB[5];
    int             fHardBinsMin[5] = {5,11,21,36,57};
    int             fHardBinsMax[5] = {11,21,36,57,200};
    bool            fIsHBins;

    //====== D pT bins ---- set up your D pT bins ======
    //const int       fPtbinsDN = 12;
    //double          fPtbinsDA[fPtbinsDN+1] = { 1,2,3,4,5,6,8,10,12,16,24,36,50 };
    const int       fPtbinsDN = 10;
    double          fPtbinsDA[fPtbinsDN+1] = { 3,4,5,6,8,10,12,16,24,36,50 };
    //const int       fPtbinsDN = 20;
    //double          fPtbinsDA[fPtbinsDN+1] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,24,30,36,50 };
  //  const int     fPtbinsDN = 9;
//    double        fPtbinsDA[fPtbinsDN+1] = { 1,2,3,4,5,6,7,8,10,24};
    
    //====== jet pT bins ---- set up your jet pT bins ======
    //const int       fPtbinsJetTrueN = 10;
    //double          fPtbinsJetTrueA[fPtbinsJetTrueN+1] = { 3,5,6,8,10,14,20,30,50,80,100 };
    const int       fPtbinsJetTrueN = 15;
    double          fPtbinsJetTrueA[fPtbinsJetTrueN+1] = { 3,5,6,8,10,12,14,16,20,25,30,40,50,60,80,100 };
    const int       fPtbinsJetMeasN = 15;
    double          fPtbinsJetMeasA[fPtbinsJetMeasN+1] = { 3,5,6,8,10,12,14,16,20,25,30,40,50,60,80,100 };
    
    const int       fPtbinsJetFinalN = 9;
    double          fPtbinsJetFinalA[fPtbinsJetFinalN+1] = { 3,4,5,6,8,10,14,20,30,50 };
    
    //====== z range ---- set up your min, max z ======
    double          fZmin = -2, fZmax = 2.; // for D-jet pT spectrum
    double          fJetptmin = 0, fJetptmax = 100;

    // ============ Prepare your config ==============
    int           fSystem = 0;            //-----! 0: pp, 1: p-Pb, Pb-Pb -- set up system
    TString       fSystemS = "pp, #sqrt{#it{s}} = 5.02 TeV";
    Int_t         fDmesonSpecie = 0;
    TString       fDmesonDsc = "Dzero";
    TString       fDmesonS = "D^{0}";
    const double  fRpar = 0.3;           //-----! jet R parameter for your studies (the one that you use in your jet finder!)
    Double_t      fJetEta = 0.6;


    //====== signal extraction config ======
    Bool_t        fUseRefl = 1;                      //-----! if to use reflections (for D0, you must have reflections files ready)
    Int_t         fBkgtype = 0;                      //-----! kExpo=0, kLin=1, kPol2=2, kNoBk=3, kPow=4, kPowEx=5
    Float_t       fSigmaSignal = 2;                  //-----! sigma for the signal region
    Float_t       fSigmaBkg[] = {-9,-4,4,9};         //-----! sigma for the SB region (both left and right side from the fit)
    Float_t       fDmass = 1.864, fDsigma = 0.010;   //-----! initial values for D mass and sigma
    Float_t       minf = 1.71, maxf = 2.1;           //-----! min/mass of the inv. mass distributions
    Int_t         fRebinMass = 2;                    //-----! rebining of the inv. mass distributions


	Int_t colors2[] = {kBlue+1,kGreen+2,kOrange+1,kViolet+1,kGray+2,kMagenta+2,kYellow+2,kCyan+1,4,6,8};
	Int_t markers2[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};
	Int_t linestyle2[] = {1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};


    Int_t fColors[] = {1,2,kViolet+5,kGreen+2,kBlue+2,kOrange-1,6,kGray+1,kCyan+1,kMagenta+2,kGreen-1,kYellow+2};
    Int_t fMarkers[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};


/* Raa_C
1.600710299340438, 0.9361111111111113
2.5680703534584826, 0.8472222222222225
3.5214780991036694, 0.5750000000000004
4.542322002367664, 0.38888888888888906
5.445416878065279, 0.255555555555556
6.478521900896329, 0.23055555555555607
7.3913411128023, 0.2250000000000001
8.973025536952477, 0.20833333333333393
11.043252156265854, 0.21111111111111125
14.149543378995428, 0.22777777777777786
20.058346017250123, 0.2694444444444448
29.98752748181971, 0.33888888888888924
43.14265178420428, 0.39722222222222237
*/
/*
1.490693739424704, 1.0557873253440786
2.4720812182741128, 1.0929844575854308
3.506627185561195, 1.042076159575767
4.494829855235947, 0.9115762787749808
5.5171084790374145, 0.7845157313889322
7.050526414739613, 0.6111070408940522
8.992855799962399, 0.5133346125444087
11.003337093438615, 0.46539617748472195
14.002021056589586, 0.44675167605723876
20.20384470765181, 0.4713357870927599
29.745111863132166, 0.5151051569650134
39.86567023876668, 0.5898517859530994
*/




    ///============== POWHEG simulations ============================
    //======= set up here names of your simulation files =======

    TString fRunB[] = {
      "AnalysisResults_FastSim_powheg+pythia6_beauty_150593961473",
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504284947",
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504259653",
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1506803374"
    };

  /*  TString fRunB[] = {
      "AnalysisResults_FastSim_powheg+pythia6_beauty_150593961473",	//central
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504284947",	//R1F05
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504259653",	//R05F1
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1506803374",	//R2F1
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504296768",	//R1F2
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504212024",	//R05F05
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504318569",	//R2F2
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504202511",	//mb5
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504197966",	//mb4.5
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504197460",	//pdf 21200
      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504199953"	//pdf 10800
    };
*/

    TString fDescB[] = {
      "central",
      "muR=1,muF=0.5" ,
      "muR=0.5,muF=1",
      "muR=2,muF=1"
    };

/*
    TString fDescB[] = {
      "central",
      "muR=1,muF=0.5" ,
      "muR=0.5,muF=1",
      "muR=2,muF=1",
      "muR=1,muF=2",
      "muR=0.5,muF=0.5",
      "muR=2,muF=2",
      "m_{b}=5",
      "m_{b}=4.5",
      "PDF 21200",
      "PDF 10800"
    };
*/
    TString fRunC[] = {
      "RAW_CHARM_POWHEG"
    };
    TString fDescC[] = {
      "central",
      "m_{c}=1.3",
      "m_{c}=1.7",
      "muR=2,muF=2",
      "muR=1,muF=2",
      "muR=2,muF=1",
      "muR=0.5,muF=0.5",
      "muR=1,muF=0.5",
      "muR=0.5,muF=1"
    };
const Int_t fBsimN = 4;
const Int_t fCsimN = 1;
