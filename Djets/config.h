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

//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------

    int fSystem = 1; // 0: pp, 1: p-Pb
    TString fSystemS = "p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV";

    const int Rpar = 3;
    const double fRpar = 0.3;

    enum DMesonSpecies {kD0toKpi=0, kDStarD0pi};
    DMesonSpecies fDmesonSpecie = 0;
    TString fDmesonDsc="Dzero";
    TString fDmesonS="D^{0}";

    //------- D pT bins
  //const int fptbinsDN = 11;
  //double fptbinsDA[fptbinsDN+1] = { 2,3,4,5,6,7,8,10,12,16,24,36 };

  const int fptbinsDN = 10;
  double fptbinsDA[fptbinsDN+1] = { 3,4,5,6,7,8,10,12,16,24,36 };

  // for Pb-Pb baseline
  //  const int fptbinsDN = 9;
  //  double fptbinsDA[fptbinsDN+1] = { 3,4,5,6,7,8,10,12,16,20 };
  

    double zmin = -2, zmax = 2.;

    //------- signal extraction config
    Bool_t  fUseRefl = 1;
    Int_t   fbkgtype = 0; // kExpo=0, kLin=1, kPol2=2, kNoBk=3, kPow=4, kPowEx=5
    Float_t fsigmaSignal = 3;
    Float_t fsigmaBkg[] = {-9,-4,4,9};
    Float_t fDmass = 1.864, fDsigma = 0.010;
    Float_t minf = 1.71, maxf = 2.1;
    Int_t   fRebinMass = 2;

    //Int_t fColors[] = {1,2,8,4,kOrange-1,6,kGray+1,kCyan+1,kMagenta+2,kViolet+5,kYellow+2};
    Int_t fColors[] = {1,2,8,4,kOrange-1,6,kGray+1,kCyan+1,kMagenta+2,kGreen+3,kViolet+5,kYellow+2};
    Int_t fMarkers[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};

    const double sigma_in = 2.09; // in bars
    const double BRDstar = 0.0257;
    const double BRDzero = 0.0393;
    const int APb = 208;

    //------- POWHEG simulations

     // New simulations, with R=0.3,0.4 (Dzero and Dstar)
    const Int_t fBsimN = 10;
    TString fRunB[] = {
    "AnalysisResults_FastSim_powheg+pythia6_beauty_1519042424",
    "AnalysisResults_FastSim_powheg+pythia6_beauty_1519293768",
    "AnalysisResults_FastSim_powheg+pythia6_beauty_1519292550",
    "AnalysisResults_FastSim_powheg+pythia6_beauty_1519275692",
    "AnalysisResults_FastSim_powheg+pythia6_beauty_1519241992",
    "AnalysisResults_FastSim_powheg+pythia6_beauty_1519127831",
    "AnalysisResults_FastSim_powheg+pythia6_beauty_1519049991",
    "AnalysisResults_FastSim_powheg+pythia6_beauty_1519079649",
    "AnalysisResults_FastSim_powheg+pythia6_beauty_1519124826",
    "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1520259013"
  };

    TString fDescB[] = {
    "central",
    "m_{b}=4.5",
    "m_{b}=5",
    "muR=2,muF=2",
    "muR=1,muF=2",
    "muR=2,muF=1",
    "muR=0.5,muF=0.5",
    "muR=1,muF=0.5" ,
    "muR=0.5,muF=1",
    "EvtGen" };

    const Int_t fCsimN = 9;
    TString fRunC[] = {
    "AnalysisResults_FastSim_powheg+pythia6_charm_1520422975",
    "AnalysisResults_FastSim_powheg+pythia6_charm_1520552390",
    "AnalysisResults_FastSim_powheg+pythia6_charm_1520532995",
    "AnalysisResults_FastSim_powheg+pythia6_charm_1520531695",
    "AnalysisResults_FastSim_powheg+pythia6_charm_1520492126",
    "AnalysisResults_FastSim_powheg+pythia6_charm_1520491792",
    "AnalysisResults_FastSim_powheg+pythia6_charm_1520423929",
    "AnalysisResults_FastSim_powheg+pythia6_charm_1520464344",
    "AnalysisResults_FastSim_powheg+pythia6_charm_1520464810"
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
    "muR=0.5,muF=1"};


//// PromptSimDirOut
/*
\subsection*{Charm simulation: jobs on grid}
\begin{table}[H]
\begin{tabular}{l l l r}
\toprule
\textbf{Sl} & \textbf{Variation type} & \textbf{Run number, grid ok}& \textbf{Merge-download}\\
\toprule
01 & central                    & 1520422975 y\\
02 & $\mu_F$ 0.5  $\mu_R$ 0.5   & 1520423929 y\\
03 & $\mu_F$ 0.5  $\mu_R$ 1.0   & 1520464344 y\\
04 & $\mu_F$ 1.0  $\mu_R$ 0.5   & 1520464810 y\\
05 & $\mu_F$ 1.0  $\mu_R$ 2.0   & 1520491792 y\\
06 & $\mu_F$ 2.0  $\mu_R$ 1.0   & 1520492126 y\\
07 & $\mu_F$ 2.0  $\mu_R$ 2.0   & 1520531695 y\\
08 & mass high  1.7                      & 1520532995 y\\
09 & mass low   1.3		              & 1520552390 y\\
\bottomrule
\end{tabular}
\caption{Run numbers for c simulation submitted on grid}
\label{tab:c-simulation}
\end{table}
*/

    // old simulations
  /*  const Int_t fBsimN = 9;
    TString fRunB[] = {
    "AnalysisResults_FastSim_powheg_beauty_1496995621",
    "AnalysisResults_FastSim_powheg_beauty_1497227464",
    "AnalysisResults_FastSim_powheg_beauty_1497228262",
    "AnalysisResults_FastSim_powheg_beauty_1497207414",
    "AnalysisResults_FastSim_powheg_beauty_1497173624",
    "AnalysisResults_FastSim_powheg_beauty_1497172782",
    "AnalysisResults_FastSim_powheg_beauty_1497132057",
    "AnalysisResults_FastSim_powheg_beauty_1497130041",
    "AnalysisResults_FastSim_powheg_beauty_1497121734"};*/

TString OUTDIRECTORY="/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_CutVarBase";
const int ND = 2;
const int NDMC = 2;
const int fptbinsJetTrueN = 7;
double fptbinsJetTrueA[fptbinsJetTrueN+1] = { 3,5,10,15,20,25,35,50 };
const int fptbinsJetMeasN = 7;
double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 3,5,10,15,20,25,35,50 };
