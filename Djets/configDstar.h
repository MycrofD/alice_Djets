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
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TDatabasePDG.h"
#include "TNtuple.h"

//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------

  enum DMesonSpecies {kD0toKpi=0, kDStarD0pi};
  DMesonSpecies fDmesonSpecie = 1;

    TString fDmesonDsc="Dstar";
    TString fDmesonS="D*";

    Bool_t  fUseRefl = 0; // no reflections for Dstar
    Int_t   fbkgtype = 5; // kExpo=0, kLin=1, kPol2=2, kNoBk=3, kPow=4, kPowEx=5
    Float_t fsigmaSignal = 3;
    Float_t fsigmaBkg[] = {-8,-5,5,13};
    Float_t fDmass = 0.1455, fDsigma = 0.00065;
    Float_t minf = 0.140, maxf = 0.155;
    Int_t   fRebinMass = 1;

    Int_t   fzMeasAxis = 0;
    Int_t   fJetPtMeasAxis = 1;
    Int_t

    //------- D pT bins
    //const int fptbinsDN = 12;
    //double fptbinsDA[fptbinsDN+1] = { 1,2,3,4,5,6,7,8,10,12,16,24,36 };
    const int fptbinsDN = 11;
    double fptbinsDA[fptbinsDN+1] = { 2,3,4,5,6,7,8,10,12,16,24,36 };

    //------- jet pT bins
    const int fptbinsJetMeasN = 10;
    double fptbinsJetMeasA[fptbinsJetMeasN+1] = { 2,3,4,5,6,8,10,14,20,30,50 };
    const int fptbinsJetTrueN = 7;
    double fptbinsJetTrueA[fptbinsJetMeasN+1] = { 4,5,6,8,10,14,20,30 };

    double zmin = -2, zmax = 2.;

    const int Rpar = 4;
