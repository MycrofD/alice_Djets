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

TString roounfoldpwd = "$HOME/ALICE_HeavyFlavour/RooUnfold-1.1.1/libRooUnfold";

enum DMesonSpecies {kD0toKpi, kDStarD0pi};

TString OUTDIRECTORY="/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_cutTight/DzeroR03_def_437_old0";
    // ========================== Prepare your config ============================================

    const int     fJetptbinsN = 5;
    double        fJetptbinsA[fJetptbinsN+1] = {2.0, 5.0, 7.0, 10.0, 15.0, 50.0};
    //this ---| is R dependent. written separately below.
    //        V
    //double fDptRangesA[] = {2,2,3,5,5};//these are Dpt mins jetpt bins, and 36 is the max in last bin
    double fDptRangesAUp[] = {5,7,10,15,36};//these are Dpt mins jetpt bins, and 36 is the max in last bin
    double fDptRespA[] = {2,3,4,5,6,7,8,10,12,16,24,36};// these are Dpt bins for response
    const int     fDptRespN = sizeof(fDptRespA)/sizeof(fDptRespA[0])-1;
    //
    int           fSystem = 0;            //-----! 0: pp, 1: p-Pb, Pb-Pb -- set up system
    TString       fSystemS = "pp, #sqrt{#it{s}} = 5.02 TeV";
    DMesonSpecies fDmesonSpecie = kD0toKpi;
    TString       fDmesonDsc = "Dzero";
    TString       fDmesonS = "D^{0}";
    const int     ND = 4;                //-----!  change these numbers based on how many D mesons you analyse in data !
    const int     NDMC = 2;              //-----!  change these numbers based on how many D mesons you analyse in MC !

    //const double  sigma_in = 0.0512;       //-----! inelastic x-section in bars
    const double  sigma_in = 0.05077;       //-----! inelastic x-section in bars
    const double  nEvScale = 1.063;      //-----! scaling factor from the number of selected events to the number of events needed for the nomrmalization, taken from the D2H normalization counter
    const double  BRDstar = 0.0257;
    const double  BRDzero = 0.0389;
    const int     APb = 208;
