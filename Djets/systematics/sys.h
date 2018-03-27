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

   
int nJetBins;
double *ptJetbins;
double ptbins[] = { 3,4,5,6,8,10,14,20,30,50 };
double ptbinstrue[] = { 5,6,8,10,14,20,30,50 };

Int_t colors[] = {1,2,4,6,kOrange-1,kGray+1,kCyan+1,kMagenta+2,kGreen+3,kViolet+5,kYellow+2,8};
Int_t markers[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};
Int_t linestyle[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};

  TPaveText *pvEn;
    TPaveText *pvJet;
    TPaveText *pvD;
    TPaveText *pvEta;
    TPaveText *pv3;

    TPaveText *pvJet2;
    TPaveText *pvD2;
    TPaveText *pvEta2;
    TPaveText *pvJetPt2;
