//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//  [Modified] A.Mohanty
//  Utrecht University
//  auro.mohanty@cern.ch
//-----------------------------------------------------------------------

#include "../DsignalExtraction/configDzero_ppz.h"

const int fptbinsZGenN=10, fJetptbinsGenN=9;
double fptbinsZGenA[fptbinsZGenN+1]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.02};
double fJetptbinsGenA[fJetptbinsGenN+1]={0,2,3,4,5,7,10,15,50,60};
const int nDim    = 5;//the five dimensions
double bins[nDim] = {fptbinsZGenN,fJetptbinsGenN,fptbinsZGenN,fJetptbinsGenN,fDptRespN};//for creating 4D ResponseMatrix
int zRec = 0, jetRec = 1, DRec = 2;
int zGen = 5, jetGen = 6;
int dim[nDim]     = {zRec,jetRec,zGen,jetGen,DRec};//for extacting 5D info from THnSparse

double jetmin = 0, jetmax = 60;

void DetRM_zjet(
bool isPrompt = 0,
TString datafile = ".root",
TString outDir = "/ResponseMatrix", //"plots",
bool postfix = 0,
TString listName = "FD",
bool isprefix=0 )
{

    gStyle->SetOptStat(0000); //Mean and RMS shown
    gStyle->SetPadRightMargin(0.1);
    gSystem->Exec(Form("mkdir %s",outDir.Data()));
    gSystem->Exec(Form("mkdir %s/plots",outDir.Data()));

    TFile *File = new TFile(datafile,"read");
    TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
    TString histName;
    if(!isprefix){
            if(fDmesonSpecie) histName = "histosDStarMBN";
            else histName = "histosD0MBN";}
    else{
            if(fDmesonSpecie) histName = "histosDStarMBN";
            else histName = "histosD0";}

    TPaveText *pv2 = new TPaveText(0.15,0.8,0.25,0.9,"brNDC");
    pv2->SetFillStyle(0);
    pv2->SetBorderSize(0);
    pv2->AddText(Form("R=0.%d",Rpar));

    THnD *hZ[NDMC];
//    TH1D* hZG[NDMC];
//    TH1D* hZR[NDMC];

return;

}
