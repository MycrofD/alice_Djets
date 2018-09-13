
#include <string>
#include <sstream>
#include <iostream>

#include "sys.h"
#include "config.h"

Int_t colors2[] = {1,2,kGreen+3,kMagenta+2,4,6,kCyan+1,8,kOrange-1,kGray+1,kViolet+5,kYellow+2};
Int_t markers2[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};
Int_t linestyle2[] = {1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};


const int nFiles = 2;
TString inDir[nFiles] = {
   "JetPtSpectrum_final"
,  "JetPtSpectrum_final"
}
/*
*/
TString desc[nFiles] = {
  "default",
  "FD Up",
  "FD Down"
};

double plotmin = 3, plotmax = 50;
const int ptbinsN = 9;
double ptbinsA[ptbinsN+1] = { 3,4,5,6,8,10,14,20,30,50 };

int nJetBins2 = 7;
double ptJetbins2[] = { 5,6,8,10,14,20,30,50 };

void sys_bfeedR(
TString inDirBase = "", 
TString inName = "", 
int truemax=50)
{

  TString inputDir = inDirBase;
//  inDir += "/";
//  inDir += input;
  gSystem->Exec(Form("mkdir %s/%s",inputDir.Data(),inName.Data()));

 yieldSys(inName,inDirBase,measmin,measmax,truemin,truemax);


return;
}

