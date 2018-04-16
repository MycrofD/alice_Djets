#ifndef SIGNALEXTRACTION_H
#define SIGNALEXTRACTION_H

/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//
//  D-jet pT spectrum extraction
//
//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------


#include "configDE.h"

    double jetmin = 0, jetmax = 50;
    double jetplotmin = 2, jetplotmax = 50;
    bool isEta = 0;
    double jetEta = 0.4;

    bool savePlots = 1;
    bool bEff = 0;
    bool isptcut = 1;
    bool isdetails = 0;

    TString fReflFilename;
    TString plotsDir;
    // mass fit params

    //------- efficiency
    double *efficiency = 0x0;

    TH1F   *hmass[fptbinsJetMeasN];
    TH1F *hmass_l[fptbinsJetMeasN];
    TH1F *hmass_u[fptbinsJetMeasN];
    TH1F *hmass_c[fptbinsJetMeasN];
    TF1  *fullfit[fptbinsJetMeasN];
    TF1  *massfit;//[fptbinsJetMeasN];
//    TF1   *bkgfit[fptbinsJetMeasN];
//    TF1  *bkgRfit[fptbinsJetMeasN];

    TH1F* hjetpt[fptbinsJetMeasN];
    TH1F *hjetpt_s[fptbinsJetMeasN];
    TH1F *hjetptsub[fptbinsJetMeasN];
    TH1F *hjetptcorr[fptbinsJetMeasN];

    TH1F *hjetptspectrum;
    TH1F *hrawjetptspectrum;
    TH1F *hjetptspectrumReb;
    TH1F *hjetptspectrumRebUnc;

    TH1F *hmean;
    TH1F *hsigma;
    TH1F *hrelErr;
    TH1F *hsign;
    TH1F *hsb;
    TH1F *hSignal;
    TH3D *hInvMassptD;
//    TH1F* hReflRS;
    TH1F *hmassjet[fptbinsJetMeasN][fptbinsDN];
    TH1F *hmassjet_scale[fptbinsJetMeasN][fptbinsDN];

    setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Width_t width);
    Bool_t rawJetSpectra(TString outdir, TString prod);
    void  saveSpectraPlots(TString outdir,TString prod);
    void  saveFitParams(TString outdir,TString prod);
    Bool_t SetReflection(AliHFInvMassFitter* &fitter, Int_t iBin, Float_t fLeftFitRange, Float_t fRightFitRange,Float_t &RS);

#endif
