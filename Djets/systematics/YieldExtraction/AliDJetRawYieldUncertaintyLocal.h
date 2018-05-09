#ifndef AliDJetRawYieldUncertaintyLocal_H
#define AliDJetRawYieldUncertaintyLocal_H

/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
//  Class to extract jet Pt spectrum yield uncertainty via multi-trial approach
//
//-----------------------------------------------------------------------
//  Author F.Colamaria
//  INFN Bari
//  fabio.colamaria@cern.ch
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
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TDatabasePDG.h"
#include "TNtuple.h"
#include "AliHFMultiTrials.h"
#include "AliAnalysisTaskDmesonJets.h"

class AliDJetRawYieldUncertaintyLocal : public TObject
{

public:

    enum DMesonSpecies {kD0toKpi, kDStarD0pi};
    enum YieldMethod {kEffScale, kSideband};

    AliDJetRawYieldUncertaintyLocal(); // default constructor
    AliDJetRawYieldUncertaintyLocal(const AliDJetRawYieldUncertaintyLocal &source);
    virtual ~AliDJetRawYieldUncertaintyLocal();

    void SetInputFilename(TString filename) {fFileNameInput=filename;}
    void SetInputDirname(TString dirname) {fDirName=dirname;}  //Dstar
    void SetInputListname(TString listname) {fListName=listname;}  //Dstar
    void SetInputObjectname(TString objname) {fObjectName=objname;}  //Dstar
    void SetInputTreename(TString treename) {fTreeName=treename;}  //Dzero
    void SetInputDBranchname(TString dname) {fDBranchName=dname;}  //Dzero
    void SetInputJetBranchname(TString jetname) {fJetBranchName=jetname;}  //Dzero

    Bool_t SetDmesonSpecie(DMesonSpecies k);
    void SetYieldMethod(YieldMethod meth) {fYieldApproach=meth;}
    void SetAllowRepetitionOfTrialExtraction(Bool_t allow) {fAllowRepetitions=allow;}
    void SetPtBinEdgesForMassPlot(Double_t ptmin, Double_t ptmax) {fpTmin=ptmin; fpTmax=ptmax;}
    void SetZedges(Double_t zmin, Double_t zmax) {fzmin=zmin; fzmax=zmax;}
    void SetMassEdgesAndBinWidthForMassPlot(Double_t mmin, Double_t mmax, Double_t mbinwidth) {fmassmin=mmin; fmassmax=mmax; fmasswidth=mbinwidth;} //Dzero
    void SetDmesonPtBins(Int_t nbins=0, Double_t* ptedges=0x0);
    void SetJetPtBins(Int_t nbins=0, Double_t* ptedges=0x0);
    void SetDmesonEfficiency(Double_t* effvalues=0x0);
    void SetRebinSpectrumIfSBApproach(Bool_t rebin) {fRebinDstarSB=rebin;}

    void SetSigmaForSignalRegion(Double_t nsigma) {fnSigmaSignReg=nsigma;}
    void SetMaxNTrialsForSidebandMethod(Int_t nmax) {fnMaxTrials=nmax;}
    void SetSigmaToFix(Double_t sigmafix) {fSigmaToFix=sigmafix;}
    void SetChi2Cut(Double_t chi2cut) {fChi2Cut=chi2cut;}
    void SetMeanSigmaVariations(Bool_t* cases);
    void SetBkgVariations(Bool_t* cases);
    void SetRebinSteps(Int_t nsteps, Int_t* cases);
    void SetMinMassSteps(Int_t nsteps, Double_t* cases);
    void SetMaxMassSteps(Int_t nsteps, Double_t* cases);
    void SetSigmaBinCounting(Int_t nsteps, Double_t* cases=0x0);
    void SetMaskOfVariations(Int_t ncases, Bool_t* cases);

    void SetFitReflections(Bool_t refl) {fFitRefl=refl;}
    void SetReflFilename(TString filename) {fReflFilenameInput=filename;}
    void SetMCSigFilename(TString filename) {fSigMCFilenameInput=filename;}
    void SetReflHistoname(TString histoname) {fReflHistoName=histoname;}
    void SetMCSigHistoname(TString histoname) {fSigMCHistoName=histoname;}
    void SetValueOfReflOverSignal(Double_t ratio, Double_t minrange=1.7, Double_t maxrange=2.1) {fFixRiflOverS=ratio; fReflRangeL=minrange; fReflRangeR=maxrange;}

    static void FitReflDistr(Int_t nPtBins, TString inputfile, TString fitType="DoubleGaus");

    Bool_t ExtractInputMassPlot();
    Bool_t ExtractInputMassPlotDzeroEffScale();
    Bool_t ExtractInputMassPlotDzeroSideband();
    Bool_t ExtractInputMassPlotDstarEffScale();
    Bool_t ExtractInputMassPlotDstarSideband();

    Bool_t RunMultiTrial();
    Bool_t CombineMultiTrialOutcomes();

    Bool_t EvaluateUncertainty();
    Bool_t EvaluateUncertainty_CoherentTrialChoice();
    Bool_t EvaluateUncertaintyDzeroEffScale();
    Bool_t EvaluateUncertaintyDzeroSideband();
    Bool_t EvaluateUncertaintyDzeroSideband_CoherentTrialChoice();
    Bool_t EvaluateUncertaintyDstarEffScale();
    Bool_t EvaluateUncertaintyDstarSideband();
    Bool_t EvaluateUncertaintyDstarSideband_CoherentTrialChoice();

    void SetDebugLevel(Int_t debug) {fDebug=debug;}
    void ClearObjects();

private:

    TFile *fFileInput;      		// file containing the task output
    TString fFileNameInput;		// name of input file
    TString fDirName; 			// name of input directory in the root file (Dstar)
    TString fListName; 			// name of input list (Dstar)
    TString fObjectName; 		// name of input container to extract the mass plot (Dstar)
    TString fTreeName; 			// name of input TTree (Dzero)
    TString fDBranchName;		// name of input branch for D meson (Dzero)
    TString fJetBranchName;		// name of input branch for jet (Dzero)

    DMesonSpecies fDmesonSpecie;	// D meson specie
    TString fDmesonLabel;		// D meson label
    YieldMethod fYieldApproach;		// method to extract jet pT spectrum
    Double_t fpTmin;	    		// pT lower edge of mass plot to evaluate variations of yields
    Double_t fpTmax;	   		// pT upper edge of mass plot to evaluate variations of yields
    Double_t fzmin;	    		// z minimum value to extract jet pT spectrum
    Double_t fzmax;	    		// z maximum value to extract jet pT spectrum
    Double_t fmassmin;	    		// mass lower edge of inv.mass plots (for Dzero, extracted from TTree)
    Double_t fmassmax;	    		// mass upper edge of inv.mass plots (for Dzero, extracted from TTree)
    Double_t fmasswidth;		// mass plots bin width (for Dzero, extracted from TTree)
    Int_t fnDbins;			// number of D-meson pT bins (for eff scaling)
    Double_t *fDbinpTedges;    		// D-meson pt bin edges values
    Int_t fnJetbins;			// number of pT-bins to be used for spectrum
    Double_t *fJetbinpTedges;		// jet pT bin edges to be used for spectrum
    Double_t *fDEffValues;    		// D-meson efficiency values

    Double_t fnSigmaSignReg;		// Number of sigma for signal region
    Double_t fSigmaToFix;		// Value of fixed sigma for MultiTrial
    Bool_t *fMeanSigmaVar;		// array of bools for mean/sigma variations
    Bool_t *fBkgVar;			// array of bools for bkg variations
    Int_t fnRebinSteps;			// number of steps for rebin
    Int_t *fRebinSteps;			// values of rebin steps
    Int_t fnMinMassSteps;		// number of steps for low mass edges
    Double_t *fMinMassSteps;		// values of low mass edges
    Int_t fnMaxMassSteps;		// number of steps for up mass edges
    Double_t *fMaxMassSteps;		// values of up mass edges
    Int_t fnSigmaBC;			// number of steps for BC (different sigma)
    Double_t *fSigmaBC;			// values of sigmas for BC
    Int_t fnMask;			// number of elements of array of variations to be kept (only for sigma/mean and bkg config)
    Bool_t *fMask;			// array of variations to be kept (only for sigma/mean and bkg config)
    Double_t fChi2Cut;			// maximum value of allowed chi2
    Int_t fnMaxTrials;			// max number of random trials for each pT(D) bin to build pT(jet) spectrum variations (sideband approach)
    Bool_t fAllowRepetitions;		// allow repetitions in the extraction of trials in a give pT(D) bin, for sideband approach
    Bool_t fRebinDstarSB;		// rebin the pt spectrum with user-defined binning, instead of using the binning from the THnSparse projection

    Bool_t fFitRefl;			// include reflection template in the mass fit
    TString fReflFilenameInput;		// name of input file for reflection template
    TString fSigMCFilenameInput;	// name of input file for MC signal
    TString fReflHistoName;		// name of reflection template histogram
    TString fSigMCHistoName;		// name of signal MC histogram
    Double_t fFixRiflOverS;		// rifl/signMC value for reflections
    Double_t fReflRangeL;		// lower range of reflection template (for rifl/signMC)
    Double_t fReflRangeR;		// upper range of reflection template (for rifl/signMC)

    TH1D* fMassPlot;		   	// mass spectra to be fitted
    TH1D* fJetYieldCentral;		// central values of the yield of jet spectrum + syst yield uncertainty
    TH1D* fJetYieldUnc;			// yield uncertainty vs jet pT bin
    TH1F** fJetSpectrSBVars;		// array of jet spectrum histograms, one per variation (sideband approach)
    TH1F** fJetPtBinYieldDistribution;  // array of histograms with yield distributions from the trials for each pT(jet)

    Int_t fDebug;			// debug level

    ClassDef(AliDJetRawYieldUncertaintyLocal,1); // class for plotting HF correlations

};

#endif
