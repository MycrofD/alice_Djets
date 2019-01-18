#ifndef AliDJetSignalExtraction_H
#define AliDJetSignalExtraction_H

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
//  Class to extract D-jet pT or z spectrum
//
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
#include "TRandom2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TDatabasePDG.h"
#include "TNtuple.h"

#include "AliHFMassFitter.h"

class AliDJetSignalExtraction : public TObject
{

public:
    
    enum DMesonSpecies {kD0toKpi, kDStarD0pi};
    enum YieldMethod {kEffScale, kSideband};

    AliDJetSignalExtraction(); // default constructor
    AliDJetSignalExtraction(const AliDJetSignalExtraction &source);
    virtual ~AliDJetSignalExtraction();

    void SetInputFilename(TString filename) {fFileNameInput=filename;}
    void SetInputDirname(TString dirname) {fDirName=dirname;} 
    void SetInputListname(TString listname) {fListName=listname;} 
    void SetInputObjectname(TString objname) {fObjectName=objname;} 

    void SetRawOutputFileName(TString filename) {fFileRawOutput=filename;}

    Bool_t SetDmesonSpecie(DMesonSpecies k);
    void SetYieldMethod(YieldMethod meth) {fYieldApproach=meth;}
    void SetPtBinEdgesForMassPlot(Double_t ptmin, Double_t ptmax) {fpTmin=ptmin; fpTmax=ptmax;}
    void SetZedges(Double_t zmin, Double_t zmax) {fzmin=zmin; fzmax=zmax;}
    void SetDmesonPtBins(Int_t nbins=0, Double_t* ptedges=0x0);
    void SetJetPtBins(Int_t nbins=0, Double_t* ptedges=0x0);
    void SetDmesonEfficiency(Double_t* effvalues=0x0);
    void SetRebinSpectrumIfSBApproach(Bool_t rebin) {fRebinDstarSB=rebin;}

  
    Bool_t ExtractInputMassPlot();
    Bool_t ExtractInputMassPlotDzeroEffScale();
    Bool_t ExtractInputMassPlotDzeroSideband();
    Bool_t ExtractInputMassPlotDstarEffScale();
    Bool_t ExtractInputMassPlotDstarSideband();


    void SetDebugLevel(Int_t debug) {fDebug=debug;}
    void ClearObjects();

private:
    
    TFile *fFileInput;      		// file containing the task output
    TString fFileNameInput;		// name of input file
    TString fDirName; 			// name of input directory in the root file (Dstar)
    TString fListName; 			// name of input list (Dstar)
    TString fObjectName; 		// name of input container to extract the mass plot (Dstar)

    TFile *fFileRawOutput;		// output file containing efficiency-corrected D-jet spectrum   


    DMesonSpecies fDmesonSpecie;	// D meson specie    
    TString fDmesonLabel;		// D meson label
    YieldMethod fYieldApproach;		// method to extract jet pT spectrum
    
    Double_t fpTmin;	    		// pT lower edge of mass plot to evaluate variations of yields
    Double_t fpTmax;	   		// pT upper edge of mass plot to evaluate variations of yields
    Double_t fzmin;	    		// z minimum value to extract jet pT spectrum    
    Double_t fzmax;	    		// z maximum value to extract jet pT spectrum
    Double_t fmassmin;	    		// mass lower edge of inv.mass plots
    Double_t fmassmax;	    		// mass upper edge of inv.mass plots
    Double_t fmasswidth;		// mass plots bin width
    Int_t fnDbins;			// number of D-meson pT bins 
    Double_t *fDbinpTedges;    		// D-meson pt bin edges values
    Int_t fnJetbins;			// number of pT-bins to be used for spectrum
    Double_t *fJetbinpTedges;		// jet pT bin edges to be used for spectrum
    Double_t *fDEffValues;    		// D-meson efficiency values

    Double_t fnSigmaSignReg;		// Number of sigma for signal region
   

    TH1D* fMassPlot;		   	// mass spectra to be fitted
    
    //TH1F** fJetSpectrSBVars;		// array of jet spectrum histograms, one per variation (sideband approach)
    //TH1F** fJetPtBinYieldDistribution;  // array of histograms with yield distributions from the trials for each pT(jet)

    Int_t fDebug;			// debug level

    ClassDef(AliDJetSignalExtraction,1); 

};

#endif

