#ifndef AliDJetsPtSpectrum_H
#define AliDJetsPtSpectrum_H

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
//  Class to extract D-jet pT spectrum
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
#include "AliHFInvMassFitter.h"

class AliDJetsPtSpectrum : public TObject
{

public:
    
    enum DMesonSpecies 	{kD0toKpi, kDStarD0pi};
    enum YieldMethod 	{kEffScale, kSideband};

    AliDJetsPtSpectrum(); // default constructor
    AliDJetsPtSpectrum(const AliDJetsPtSpectrum &source);
    virtual ~AliDJetsPtSpectrum();

    void SetOutputFilename(TString filename) 							{fFileNameInput=filename;}
    void SetOutputDirectory(TString filename) 							{fOutputDir=filename;}

    Bool_t SetDmesonSpecie(DMesonSpecies k);
    void SetYieldMethod(YieldMethod meth) 								{fYieldApproach=meth;}
    void SetPtBinEdgesForMassPlot(Double_t ptmin, Double_t ptmax) {fpTmin=ptmin; fpTmax=ptmax;}
    void SetZedges(Double_t zmin, Double_t zmax) 						{fzmin=zmin; fzmax=zmax;}
    void SetDmesonPtBins(Int_t nbins=0, Double_t* ptedges=0x0);
    void SetJetPtBins(Int_t nbins=0, Double_t* ptedges=0x0);
    void SetJetEta(Double_t etamin, Double_t etamax)					{fjetetamin=etamin; fjetetamax=etamax;}
    void SetEffCorrection(Bool_t eff=1)									{fEffCorr=eff;}	 
	 
	 void SetSigmaSign(Double_t sigma);										{fnSigmaSign=sigma;}
	 Bool_t SetSigmaBkg(Int_t nbins, Double_t* sigma=0x0);
	 void SetBkgType(Int bkgType=5)											{fbkgtype=bkgType;}
	 void SetDMass(Double_t mass)												{fDmass=mass;}
	 void SetDSigma(Double_t sigma)											{fDsigma=sigma;}
  
    Bool_t ExtractInputMassPlot(const Int_t nFiles, TString *datafiles, Bool_t postfix=0, TString listName="Cut");
    Bool_t ExtractInputMassPlotEffScale();
    Bool_t ExtractInputMassPlotSideband();


    void SetDebugLevel(Int_t debug) {fDebug=debug;}

    void setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Size_t size=0.9, Width_t width=2, int scale=0);
    void SaveCanvas(TCanvas *c, TString name="tmp");

private:
    
    TFile*			fFileOutput;			// output file containing efficiency-corrected D-jet spectrum
    TString 		fFileNameOutput;		// name of output file
    TString			fOutputDir;				// name of output directory

    DMesonSpecies fDmesonSpecie;			// D meson specie    
    TString 		fDmesonLabel;			// D meson label
    YieldMethod 	fYieldApproach;		// method to extract jet pT spectrum
    Int_t 			fNDdata;					// number of D mesons candidates per event analyzed in data
    Int_t 			fNDMC;					// number of D mesons candidates per event analyzed in MC
    
    Double_t 		fpTmin;	    			// pT lower edge of mass plot to evaluate variations of yields
    Double_t 		fpTmax;	   			// pT upper edge of mass plot to evaluate variations of yields
    Double_t 		fzmin;	    			// z minimum value to extract jet pT spectrum    
    Double_t 		fzmax;	    			// z maximum value to extract jet pT spectrum
    Double_t 		fmassmin;	    		// mass lower edge of inv.mass plots
    Double_t 		fmassmax;	    		// mass upper edge of inv.mass plots
    Double_t		fjetetamin;
    Double_t		fjetetamax;
   
    Int_t 			fnDbins;					// number of D-meson pT bins 
    Double_t* 		fDbinpTedges;   		// D-meson pt bin edges values
    Int_t 			fnJetbins;				// number of pT-bins to be used for spectrum
    Double_t* 		fJetbinpTedges;		// jet pT bin edges to be used for spectrum
    Bool_t			fEffCorr;				// if to apply efficiency correction
    Double_t* 		fDEffValues;    		// D-meson efficiency values

    Double_t 		fnSigmaSign;			// Sigma for signal region
    Double_t 		fnSigmaBkg[4];			// Sigmas for the SB regions
    Int_t 			fbkgtype;				// Bkg type for inv. mass fitting
    Double_t 		fDmass;					// initial D mass for fitting
    Double_t 		fDsigma;					// initial D sigma for fitting
    
    TH3D* 			fInvMass3D;				// raw inv. mass plot

    TH1D* 			fRawPtSpectrum;		// raw D-jet pT spectrum
    TH1D* 			fEffCorrPtSpectrum;	// efficiency-corrected D-jet pT spectrum
    
  
    Int_t fDebug;			// debug level

    ClassDef(AliDJetsPtSpectrum,1); 

};

#endif

