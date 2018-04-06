#ifndef ALIANALYSISTASKFLAVOURJETCORRELATIONSTEST_H
#define ALIANALYSISTASKFLAVOURJETCORRELATIONSTEST_H
/**************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------
// Author : S. Antônio (University of São Paulo) antonio.silva@cern.ch
//          A. Grelli,  Utrecht University
//          C. Bianchin, Utrecht University
//          X. Zhang, LBNL
//	        B. Trzeciak, Utrecht Univeristy
//-----------------------------------------------------------------------


#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include "AliAODEvent.h"
#include "AliAnalysisTaskEmcalJet.h"

class TParticle;
class TClonesArray;
class THnSparse;
class AliMCParticle;
class AliAODMCParticle;
class AliRDHFCuts;
class AliEmcalJet;
class AliAODRecoDecayHF;
class AliAODRecoCascadeHF;
class AliAODEvent;
class AliParticleContainer;
class AliClusterContainer;
class AliJetContainer;

class AliAnalysisTaskFlavourJetCorrelationsTest : public AliAnalysisTaskEmcalJet
{

public:

   enum ECandidateType{ kD0toKpi, kDstartoKpipi };
   enum ECorrelationMethod{ kConstituent, kAngular, kResponseMatrix };

   AliAnalysisTaskFlavourJetCorrelationsTest();
   AliAnalysisTaskFlavourJetCorrelationsTest(const Char_t* name,AliRDHFCuts* cuts, ECandidateType candtype);
   virtual ~AliAnalysisTaskFlavourJetCorrelationsTest();

   virtual void     UserCreateOutputObjects();
   virtual Bool_t   Run();
   virtual void     Terminate(Option_t *);
   virtual void     Init();
   virtual void     LocalInit() {Init();}

   // inizializations
   Bool_t DefineHistoForAnalysis();

   // set MC usage
   void   SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
   Bool_t GetMC() const {return fUseMCInfo;}
   // set use Pythia info only for MC
   void   SetUsePythia(Bool_t theUsePythia) {fUsePythia = theUsePythia;}
   Bool_t GetUsePythia() const {return fUsePythia;}
   // set usage of reconstructed tracks
   void   SetUseReco(Bool_t reco) {fUseReco=reco;}
   Bool_t GetUseReco() const {return fUseReco;}

   void SetUseHFJet(Bool_t useHFJet=kTRUE) { fUseHFJet = useHFJet;}
   Bool_t GetUseHFJet() const {return fUseHFJet;}

   void SetHFLeadJet(Bool_t HFLeadJet=kFALSE) { fHFLeadJet = HFLeadJet;}
   Bool_t GetHFLeadJet() const  {return fHFLeadJet;}

   void SetMassLimits(Double_t range, Int_t pdg);
   void SetMassLimits(Double_t lowlimit, Double_t uplimit);

   void SetCorrelationMethod(Int_t c) {fCorrelationMethod=c;}
   Int_t GetCorrelationMethod() const {return fCorrelationMethod;}

   void SetUseCandArray(Bool_t b) {fUseCandArray=b;}
   Bool_t GetUseCandArray() const {return fUseCandArray;}

   void SetUseSBArray(Bool_t b) {fUseSBArray=b;}
   Bool_t GetUseSBArray() const {return fUseSBArray;}

   // Array of D0 width for the Dstar
   Bool_t SetD0WidthForDStar(Int_t nptbins,Float_t* width);
   void ConstituentCorrelationMethod(Bool_t IsBkg, AliAODEvent* aodEvent);
   void AngularCorrelationMethod(Bool_t IsBkg, AliAODEvent* aodEvent);
   void CreateResponseMatrix(AliEmcalJet* jet);
   void CreateMCResponseMatrix(AliEmcalJet* MCjet, AliAODEvent* aodEvent);
   void FillDJetHistograms(AliEmcalJet* jet, Double_t rho, Bool_t IsBkg, AliAODEvent* aodEvent);
   void GetHFJet(AliEmcalJet*& jet, Bool_t IsBkg);
   void FillHistogramsD0JetCorr(AliAODRecoDecayHF* candidate, Double_t z, Double_t ptD, Double_t ptj, Double_t jetTracks, Bool_t IsBkg, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc, AliAODEvent* aodEvent);
   void FillHistogramsDstarJetCorr(AliAODRecoCascadeHF* dstar, Double_t z, Double_t ptD, Double_t ptj, Double_t jetTracks, Bool_t IsBkg, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc);
   void FillHistogramsMCGenDJetCorr(Double_t z,Double_t ptD,Double_t ptjet, Double_t jetTracks, Double_t yD, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc);
   void FindMCJet(AliEmcalJet*& mcjet);
   Int_t IsDzeroSideBand(AliAODRecoCascadeHF *candDstar);
   Bool_t InEMCalAcceptance(AliVParticle *vpart);
   AliEmcalJet* GetJetCone(Double_t radius, Double_t eta, Double_t phi);

   void SetAnalyseDBackground(Bool_t b){ fAnalyseDBkg=b; }
   Bool_t GetAnalyseDBackground() const {return fAnalyseDBkg;}
   void SetBuildResponseMatrix(Bool_t b){ fBuildRM=b; }
   Bool_t GetBuildResponseMatrix() const {return fBuildRM;}
   void SetBuildResponseMatrixEff(Bool_t b){ fBuildRMEff=b; }
   Bool_t GetBuildResponseMatrixEff() const {return fBuildRMEff;}


private:

   AliAnalysisTaskFlavourJetCorrelationsTest(const AliAnalysisTaskFlavourJetCorrelationsTest &source);
   AliAnalysisTaskFlavourJetCorrelationsTest& operator=(const AliAnalysisTaskFlavourJetCorrelationsTest& source);

   Double_t Z(AliVParticle* part, AliEmcalJet* jet, Double_t rho) const;
   Double_t Z(AliVParticle* part, AliEmcalJet* jet) const;
   Double_t Z(Double_t* p, Double_t *pj) const;
   Double_t ZT(Double_t* p, Double_t *pj) const;
   Float_t DeltaR(AliEmcalJet *p1, AliVParticle *p2, Double_t rho) const;
   Float_t CheckDeltaR(AliEmcalJet *p1, AliVParticle *p2) const;


   Bool_t fUseMCInfo;               // Use MC info
   Bool_t fUseReco;                 // use reconstructed tracks when running on MC
   Bool_t fUsePythia;		    // Use Pythia info only for MC
   Bool_t fBuildRM;                 // flag to switch on/off the Response Matrix (Needs MC)
   Bool_t fBuildRMEff;              // flag to switch on/off the Response Matrix with efficiencies (Needs MC)
   Bool_t fUseHFJet;                 // if find HF jet
   Bool_t fHFLeadJet;                 // if HF jet should be the leading jet

   Int_t  fCandidateType;           // Dstar or D0
   Int_t  fCorrelationMethod;       // Method to correlate D mesons and jets
   Int_t  fPDGmother;               // PDG code of D meson
   Int_t  fNProngs;                 // number of prong of the decay channel
   Int_t  fPDGdaughters[4];         // PDG codes of daughters
   Float_t fSigmaD0[30];            // Sigma of D0 for D*
   TString fBranchName;             // AOD branch name
   AliRDHFCuts *fCuts;              // Cuts

   Double_t fMinMass;               // mass lower limit histogram
   Double_t fMaxMass;               // mass upper limit histogram

   Int_t fLeadingJetIds[2];	// Ids of two leading jets

   TClonesArray *fCandidateArray;   //! contains candidates selected by AliRDHFCuts
   TClonesArray *fSideBandArray;    //! contains candidates selected by AliRDHFCuts::IsSelected(kTracks), to be used for side bands (DStar case only!!)
   Bool_t fAnalyseDBkg;             // flag to switch off/on the SB analysis (default is off)


   Int_t fNAxesBigSparse;           // number of axis
   Bool_t fUseCandArray;            //! Use D meson candidates array
   Bool_t fUseSBArray;              //! Use D meson SB array

   // Histograms
   TH1I* fhstat;                    //!
   TH1F* fhCentDjet;                //!
   //generic jet and jet track distributions
   TH1I *fhNjets;			//!
   TH1I *fhNaccjets;			//!
   TH1F* fhPtJetTrks;		    //!
   TH1F* fhPhiJetTrks;              //!
   TH1F* fhEtaJetTrks;              //!
   TH1F* fhPtJet;                   //!
   TH1F* fhPtJetRaw;                //!
   TH1F* fhPhiJet;                  //!
   TH1F* fhEtaJet;                  //!
   TH2F *fhPtJetArea;			//!
   TH2F *fhJetAreaCent;			//!
   TH2D* fhRMRes;			//!

   //D mesons
   TH2F* fhInvMassptD;              //!
   TH2F* fhDiffSideBand;            //!
   TH2F* fhInvMassptDbg;            //!
   TH1F* fhPtPion;                  //!
                                    //!
   TH1I *fhDleadStat;			//!
   TH2F *fhDleadStatJetPt;		//!

   TH1F *fhLeadPt;				//!
   TH1F *fhLeadPt2;				//!
   TH1F *fhLeadEta;				//!
   TH1F *fhLeadEta2;				//!
   TH1F *fhLeadPhi;				//!
   TH1F *fhLeadPhi2;				//!
   TH1F *fhLeadArea;				//!
   TH1F *fhLeadArea2;				//!
   TH2F *fhLeadDeltaEtaDeltaPhi;				//!

   TH2F *fhRhoMult;			//!
   TH2F *fhRhoLeadPt;			//!
   TH2F *fhRhoLeadArea; 		//!
   TH2F *fhRhoLeadPtRaw;		//!
   TH3F *fhRhoLeadPtMult;		//!
   TH3F *fhRhoLeadPtMultRaw;		//!


   TH1F *fDeltaPT;			//!
   TH1F *fDeltaPT_excl_lead;		//!
   TH1F *fDeltaPT_trans;		//!
   TH2F *fDeltaPTCent;			//!
   TH2F *fDeltaPTCentPtCone;		//!
   TH2F *fDeltaPTCent_excl_lead;	//!
   TH2F *fDeltaPTCent_trans;	//!
   TH2F *fDeltaPTCentPtCone_excl_lead;	//!
   TH2F *fDeltaPTCentPtCone_trans;	//!
   TH2F *fDeltaPTRho;			//!
   TH2F *fDeltaPTRhoPtCone;		//!
   TH2F *fDeltaPTRho_excl_lead;		//!
   TH2F *fDeltaPTRho_trans;		//!
   TH2F *fDeltaPTRhoPtCone_excl_lead;	//!
   TH2F *fDeltaPTRhoPtCone_trans;	//!
   TH2F *fDeltaPTLeadPt;		//!
   TH2F *fDeltaPTLeadPtRaw;		//!
   TH2F *fDeltaPTLeadPtPtCone;		//!
   TH2F *fDeltaPTLeadPtRawPtCone;	//!
   TH2F *fDeltaPTLeadPt_excl_lead;	//!
   TH2F *fDeltaPTLeadPt_trans;	//!
   TH2F *fDeltaPTLeadPtRaw_excl_lead;	//!
   TH2F *fDeltaPTLeadPtRaw_trans;	//!
   TH2F *fDeltaPTLeadPtPtCone_excl_lead;	//!
   TH2F *fDeltaPTLeadPtPtCone_trans;	//!
   TH2F *fDeltaPTLeadPtRawPtCone_excl_lead;	//!
   TH2F *fDeltaPTLeadPtRawPtCone_trans;	//!
   TH2F *fRandConeEtaPhi;		//!
   TH2F *fRandConeEtaPhi_excl_lead;	//!
   TH2F *fRandConeEtaPhi_trans;	//!
   TH2F *fhRandomConeDeltaEtaDeltaPhi;	//!
   TH2F *fhRandomConeDeltaEtaDeltaPhi_trans;  //!
   TH2F *fEtaRandConeOverlap;		//!
   TH2F *fPhiRandConeOverlap;		//!

   //main histograms
   THnSparse* fhsDphiz;             //!
   THnSparse* fResponseMatrix;      //!


   ClassDef(AliAnalysisTaskFlavourJetCorrelationsTest,8); // class for charm-jet CorrelationsExch
};

#endif
