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
//
//
//  Analysis Taks for reconstructed particle correlation
//  (first implementation done for D mesons) with jets
//  (use the so called Emcal framework)
//
//-----------------------------------------------------------------------
// Authors:
// C. Bianchin (Utrecht University) chiara.bianchin@cern.ch
// S. Antônio (University of São Paulo) antonio.silva@cern.ch
// A. Grelli (Utrecht University) a.grelli@uu.nl
// X. Zhang (LBNL)  XMZhang@lbl.gov
// B. Trzeciak (Utrecht University) barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include "TROOT.h"
#include <THnSparse.h>
#include <TSystem.h>
#include <TObjectTable.h>
#include "AliMultSelection.h"

#include "AliAnalysisTaskFlavourJetCorrelationsRC.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliParticleContainer.h"
#include "AliEmcalParticle.h"
#include "AliLocalRhoParameter.h"
#include "AliAnalysisTaskLocalRho.h"

ClassImp(AliAnalysisTaskFlavourJetCorrelationsRC)


//_______________________________________________________________________________

AliAnalysisTaskFlavourJetCorrelationsRC::AliAnalysisTaskFlavourJetCorrelationsRC() :
AliAnalysisTaskEmcalJet("AliAnalysisTaskFlavourJetCorrelationsRC",kTRUE),
fUseMCInfo(kTRUE),
fUseReco(kTRUE),
fUsePythia(kFALSE),
fCandidateType(),
fCorrelationMethod(),
fPDGmother(),
fNProngs(),
fPDGdaughters(),
fBranchName(),
fCuts(0),
fMinMass(),
fMaxMass(),
fLeadingJetIds(),
fCandidateArray(0),
fSideBandArray(0),
fAnalyseDBkg(kFALSE),
fBuildRM(kFALSE),
fBuildRMEff(kFALSE),
fUseHFJet(kTRUE),
fNAxesBigSparse(9),
fUseCandArray(kFALSE),
fUseSBArray(kFALSE),
fhstat(),
fhCentDjet(),
fhPtJetTrks(),
fhPhiJetTrks(),
fhEtaJetTrks(),
fhPtJet(),
fhPtJetRaw(),
fhPhiJet(),
fhEtaJet(),
fhNjets(),
fhNaccjets(),
fhPtJetArea(),
fhJetAreaCent(),
fhRMRes(),
fhDleadStat(),
fhDleadStatJetPt(),
fhLeadPt(),
fhLeadPt2(),
fhLeadEta(),
fhLeadEta2(),
fhLeadPhi(),
fhLeadPhi2(),
fhLeadArea(),
fhLeadArea2(),
fhLeadDeltaEtaDeltaPhi(),
fhInvMassptD(),
fhDiffSideBand(),
fhInvMassptDbg(),
fhPtPion(),
fhsDphiz(),
fResponseMatrix(),
fhRhoMult(),
fhRhoLeadPt(),
fhRhoLeadArea(),
fhRhoLeadPtRaw(),
fhRhoLeadPtMult(),
fhRhoLeadPtMultRaw(),
fDeltaPT(),
fDeltaPT_excl_lead(),
fDeltaPT_trans(),
fDeltaPTCent(),
fDeltaPTCentPtCone(),
fDeltaPTCent_excl_lead(),
fDeltaPTCent_trans(),
fDeltaPTCentPtCone_excl_lead(),
fDeltaPTCentPtCone_trans(),
fDeltaPTRho(),
fDeltaPTRhoPtCone(),
fDeltaPTRho_excl_lead(),
fDeltaPTRho_trans(),
fDeltaPTRhoPtCone_excl_lead(),
fDeltaPTRhoPtCone_trans(),
fDeltaPTLeadPt(),
fDeltaPTLeadPtRaw(),
fDeltaPTLeadPtPtCone(),
fDeltaPTLeadPtRawPtCone(),
fDeltaPTLeadPt_excl_lead(),
fDeltaPTLeadPt_trans(),
fDeltaPTLeadPtRaw_excl_lead(),
fDeltaPTLeadPtRaw_trans(),
fDeltaPTLeadPtPtCone_excl_lead(),
fDeltaPTLeadPtPtCone_trans(),
fDeltaPTLeadPtRawPtCone_excl_lead(),
fDeltaPTLeadPtRawPtCone_trans(),
fRandConeEtaPhi(),
fRandConeEtaPhi_excl_lead(),
fRandConeEtaPhi_trans(),
fhRandomConeDeltaEtaDeltaPhi(),
fhRandomConeDeltaEtaDeltaPhi_trans(),
fEtaRandConeOverlap(),
fPhiRandConeOverlap()

{
   //
   // Default ctor
}

//_______________________________________________________________________________

AliAnalysisTaskFlavourJetCorrelationsRC::AliAnalysisTaskFlavourJetCorrelationsRC(const Char_t* name, AliRDHFCuts* cuts,ECandidateType candtype) :
AliAnalysisTaskEmcalJet(name,kTRUE),
fUseMCInfo(kTRUE),
fUseReco(kTRUE),
fUsePythia(kFALSE),
fCandidateType(),
fCorrelationMethod(),
fPDGmother(),
fNProngs(),
fPDGdaughters(),
fBranchName(),
fCuts(0),
fMinMass(),
fMaxMass(),
fLeadingJetIds(),
fCandidateArray(0),
fSideBandArray(0),
fAnalyseDBkg(kFALSE),
fBuildRM(kFALSE),
fBuildRMEff(kFALSE),
fUseHFJet(kTRUE),
fNAxesBigSparse(9),
fUseCandArray(kFALSE),
fUseSBArray(kFALSE),
fhstat(),
fhCentDjet(),
fhPtJetTrks(),
fhPhiJetTrks(),
fhEtaJetTrks(),
fhPtJet(),
fhPtJetRaw(),
fhPhiJet(),
fhEtaJet(),
fhNjets(),
fhNaccjets(),
fhPtJetArea(),
fhJetAreaCent(),
fhRMRes(),
fhDleadStat(),
fhDleadStatJetPt(),
fhLeadPt(),
fhLeadPt2(),
fhLeadEta(),
fhLeadEta2(),
fhLeadPhi(),
fhLeadPhi2(),
fhLeadArea(),
fhLeadArea2(),
fhLeadDeltaEtaDeltaPhi(),
fhInvMassptD(),
fhDiffSideBand(),
fhInvMassptDbg(),
fhPtPion(),
fhsDphiz(),
fResponseMatrix(),
fhRhoMult(),
fhRhoLeadPt(),
fhRhoLeadArea(),
fhRhoLeadPtRaw(),
fhRhoLeadPtMult(),
fhRhoLeadPtMultRaw(),
fDeltaPT(),
fDeltaPT_excl_lead(),
fDeltaPT_trans(),
fDeltaPTCent(),
fDeltaPTCentPtCone(),
fDeltaPTCent_excl_lead(),
fDeltaPTCent_trans(),
fDeltaPTCentPtCone_excl_lead(),
fDeltaPTCentPtCone_trans(),
fDeltaPTRho(),
fDeltaPTRhoPtCone(),
fDeltaPTRho_excl_lead(),
fDeltaPTRho_trans(),
fDeltaPTRhoPtCone_excl_lead(),
fDeltaPTRhoPtCone_trans(),
fDeltaPTLeadPt(),
fDeltaPTLeadPtRaw(),
fDeltaPTLeadPtPtCone(),
fDeltaPTLeadPtRawPtCone(),
fDeltaPTLeadPt_excl_lead(),
fDeltaPTLeadPt_trans(),
fDeltaPTLeadPtRaw_excl_lead(),
fDeltaPTLeadPtRaw_trans(),
fDeltaPTLeadPtPtCone_excl_lead(),
fDeltaPTLeadPtPtCone_trans(),
fDeltaPTLeadPtRawPtCone_excl_lead(),
fDeltaPTLeadPtRawPtCone_trans(),
fRandConeEtaPhi(),
fRandConeEtaPhi_excl_lead(),
fRandConeEtaPhi_trans(),
fhRandomConeDeltaEtaDeltaPhi(),
fhRandomConeDeltaEtaDeltaPhi_trans(),
fEtaRandConeOverlap(),
fPhiRandConeOverlap()
{
   //
   // Constructor. Initialization of Inputs and Outputs
   //

   Info("AliAnalysisTaskFlavourJetCorrelationsRC","Calling Constructor");
   fCuts=cuts;
   fCandidateType=candtype;
   const Int_t nptbins=fCuts->GetNPtBins();
   Float_t defaultSigmaD013[20]={0.012, 0.012, 0.012, 0.015, 0.015,0.018,0.018,0.020,0.020,0.030,0.030,0.037,0.040,0.040,0.040,0.040,0.040,0.040,0.040,0.040};
   fLeadingJetIds[0]=-1;
   fLeadingJetIds[1]=-1;

   switch(fCandidateType){
   case 0:
      fPDGmother=421;
      fNProngs=2;
      fPDGdaughters[0]=211;//pi
      fPDGdaughters[1]=321;//K
      fPDGdaughters[2]=0; //empty
      fPDGdaughters[3]=0; //empty
      fBranchName="D0toKpi";
      break;
   case 1:
      fPDGmother=413;
      fNProngs=3;
      fPDGdaughters[1]=211;//pi soft
      fPDGdaughters[0]=421;//D0
      fPDGdaughters[2]=211;//pi fromD0
      fPDGdaughters[3]=321; // K from D0
      fBranchName="Dstar";

      if(nptbins<20){
      	 for(Int_t ipt=0;ipt<nptbins;ipt++) fSigmaD0[ipt]= defaultSigmaD013[ipt];
      } else {
      	 AliFatal(Form("Default sigma D0 not enough for %d pt bins, use SetSigmaD0ForDStar to set them",nptbins));
      }
      break;
   default:
      printf("%d not accepted!!\n",fCandidateType);
      break;
   }

   if(fCandidateType==kD0toKpi)SetMassLimits(0.15,fPDGmother);
   if(fCandidateType==kDstartoKpipi) SetMassLimits(0.015, fPDGmother);
   if(fUseCandArray) DefineInput(1, TClonesArray::Class());
   if(fUseSBArray) DefineInput(2, TClonesArray::Class());

      DefineOutput(1,TList::Class()); // histos
      DefineOutput(2,AliRDHFCuts::Class()); // my cuts

}

//_______________________________________________________________________________

AliAnalysisTaskFlavourJetCorrelationsRC::~AliAnalysisTaskFlavourJetCorrelationsRC() {
   //
   // destructor
   //

   Info("~AliAnalysisTaskFlavourJetCorrelationsRC","Calling Destructor");

   delete fCuts;

}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelationsRC::Init(){
   //
   // Initialization
   //

   if(fDebug > 1) printf("AnalysisTaskRecoJetCorrelations::Init() \n");

   switch(fCandidateType){
   case 0:
      {
      	 AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fCuts)));
      	 copyfCuts->SetName("AnalysisCutsDzero");
      	 // Post the data
      	 PostData(2,copyfCuts);
      }
      break;
   case 1:
      {
      	 AliRDHFCutsDStartoKpipi* copyfCuts=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
      	 copyfCuts->SetName("AnalysisCutsDStar");
      	 // Post the cuts
      	 PostData(2,copyfCuts);
      }
      break;
   default:
      return;
   }

   return;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelationsRC::UserCreateOutputObjects() {
   // output
   Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
   AliAnalysisTaskEmcal::UserCreateOutputObjects();

   // define histograms
   // the TList fOutput is already defined in  AliAnalysisTaskEmcal::UserCreateOutputObjects()
   DefineHistoForAnalysis();
   PostData(1,fOutput);

   return;
}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskFlavourJetCorrelationsRC::Run()
{
   // user exec from AliAnalysisTaskEmcal is used

   // Load the event
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  if(!aodEvent) { return kFALSE; }

  Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
  if (matchingAODdeltaAODlevel<=0) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      return kFALSE;
  }

   TClonesArray* mcArray = 0x0;
   if (fUseMCInfo) { //not used at the moment,uncomment return if you use
      mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcArray) {
      	 printf("AliAnalysisTaskFlavourJetCorrelationsRC::UserExec: MC particles not found!\n");
      }
   }

    //D meson candidates. Also background if is MC

    if(fHFLeadJet && !fUseHFJet) { printf("AliAnalysisTaskFlavourJetCorrelationsRC::UserExec: HF set as the leading jet, but HFjet FALSE!\n"); return kFALSE; }

    if(fUseHFJet) {
      if(fUseCandArray)
      {
          fCandidateArray = dynamic_cast<TClonesArray*>(GetInputData(1));
          if (!fCandidateArray) return kFALSE;
          for(Int_t icand=0; icand<fCandidateArray->GetEntriesFast(); icand++)
          {
              fhstat->Fill(2);
          }
      }
      if(fUseSBArray)
      {
          fSideBandArray = dynamic_cast<TClonesArray*>(GetInputData(2));
          if (!fSideBandArray) return kFALSE;
      }
    }

    fhstat->Fill(0);

   // fix for temporary bug in ESDfilter
   // the AODs with null vertex pointer didn't pass the PhysSel
   if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return kFALSE;

   //Event selection
   Bool_t iseventselected=fCuts->IsEventSelected(aodEvent);
   TString firedTriggerClasses=((AliAODEvent*)aodEvent)->GetFiredTriggerClasses();
   if(!iseventselected) return kFALSE;

   fhstat->Fill(1);
    Float_t lPercentile = 300;
  /*  if(aodEvent->GetRunNumber()>200000)
    {
      AliMultSelection *MultSelection = 0x0;
      MultSelection = (AliMultSelection*)aodEvent->FindListObject("MultSelection");
      if(!MultSelection) {
          //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
          AliWarning("AliMultSelection object not found!");
      }else{
          lPercentile = MultSelection->GetMultiplicityPercentile("V0A");
      }
      fhCentDjet->Fill(lPercentile);
    }
    else fhCentDjet->Fill(fCent);*/

// for MC response matrix of efficiency studies, fMultCand option only
if(fUseMCInfo && fBuildRMEff){

    AliJetContainer* mcjets = 0;
    if(!fAnalyseDBkg) mcjets = GetJetContainer(1);
    else mcjets = GetJetContainer(2);
    if(!mcjets) return kFALSE;

    AliParticleContainer *MCParticlesCont = mcjets->GetParticleContainer();

    mcjets->ResetCurrentID();
    AliEmcalJet* jet=0;

    while ((jet = mcjets->GetNextJet()))
    {
        UInt_t rejectionReason = 0;
        Bool_t OKjet = mcjets->AcceptJet(jet, rejectionReason);
        if(!OKjet) {
            fhstat->Fill(5);
            continue;
        }

        fhstat->Fill(3); //Jet accepted
        fhPhiJet->Fill(jet->Phi());
        fhEtaJet->Fill(jet->Eta());
        fhPtJet->Fill(jet->Pt());

        Int_t ntrjet=  jet->GetNumberOfTracks();

        for(Int_t itrk=0;itrk<ntrjet;itrk++)
        {
            AliAODMCParticle* jetTrk=(AliAODMCParticle*)jet->TrackAt(itrk,MCParticlesCont->GetArray());
            if (!jetTrk) continue;
            fhPtJetTrks->Fill(jetTrk->Pt());
            fhPhiJetTrks->Fill(jetTrk->Phi());
            fhEtaJetTrks->Fill(jetTrk->Eta());
        } //end loop on jet tracks

    } // end loop on mc jets

    // Get HF accepted MC jet
    AliEmcalJet* MCjet = 0;
    FindMCJet(MCjet);
    if(!MCjet) return kFALSE;
    //if( TMath::Abs(MCjet->Eta()) > (0.9 - mcjets->GetJetRadius()) ) return kFALSE;

    if(fCorrelationMethod==kConstituent)
    {
            if(fBuildRMEff==kTRUE) CreateMCResponseMatrix(MCjet, aodEvent);
    }
    /* the other method not enabled for now
     * else if(fCorrelationMethod==kAngular)
    {
        if(fCandidateArray->GetEntriesFast()>0) AngularCorrelationMethod(kFALSE,aodEvent);
        if(fAnalyseDBkg==kTRUE && fSideBandArray->GetEntriesFast()>0) AngularCorrelationMethod(kTRUE,aodEvent);
    }*/

}

else {

    AliJetContainer* JetCont = GetJetContainer(0);
    if(!JetCont) return kFALSE;
    AliParticleContainer *ParticlesCont = JetCont->GetParticleContainer();

    AliJetContainer* JetContSB = 0;
    AliParticleContainer *ParticlesContSB = 0;
    if(fAnalyseDBkg)
    {
        JetContSB = GetJetContainer(1);
        if(!JetContSB) return kFALSE;
        ParticlesContSB = JetContSB->GetParticleContainer();
    }

    AliEmcalJet* HFjet = NULL;
    if(fUseHFJet){
      GetHFJet(HFjet,kFALSE);
      if(!HFjet) return kFALSE;
    }

    JetCont->ResetCurrentID();
    AliEmcalJet* jet=0;

    Int_t maxJetIds[]   = {-1, -1};
    Float_t maxJetPts[] = {0,  0};
    Float_t maxJetPtsCorr[] = {0,  0};
    Float_t maxJetAreas[] = {0, 0};
    Float_t maxJetEta[] = {0, 0};
    Float_t maxJetPhi[] = {0, 0};
    //const Int_t Njets = JetCont->GetEntries();
    //if(!Njets) return kFALSE;
    Int_t counter = 0;
    Int_t accjet = 0;
    AliEmcalJet *leadingJet = 0;
    while ((jet = JetCont->GetNextJet()))
    //for (Int_t ij = 0; ij < Njets; ++ij)
    {
        //jet = static_cast<AliEmcalJet*>(JetCont->At(ij));
        counter++;
        if (!jet) {
            //AliError(Form("%s: Could not receive jet %d", GetName(), ij));
            continue;
        }

        UInt_t rejectionReason = 0;
        Bool_t OKjet = JetCont->AcceptJet(jet, rejectionReason);
        if(!OKjet) {
            fhstat->Fill(5);
            continue;
        }
        accjet++;

         if (jet->Pt() > maxJetPts[0]) {
            maxJetPts[1] = maxJetPts[0];
            maxJetPtsCorr[1] = maxJetPtsCorr[0];
            maxJetIds[1] = maxJetIds[0];
            maxJetAreas[1] = maxJetAreas[0];
            maxJetEta[1] = maxJetEta[0];
            maxJetPhi[1] = maxJetPhi[0];
            maxJetPts[0] = jet->Pt();
            maxJetPtsCorr[0] = jet->Pt() - jet->Area()*JetCont->GetRhoVal();
            maxJetAreas[0] = jet->Area();
            maxJetEta[0] = jet->Eta();
            maxJetPhi[0] = jet->Phi();
            maxJetIds[0] = counter-1;
            leadingJet = static_cast<AliEmcalJet*>(jet);
        } else if (jet->Pt() > maxJetPts[1]) {
            maxJetPts[1] = jet->Pt();
            maxJetPtsCorr[1] = jet->Pt() - jet->Area()*JetCont->GetRhoVal();
            maxJetAreas[1] = jet->Area();
            maxJetEta[1] = jet->Eta();
            maxJetPhi[1] = jet->Phi();
            maxJetIds[1] = counter-1;
        }

        Double_t JetPtCorr = 0;
        if(fUseMCInfo && fUsePythia){
            JetPtCorr = jet->Pt();
        }
        else {
            JetPtCorr = jet->Pt() - jet->Area()*JetCont->GetRhoVal(); //background subtraction
            if(fLocalRho)
            {
                JetPtCorr = jet->Pt() - jet->Area()*fLocalRho->GetLocalVal(jet->Phi(),JetCont->GetJetRadius(),JetCont->GetRhoVal()); //
            }
        }
        fhstat->Fill(3); //Jet accepted
        fhPhiJet->Fill(jet->Phi());
        fhEtaJet->Fill(jet->Eta());
        fhPtJet->Fill(JetPtCorr);
        fhPtJetRaw->Fill(jet->Pt());
        fhPtJetArea->Fill(JetPtCorr,jet->Area());
        fhJetAreaCent->Fill(lPercentile,jet->Area());
    }

    fhNjets->Fill(counter);
    fhNaccjets->Fill(accjet);

    if(!accjet) return kTRUE;

    //Distribution of all particles in the event
    Int_t ntrarr=ParticlesCont->GetNParticles();
    for(Int_t i=0;i<ntrarr;i++)
    {
        AliVParticle* jetTrk= ParticlesCont->GetParticle(i);
        if (!jetTrk) continue;
        AliEmcalParticle* emcpart = dynamic_cast<AliEmcalParticle*>(jetTrk);
        if (emcpart) jetTrk = emcpart->GetTrack();
        fhPtJetTrks->Fill(jetTrk->Pt());
        fhPhiJetTrks->Fill(jetTrk->Phi());
        fhEtaJetTrks->Fill(jetTrk->Eta());
    }


    fLeadingJetIds[0]=maxJetIds[0];
    fLeadingJetIds[1]=maxJetIds[1];
    Double_t rho = JetCont->GetRhoVal();

/*
      fhRhoMult->Fill(lPercentile,rho);

      if(fUseHFJet && fHFLeadJet){
        fhRhoLeadPt->Fill(HFjet->Pt() - HFjet->Area()*rho,rho);
        fhRhoLeadPtMult->Fill(HFjet->Pt() - HFjet->Area()*rho,rho);
        fhRhoLeadPtRaw->Fill(HFjet->Pt(),rho);
        fhRhoLeadPtMultRaw->Fill(HFjet->Pt(),lPercentile,rho);
        fhRhoLeadArea->Fill(HFjet->Area(),rho);
        fhLeadPt->Fill(HFjet->Pt());
        fhLeadEta->Fill(HFjet->Eta());
        fhLeadPhi->Fill(HFjet->Phi());
        fhLeadArea->Fill(HFjet->Area());

        if( (maxJetIds[0] != -1) && ((Float_t)HFjet->Eta() != maxJetEta[0] && (Float_t)HFjet->Phi() != maxJetPhi[0]) ){
            fhLeadPt2->Fill(maxJetPtsCorr[0]);
            fhLeadEta2->Fill(maxJetEta[0]);
            fhLeadPhi2->Fill(maxJetPhi[0]);
            fhLeadArea2->Fill(maxJetAreas[0]);
        }
        else if( (maxJetIds[1] != -1) && ((Float_t)HFjet->Eta() != maxJetEta[1] && (Float_t)HFjet->Phi() != maxJetPhi[1]) ){
            fhLeadPt2->Fill(maxJetPtsCorr[1]);
            fhLeadEta2->Fill(maxJetEta[1]);
            fhLeadPhi2->Fill(maxJetPhi[1]);
            fhLeadArea2->Fill(maxJetAreas[1]);
        }

      }
      else {
        if(maxJetIds[0] != -1)  {
            fhRhoLeadPt->Fill(maxJetPtsCorr[0],rho);
            fhRhoLeadPtMult->Fill(maxJetPtsCorr[0],lPercentile,rho);
            fhRhoLeadPtRaw->Fill(maxJetPts[0],rho);
            fhRhoLeadPtMultRaw->Fill(maxJetPts[0],lPercentile,rho);
            fhRhoLeadArea->Fill(maxJetAreas[0],rho);
            fhLeadPt->Fill(maxJetPts[0]);
            fhLeadEta->Fill(maxJetEta[0]);
            fhLeadPhi->Fill(maxJetPhi[0]);
            fhLeadArea->Fill(maxJetAreas[0]);
        }
        if(maxJetIds[1] != -1){
            fhLeadPt2->Fill(maxJetPtsCorr[1]);
            fhLeadEta2->Fill(maxJetEta[1]);
            fhLeadPhi2->Fill(maxJetPhi[1]);
            fhLeadArea2->Fill(maxJetAreas[1]);
        }

        if(maxJetIds[0] != -1 && maxJetIds[1] != -1) {
            Double_t deltaEta = maxJetEta[0]-maxJetEta[1];
            Double_t deltaPhi = TMath::Abs(maxJetPhi[0]-maxJetPhi[1]);
            if(deltaPhi>TMath::Pi()) deltaPhi = (2*TMath::Pi())-deltaPhi;
            fhLeadDeltaEtaDeltaPhi->Fill(deltaPhi,deltaEta);
        }

        if(fUseHFJet){
            if( (maxJetIds[0] != -1) && ((Float_t)HFjet->Eta() == maxJetEta[0] && (Float_t)HFjet->Phi() == maxJetPhi[0]) ){
              fhDleadStatJetPt->Fill(1.,HFjet->Pt()-HFjet->Area()*rho);
              fhDleadStat->Fill(1);
            }
            else if( (maxJetIds[1] != -1) && ((Float_t)HFjet->Eta() == maxJetEta[1] && (Float_t)HFjet->Phi() == maxJetPhi[1]) ){
              fhDleadStatJetPt->Fill(2.,HFjet->Pt()-HFjet->Area()*rho);
              fhDleadStat->Fill(2);
            }
            else {
              fhDleadStatJetPt->Fill(0.,HFjet->Pt()-HFjet->Area()*rho);
              fhDleadStat->Fill(0);
            }

        }

      }
*/

//*************************************
// RANDOM CONES ====================
//*************************************

TRandom3 rphi(0);
TRandom3 reta(0);

Double_t conephi, coneeta, ptcone;
Float_t Bkgpt = 0, deltaPt = 0;
AliVTrack* track = 0;

//*************************************
// RANDOM CONES, no lead. jet exclution
//*************************************
/*    conephi = rphi.Uniform(0,TMath::TwoPi());
    coneeta = rphi.Uniform(-(0.9-JetCont->GetJetRadius()),0.9-JetCont->GetJetRadius());

    ptcone = 0.;
    JetCont->ResetCurrentID();
    ParticlesCont->ResetCurrentID();
    track = 0;
    while ((track = static_cast<AliVTrack*>(ParticlesCont->GetNextAcceptParticle()))) {
        Double_t DeltaPhi = TMath::Abs(track->Phi()-conephi);
        if(DeltaPhi>TMath::Pi()) DeltaPhi = TMath::Abs(TMath::TwoPi() - DeltaPhi);
        Double_t DeltaEta = TMath::Abs(track->Eta()-coneeta);
        Double_t DeltaRCone = TMath::Sqrt((DeltaPhi*DeltaPhi)+(DeltaEta*DeltaEta));

        if(DeltaRCone<JetCont->GetJetRadius())
        {
            ptcone = ptcone + track->Pt();
        }
    }
    Bkgpt = JetCont->GetRhoVal()*TMath::Pi()*JetCont->GetJetRadius()*JetCont->GetJetRadius();
    deltaPt = ptcone - Bkgpt;

    fDeltaPT->Fill(deltaPt);
    fDeltaPTCent->Fill(lPercentile,deltaPt);
    fDeltaPTCentPtCone->Fill(lPercentile,ptcone);
    fDeltaPTRho->Fill(rho,deltaPt);
    fDeltaPTRhoPtCone->Fill(rho,ptcone);
    fRandConeEtaPhi->Fill(coneeta, conephi);

    if (leadingJet) {
        Float_t leadJetPt = leadingJet->Pt();
        Float_t leadJetPtCorr = leadingJet->Pt() - leadingJet->Area()*rho;
        fDeltaPTLeadPt->Fill(leadJetPtCorr,deltaPt);
        fDeltaPTLeadPtPtCone->Fill(leadJetPtCorr,ptcone);
        fDeltaPTLeadPtRaw->Fill(leadJetPt,deltaPt);
        fDeltaPTLeadPtRawPtCone->Fill(leadJetPt,ptcone);
    }
*/

//*************************************
// RANDOM CONES, EXCL LEADING =========
//*************************************

    if(!fHFLeadJet && !leadingJet) return kTRUE;
    if(fUseHFJet && fHFLeadJet) leadingJet = HFjet;

    JetCont->ResetCurrentID();
    ParticlesCont->ResetCurrentID();
    const int MAXCONE = 20;
    int isOverlap =0, isRanCone = 0;

    for(int icone=0; icone<MAXCONE; icone++){
        conephi = rphi.Uniform(0,TMath::TwoPi());
        coneeta = reta.Uniform(-(0.9-JetCont->GetJetRadius()),0.9-JetCont->GetJetRadius());
        ptcone = 0;
        int itracks = 0;
        std::vector< AliVParticle* > coneParticles;
        std::vector< Int_t > coneParticlesId;

        track = 0;
        for(Int_t counter=0;counter<ntrarr;counter++)
        {
            AliVParticle* track= ParticlesCont->GetParticle(counter);
            if (!track) continue;
            AliEmcalParticle* emcpart = dynamic_cast<AliEmcalParticle*>(track);
            if (emcpart) track = emcpart->GetTrack();

            Double_t DeltaPhi = TMath::Abs(track->Phi()-conephi);
            if(DeltaPhi>TMath::Pi()) DeltaPhi = TMath::Abs(TMath::TwoPi() - DeltaPhi);
            Double_t DeltaEta = TMath::Abs(track->Eta()-coneeta);
            Double_t DeltaRCone = TMath::Sqrt((DeltaPhi*DeltaPhi)+(DeltaEta*DeltaEta));

            if(DeltaRCone<JetCont->GetJetRadius())
            {
                coneParticles.push_back(track);
                coneParticlesId.push_back(counter);
                ptcone = ptcone + track->Pt();
            }
        }


        isOverlap = 0; isRanCone = 0;
        for (int counter=0; counter<coneParticles.size(); counter++){
                AliVParticle* part = coneParticles.at(counter);
                if (!part) continue;
                AliEmcalParticle* emcpart = dynamic_cast<AliEmcalParticle*>(part);
                if (emcpart) part = emcpart->GetTrack();
                Float_t parteta = part->Eta();
                Float_t partphi = part->Phi();

                    if ( CheckDeltaR(leadingJet,part)<JetCont->GetJetRadius()) {
                        fEtaRandConeOverlap->Fill(parteta, coneeta);
                        fPhiRandConeOverlap->Fill(partphi, conephi);
                        isOverlap = 1;
                        break;
                    }
                if(isOverlap) {   break;  }
        }

        if(!isOverlap) isRanCone = 1;
        if(isRanCone) break;
    }

   if(isRanCone){
      Bkgpt = JetCont->GetRhoVal()*TMath::Pi()*JetCont->GetJetRadius()*JetCont->GetJetRadius();
      deltaPt = ptcone - Bkgpt;

      fDeltaPT_excl_lead->Fill(deltaPt);
      //fDeltaPTCent_excl_lead->Fill(lPercentile,deltaPt);
      //fDeltaPTCentPtCone_excl_lead->Fill(lPercentile,ptcone);
      fDeltaPTRho_excl_lead->Fill(rho,deltaPt);
      fDeltaPTRhoPtCone_excl_lead->Fill(rho,ptcone);
      fRandConeEtaPhi_excl_lead->Fill(coneeta, conephi);

      Float_t leadJetPt = leadingJet->Pt();
      Float_t leadJetPtCorr = leadingJet->Pt() - leadingJet->Area()*rho;
      fDeltaPTLeadPt_excl_lead->Fill(leadJetPtCorr,deltaPt);
      fDeltaPTLeadPtPtCone_excl_lead->Fill(leadJetPtCorr,ptcone);
      fDeltaPTLeadPtRaw_excl_lead->Fill(leadJetPt,deltaPt);
      fDeltaPTLeadPtRawPtCone_excl_lead->Fill(leadJetPt,ptcone);
      //Double_t deltaEta = leadingJet->Eta() - coneeta;
     // Double_t deltaPhi = TMath::Abs(leadingJet->Phi() - conephi);
      //if(deltaPhi>TMath::Pi()) deltaPhi = (2*TMath::Pi())-deltaPhi;
      //fhRandomConeDeltaEtaDeltaPhi->Fill(deltaPhi,deltaEta);

  }

coneeta = 0; conephi = 0;
ptcone = 0.; Bkgpt = 0; deltaPt = 0;
//*************************************
// RANDOM CONES, TRANS PLANE =========
//*************************************

  coneeta = reta.Uniform(-(0.9-JetCont->GetJetRadius()),0.9-JetCont->GetJetRadius());
  Double_t double_radius = JetCont->GetJetRadius()*2;
  Double_t minPhi = leadingJet->Phi() + double_radius;
  Double_t maxPhi = leadingJet->Phi() + TMath::Pi() - double_radius;
  Double_t sign = rphi.Rndm();
  if (sign > 0.5) {
    minPhi += TMath::Pi();
    maxPhi += TMath::Pi();
  }
  conephi = TVector2::Phi_0_2pi(rphi.Uniform(minPhi, maxPhi));

  JetCont->ResetCurrentID();
  ParticlesCont->ResetCurrentID();
  track = 0;
  while ((track = static_cast<AliVTrack*>(ParticlesCont->GetNextAcceptParticle()))) {
       Double_t DeltaPhi = TMath::Abs(track->Phi()-conephi);
       if(DeltaPhi>TMath::Pi()) DeltaPhi = TMath::Abs(TMath::TwoPi() - DeltaPhi);
       Double_t DeltaEta = TMath::Abs(track->Eta()-coneeta);
       Double_t DeltaRCone = TMath::Sqrt((DeltaPhi*DeltaPhi)+(DeltaEta*DeltaEta));

       if(DeltaRCone<JetCont->GetJetRadius())
       {
           ptcone = ptcone + track->Pt();
       }
   }
   Bkgpt = JetCont->GetRhoVal()*TMath::Pi()*JetCont->GetJetRadius()*JetCont->GetJetRadius();
   deltaPt = ptcone - Bkgpt;

   fDeltaPT_trans->Fill(deltaPt);
   //fDeltaPTCent_trans->Fill(lPercentile,deltaPt);
   //fDeltaPTCentPtCone_trans->Fill(lPercentile,ptcone);
   fDeltaPTRho_trans->Fill(rho,deltaPt);
   fDeltaPTRhoPtCone_trans->Fill(rho,ptcone);
   fRandConeEtaPhi_trans->Fill(coneeta, conephi);

   Float_t leadJetPt = leadingJet->Pt();
   Float_t leadJetPtCorr = leadingJet->Pt() - leadingJet->Area()*rho;
   fDeltaPTLeadPt_trans->Fill(leadJetPtCorr,deltaPt);
   fDeltaPTLeadPtPtCone_trans->Fill(leadJetPtCorr,ptcone);
   fDeltaPTLeadPtRaw_trans->Fill(leadJetPt,deltaPt);
   fDeltaPTLeadPtRawPtCone_trans->Fill(leadJetPt,ptcone);
   //Double_t deltaEta = leadingJet->Eta() - coneeta;
   //Double_t deltaPhi = TMath::Abs(leadingJet->Phi() - conephi);
   //if(deltaPhi>TMath::Pi()) deltaPhi = (2*TMath::Pi())-deltaPhi;
   //fhRandomConeDeltaEtaDeltaPhi_trans->Fill(deltaPhi,deltaEta);


//*************************************
//*************************************

     /*
    if(ParticlesCont->GetNParticles()>0) fhstat->Fill(2);

    if(fCorrelationMethod==kConstituent)
    {
        if(ParticlesCont->GetNParticles()>0) ConstituentCorrelationMethod(kFALSE,aodEvent);
        if(fAnalyseDBkg==kTRUE && ParticlesContSB->GetNParticles()>0) ConstituentCorrelationMethod(kTRUE,aodEvent);
    }

    else if(fCorrelationMethod==kAngular)
    {
        if(fCandidateArray->GetEntriesFast()>0) AngularCorrelationMethod(kFALSE,aodEvent);
        if(fAnalyseDBkg==kTRUE && fSideBandArray->GetEntriesFast()>0) AngularCorrelationMethod(kTRUE,aodEvent);
    }
    */
}


   PostData(1,fOutput);
   return kTRUE;
}
void AliAnalysisTaskFlavourJetCorrelationsRC::ConstituentCorrelationMethod(Bool_t IsBkg, AliAODEvent* aodEvent)
{

    AliJetContainer* JetCont = 0;

    if(!IsBkg) JetCont = GetJetContainer(0);
    else JetCont = GetJetContainer(1);

    Double_t rho = 0;
    if(!JetCont->GetRhoName().IsNull()) rho = JetCont->GetRhoVal();

    AliEmcalJet *jet;

    GetHFJet(jet,IsBkg);

    if(jet)
    {
        if(fLocalRho) rho = fLocalRho->GetLocalVal(jet->Phi(),JetCont->GetJetRadius(),JetCont->GetRhoVal());
        FillDJetHistograms(jet,rho,IsBkg,aodEvent);
    }

}
void AliAnalysisTaskFlavourJetCorrelationsRC::AngularCorrelationMethod(Bool_t IsBkg, AliAODEvent* aodEvent)
{
    AliJetContainer* JetCont = GetJetContainer(0);

    Int_t ncand = 0;
    if(!IsBkg) ncand = fCandidateArray->GetEntriesFast();
    else ncand = fSideBandArray->GetEntriesFast();

    Double_t rho = 0;
    if(!JetCont->GetRhoName().IsNull()) rho = JetCont->GetRhoVal();

    for(Int_t icand = 0; icand<ncand; icand++)
    {
        AliVParticle* charm=0x0;
        if(!IsBkg) charm = (AliVParticle*)fCandidateArray->At(icand);
        else charm = (AliVParticle*)fSideBandArray->At(icand);
        if(!charm) continue;

        Int_t JetTag = AliEmcalJet::kD0;
        if (fCandidateType == kDstartoKpipi) JetTag = AliEmcalJet::kDStar;
        //loop over jets
        JetCont->ResetCurrentID();
        AliEmcalJet* jet=0;
        while ((jet = JetCont->GetNextJet()))
        {
            UInt_t rejectionReason = 0;
            Bool_t OKjet = JetCont->AcceptJet(jet, rejectionReason);
            if(!OKjet) continue;

            if(fLocalRho) rho = fLocalRho->GetLocalVal(jet->Phi(),JetCont->GetJetRadius(),JetCont->GetRhoVal());

            if(DeltaR(jet,charm,rho)<JetCont->GetJetRadius() && CheckDeltaR(jet,charm)<JetCont->GetJetRadius())
            {
                if(!IsBkg) jet->AddFlavourTag(JetTag);
                jet->AddFlavourTrack(charm);
                FillDJetHistograms(jet,rho,IsBkg,aodEvent);
            }
        }
    }

}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelationsRC::CreateResponseMatrix(AliEmcalJet* jet)
{
    AliJetContainer* JetContRec = GetJetContainer(0);

    Double_t rho = 0;
    if(!fUsePythia) rho = JetContRec->GetRhoVal();

    AliEmcalJet* MCjet = 0;

    FindMCJet(MCjet);

    if(!jet) AliDebug(2, "No Reconstructed Level Jet Found!");
    else if(!MCjet) AliDebug(2, "No Generated Level Jet Found!");
    else
    {
        if(fLocalRho) rho = fLocalRho->GetLocalVal(jet->Phi(),JetContRec->GetJetRadius(),JetContRec->GetRhoVal());
        AliVParticle *Drec = jet->GetFlavourTrack();
        AliVParticle *Dgen = MCjet->GetFlavourTrack();
        Double_t zRec = Z(Drec,jet,rho);
        Double_t zGen = Z(Dgen,MCjet,0);
        Double_t JetPtRec = jet->Pt() - rho*jet->Area();
        Double_t JetPtGen = MCjet->Pt();
        Double_t etaRec = jet->Eta();
        Double_t etaGen = MCjet->Eta();
        Double_t nTracksRec = jet->GetNumberOfTracks();
        Double_t nTracksGen = MCjet->GetNumberOfTracks();

        Double_t res = (JetPtRec-JetPtGen)/JetPtGen;
        fhRMRes -> Fill(res,JetPtGen);

        Double_t fillRM[10] = {zRec,JetPtRec,Drec->Pt(),etaRec,nTracksRec,zGen,JetPtGen,Dgen->Pt(),etaGen,nTracksGen};
        fResponseMatrix->Fill(fillRM,1.);
    }
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelationsRC::CreateMCResponseMatrix(AliEmcalJet* MCjet, AliAODEvent* aodEvent)
{

    if(!MCjet) AliDebug(2, "No Generated Level Jet Found!");

    Int_t pdgMeson = 413;
    if (fCandidateType == kD0toKpi) pdgMeson = 421;

    AliVParticle *Dgen = MCjet->GetFlavourTrack();
    Double_t zGen = Z(Dgen,MCjet,0);
    Double_t JetPtGen = MCjet->Pt();
    Double_t JetEtaGen = MCjet->Eta();
    Double_t DPtGen = Dgen->Pt();
    Double_t DEtaGen = Dgen->Eta();
    Double_t DYGen = Dgen->Y();
    Double_t nTracksGen = MCjet->GetNumberOfTracks();


    Double_t zRec = -999;
    Double_t JetPtRec = -999;
    Double_t DPtRec = -999;
    Double_t JetEtaRec = -999;
    Double_t DYRec = -999;
    Double_t nTracksRec = -999;

    AliJetContainer* JetContRec = GetJetContainer(0);
    if(JetContRec){
        AliEmcalJet* jet;
        GetHFJet(jet,kFALSE);

        if (jet){

            AliVParticle *Drec = jet->GetFlavourTrack();

            Double_t rho = 0;
            if(!fUsePythia){
                rho = JetContRec->GetRhoVal();
                if(fLocalRho) rho = fLocalRho->GetLocalVal(jet->Phi(),JetContRec->GetJetRadius(),JetContRec->GetRhoVal());
            }
            zRec = 0;
            if(rho>0) zRec = Z(Drec,jet,rho);
            else zRec = Z(Drec,jet);
            JetPtRec = jet->Pt() - rho*jet->Area();
            JetEtaRec = jet->Eta();
            nTracksRec = jet->GetNumberOfTracks();

            DPtRec = Drec->Pt();
            DYRec = Drec->Y();

            Bool_t bDInEMCalAcc=InEMCalAcceptance(Drec);
            Bool_t bJetInEMCalAcc=InEMCalAcceptance(jet);

            Double_t res = (JetPtRec-JetPtGen)/JetPtGen;
            fhRMRes -> Fill(res,JetPtGen);

            if(fCandidateType==kD0toKpi)
            {
                AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)Drec;
                FillHistogramsD0JetCorr(dzero,zRec,Drec->Pt(),JetPtRec, nTracksRec, kFALSE,bDInEMCalAcc,bJetInEMCalAcc,aodEvent);
            }
            else if(fCandidateType==kDstartoKpipi)
            {
                AliAODRecoCascadeHF* dstar=(AliAODRecoCascadeHF*)Drec;
                FillHistogramsDstarJetCorr(dstar,zRec,Drec->Pt(),JetPtRec, nTracksRec, kFALSE,bDInEMCalAcc,bJetInEMCalAcc);
            }

        } // if HF reco jet
    } // if jet cont reco

    //Double_t fillRM[8] = {zRec,JetPtRec,DPtRec,JetEtaRec,zGen,JetPtGen,DPtGen,JetEtaGen};
    Double_t fillRM[10] = {zRec,JetPtRec,DPtRec,DYRec,nTracksRec,zGen,JetPtGen,DPtGen,DYGen,nTracksGen};
    fResponseMatrix->Fill(fillRM,1.);

}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelationsRC::FillDJetHistograms(AliEmcalJet* jet, Double_t rho, Bool_t IsBkg, AliAODEvent* aodEvent)
{

    AliVParticle *Dmeson = jet->GetFlavourTrack(0);
    Double_t JetPtCorr = jet->Pt() - rho*jet->Area();
    Double_t z = 0;
    if(rho>0) z = Z(Dmeson,jet,rho);
    else z = Z(Dmeson,jet);

    if(IsBkg==kFALSE && fBuildRM==kTRUE) CreateResponseMatrix(jet);

    Bool_t bDInEMCalAcc=InEMCalAcceptance(Dmeson);
    Bool_t bJetInEMCalAcc=InEMCalAcceptance(jet);

    if(fCandidateType==kD0toKpi)
    {
        AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)Dmeson;
        FillHistogramsD0JetCorr(dzero,z,Dmeson->Pt(),JetPtCorr, jet->GetNumberOfTracks(), IsBkg,bDInEMCalAcc,bJetInEMCalAcc,aodEvent);
    }
    else if(fCandidateType==kDstartoKpipi)
    {
        AliAODRecoCascadeHF* dstar=(AliAODRecoCascadeHF*)Dmeson;
        FillHistogramsDstarJetCorr(dstar,z,Dmeson->Pt(),JetPtCorr, jet->GetNumberOfTracks(), IsBkg,bDInEMCalAcc,bJetInEMCalAcc);
    }


}

AliEmcalJet* AliAnalysisTaskFlavourJetCorrelationsRC::GetJetCone(Double_t radius, Double_t eta, Double_t phi)
{
  Double_t pt = 0;
  Double_t maxD2 = radius * radius;

  const Int_t MAX_CONSTITUENTS = 100;

  static std::array<Int_t, MAX_CONSTITUENTS> tracks;
  static std::array<Int_t, MAX_CONSTITUENTS> clusters;

  Int_t nParticles = 0; // placeholder
  //Int_t nParticles = SumParticles<AliParticleContainer, MAX_CONSTITUENTS>(pt, eta, phi, maxD2, fParticleCollArray, tracks);
 // Int_t nClusters = SumParticles<AliClusterContainer, MAX_CONSTITUENTS>(pt, eta, phi, maxD2, fClusterCollArray, clusters);

  AliEmcalJet* jet = new AliEmcalJet(pt, eta, phi, 0);
  jet->SetArea(maxD2*TMath::Pi());

  jet->SetNumberOfTracks(nParticles);
  //jet->SetNumberOfClusters(nClusters);

  for (Int_t i = 0; i < nParticles; i++) jet->AddTrackAt(tracks[i], i);
  //for (Int_t i = 0; i < nClusters; i++) jet->AddClusterAt(clusters[i], i);

  return jet;
}

/*
template <class T, Int_t MAX_CONSTITUENTS>
Int_t AliAnalysisTaskJetUEStudies::SumParticles(Double_t& pt, Double_t eta, Double_t phi, Double_t maxD2, std::map<std::string, T*>& CollArray, std::array<Int_t, MAX_CONSTITUENTS>& ConstList)
{
  auto IndexMap = T::GetEmcalContainerIndexMap();
  Int_t N = 0;
  for (auto coll : CollArray) {
    auto itcont = coll.second->accepted_momentum();
    for (auto itmom = itcont.begin(); itmom != itcont.end(); itmom++) {
      Double_t etaDiff = eta - itmom->first.Eta();
      Double_t phiDiff = AliEmcalContainer::RelativePhi(phi, itmom->first.Phi());
      Double_t d2 = etaDiff*etaDiff + phiDiff*phiDiff;
      if (d2 > maxD2) continue;
      ConstList[N] = IndexMap.GlobalIndexFromLocalIndex(coll.second->GetArray(), itmom.current_index());
      pt += itmom->first.Pt();
      N++;
      if (N >= MAX_CONSTITUENTS) {
        AliError(Form("Reached the maximum number of constituents = %d", MAX_CONSTITUENTS));
        return N;
      }
    }
  }
  return N;
}
*/

void AliAnalysisTaskFlavourJetCorrelationsRC::GetHFJet(AliEmcalJet*& jet, Bool_t IsBkg)
{
    AliJetContainer* JetCont = 0;

    if(!IsBkg) JetCont = GetJetContainer(0);
    else JetCont = GetJetContainer(1);

    AliParticleContainer *ParticlesCont = JetCont->GetParticleContainer();

    Bool_t JetIsHF = kFALSE;

    JetCont->ResetCurrentID();
    int counter = 0;
    while ((jet = JetCont->GetNextJet()))
    {
        counter++;
        UInt_t rejectionReason = 0;
        Bool_t OKjet = JetCont->AcceptJet(jet, rejectionReason);
        if(!OKjet) continue;

        Int_t JetTag = AliEmcalJet::kD0;
        TString recoDecayClassName("AliAODRecoDecayHF2Prong");
        if (fCandidateType == kDstartoKpipi)
        {
            JetTag = AliEmcalJet::kDStar;
            recoDecayClassName = "AliAODRecoCascadeHF";
        }
        //loop on jet particles
        Int_t ntrjet=  jet->GetNumberOfTracks();
        for(Int_t itrk=0;itrk<ntrjet;itrk++)
        {
            AliVParticle* jetTrk=jet->TrackAt(itrk,ParticlesCont->GetArray());
            if(!jetTrk) continue;
            AliEmcalParticle* emcpart = dynamic_cast<AliEmcalParticle*>(jetTrk);
            if(emcpart) jetTrk = emcpart->GetTrack();
            if(strcmp(jetTrk->ClassName(),recoDecayClassName)==0)
            {
                JetIsHF = kTRUE;
                if(!IsBkg) jet->AddFlavourTag(JetTag);
                jet->AddFlavourTrack(jetTrk);
                break;
            }
        } //end loop on jet tracks
        if(JetIsHF) break;
    } //end of jet loop

    if(!JetIsHF) jet = 0;

}
//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelationsRC::FindMCJet(AliEmcalJet*& mcjet)
{
    Bool_t HFMCjet = kFALSE;

    AliJetContainer* mcjets = 0;

    if(!fAnalyseDBkg) mcjets = GetJetContainer(1);
    else mcjets = GetJetContainer(2);

    AliParticleContainer *ParticlesCont = mcjets->GetParticleContainer();
    mcjets->ResetCurrentID();

    Int_t njet=0;

    while ((mcjet = mcjets->GetNextAcceptJet()))
    {
        njet++;
        //loop on jet particles
        Int_t ntrjet=  mcjet->GetNumberOfTracks();

        for(Int_t itrk=0;itrk<ntrjet;itrk++)
        {
            AliAODMCParticle* jetTrk=(AliAODMCParticle*)mcjet->TrackAt(itrk,ParticlesCont->GetArray());

            if(TMath::Abs(jetTrk->GetPdgCode())==fPDGmother)
            {
                HFMCjet=kTRUE;
                mcjet->AddFlavourTrack(jetTrk);
                break;
            }
        } //end loop on jet tracks
        if(HFMCjet==kTRUE) break;
    }

    if(!HFMCjet) mcjet=0;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelationsRC::Terminate(Option_t*)
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   Info("Terminate"," terminate");
   AliAnalysisTaskSE::Terminate();

   fOutput = dynamic_cast<AliEmcalList*> (GetOutputData(1));
   if (!fOutput) {
      printf("ERROR: fOutput not available\n");
      return;
   }
}

//_______________________________________________________________________________

void  AliAnalysisTaskFlavourJetCorrelationsRC::SetMassLimits(Double_t range, Int_t pdg){
   Float_t mass=0;
   Int_t abspdg=TMath::Abs(pdg);

   mass=TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
   // compute the Delta mass for the D*
   if(fCandidateType==kDstartoKpipi){
      Float_t mass1=0;
      mass1=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      mass = mass-mass1;
   }

   fMinMass = mass-range;
   fMaxMass = mass+range;

   AliInfo(Form("Setting mass limits to %f, %f",fMinMass,fMaxMass));
   if (fMinMass<0 || fMaxMass<=0 || fMaxMass<fMinMass) AliFatal("Wrong mass limits!\n");
}

//_______________________________________________________________________________

void  AliAnalysisTaskFlavourJetCorrelationsRC::SetMassLimits(Double_t lowlimit, Double_t uplimit){
   if(uplimit>lowlimit)
   {
      fMinMass = lowlimit;
      fMaxMass = uplimit;
   }
   else{
      printf("Error! Lower limit larger than upper limit!\n");
      fMinMass = uplimit - uplimit*0.2;
      fMaxMass = uplimit;
   }
}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskFlavourJetCorrelationsRC::SetD0WidthForDStar(Int_t nptbins,Float_t *width){
   if(nptbins>30) {
      AliInfo("Maximum number of bins allowed is 30!");
      return kFALSE;
   }
   if(!width) return kFALSE;
   for(Int_t ipt=0;ipt<nptbins;ipt++) fSigmaD0[ipt]=width[ipt];
   return kTRUE;
}

//_______________________________________________________________________________

Double_t AliAnalysisTaskFlavourJetCorrelationsRC::Z(AliVParticle* part,AliEmcalJet* jet, Double_t rho) const{

    Double_t p[3],pj[3];
    Bool_t okpp=part->PxPyPz(p);
    Bool_t okpj=jet->PxPyPz(pj);

    if(!okpp || !okpj)
    {
        printf("Problems getting momenta\n");
        return -999;
    }

    //Background Subtraction
    //It corrects the each component of the jet momentum for Z calculation

    pj[0] = jet->Px() - jet->Area()*(rho*TMath::Cos(jet->AreaPhi()));
    pj[1] = jet->Py() - jet->Area()*(rho*TMath::Sin(jet->AreaPhi()));
    pj[2] = jet->Pz() - jet->Area()*(rho*TMath::SinH(jet->AreaEta()));

    return Z(p,pj);
}
//_______________________________________________________________________________

Double_t AliAnalysisTaskFlavourJetCorrelationsRC::Z(AliVParticle* part,AliEmcalJet* jet) const{

    Double_t p[3],pj[3];
    Bool_t okpp=part->PxPyPz(p);
    Bool_t okpj=jet->PxPyPz(pj);

    if(!okpp || !okpj)
    {
        printf("Problems getting momenta\n");
        return -999;
    }

    return Z(p,pj);
}
//_______________________________________________________________________________
Double_t AliAnalysisTaskFlavourJetCorrelationsRC::Z(Double_t* p, Double_t *pj) const{

   Double_t pjet2=pj[0]*pj[0]+pj[1]*pj[1]+pj[2]*pj[2];
   Double_t z=(p[0]*pj[0]+p[1]*pj[1]+p[2]*pj[2])/(pjet2);
   return z;
}


//_______________________________________________________________________________
Double_t AliAnalysisTaskFlavourJetCorrelationsRC::ZT(Double_t* p, Double_t *pj) const{

   Double_t pjet2=pj[0]*pj[0]+pj[1]*pj[1];
   Double_t z=(p[0]*pj[0]+p[1]*pj[1])/(pjet2);
   return z;
}

//_______________________________________________________________________________

Bool_t  AliAnalysisTaskFlavourJetCorrelationsRC::DefineHistoForAnalysis(){

   // Statistics
   Int_t nbins=8;
   if(fUseMCInfo) nbins+=2;

   fhstat=new TH1I("hstat","Statistics",nbins,-0.5,nbins-0.5);
   fhstat->GetXaxis()->SetBinLabel(1,"N ev anal");
   fhstat->GetXaxis()->SetBinLabel(2,"N ev sel");
   fhstat->GetXaxis()->SetBinLabel(3,"N cand sel");
   fhstat->GetXaxis()->SetBinLabel(4,"N jets");
   fhstat->GetXaxis()->SetBinLabel(5,"N cand in jet");
   fhstat->GetXaxis()->SetBinLabel(6,"N jet rej");
   fhstat->GetXaxis()->SetBinLabel(7,"N cand sel & !jet");
   fhstat->GetXaxis()->SetBinLabel(8,"N jets & cand");
   if(fUseMCInfo) {
    fhstat->GetXaxis()->SetBinLabel(3,"N Signal sel & jet");
    fhstat->GetXaxis()->SetBinLabel(5,"N Signal in jet");
    fhstat->GetXaxis()->SetBinLabel(9,"N Bkg sel & jet");
    fhstat->GetXaxis()->SetBinLabel(10,"N Bkg in jet");
   }
   fhstat->SetNdivisions(1);
   fOutput->Add(fhstat);

   fhDleadStat = new TH1I("fhDleadStat","leading D-jet",3,-0.5,2.5);
   fhDleadStat->GetXaxis()->SetBinLabel(1,"Not lead. jet");
   fhDleadStat->GetXaxis()->SetBinLabel(2,"Lead. jet");
   fhDleadStat->GetXaxis()->SetBinLabel(3,"2nd Lead. jet");
   fhDleadStat->SetNdivisions(1);
   fOutput->Add(fhDleadStat);

   fhDleadStatJetPt = new TH2F("fhDleadStatJetPt","leading D-jet vs jet p_{T}",3,-0.5,2.5,400,0,200);
   fhDleadStatJetPt->GetXaxis()->SetBinLabel(1,"Not lead. jet");
   fhDleadStatJetPt->GetXaxis()->SetBinLabel(2,"Lead. jet");
   fhDleadStatJetPt->GetXaxis()->SetBinLabel(3,"2nd Lead. jet");
   fhDleadStatJetPt->GetYaxis()->SetTitle("p_{T,jet}^{lead}");
   //fhDleadStat->SetNdivisions(1);
   fOutput->Add(fhDleadStatJetPt);


   fhCentDjet=new TH1F("hCentDjet","Centrality",100,0,100);
   fOutput->Add(fhCentDjet);

   const Int_t nbinsmass=300;
   const Int_t nbinspttrack=400;
   const Int_t nbinsptjet=300;
   const Int_t nbinsptD=100;
   const Int_t nbinsphi=200;
   const Int_t nbinseta=100;

   //binning for THnSparse
   const Int_t nbinsSpsmass=120;
   const Int_t nbinsSpsptjet=400;
   const Int_t nbinsSpsptD=100;
   const Int_t nbinsSpsz=160;
   const Int_t nbinsSpsy=150;

   const Float_t pttracklims[2]={0.,100.};
   const Float_t ptjetlims[2]={-30.,120.};
   const Float_t ptDlims[2]={0.,50.};
   const Float_t zlims[2]={-1.2,2};
   const Float_t philims[2]={0.,6.3};
   const Float_t etalims[2]={-1.5,1.5};

   // jet related fistograms

   fhNjets = new TH1I("fhNjets","number of jets; nJets",40,-0.5,39.5);
   fhNjets->Sumw2();
   fhNaccjets = new TH1I("fhNaccjets","number of accepted jets; nAccJets",40,-0.5,39.5);
   fhNaccjets->Sumw2();
   fhPhiJetTrks    = new TH1F("hPhiJetTrks","Jet tracks #phi distribution; #phi (rad)",  nbinsphi,philims[0],philims[1]);
   fhPhiJetTrks->Sumw2();
   fhEtaJetTrks    = new TH1F("hEtaJetTrks","Jet tracks #eta distribution; #eta",  nbinseta,etalims[0],etalims[1]);
   fhEtaJetTrks->Sumw2();
   fhPtJetTrks     = new TH1F("hPtJetTrks",  "Jet tracks Pt distribution; p_{T} (GeV/c)",nbinspttrack,pttracklims[0],pttracklims[1]);
   fhPtJetTrks->Sumw2();

   fhPhiJet    = new TH1F("hPhiJet","Jet #phi distribution; #phi (rad)",  nbinsphi,philims[0],philims[1]);
   fhPhiJet->Sumw2();
   fhEtaJet    = new TH1F("hEtaJet","Jet #eta distribution; #eta", nbinseta,etalims[0],etalims[1]);
   fhEtaJet->Sumw2();
   fhPtJet      = new TH1F("hPtJet",  "Jet Pt distribution; p_{T}^{jet}-A*#rho(GeV/c)",nbinsptjet,ptjetlims[0],ptjetlims[1]);
   fhPtJet->Sumw2();
   fhPtJetRaw      = new TH1F("hPtJetRaw",  "Jet Pt distribution; p_{T} (GeV/c)",200,0,100);
   fhPtJetRaw->Sumw2();

   fhPtJetArea = new TH2F("fhPtJetArea","Jet Area vs Jet p_{T}",nbinsptjet,ptjetlims[0],ptjetlims[1],100,0,1);
   fhPtJetArea->GetXaxis()->SetTitle("p_{R}^{jet}-A*#rho (GeV/c)");
   fhPtJetArea->GetYaxis()->SetTitle("Area");
   fhPtJetArea->Sumw2();

   fhJetAreaCent = new TH2F("fhJetAreaCent","Jet Area vs centrality",100,0,100,100,0,1);
   fhJetAreaCent->GetXaxis()->SetTitle("centrality");
   fhJetAreaCent->GetYaxis()->SetTitle("Area");
   fhJetAreaCent->Sumw2();

   fOutput->Add(fhNjets);
   fOutput->Add(fhNaccjets);
   fOutput->Add(fhPhiJetTrks);
   fOutput->Add(fhEtaJetTrks);
   fOutput->Add(fhPtJetTrks);
   fOutput->Add(fhPhiJet);
   fOutput->Add(fhEtaJet);
   fOutput->Add(fhPtJet);
   fOutput->Add(fhPtJetRaw);
   fOutput->Add(fhPtJetArea);



   fhLeadPt = new TH1F("fhLeadPt","leading jet p_{T}; p_{T} (GeV/c)",400,0,200);
   fhLeadPt->Sumw2();
   fhLeadPt2 = new TH1F("fhLeadPt2","2nd leading jet p_{T}; p_{T} (GeV/c)",400,0,200);
   fhLeadPt2->Sumw2();
   fhLeadEta = new TH1F("fhLeadEta","leading jet #eta",300,-1.5,1.5);
   fhLeadEta->Sumw2();
   fhLeadEta2 = new TH1F("fhLeadEta2","2nd leading jet #eta",300,-1.5,1.5);
   fhLeadEta2->Sumw2();
   fhLeadPhi = new TH1F("fhLeadPhi","leading jet #phi",130,0,6.5);
   fhLeadPhi->Sumw2();
   fhLeadPhi2 = new TH1F("fhLeadPhi2","2nd leading jet #phi",130,0,6.5);
   fhLeadPhi2->Sumw2();
   fhLeadArea = new TH1F("fhLeadArea","leading jet area",120,0,1.2);
   fhLeadArea->Sumw2();
   fhLeadArea2 = new TH1F("fhLeadArea2","2nd leading jet area",120,0,1.2);
   fhLeadArea2->Sumw2();
   fhLeadDeltaEtaDeltaPhi = new TH2F("fhLeadDeltaEtaDeltaPhi","leading jets #Delta #eta #Delta #phi",320,0,3.2,240,-1.2,1.2);
   fhLeadDeltaEtaDeltaPhi->GetXaxis()->SetTitle("#Delta #phi");
   fhLeadDeltaEtaDeltaPhi->GetYaxis()->SetTitle("#Delta #eta");
   fhLeadDeltaEtaDeltaPhi->Sumw2();

   fOutput->Add(fhLeadPt);
   fOutput->Add(fhLeadPt2);
   fOutput->Add(fhLeadEta);
   fOutput->Add(fhLeadEta2);
   fOutput->Add(fhLeadPhi);
   fOutput->Add(fhLeadPhi2);
   fOutput->Add(fhLeadArea);
   fOutput->Add(fhLeadArea2);
   fOutput->Add(fhLeadDeltaEtaDeltaPhi);


   //============== UE ===================
   fhRhoMult = new TH2F("fhRhoMult","#rho vs mult",100,0,100,200,0,20);
   fhRhoMult->GetXaxis()->SetTitle("centrality");
   fhRhoMult->GetYaxis()->SetTitle("#rho (GeV/c*rad^{-1})");
   fhRhoMult->Sumw2();

   fhRhoLeadArea = new TH2F("fhRhoLeadArea","#rho vs lead area",120,0,1.2,200,0,20);
   fhRhoLeadArea->GetXaxis()->SetTitle("Area");
   fhRhoLeadArea->GetYaxis()->SetTitle("#rho (GeV/c*rad^{-1})");
   fhRhoLeadArea->Sumw2();

   fhRhoLeadPt = new TH2F("fhRhoLeadPt","#rho vs lead jet p_{T}",400,-20,180,200,0,20);
   fhRhoLeadPt->GetXaxis()->SetTitle("p_{T,lead}^{jet}-A*#rho (GeV/c)");
   fhRhoLeadPt->GetYaxis()->SetTitle("#rho (GeV/c*rad^{-1})");
   fhRhoLeadPt->Sumw2();

   fhRhoLeadPtRaw = new TH2F("fhRhoLeadPtRaw","#rho vs lead jet p_{T}",400,-20,180,200,0,20);
   fhRhoLeadPtRaw->GetXaxis()->SetTitle("p_{T,lead}^{jet} (GeV/c)");
   fhRhoLeadPtRaw->GetYaxis()->SetTitle("#rho (GeV/c*rad^{-1})");
   fhRhoLeadPtRaw->Sumw2();

   fhRhoLeadPtMult = new TH3F("fhRhoLeadPtMult","#rho vs mult vs lead jet p_{T}",400,-20,180,100,0,100,200,0,20);
   fhRhoLeadPtMult->GetXaxis()->SetTitle("p_{T,lead}^{jet} (GeV/c)");
   fhRhoLeadPtMult->GetYaxis()->SetTitle("centrality");
   fhRhoLeadPtMult->GetZaxis()->SetTitle("#rho (GeV/c*rad^{-1})");
   fhRhoLeadPtMult->Sumw2();

   fhRhoLeadPtMultRaw = new TH3F("fhRhoLeadPtMultRaw","#rho vs mult vs lead jet p_{T}",400,-20,180,100,0,100,200,0,20);
   fhRhoLeadPtMultRaw->GetXaxis()->SetTitle("p_{T,lead}^{jet}-A*#rho (GeV/c)");
   fhRhoLeadPtMultRaw->GetYaxis()->SetTitle("centrality");
   fhRhoLeadPtMultRaw->GetZaxis()->SetTitle("#rho (GeV/c*rad^{-1})");
   fhRhoLeadPtMultRaw->Sumw2();

   fOutput->Add(fhRhoMult);
   fOutput->Add(fhRhoLeadArea);
   fOutput->Add(fhRhoLeadPt);
   fOutput->Add(fhRhoLeadPtRaw);
   fOutput->Add(fhRhoLeadPtMult);
   fOutput->Add(fhRhoLeadPtMultRaw);


   //============= RANDOM CONES =====================================
   fDeltaPT = new TH1F("fDeltaPT", "#Delta p_{T};#Delta p_{T} (GeV/c);Entries", 600, -50., 100.);
   fDeltaPT->Sumw2();
   fOutput->Add(fDeltaPT);

   fDeltaPT_excl_lead = new TH1F("fDeltaPT_excl_lead", "#Delta p_{T}, excl. lead;#Delta p_{T} (GeV/c);Entries", 600, -50., 100.);
   fDeltaPT_excl_lead->Sumw2();
   fOutput->Add(fDeltaPT_excl_lead);

   fDeltaPT_trans = new TH1F("fDeltaPT_trans", "#Delta p_{T}, transverse plane;#Delta p_{T} (GeV/c);Entries", 600, -50., 100.);
   fDeltaPT_trans->Sumw2();
   fOutput->Add(fDeltaPT_trans);

   fDeltaPTCent = new TH2F("fDeltaPTCent","#Delta p_{T} vs cent",100,0,100,600,-50,100);
   fDeltaPTCent->Sumw2();
   fDeltaPTCent->GetXaxis()->SetTitle("centrality");
   fDeltaPTCent->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTCent);

   fDeltaPTCentPtCone = new TH2F("fDeltaPTCentPtCone","#Delta p_{T} vs cent",100,0,100,600,-50,100);
   fDeltaPTCentPtCone->Sumw2();
   fDeltaPTCentPtCone->GetXaxis()->SetTitle("centrality");
   fDeltaPTCentPtCone->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTCentPtCone);

   fDeltaPTCent_excl_lead = new TH2F("fDeltaPTCent_excl_lead","#Delta p_{T} vs cent, excl. lead",100,0,100,600,-50,100);
   fDeltaPTCent_excl_lead->Sumw2();
   fDeltaPTCent_excl_lead->GetXaxis()->SetTitle("centrality");
   fDeltaPTCent_excl_lead->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTCent_excl_lead);

   fDeltaPTCent_trans = new TH2F("fDeltaPTCent_trans","#Delta p_{T} vs cent, transverse plane",100,0,100,600,-50,100);
   fDeltaPTCent_trans->Sumw2();
   fDeltaPTCent_trans->GetXaxis()->SetTitle("centrality");
   fDeltaPTCent_trans->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTCent_trans);

   fDeltaPTCentPtCone_excl_lead = new TH2F("fDeltaPTCentPtCone_excl_lead","#Delta p_{T} vs cent, excl. lead",100,0,100,600,-50,100);
   fDeltaPTCentPtCone_excl_lead->Sumw2();
   fDeltaPTCentPtCone_excl_lead->GetXaxis()->SetTitle("centrality");
   fDeltaPTCentPtCone_excl_lead->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTCentPtCone_excl_lead);

   fDeltaPTCentPtCone_trans = new TH2F("fDeltaPTCentPtCone_trans","#Delta p_{T} vs cent, transverse plane",100,0,100,600,-50,100);
   fDeltaPTCentPtCone_trans->Sumw2();
   fDeltaPTCentPtCone_trans->GetXaxis()->SetTitle("centrality");
   fDeltaPTCentPtCone_trans->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTCentPtCone_trans);

   fDeltaPTRho = new TH2F("fDeltaPTRho","#Delta p_{T} vs #rho",100,0,10,600,-50,100);
   fDeltaPTRho->Sumw2();
   fDeltaPTRho->GetXaxis()->SetTitle("#rho (GeV/c*rad^{-1})");
   fDeltaPTRho->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTRho);

   fDeltaPTRhoPtCone = new TH2F("fDeltaPTRhoPtCone","#Delta p_{T} vs #rho",100,0,10,600,-50,100);
   fDeltaPTRhoPtCone->Sumw2();
   fDeltaPTRhoPtCone->GetXaxis()->SetTitle("#rho (GeV/c*rad^{-1})");
   fDeltaPTRhoPtCone->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTRhoPtCone);

   fDeltaPTRho_excl_lead = new TH2F("fDeltaPTRho_excl_lead","#Delta p_{T} vs #rho, excl. lead",100,0,10,600,-50,100);
   fDeltaPTRho_excl_lead->Sumw2();
   fDeltaPTRho_excl_lead->GetXaxis()->SetTitle("#rho (GeV/c*rad^{-1})");
   fDeltaPTRho_excl_lead->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTRho_excl_lead);

   fDeltaPTRho_trans = new TH2F("fDeltaPTRho_trans","#Delta p_{T} vs #rho, transverse plane",100,0,10,600,-50,100);
   fDeltaPTRho_trans->Sumw2();
   fDeltaPTRho_trans->GetXaxis()->SetTitle("#rho (GeV/c*rad^{-1})");
   fDeltaPTRho_trans->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTRho_trans);

   fDeltaPTRhoPtCone_excl_lead = new TH2F("fDeltaPTRhoPtCone_excl_lead","#Delta p_{T} vs #rho, excl. lead",100,0,10,600,-50,100);
   fDeltaPTRhoPtCone_excl_lead->Sumw2();
   fDeltaPTRhoPtCone_excl_lead->GetXaxis()->SetTitle("#rho (GeV/c*rad^{-1})");
   fDeltaPTRhoPtCone_excl_lead->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTRhoPtCone_excl_lead);

   fDeltaPTRhoPtCone_trans = new TH2F("fDeltaPTRhoPtCone_trans","#Delta p_{T} vs #rho, transverse plane",100,0,10,600,-50,100);
   fDeltaPTRhoPtCone_trans->Sumw2();
   fDeltaPTRhoPtCone_trans->GetXaxis()->SetTitle("#rho (GeV/c*rad^{-1})");
   fDeltaPTRhoPtCone_trans->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTRhoPtCone_trans);

   fDeltaPTLeadPt = new TH2F("fDeltaPTLeadPt","#Delta p_{T} lead jet",800,-20,180,600,-50,100);
   fDeltaPTLeadPt->Sumw2();
   fDeltaPTLeadPt->GetXaxis()->SetTitle("p_{T}^{lead}-A*#rho (GeV/c)");
   fDeltaPTLeadPt->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTLeadPt);

   fDeltaPTLeadPtRaw = new TH2F("fDeltaPTLeadPtRaw","#Delta p_{T} lead jet",800,-20,180,600,-50,100);
   fDeltaPTLeadPtRaw->Sumw2();
   fDeltaPTLeadPtRaw->GetXaxis()->SetTitle("p_{T}^{lead}(GeV/c)");
   fDeltaPTLeadPtRaw->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTLeadPtRaw);

   fDeltaPTLeadPtPtCone = new TH2F("fDeltaPTLeadPtPtCone","#Delta p_{T} lead jet",800,-20,180,600,-50,100);
   fDeltaPTLeadPtPtCone->Sumw2();
   fDeltaPTLeadPtPtCone->GetXaxis()->SetTitle("p_{T}^{lead}-A*#rho (GeV/c)");
   fDeltaPTLeadPtPtCone->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTLeadPtPtCone);

   fDeltaPTLeadPtRawPtCone = new TH2F("fDeltaPTLeadPtRawPtCone","#Delta p_{T} lead jet",800,-20,180,600,-50,100);
   fDeltaPTLeadPtRawPtCone->Sumw2();
   fDeltaPTLeadPtRawPtCone->GetXaxis()->SetTitle("p_{T}^{lead} (GeV/c)");
   fDeltaPTLeadPtRawPtCone->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTLeadPtRawPtCone);

   fDeltaPTLeadPt_excl_lead = new TH2F("fDeltaPTLeadPt_excl_lead","#Delta p_{T} lead jet, excl. lead",800,-20,180,600,-50,100);
   fDeltaPTLeadPt_excl_lead->Sumw2();
   fDeltaPTLeadPt_excl_lead->GetXaxis()->SetTitle("p_{T}^{lead}-A*#rho (GeV/c)");
   fDeltaPTLeadPt_excl_lead->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTLeadPt_excl_lead);

   fDeltaPTLeadPt_trans = new TH2F("fDeltaPTLeadPt_trans","#Delta p_{T} lead jet, transverse plane",800,-20,180,600,-50,100);
   fDeltaPTLeadPt_trans->Sumw2();
   fDeltaPTLeadPt_trans->GetXaxis()->SetTitle("p_{T}^{lead}-A*#rho (GeV/c)");
   fDeltaPTLeadPt_trans->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTLeadPt_trans);

   fDeltaPTLeadPtRaw_excl_lead = new TH2F("fDeltaPTLeadPtRaw_excl_lead","#Delta p_{T} lead jet, excl. lead",800,-20,180,600,-50,100);
   fDeltaPTLeadPtRaw_excl_lead->Sumw2();
   fDeltaPTLeadPtRaw_excl_lead->GetXaxis()->SetTitle("p_{T}^{lead} (GeV/c)");
   fDeltaPTLeadPtRaw_excl_lead->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTLeadPtRaw_excl_lead);

   fDeltaPTLeadPtRaw_trans = new TH2F("fDeltaPTLeadPtRaw_trans","#Delta p_{T} lead jet, transverse plane",800,-20,180,600,-50,100);
   fDeltaPTLeadPtRaw_trans->Sumw2();
   fDeltaPTLeadPtRaw_trans->GetXaxis()->SetTitle("p_{T}^{lead} (GeV/c)");
   fDeltaPTLeadPtRaw_trans->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTLeadPtRaw_trans);

   fDeltaPTLeadPtPtCone_excl_lead = new TH2F("fDeltaPTLeadPtPtCone_excl_lead","#Delta p_{T} lead jet, excl. lead",800,-20,180,600,-50,100);
   fDeltaPTLeadPtPtCone_excl_lead->Sumw2();
   fDeltaPTLeadPtPtCone_excl_lead->GetXaxis()->SetTitle("p_{T}^{lead}-A*#rho (GeV/c)");
   fDeltaPTLeadPtPtCone_excl_lead->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTLeadPtPtCone_excl_lead);

   fDeltaPTLeadPtPtCone_trans = new TH2F("fDeltaPTLeadPtPtCone_trans","#Delta p_{T} lead jet, transverse plane",800,-20,180,600,-50,100);
   fDeltaPTLeadPtPtCone_trans->Sumw2();
   fDeltaPTLeadPtPtCone_trans->GetXaxis()->SetTitle("p_{T}^{lead}-A*#rho (GeV/c)");
   fDeltaPTLeadPtPtCone_trans->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTLeadPtPtCone_trans);

   fDeltaPTLeadPtRawPtCone_excl_lead = new TH2F("fDeltaPTLeadPtRawPtCone_excl_lead","#Delta p_{T} lead jet, excl. lead",800,-20,180,600,-50,100);
   fDeltaPTLeadPtRawPtCone_excl_lead->Sumw2();
   fDeltaPTLeadPtRawPtCone_excl_lead->GetXaxis()->SetTitle("p_{T}^{lead} (GeV/c)");
   fDeltaPTLeadPtRawPtCone_excl_lead->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTLeadPtRawPtCone_excl_lead);

   fDeltaPTLeadPtRawPtCone_trans = new TH2F("fDeltaPTLeadPtRawPtCone_trans","#Delta p_{T} lead jet, transverse plane",800,-20,180,600,-50,100);
   fDeltaPTLeadPtRawPtCone_trans->Sumw2();
   fDeltaPTLeadPtRawPtCone_trans->GetXaxis()->SetTitle("p_{T}^{lead} (GeV/c)");
   fDeltaPTLeadPtRawPtCone_trans->GetYaxis()->SetTitle("#Delta p_{T} (GeV/c)");
   fOutput->Add(fDeltaPTLeadPtRawPtCone_trans);

   fRandConeEtaPhi = new TH2F("fRandConeEtaPhi","rand cone jet #eta vs jet #phi",300,-1.5,1.5,130,0,6.5);
   fRandConeEtaPhi->Sumw2();
   fRandConeEtaPhi->GetXaxis()->SetTitle("jet #eta");
   fRandConeEtaPhi->GetYaxis()->SetTitle("jet #phi");
   fOutput->Add(fRandConeEtaPhi);

   fRandConeEtaPhi_excl_lead = new TH2F("fRandConeEtaPhi_excl_lead","rand cone jet #eta vs jet #phi, excl. lead",300,-1.5,1.5,130,0,6.5);
   fRandConeEtaPhi_excl_lead->Sumw2();
   fRandConeEtaPhi_excl_lead->GetXaxis()->SetTitle("jet #eta");
   fRandConeEtaPhi_excl_lead->GetYaxis()->SetTitle("jet #phi");
   fOutput->Add(fRandConeEtaPhi_excl_lead);

   fRandConeEtaPhi_trans = new TH2F("fRandConeEtaPhi_trans","rand cone jet #eta vs jet #phi, transverse plane",300,-1.5,1.5,130,0,6.5);
   fRandConeEtaPhi_trans->Sumw2();
   fRandConeEtaPhi_trans->GetXaxis()->SetTitle("jet #eta");
   fRandConeEtaPhi_trans->GetYaxis()->SetTitle("jet #phi");
   fOutput->Add(fRandConeEtaPhi_trans);

   fhRandomConeDeltaEtaDeltaPhi = new TH2F("fhRandomConeDeltaEtaDeltaPhi","leading jet vs random cone #Delta #eta #Delta #phi",320,0,3.2,240,-1.2,1.2);
   fhRandomConeDeltaEtaDeltaPhi->GetXaxis()->SetTitle("#Delta #phi");
   fhRandomConeDeltaEtaDeltaPhi->GetYaxis()->SetTitle("#Delta #eta");
   fhRandomConeDeltaEtaDeltaPhi->Sumw2();
   fOutput->Add(fhRandomConeDeltaEtaDeltaPhi);

   fhRandomConeDeltaEtaDeltaPhi_trans = new TH2F("fhRandomConeDeltaEtaDeltaPhi_trans","leading jet vs random cone #Delta #eta #Delta #phi, trans plane",320,0,3.2,240,-1.2,1.2);
   fhRandomConeDeltaEtaDeltaPhi_trans->GetXaxis()->SetTitle("#Delta #phi");
   fhRandomConeDeltaEtaDeltaPhi_trans->GetYaxis()->SetTitle("#Delta #eta");
   fhRandomConeDeltaEtaDeltaPhi_trans->Sumw2();
   fOutput->Add(fhRandomConeDeltaEtaDeltaPhi_trans);

   fEtaRandConeOverlap = new TH2F("fEtaRandConeOverlap","#eta of ran cone overlapign with leading jet",240,-1.2,1.2,240,-1.2,1.2);
   fEtaRandConeOverlap->GetXaxis()->SetTitle("#eta_{cone}");
   fEtaRandConeOverlap->GetYaxis()->SetTitle("#eta_{jet trk}");
   fEtaRandConeOverlap->Sumw2();
   fOutput->Add(fEtaRandConeOverlap);

   fPhiRandConeOverlap = new TH2F("fPhiRandConeOverlap","#phi of ran cone overlapign with leading jet",130,0,6.5,130,0,6.5);
   fPhiRandConeOverlap->GetXaxis()->SetTitle("#phi_{cone}");
   fPhiRandConeOverlap->GetYaxis()->SetTitle("#phi_{jet trk}");
   fPhiRandConeOverlap->Sumw2();
   fOutput->Add(fPhiRandConeOverlap);

      if(fCandidateType==kDstartoKpipi)
      {
	 if(fAnalyseDBkg){
      	    fhDiffSideBand = new TH2F("hDiffSideBand","M(kpipi)-M(kpi) Side Band Background",nbinsmass,fMinMass,fMaxMass,nbinsptD, ptDlims[0],ptDlims[1]);
      	    fhDiffSideBand->SetStats(kTRUE);
      	    fhDiffSideBand->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV");
      	    fhDiffSideBand->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
      	    fhDiffSideBand->Sumw2();
      	    fOutput->Add(fhDiffSideBand);
      	 }

      	 fhPtPion = new TH1F("hPtPion","Primary pions candidates pt ",500,0,10);
      	 fhPtPion->SetStats(kTRUE);
      	 fhPtPion->GetXaxis()->SetTitle("GeV/c");
      	 fhPtPion->GetYaxis()->SetTitle("Entries");
      	 fhPtPion->Sumw2();
      	 fOutput->Add(fhPtPion);

      }


      fhInvMassptD = new TH2F("hInvMassptD","D (Delta R < Rjet) invariant mass distribution p_{T}^{j} > threshold",nbinsmass,fMinMass,fMaxMass,nbinsptD,ptDlims[0],ptDlims[1]);
      fhInvMassptD->SetStats(kTRUE);
      fhInvMassptD->GetXaxis()->SetTitle("mass (GeV)");
      fhInvMassptD->GetYaxis()->SetTitle("#it{p}_{t}^{D} (GeV/c)");
      fhInvMassptD->Sumw2();

      fOutput->Add(fhInvMassptD);

      if(fUseMCInfo){
	 fhInvMassptDbg = new TH2F("hInvMassptDbg","Background D (Delta R < Rjet) invariant mass distribution p_{T}^{j} > threshold",nbinsmass,fMinMass,fMaxMass,nbinsptD,ptDlims[0],ptDlims[1]);
      	 fhInvMassptDbg->GetXaxis()->SetTitle(Form("%s (GeV)",(fCandidateType==kDstartoKpipi) ? "M(Kpipi)" : "M(Kpi)"));
      	 fhInvMassptDbg->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
      	 fhInvMassptDbg->Sumw2();
      	 fOutput->Add(fhInvMassptDbg);

      }

    fhsDphiz=0x0; //definition below according to the switches

    if(fUseMCInfo){
        AliInfo("Creating a 9 axes container (MB background candidates)");
        const Int_t nAxis=9;
        const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsmass,nbinsSpsy, 50, 2, 2, 2};
        const Double_t minSparse[nAxis]={zlims[0],ptjetlims[0],ptDlims[0],fMinMass,etalims[0], 0,  -0.5,-0.5,-0.5};
        const Double_t maxSparse[nAxis]={zlims[1],ptjetlims[1],ptDlims[1],fMaxMass,etalims[1], 50, 1.5, 1.5 , 1.5};
        fNAxesBigSparse=nAxis;
        fhsDphiz=new THnSparseF("hsDphiz","Z, p_{T}^{jet}, p_{T}^{D}, mass, y^{D}, num.of.tracks, Bkg?, D in EMCal acc?, jet in EMCal acc?", nAxis, nbinsSparse, minSparse, maxSparse);

    } else{
        AliInfo("Creating a 6 axes container");
        const Int_t nAxis=6;
        const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsmass,nbinsSpsy,50};
        const Double_t minSparse[nAxis]={zlims[0],ptjetlims[0],ptDlims[0],fMinMass,etalims[0],0};
        const Double_t maxSparse[nAxis]={zlims[1],ptjetlims[1],ptDlims[1],fMaxMass,etalims[1],50};
        fNAxesBigSparse=nAxis;

        fhsDphiz=new THnSparseF("hsDphiz","Z, p_{T}^{jet}, p_{T}^{D}, mass, y^{D}, num.of.tracks", nAxis, nbinsSparse, minSparse, maxSparse);
    }

    if(!fhsDphiz) AliFatal("No THnSparse created");
    fhsDphiz->Sumw2();

    fOutput->Add(fhsDphiz);

    if(fBuildRM == kTRUE || fBuildRMEff == kTRUE)
    {
        const Int_t nRMAxis=10;
        const Int_t RMnbins[nRMAxis] = {nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsy,50,nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsy,50};
        const Double_t RMmin[nRMAxis]={zlims[0],ptjetlims[0],ptDlims[0],etalims[0],0,zlims[0],ptjetlims[0],ptDlims[0],etalims[0],0};
        const Double_t RMmax[nRMAxis]={zlims[1],ptjetlims[1],ptDlims[1],etalims[1],50,zlims[1],ptjetlims[1],ptDlims[1],etalims[1],50};
        fResponseMatrix = new THnSparseF("ResponseMatrix","z, jet pt, D pt, #eta^{jet}, num.of.tracks: Rec and Gen",nRMAxis,RMnbins,RMmin,RMmax);
        fResponseMatrix->Sumw2();
        fOutput->Add(fResponseMatrix);

        fhRMRes      = new TH2D("fhRMRes", "RM resolution",200,-1,1,nbinsptjet,ptjetlims[0],ptjetlims[1]);
        fhRMRes->Sumw2();
        fOutput->Add(fhRMRes);


    }

   PostData(1, fOutput);

   return kTRUE;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelationsRC::FillHistogramsD0JetCorr(AliAODRecoDecayHF* candidate, Double_t z, Double_t ptD, Double_t ptj, Double_t jetTracks, Bool_t IsBkg, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc, AliAODEvent* aodEvent){


   Double_t masses[2]={0.,0.};
   Int_t pdgdaughtersD0[2]={211,321};//pi,K
   Int_t pdgdaughtersD0bar[2]={321,211};//K,pi
   Int_t pdgMeson = 413;
   if (fCandidateType == kD0toKpi) pdgMeson = 421;

   masses[0]=candidate->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0); //D0
   masses[1]=candidate->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0bar); //D0bar

   Double_t *point=0x0;

   if(!fAnalyseDBkg)
   {
      point=new Double_t[6];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=masses[0];
      point[4]=candidate->Y(pdgMeson);
      point[5]=jetTracks;

   }
   else
   {
      point=new Double_t[9];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=masses[0];
      point[4]=candidate->Y(pdgMeson);
      point[5]=jetTracks;
      point[6]=static_cast<Double_t>(IsBkg ? 1 : 0);
      point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);

   }
    Int_t isselected=fCuts->IsSelected(candidate,AliRDHFCuts::kAll,aodEvent);
    if(isselected==1 || isselected==3)
    {

        fhInvMassptD->Fill(masses[0],ptD);
        fhsDphiz->Fill(point,1.);

    }
    if(isselected>=2)
    {
        fhInvMassptD->Fill(masses[1],ptD);
        point[3]=masses[1];
        fhsDphiz->Fill(point,1.);

    }
    delete[] point;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelationsRC::FillHistogramsDstarJetCorr(AliAODRecoCascadeHF* dstar,  Double_t z, Double_t ptD, Double_t ptj, Double_t jetTracks, Bool_t IsBkg, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc){

    AliAODTrack *softpi = (AliAODTrack*)dstar->GetBachelor();
    Double_t deltamass= dstar->DeltaInvMass();
    Int_t pdgMeson = 413;
    if (fCandidateType == kD0toKpi) pdgMeson = 421;

    fhPtPion->Fill(softpi->Pt());

    Double_t *point=0x0;
    if(!fAnalyseDBkg)
    {
        point=new Double_t[6];
        point[0]=z;
        point[1]=ptj;
        point[2]=ptD;
        point[3]=deltamass;
        point[4]=dstar->Y(pdgMeson);
        point[5]=jetTracks;
    }
    else
    {
        point=new Double_t[9];
        point[0]=z;
        point[1]=ptj;
        point[2]=ptD;
        point[3]=deltamass;
        point[4]=dstar->Y(pdgMeson);
        point[5]=jetTracks;
        point[6]=static_cast<Double_t>(IsBkg ? 1 : 0);
        point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
        point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
    }

    if(!point){
        AliError(Form("Numer of THnSparse entries %d not valid", fNAxesBigSparse));
        return;
    }


    fhInvMassptD->Fill(deltamass,ptD);
    fhsDphiz->Fill(point,1.);
    delete[] point;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelationsRC::FillHistogramsMCGenDJetCorr(Double_t z,Double_t ptD,Double_t ptjet, Double_t jetTracks, Double_t yD, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc){

    Double_t pdgmass=0;
    if(fCandidateType==kD0toKpi) pdgmass=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    if(fCandidateType==kDstartoKpipi) pdgmass=TDatabasePDG::Instance()->GetParticle(413)->Mass()-TDatabasePDG::Instance()->GetParticle(421)->Mass();

    Double_t *point=0x0;

    if(fNAxesBigSparse==6){
        point=new Double_t[6];
        point[0]=z;
        point[1]=ptjet;
        point[2]=ptD;
        point[3]=pdgmass;
        point[4]=yD;
        point[5]=jetTracks;
    }
    if(fNAxesBigSparse==9){
        point=new Double_t[9];
        point[0]=z;
        point[1]=ptjet;
        point[2]=ptD;
        point[3]=pdgmass;
        point[4]=yD;
        point[5]=jetTracks;
        point[6]=1;
        point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
        point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
    }

    if(!point){
        AliError(Form("Numer of THnSparse entries %d not valid", fNAxesBigSparse));
        return;
    }


    point[3]=pdgmass;
    fhsDphiz->Fill(point,1.);

    delete[] point;
}

//_______________________________________________________________________________

Float_t AliAnalysisTaskFlavourJetCorrelationsRC::DeltaR(AliEmcalJet *p1, AliVParticle *p2, Double_t rho) const {
   //Calculate DeltaR between p1 and p2: DeltaR=sqrt(Delataphi^2+DeltaEta^2)
   //It recalculates the eta-phi values if it was asked for background subtraction of the jets
   if(!p1 || !p2) return -1;

    Double_t phi1=p1->Phi(),eta1=p1->Eta();

    if (rho>0)
    {
            Double_t pj[3];
            Bool_t okpj=p1->PxPyPz(pj);
            if(!okpj){
                printf("Problems getting momenta\n");
                return -999;
            }
            pj[0] = p1->Px() - p1->Area()*(rho*TMath::Cos(p1->AreaPhi()));
            pj[1] = p1->Py() - p1->Area()*(rho*TMath::Sin(p1->AreaPhi()));
            pj[2] = p1->Pz() - p1->Area()*(rho*TMath::SinH(p1->AreaEta()));
            //Image of the function Arccos(px/pt) where pt = sqrt(px*px+py*py) is:
            //[0,pi]    if py > 0 and
            //[pi,2pi]  if py < 0
            phi1 = TMath::ACos(pj[0]/TMath::Sqrt(pj[0]*pj[0]+pj[1]*pj[1]));
            if(pj[1]<0) phi1 = 2*TMath::Pi() - phi1;
            eta1 = TMath::ASinH(pj[2]/TMath::Sqrt(pj[0]*pj[0]+pj[1]*pj[1]));
    }

   Double_t phi2 = p2->Phi(),eta2 = p2->Eta() ;

   Double_t dPhi=TMath::Abs(phi1-phi2);
   if(dPhi>TMath::Pi()) dPhi = (2*TMath::Pi())-dPhi;


   Double_t dEta=eta1-eta2;
   Double_t deltaR=TMath::Sqrt(dEta*dEta + dPhi*dPhi );
   return deltaR;

}
//_______________________________________________________________________________
Float_t AliAnalysisTaskFlavourJetCorrelationsRC::CheckDeltaR(AliEmcalJet *p1, AliVParticle *p2) const {
    //Calculate DeltaR between p1 and p2: DeltaR=sqrt(Delataphi^2+DeltaEta^2)
    if(!p1 || !p2) return -1;

    Double_t phi1=p1->Phi(),eta1=p1->Eta();

    Double_t phi2 = p2->Phi(),eta2 = p2->Eta() ;

    Double_t dPhi=TMath::Abs(phi1-phi2);
    if(dPhi>TMath::Pi()) dPhi = (2*TMath::Pi())-dPhi;


    Double_t dEta=eta1-eta2;
    Double_t deltaR=TMath::Sqrt(dEta*dEta + dPhi*dPhi );
    return deltaR;

}

//_______________________________________________________________________________

Int_t AliAnalysisTaskFlavourJetCorrelationsRC::IsDzeroSideBand(AliAODRecoCascadeHF *candDstar){

   Double_t ptD=candDstar->Pt();
   Int_t bin = fCuts->PtBin(ptD);
   if (bin < 0){
      // /PWGHF/vertexingHF/AliRDHFCuts::PtBin(Double_t) const may return values below zero depending on config.
      bin = 9999; // void result code for coverity (bin later in the code non-zero) - will coverity pick up on this?
      return -1;  // out of bounds
   }

   Double_t invM=candDstar->InvMassD0();
   Double_t mPDGD0=TDatabasePDG::Instance()->GetParticle(421)->Mass();

   Float_t fourSigmal= mPDGD0-4.*fSigmaD0[bin] , sixSigmal= mPDGD0-8.*fSigmaD0[bin];
   Float_t fourSigmar= mPDGD0+4.*fSigmaD0[bin] , sixSigmar= mPDGD0+8.*fSigmaD0[bin];

   if((invM>=sixSigmal && invM<fourSigmal) || (invM>fourSigmar && invM<=sixSigmar)) return 1;
   else return 0;

}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskFlavourJetCorrelationsRC::InEMCalAcceptance(AliVParticle *vpart){
   //check eta phi of a VParticle: return true if it is in the EMCal acceptance, false otherwise

   Double_t phiEMCal[2]={1.405,3.135},etaEMCal[2]={-0.7,0.7};
   Bool_t binEMCal=kTRUE;
   Double_t phi=vpart->Phi(), eta=vpart->Eta();
   if(phi<phiEMCal[0] || phi>phiEMCal[1]) binEMCal=kFALSE;
   if(eta<etaEMCal[0] || eta>etaEMCal[1]) binEMCal=kFALSE;
   return binEMCal;


}
