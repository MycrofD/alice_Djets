void AddTasksFlavourJetMyData(const Int_t iCandType = 1 /*0 = D0, 1=Dstar...*/,
   const TString sCutFile = "cutsHF/D0toKpiCutsppRecVtxNoPileupRejNoEMCAL.root",
   const Double_t dJetPtCut   = 0.,
   const Double_t dJetAreaCut = 0.,
   const char *acctype = "TPC",
   const TString sRunPeriod = "LHC10b",
   const Int_t    uBeamType = 0,
   const UInt_t uTriggerMask = AliVEvent::kAny, /*for jets; the D mesons trigger is defined in the cut object*/
   const UInt_t bUseHFJet = kTRUE,
   const UInt_t bLeadHFJet = kFALSE,
   const Bool_t bIsMC = kFALSE,
   const Bool_t bIsReco = kFALSE,
   const Bool_t bIsMap = kFALSE,
   const Bool_t bRM = kFALSE,
   const Bool_t bRMEff = kFALSE,
   const Bool_t bPythia = kFALSE,
   const Bool_t isPrompt = kTRUE,
   TString sText="",/*completes the name of the candidate task lists*/
   Bool_t doBkg = kFALSE
   )
{
   const TString sInputTrkMC  = "MCParticlesSelected";
   const TString sInputTrkRec  = "tracks";
   const TString sUsedTrks  = "PicoTracks";
    //PicoTracks
   const TString sUsedClus  = "CaloClustersCorr";
   TString rhoName = "";
   TString rhoNameBkg = "";
   TString rhoNameMC = "";
   TString sInputTrk = bIsReco ? sInputTrkRec : sInputTrkMC;
   const Int_t iJetAlgo = 1;
   const Int_t iJetType = 1; /*0=Full, 1=Charged, 2=Neutral*/
   /*
   const Int_t    nRadius = 3;
   const Double_t aRadius[] = {  0.2,   0.4,   0.6  };
   const TString  sRadius[] = { "R02", "R04", "R06" };
   */
   const Int_t    nRadius = 1;
   const Double_t aRadius[] = {  0.3  };
   const TString  sRadius[] = { "R03" };


   //=============================================================================

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

   if (!mgr) {
      ::Error("AddTasksFlavourJet.C::AddTasksFlavourJet", "No analysis manager to connect to.");
      return;
   }

   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD") && !type.Contains("AOD")) {
      ::Error("AddTasksFlavourJet.C::AddTasksFlavourJet", "Task manager to have an ESD or AOD input handler.");
      return;
   }

   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTasksFlavourJet.C::AddTasksFlavourJet", "This task requires an input event handler");
      return;
   }
   //=============================================================================

    AliAnalysisTaskEmcal::TriggerType trType=AliAnalysisTaskEmcal::kND;

   UInt_t uAnaType = (((iJetType==0) ||     (iJetType==2)) ? 1 : 0);
   Int_t  iLeading =  ((iJetType==0) ? 3 : ((iJetType==1)  ? 0 : 1));
   Int_t leadHadType=0; /* 0=charged, 1=neutral, 2=both*/


    //Centrality Selection
   gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
   AliMultSelectionTask *taskMult = AddTaskMultSelection();

   //D mesons -- PID
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
   AliAnalysisTaskSE *taskRespPID = AddTaskPIDResponse(bIsMC,kFALSE,kTRUE,"1");

   gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
   // -- Physics selection task
   if(!bIsMC){
        gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
        //AliPhysicsSelectionTask *physSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kTRUE, uTriggerMask, 5, 5, 10, kTRUE, -1, -1, -1, -1);
        AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(kFALSE,kTRUE);
        if (!physSelTask) {
                cout << "no physSelTask";
                return;
        }
    }
    else {
        AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(kTRUE,kFALSE);
        if (!physSelTask) {
                cout << "no physSelTask";
                return;
        }
   }
   // --

 //   gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskImproveITS.C");
  // AliAnalysisTaskSEImproveITS *taskImprover = AddTaskImproveITS(kFALSE,"alien:///alice/cern.ch/user/a/afestant/filesForImprover/pPb2016/ITSgraphs_Current.root","alien:///alice/cern.ch/user/a/afestant/filesForImprover/pPb2016/ITSgraphs_NewAll-X0.3-Res4um.root",0 );

    //D meson filtering task
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskSEDmesonsFilterCJ.C");
    //D-jet correlation task
   //gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskDFilterAndCorrelations.C");

    //D meson filtering task
    gROOT->LoadMacro("AddTaskSEDmesonsFilterCJ.C");
    //D-jet correlation task
    gROOT->LoadMacro("AddTaskDFilterAndCorrelations.C");

    //Jet task
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
    //Rho task
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskLocalRho.C");


    //In Pb-Pb there are no events with more than 5 candidates. For pp or p-Pb this number is probably smaller
    for(Int_t i=0; i<1  ; i++)
    {
        TString TaskText = sText;
        TaskText += Form("N%d",i);

        //Filtering task. It will create a set of particles with the D meson instead of the daughters
        AliAnalysisTaskSEDmesonsFilterCJTest *filter = AddTaskSEDmesonsFilterCJ(iCandType,sCutFile,bIsMC,bIsReco,TaskText);
        filter->SetCombineDmesons(kTRUE);
        filter->SetMultipleCandidates(kTRUE); //Analyse one candidate per event
        filter->SetAnalysedCandidate(i); //Number of the candidate that will be analysed (0 = first candidate)
        filter->SetUseHFJet(bUseHFJet);

        // set filter bits if needed
        //This is the particle container with the tracks of the event
        AliTrackContainer* trackCont1 = filter->AddTrackContainer("tracks");
        trackCont1->SetClassName("AliAODTrack");
        trackCont1->SetTrackFilterType(AliEmcalTrackSelection::kCustomTrackFilter);
        trackCont1->SetAODFilterBits((1<<4)|(1<<9));


        TString candArrName = "Dcandidates";
        TString sbArrName = "DSBcandidates";
        TString DcandAndTracks = "DcandidatesAndTracks";
        TString DBkgAndTracks = "DSBcandidatesAndTracks";
        TString MCDcandAndTracks = "MCDcandidatesAndTracks";
        TString candName;
        if(iCandType==0) candName = "D0";
        if(iCandType==1) candName = "DStar";

        candArrName += candName;
        candArrName += TaskText;
        sbArrName += candName;
        sbArrName += TaskText;
        DcandAndTracks += candName;
        DcandAndTracks += TaskText;
        DBkgAndTracks += candName;
        DBkgAndTracks += TaskText;
        MCDcandAndTracks += candName;
        MCDcandAndTracks += TaskText;

        if (bIsMC) {
            candArrName += "MC";
            sbArrName  += "MC";
            DcandAndTracks += "MC";
            DBkgAndTracks += "MC";
            MCDcandAndTracks += "MC";
            if (bIsReco)
            {
                candArrName += "rec";
                sbArrName += "rec";
                DcandAndTracks += "rec";
                DBkgAndTracks += "rec";
                MCDcandAndTracks += "rec";

            }
        }

        TString AKTJet = "AKTJet";
        AKTJet += TaskText;
        TString KTJet = "KTJet";
        KTJet += TaskText;

        AliEmcalJetTask *taskFJDandTracks = AddTaskEmcalJet(DcandAndTracks,"",1,aRadius[0],AliJetContainer::kFullJet,0.15,0.30,0.005,1,AKTJet,0.,kFALSE,kFALSE);
        taskFJDandTracks->SelectCollisionCandidates(uTriggerMask);

        rhoName = "Rho";
        rhoName += TaskText;
        AliEmcalJetTask *taskFJ2 = AddTaskEmcalJet(DcandAndTracks,"",0,aRadius[0],AliJetContainer::kFullJet,0.15,0.30,0.005,1,KTJet,0.,kFALSE,kFALSE);
        taskFJ2->SelectCollisionCandidates(uTriggerMask);

        // for pPb use rhoSparse
        rhotask = (AliAnalysisTaskRhoSparse*) AddTaskRhoSparse(taskFJ2->GetName(), taskFJDandTracks->GetName(), DcandAndTracks,"", rhoName, aRadius[0], "TPCFID", 0.01, 0., 0, 0, 2, kTRUE);
        //rhotask = (AliAnalysisTaskRhoSparse*) AddTaskRhoSparse(taskFJ2->GetName(), taskFJDandTracks->GetName(), DcandAndTracks,"", rhoName, aRadius[0], "TPCFID", 0.01, 0., 0, 0, 1, kTRUE);
        rhotask->SelectCollisionCandidates(uTriggerMask);
        rhotask->SetVzRange(-10,10);

        //For Data. Comment this part if you run on Monte Carlo
        AliAnalysisTaskFlavourJetCorrelationsTest *CorrTask = AddTaskDFilterAndCorrelations(
                                                                                         iCandType,
                                                                                         sCutFile,
                                                                                         bIsMC,
                                                                                         bIsReco,
                                                                                         TaskText,
                                                                                         taskFJDandTracks->GetName(),
                                                                                         DcandAndTracks,
                                                                                         "",
                                                                                         rhoName,
                                                                                         "",
                                                                                         "",
                                                                                         "",
                                                                                         "",
                                                                                         aRadius[0],
                                                                                         dJetPtCut,
                                                                                         acctype,
                                                                                         dJetAreaCut,
                                                                                         AliAnalysisTaskFlavourJetCorrelationsTest::kConstituent);


      CorrTask->SetUseHFJet(bUseHFJet);
      CorrTask->SetHFLeadJet(bLeadHFJet);

    }

   return;
}
