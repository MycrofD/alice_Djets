void AddTasksUE( const UInt_t uTriggerMask = AliVEvent::kAny)
{
  
   const Bool_t bIsMC = kFALSE;
    
   //=============================================================================
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   
   if (!mgr) {
      ::Error("AddTasksUE.C::AddTasksUE", "No analysis manager to connect to.");
      return;
   }
   
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD") && !type.Contains("AOD")) {
      ::Error("AddTasksUE.C::AddTasksUE", "Task manager to have an ESD or AOD input handler.");
      return;
   }
   
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTasksUE.C::AddTasksUE", "This task requires an input event handler");
      return;
   }
   //=============================================================================
   
   
   
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

    
    //Jet task
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
    //Rho task
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");
    
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetUEStudies.C");
    gROOT->LoadMacro("AddTaskJetUEStudies.C");
    
    
    
        AliEmcalJetTask *jetAKT = AddTaskEmcalJet("usedefault","", AliJetContainer::antikt_algorithm,0.4,AliJetContainer::kChargedJet,0.15,0.30,0.005,AliJetContainer::pt_scheme,"Jet",0.,kTRUE,kFALSE);
        jetAKT->SelectCollisionCandidates(AliVEvent::kINT7); 
        jetAKT->SetNeedEmcalGeom(kFALSE);
        //jetAKT->SetForceBeamType(AliAnalysisTaskEmcal::kpA);
        jetAKT->GetTrackContainer(0)->SetTrackFilterType(AliEmcalTrackSelection::kCustomTrackFilter);
        jetAKT->GetTrackContainer(0)->SetAODFilterBits((1<<4)|(1<<9));

      AliEmcalJetTask *jetKT = AddTaskEmcalJet("usedefault","", AliJetContainer::kt_algorithm,0.4,AliJetContainer::kChargedJet,0.15,0.30,0.005,AliJetContainer::pt_scheme,"Jet",0.,kTRUE,kFALSE);
      jetKT->SelectCollisionCandidates(AliVEvent::kINT7); 
      jetKT->SetNeedEmcalGeom(kFALSE);
      jetKT->GetTrackContainer(0)->SetTrackFilterType(AliEmcalTrackSelection::kCustomTrackFilter);
      jetKT->GetTrackContainer(0)->SetAODFilterBits((1<<4)|(1<<9));

    // for pPb use rhoSparse
    AliAnalysisTaskRhoSparse  *rhotask = (AliAnalysisTaskRhoSparse*) AddTaskRhoSparse("Jet_KTChargedR040_tracks_pT0150_pt_scheme", "Jet_AKTChargedR040_tracks_pT0150_pt_scheme", "tracks","", "Rho", 0.4, "TPCFID", 0.01, 0., 0, 0x0, 2, kTRUE);
    rhotask->SelectCollisionCandidates(AliVEvent::kINT7);
    rhotask->SetVzRange(-10,10);
    rhotask->SetNeedEmcalGeom(kFALSE);
    
   /* AliJetContainer *jetCont = new AliJetContainer(AliJetContainer::kChargedJet, AliJetContainer::kt_algorithm, AliJetContainer::pt_scheme, 0.4, rhotask-       GetParticleContainer("tracks"), 0);
    jetCont->SetJetPtCut(1);  
    jetCont->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
    jetCont->SetName("Signal");
    rhotask->AdoptJetContainer(jetCont);*/


    //AliAnalysisTaskJetUEStudiesMy *JetUEtask = (AliAnalysisTaskJetUEStudiesMy*)AddTaskJetUEStudies("usedefault", "", 0.15, 0.30); 
    AliAnalysisTaskJetUEStudies *JetUEtask = (AliAnalysisTaskJetUEStudies*)AddTaskJetUEStudies("usedefault", "", 0.15, 0.30); 
    JetUEtask->SelectCollisionCandidates(AliVEvent::kAnyINT);
    JetUEtask->SetVzRange(-10, 10);
    AliJetContainer* jetCont2 = JetUEtask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 0.4, AliJetContainer::kTPCfid, "tracks", "");
    jetCont2->SetRhoName("Rho");

   
   
   return;
}
