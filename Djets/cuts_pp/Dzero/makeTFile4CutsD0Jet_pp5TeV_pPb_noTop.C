#include <Riostream.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>

//Use:
//Set hard coded commentet with //set this!!
// root[] .L makeInputD0tasks.C++
// root[] makeInputAliAnalysisTaskSED0Correlations()
//similar macros for the other D mesons

//Author: Fabio Colamaria, fabio.colamaria@ba.infn.it

//macro to make a .root file which contains an AliRDHFCutsD0toKpi for AliAnalysisTaskSED0Mass task

/***************************************************************************/
/***************************************************************************/
//  
//  Selections for pPb 2016 cent-integrated preliminaries (SQM2017)
//
//  *** D0/Event CUTS ***
//  - Optimized cuts (std2013 pPb + tighter cosThPoint + NormLxy + topomatic)
//  - Pileup rejection with 5 contributors, <0.8 cm (F.Prino)
//  - ZNA centrality estimator (though 0-100 centrality used)
//  - Reference pT bins (0-0.5-1-2-3-4-5-6-7-8-10-12-16-24-36)
//  
//  ----> For MC: no pileup selection; don't use centrality
//
//  *** ASSOC. TRACK CUTS ***
//  - Filterbit0
//  - 2 ITS cls, no ITSrefit
//  - TPC cuts using crossed rows
//  - DCAxy = DCAz = 1 cm (not 0.25cm)
//  - SB ranges in reference pT bins
//  - Optimized pools: -10,-1.5,3.5,10 cm; 0,35,55,500 trkl
//  
/***************************************************************************/
/***************************************************************************/

void makeInputAliAnalysisTaskSED0Correlations_pPb(Bool_t isOnData = kTRUE){

//____________________________________________________

  //Set Centrality
  Float_t minc=0,maxc=100;

// Cuts for D0 cuts

  AliRDHFCutsD0toKpi* RDHFD0Corr=new AliRDHFCutsD0toKpi();
  RDHFD0Corr->SetName("D0toKpiCuts");
  RDHFD0Corr->SetTitle("Cuts for D0 analysis");

  // PILE UP REJECTION
  if(isOnData) {
    RDHFD0Corr->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);  	   //per DATI (spegni per MC)		
    RDHFD0Corr->ConfigurePileupCuts(5,0.8);  //-> Overloaded by mult dep pileup (F.Prino) //per DATI (spegni per MC) 
    //RDHFD0Corr->SetUseMultDepPileupCut(kTRUE);  	   //NON USATO PIU! //per DATI (spegni per MC)		 
  }

  //Event cuts
  RDHFD0Corr->SetMinVtxContr(1);
  RDHFD0Corr->SetMaxVtxZ(10.);
  RDHFD0Corr->SetCutOnzVertexSPD(2); 

  RDHFD0Corr->SetSelectCandTrackSPDFirst(kTRUE, 4); //Cristina 2016...

  //Trigger selection
  RDHFD0Corr->SetTriggerClass("");
  if(isOnData) RDHFD0Corr->SetTriggerMask(AliVEvent::kINT7);
  else RDHFD0Corr->SetTriggerMask(AliVEvent::kINT7); //A&E said to prefer kINT7 to kMB also in pPb MC

  //Quality tracks for daughters
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4); // default is 5
  //esdTrackCuts->SetMinNClustersTPC(120);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.3,1.e10);

  RDHFD0Corr->AddTrackCuts(esdTrackCuts);

  //D0 selection topological cuts
  const Int_t nptbins =14;
  const Double_t ptmax = 9999.;
  const Int_t nvars=11;
  Float_t ptbins[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=0.5;	
  ptbins[2]=1.;
  ptbins[3]=2.;
  ptbins[4]=3.;
  ptbins[5]=4.;
  ptbins[6]=5.;
  ptbins[7]=6.;
  ptbins[8]=7.;
  ptbins[9]=8.;
  ptbins[10]=10.;
  ptbins[11]=12.;
  ptbins[12]=16.;
  ptbins[13]=24.;
  ptbins[14]=ptmax;

  RDHFD0Corr->SetGlobalIndex(nvars,nptbins);
  RDHFD0Corr->SetPtBins(nptbins+1,ptbins);
												 //        dca    ThSt ptk ptp     d0k       d0p         d02      thp
  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.400,300.*1E-4,0.8,0.5,0.5,1000.*1E-4,1000.*1E-4,-50000.*1E-8,0.95,0.,5.},   /*BIN 0*/ /* pt<0.5*/ 
                                                  {0.400,300.*1E-4,0.8,0.5,0.5,1000.*1E-4,1000.*1E-4,-50000.*1E-8,0.95,0.,5.},   /*BIN 1*/ /* 0.5<pt<1*/
                                                  {0.400,300.*1E-4,0.8,0.5,0.5,1000.*1E-4,1000.*1E-4,-35000.*1E-8,0.95,0.,5.}, /*BIN 2*/ /* 1<pt<2 */
                                                  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-30000.*1E-8,0.95,0.,5.}, /*BIN 3*/ /* 2<pt<3 */
                                                  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-30000.*1E-8,0.95,0.,5.}, /*BIN 4*/ /* 3<pt<4 */
                                                  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-15000.*1E-8,0.95,0.,5.}, /*BIN 5*/ /* 4<pt<5 */
                                                  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-10000.*1E-8,0.95,0.,4.}, /*BIN 6*/ /* 5<pt<6 */
                                                  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.95,0.,4.},  /*BIN 7*/ /* 6<pt<7 */
                                                  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.95,0.,4.},  /*BIN 8*/ /* 7<pt<8 */
                                                  {0.400,300.*1E-4,0.9,0.7,0.7,1000.*1E-4,1000.*1E-4,-5000.*1E-8,0.95,0.,3.},  /*BIN 9*/ /* 8<pt<10 */
                                                  {0.400,300.*1E-4,0.9,0.7,0.7,1000.*1E-4,1000.*1E-4,-5000.*1E-8,0.95,0.,3.},  /*BIN 10*/ /* 10<pt<12 */
                                                  {0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4, 10000.*1E-8,0.95,0.,3.},  /*BIN 11*/ /* 12<pt<16 */
                                                  {0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4, 10000.*1E-8,0.90,0.,3.},  /*BIN 12*/ /* 16<pt<24 */
                                                  {0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4, 10000.*1E-8,0.90,0.,3.}}; /*BIN 13*/ /* pt>24 */
  
//topomatic cut
  Float_t d0Topomatic[nptbins] = {3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.};
  RDHFD0Corr->Setd0MeasMinusExpCut(nptbins,d0Topomatic);
  
  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
 
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];      
    }
  }
 
  RDHFD0Corr->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);
  RDHFD0Corr->SetUseSpecialCuts(kTRUE);
  RDHFD0Corr->SetRemoveDaughtersFromPrim(kTRUE);
  
  for(Int_t iv=0;iv<nvars;iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand=NULL;
 
  //D0 pid settings
  Bool_t pidflag=kTRUE;
  RDHFD0Corr->SetUsePID(pidflag);
  if(pidflag) cout<<"PID is used"<<endl;
  else cout<<"PID is not used"<<endl;

  // PID SETTINGS
  AliAODPidHF* pidObj=new AliAODPidHF();
  //pidObj->SetName("pid4D0");
  Int_t mode=1;
  const Int_t nlims=2;
  Double_t plims[nlims]={0.6,0.8}; //TPC limits in momentum [GeV/c]
  Bool_t compat=kTRUE; //effective only for this mode
  Bool_t asym=kTRUE;
  Double_t sigmas[5]={2.,1.,0.,3.,0.}; //to be checked and to be modified with new implementation of setters by Rossella
  pidObj->SetAsym(asym);// if you want to use the asymmetric bands in TPC
  pidObj->SetMatch(mode);
  pidObj->SetPLimit(plims,nlims);
  pidObj->SetSigma(sigmas);
  pidObj->SetCompat(compat);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  pidObj->SetPCompatTOF(2.);
  pidObj->SetSigmaForTPCCompat(3.);
  pidObj->SetSigmaForTOFCompat(3.);
  pidObj->SetOldPid(kFALSE);
  RDHFD0Corr->SetPidHF(pidObj);
  RDHFD0Corr->SetUsePID(kTRUE);
  RDHFD0Corr->SetUseDefaultPID(kFALSE); //to use the AliAODPidHF
  RDHFD0Corr->SetLowPt(kFALSE);
  //RDHFD0Corr->SetMaximumPforPID(4.);

  //activate pileup rejection (for pp)
//  RDHFD0Corr->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  //centrality selection
  TString cent="";
  RDHFD0Corr->SetUseCentrality(AliRDHFCuts::kCentOff);
	cent=Form("%.0f%.0f",minc,maxc);
  
 /* if(isOnData) {
    RDHFD0Corr->SetMinCentrality(minc);
    RDHFD0Corr->SetMaxCentrality(maxc);
    cent=Form("%.0f%.0f",minc,maxc);
    RDHFD0Corr->SetUseCentrality(AliRDHFCuts::kCentZNA); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
  } else {
	RDHFD0Corr->SetUseCentrality(AliRDHFCuts::kCentOff);
	cent=Form("%.0f%.0f",minc,maxc);
  }*/

  //temporary
  //RDHFD0Corr->SetFixRefs();

  cout<<"This is the object I'm going to save:"<<endl;
  RDHFD0Corr->PrintAll();

  TFile* fout;
  if(isOnData) fout=new TFile("D0toKpiCuts_pp5TeV_pPbCuts_topo3.root","recreate");   //set this!! 
  else fout=new TFile("D0toKpiCuts_pp5TeV_pPbCuts_topo3.root","recreate");   //set this!! 

  fout->cd();
  RDHFD0Corr->Write();
  fout->Close();


}

