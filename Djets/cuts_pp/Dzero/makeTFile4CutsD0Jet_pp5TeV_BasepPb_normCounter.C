#include <Riostream.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>
#include <AliRDHFCutsD0toKpi.h>

/***************************************************************************/
/***************************************************************************/
// Cuts that p-Pb 2016 analysis
// Event selection changed to what is done in D2H
/***************************************************************************/
/***************************************************************************/


//void makeInputAliAnalysisTaskSED0Correlations_pPb(Bool_t isOnData = kTRUE){
void makeTFile4CutsD0Jet_pp5TeV_BasepPb_normCounter(Bool_t isOnData = kTRUE){

//____________________________________________________

AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
RDHFD0toKpi->SetName("D0toKpiCuts");
RDHFD0toKpi->SetTitle("Cuts for D0 analysis");

RDHFD0toKpi->SetSelectCandTrackSPDFirst(kTRUE,3.);
// PILE UP REJECTION
RDHFD0toKpi->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);

// EVENT CUTS
RDHFD0toKpi->SetMinVtxContr(1);
// MAX Z-VERTEX CUT
RDHFD0toKpi->SetMaxVtxZ(10.);
RDHFD0toKpi->SetCutOnzVertexSPD(0);

RDHFD0toKpi->SetTriggerMask(0);
RDHFD0toKpi->SetTriggerMask(AliVEvent::kINT7);
RDHFD0toKpi->SetTriggerClass("");

  //Quality tracks for daughters
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4); // default is 5
  //esdTrackCuts->SetMinNClustersTPC(120);
  //esdTrackCuts->SetMinNClustersTPC(70);
  //esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.1);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.3,1.e10);

  RDHFD0toKpi->AddTrackCuts(esdTrackCuts);

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

  RDHFD0toKpi->SetGlobalIndex(nvars,nptbins);
  RDHFD0toKpi->SetPtBins(nptbins+1,ptbins);
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
  Float_t d0Topomatic[nptbins] = {2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.};
  RDHFD0toKpi->Setd0MeasMinusExpCut(nptbins,d0Topomatic);

  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];

  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];
    }
  }

  RDHFD0toKpi->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);
  RDHFD0toKpi->SetUseSpecialCuts(kTRUE);
  RDHFD0toKpi->SetRemoveDaughtersFromPrim(kTRUE);

  for(Int_t iv=0;iv<nvars;iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand=NULL;

  //D0 pid settings
  Bool_t pidflag=kTRUE;
  RDHFD0toKpi->SetUsePID(pidflag);
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
  RDHFD0toKpi->SetPidHF(pidObj);
  RDHFD0toKpi->SetUsePID(kTRUE);
  RDHFD0toKpi->SetUseDefaultPID(kFALSE); //to use the AliAODPidHF
  RDHFD0toKpi->SetLowPt(kFALSE);
  //RDHFD0toKpi->SetMaximumPforPID(4.);

  RDHFD0toKpi->SetRemoveDaughtersFromPrim(kTRUE); //activate for pp
  RDHFD0toKpi->SetPidHF(pidObj);
  RDHFD0toKpi->SetUseDefaultPID(kFALSE);
  RDHFD0toKpi->SetLowPt(kFALSE);
  TString cent="";
  RDHFD0toKpi->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid

  cout<<"This is the object I'm going to save:"<<endl;
  RDHFD0toKpi->PrintAll();

  TFile* fout;
  if(isOnData) fout=new TFile("D0toKpiCuts_pp5TeV_pPbBase_newNorm.root","recreate");   //set this!!
  else fout=new TFile("D0toKpiCuts_pp5TeV_pPbBase_newNorm.root","recreate");   //set this!!

  fout->cd();
  RDHFD0toKpi->Write();
  fout->Close();


}
