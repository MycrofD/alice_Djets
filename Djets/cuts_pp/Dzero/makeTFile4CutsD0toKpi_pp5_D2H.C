#include <Riostream.h>
#include <TFile.h>
#include <AliRDHFCutsD0toKpi.h>
#include <TClonesArray.h>
#include <TParameter.h>


void makeInputAliAnalysisTaskSED0Mass(TString name){

  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  RDHFD0toKpi->SetName("D0toKpiCuts");
  RDHFD0toKpi->SetTitle("Cuts for D0 analysis");

 RDHFD0toKpi->SetSelectCandTrackSPDFirst(kTRUE,3.);
   // PILE UP REJECTION
   //
   RDHFD0toKpi->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
  //

  // EVENT CUTS
  RDHFD0toKpi->SetMinVtxContr(1);
 // MAX Z-VERTEX CUT
  RDHFD0toKpi->SetMaxVtxZ(10.);
//  RDHFD0toKpi->SetCutOnzVertexSPD(2);


 AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  //Trigger mask                                
RDHFD0toKpi->SetTriggerMask(0);               
RDHFD0toKpi->SetTriggerMask(AliVEvent::kINT7);
RDHFD0toKpi->SetTriggerClass("");             
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  //esdTrackCuts->SetMinNClustersITS(4); // default is 5
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny); 
 // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  RDHFD0toKpi->AddTrackCuts(esdTrackCuts);

  const Int_t nvars=11;

  const Int_t nptbins=17;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
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
  ptbins[13]=20.;
  ptbins[14]=24.;
  ptbins[15]=36.;
  ptbins[16]=50.;
  ptbins[17]=9999.;

   RDHFD0toKpi->SetPtBins(nptbins+1,ptbins);
   RDHFD0toKpi->SetGlobalIndex(nvars,nptbins);
   RDHFD0toKpi->SetPtBins(nptbins+1,ptbins);
                                                  //m    dca      cost*  ptk ptpi  d0k          d0pi       d0d0          cosp  cosxy normdxy 
Float_t cutsMatrixD0toKpiStand[nptbins][nvars]=  {{0.400,350.*1E-4, 0.8, 0.5, 0.5, 1000.*1E-4, 1000.*1E-4, -2000. *1E-8, 0.7,  0.,0.},/* pt<0.5*/
                                                  {0.400,350.*1E-4, 0.8, 0.5, 0.5, 1000.*1E-4, 1000.*1E-4, -2000. *1E-8, 0.7,  0.,0.},/* 0.5<pt<1*/
                                                  {0.400,300.*1E-4, 0.8, 0.4, 0.4, 1000.*1E-4, 1000.*1E-4, -25000.*1E-8, 0.8,  0.,0.},/* 1<pt<2 */
                                                  {0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, -20000.*1E-8, 0.9,   0.,0.},/* 2<pt<3 *///d0d0 e cosp
                                                  {0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, -12000.*1E-8, 0.9,   0.,0.},/* 3<pt<4 */
                                                  {0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, -8000. *1E-8, 0.85,  0.,0.},/* 4<pt<5 */
                                                  {0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, -8000. *1E-8, 0.85,  0.,0.},/* 5<pt<6 */
                                                  {0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, -8000. *1E-8, 0.85,  0.,0.},/* 6<pt<7 */
                                                  {0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, -7000. *1E-8, 0.85,  0.,0.},/* 7<pt<8 */
                                                  {0.400,300.*1E-4, 0.9, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, -5000. *1E-8, 0.85,  0.,0.},/* 8<pt<10 */
                                                  {0.400,300.*1E-4, 0.9, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, -5000. *1E-8, 0.85,  0.,0.},/* 10<pt<12 */
                                                  {0.400,300.*1E-4, 1.0, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 10000. *1E-8, 0.85,  0.,0.},/* 12<pt<16 */
                                                  {0.400,300.*1E-4, 1.0, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999.*1E-8, 0.85,  0.,0.},/* 16<pt<20 */
                                                  {0.400,300.*1E-4, 1.0, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999.*1E-8, 0.85,  0.,0.},/* 20<pt<24 */
                                                  {0.400,300.*1E-4, 1.0, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999.*1E-8, 0.85,  0.,0.},/* 24<pt<36 */
                                                  {0.400,300.*1E-4, 1.0, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999.*1E-8, 0.85,  0.,0.},/* 36<pt<50 */
                                                  {0.400,300.*1E-4, 1.0, 0.6, 0.6, 1000.*1E-4, 1000.*1E-4, 999999.*1E-8, 0.8,   0.,0.}};/* pt>50 */



  //setting my cut values
    //cuts order
    //       printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
    //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
    //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
    //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
    //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
    //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
    //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
    //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
    //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);
    //     printf("    |cosThetaPointXY| < %f\n",fD0toKpiCuts[9]);
    //     printf("    NormDecayLenghtXY    > %f\n",fD0toKpiCuts[10]);

 //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];

  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];
    }
  }


  //here add what you need! Find examples below
  //Event selection
  RDHFD0toKpi->SetUsePhysicsSelection(kTRUE);
  
   RDHFD0toKpi->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);
   RDHFD0toKpi->SetUseSpecialCuts(kTRUE);
   RDHFD0toKpi->SetRemoveDaughtersFromPrim(kTRUE);
  for(Int_t iv=0;iv<nvars;iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand=NULL;

  Bool_t pidflag=kTRUE;
  RDHFD0toKpi->SetUsePID(pidflag);
  if(pidflag) cout<<"PID is used"<<endl;
  else cout<<"PID is not used"<<endl;

    //pid settings
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
  pidObj->SetPCompatTOF(1.5);
  pidObj->SetSigmaForTPCCompat(3.);
  pidObj->SetSigmaForTOFCompat(3.);
  pidObj->SetOldPid(kFALSE);
  RDHFD0toKpi->SetPidHF(pidObj);


  //activate pileup rejection (for pp)
  //RDHFD0toKpi->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);
  //Do not recalculate the vertex
   RDHFD0toKpi->SetRemoveDaughtersFromPrim(kTRUE); //activate for pp
   RDHFD0toKpi->SetPidHF(pidObj);
   RDHFD0toKpi->SetUseDefaultPID(kFALSE);
   RDHFD0toKpi->SetLowPt(kFALSE);
  TString cent="";
  RDHFD0toKpi->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid

  //temporary
  //RDHFD0toKpi->SetFixRefs();

  RDHFD0toKpi->PrintAll();
  TFile* fout=new TFile(Form("D0toKpiCutsNewCut_%s.root",name.Data()), "recreate");   //set this!! 

  fout->cd();
  RDHFD0toKpi->Write();
  fout->Close();

  delete pidObj;
  pidObj=NULL;

  return;

}
 
