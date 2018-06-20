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

#include "AliDJetRawYieldUncertaintyLocal.h"

//___________________________________________________________________________________________
AliDJetRawYieldUncertaintyLocal::AliDJetRawYieldUncertaintyLocal():
fFileInput(0x0),
fDmesonSpecie(kD0toKpi),
fDmesonLabel("Dzero"),
fYieldApproach(kEffScale),
fpTmin(0),
fpTmax(99),
fnDbins(0),
fDbinpTedges(0x0),
fDEffValues(0x0),
fMassPlot(0x0)
{


}

//___________________________________________________________________________________________
AliDJetRawYieldUncertaintyLocal::AliDJetRawYieldUncertaintyLocal(const AliDJetRawYieldUncertaintyLocal &source):
// copy constructor
fFileInput(source.fFileInput),
fDmesonSpecie(source.fDmesonSpecie),
fDmesonLabel(source.fDmesonLabel),
fYieldApproach(source.fYieldApproach),
fpTmin(source.fpTmin),
fpTmax(source.fpTmax),
fnDbins(source.fnDbins),
fDbinpTedges(source.fDbinpTedges),
fDEffValues(source.fDEffValues),
fMassPlot(source.fMassPlot)
{

}

//___________________________________________________________________________________________
AliDJetRawYieldUncertaintyLocal::~AliDJetRawYieldUncertaintyLocal() {
//destructor

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::SetDmesonSpecie(DMesonSpecies k){

  if(k<0 || k>1) {
    printf("Error! D meson specie not correctly set!\n");
    return kFALSE;
  } else if(k==0) fDmesonLabel="Dzero";
  else fDmesonLabel="Dstar";

  fDmesonSpecie=k;
  return kTRUE;

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::ExtractInputMassPlot(){

  std::cout << "Configuration:\nD meson: " << fDmesonLabel << "\nMethod (eff.scale/sideband): " << fYieldApproach << std::endl;

  fFileInput = TFile::Open(fFileNameInput.Data(),"read");
  if(!fFileInput){
    std::cout << "File " << fFileInput << " cannot be opened! check your file path!" << std::endl; return kFALSE;
  }

  Bool_t success = kTRUE;
  if(fDmesonSpecie==kD0toKpi && fYieldApproach==kEffScale) success = ExtractInputMassPlotDzeroEffScale();
  if(fDmesonSpecie==kD0toKpi && fYieldApproach==kSideband) success = ExtractInputMassPlotDzeroSideband();
  if(fDmesonSpecie==kDStarD0pi && fYieldApproach==kEffScale) success = ExtractInputMassPlotDstarEffScale();
  if(fDmesonSpecie==kDStarD0pi && fYieldApproach==kSideband) success = ExtractInputMassPlotDstarSideband();

  if(success) std::cout << "Extracted mass spectrum for fit variations" << std::endl;

  //std::cout << "Mass spectrum entries: " << fMassPlot->GetEntries() << std::endl;
  //fMassPlot->SaveAs("spettrotemp.root"); //DEBUG TEMP

  //fFileInput->Close();

  return success;

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::ExtractInputMassPlotDzeroEffScale() {
/*
    std::cout << "Extracting input mass plot: " << fpTmin << " to " << fpTmax << std::endl;

    TTree *tree = (TTree*)fFileInput->Get(fTreeName);
    AliAnalysisTaskDmesonJets::AliD0InfoSummary *brD = 0;
    AliAnalysisTaskDmesonJets::AliJetInfoSummary *brJet = 0;
    tree->SetBranchAddress(fDBranchName,&brD);
    tree->SetBranchAddress(fJetBranchName,&brJet);

    if(!tree || !brD || !brJet) {
      std::cout << "Error in setting the tree/branch names! Exiting..." << std::endl;
      return kFALSE;
    }

    TH1F* hmassjet[fnDbins];
    TH1F* hmassjet_scale[fnDbins];
    TH1F *hmass;

    for(int j=0; j<fnDbins; j++) {

      hmassjet[j] = new TH1F(Form("hmassjet_%d",j),"hmassjet",(fmassmax-fmassmin)/fmasswidth,fmassmin,fmassmax);
      hmassjet[j]->Sumw2();

      for(int k=0; k<tree->GetEntries(); k++) {
        tree->GetEntry(k);
        if(TMath::Abs(brJet->fEta)>0.5) continue;
        if(brJet->fPt >= fpTmin && brJet->fPt < fpTmax && brD->fPt >= fDbinpTedges[j] && brD->fPt < fDbinpTedges[j+1]) hmassjet[j]->Fill(brD->fInvMass);
      }

      hmassjet_scale[j] = (TH1F*)hmassjet[j]->Clone(Form("hmassjet%d_scale",j));
      hmassjet_scale[j]->Scale(1./fDEffValues[j]);

      if(!j) hmass = (TH1F*)hmassjet_scale[j]->Clone("hmass");
      else hmass->Add(hmassjet_scale[j]);

    }  //end of D-meson pT bin loop

    fMassPlot = (TH1D*)hmass->Clone("inputSpectrum");

    if(!fMassPlot) {
      std::cout << "Error in extracting the mass plot! Exiting..." << std::endl;
      return kFALSE;
    }
*/
    return kTRUE;

}
/*
//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::ExtractInputMassPlotDzeroSideband() {

    double jetmin = 5, jetmax = 24;

    std::cout << "Extracting input mass plot: " << fpTmin << " to " << fpTmax << std::endl;

    TTree *tree = (TTree*)fFileInput->Get(fTreeName);
    AliAnalysisTaskDmesonJets::AliD0InfoSummary *brD = 0;
    AliAnalysisTaskDmesonJets::AliJetInfoSummary *brJet = 0;
    tree->SetBranchAddress(fDBranchName,&brD);
    tree->SetBranchAddress(fJetBranchName,&brJet);

    if(!tree || !brD || !brJet) {
      std::cout << "Error in setting the tree/branch names! Exiting..." << std::endl;
      return kFALSE;
    }

    fMassPlot = new TH1D("inputSpectrum","hmassjet",(fmassmax-fmassmin)/fmasswidth,fmassmin,fmassmax);
    fMassPlot->Sumw2();

    for(int k=0; k<tree->GetEntries(); k++) {
      tree->GetEntry(k);
      if(TMath::Abs(brJet->fEta)>0.5) continue;
      if(brJet->fPt > jetmin && brJet->fPt <= jetmax && brD->fPt > fpTmin && brD->fPt <= fpTmax) fMassPlot->Fill(brD->fInvMass);
    }

    if(!fMassPlot) {
      std::cout << "Error in extracting the mass plot! Exiting..." << std::endl;
      return kFALSE;
    }

    return kTRUE;

}
*/

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::ExtractInputMassPlotDzeroSideband() {

    std::cout << "Extracting input mass plot: " << fpTmin << " to " << fpTmax << std::endl;

    double jetmin = 0, jetmax = 50;
    const int ND = 4;
    TDirectoryFile* dir=(TDirectoryFile*)fFileInput->Get(fDirName.Data());
    TList *histList;
    THnSparseF *sparse;
    TH3D* hInvMassptD;
    for(int i=0;i<ND; i++){
        histList =  (TList*)dir->Get(Form("%s%d",fListName.Data(),i));
        sparse = (THnSparseF*)histList->FindObject(fObjectName.Data());
        sparse->GetAxis(0)->SetRangeUser(fzmin,fzmax);
        //sparse->GetAxis(1)->SetRangeUser(jetmin,jetmax);
        if(i==0) hInvMassptD=(TH3D*)sparse->Projection(3,1,2);
        else hInvMassptD->Add((TH3D*)sparse->Projection(3,1,2));
    }
    if(!hInvMassptD) return kFALSE;

	// TCanvas *c = new TCanvas;
   //hInvMassptD->Draw();

    fMassPlot = (TH1D*)(hInvMassptD->ProjectionX("projX",hInvMassptD->GetYaxis()->FindBin(jetmin),hInvMassptD->GetYaxis()->FindBin(jetmax)-1,hInvMassptD->GetZaxis()->FindBin(fpTmin), hInvMassptD->GetZaxis()->FindBin(fpTmax)-1))->Clone("inputSpectrum");

    if(!fMassPlot) {
      std::cout << "Error in extracting the mass plot! Exiting..." << std::endl;
      return kFALSE;
    }

    return kTRUE;

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::ExtractInputMassPlotDstarEffScale() {

    std::cout << "Extracting input mass plot: " << fpTmin << " to " << fpTmax << std::endl;

    TDirectoryFile* dir=(TDirectoryFile*)fFileInput->Get(fDirName.Data());

    TList *histList =  (TList*)dir->Get(Form("%s0",fListName.Data()));
    THnSparseF *sparse = (THnSparseF*)histList->FindObject(fObjectName.Data());
    sparse->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD=(TH3D*)sparse->Projection(3,1,2);
    hInvMassptD->SetName("hInvMassptD");

    TList *histList_1 =  (TList*)dir->Get(Form("%s1",fListName.Data()));
    THnSparseF *sparse_1 = (THnSparseF*)histList_1->FindObject("hsDphiz");
    sparse_1->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD_1=(TH3D*)sparse_1->Projection(3,1,2);
    hInvMassptD_1->SetName("hInvMassptD_1");
    hInvMassptD->Add(hInvMassptD_1);

    TList *histList_2 =  (TList*)dir->Get(Form("%s2",fListName.Data()));
    THnSparseF *sparse_2 = (THnSparseF*)histList_2->FindObject("hsDphiz");
    sparse_2->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD_2=(TH3D*)sparse_2->Projection(3,1,2);
    hInvMassptD_2->SetName("hInvMassptD_2");
    hInvMassptD->Add(hInvMassptD_2);

    TList *histList_3 =  (TList*)dir->Get(Form("%s3",fListName.Data()));
    THnSparseF *sparse_3 = (THnSparseF*)histList_3->FindObject("hsDphiz");
    sparse_3->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD_3=(TH3D*)sparse_3->Projection(3,1,2);
    hInvMassptD_3->SetName("hInvMassptD_3");
    hInvMassptD->Add(hInvMassptD_3);

   /* TList *histList_4 =  (TList*)dir->Get(Form("%s4",fListName.Data()));
    THnSparseF *sparse_4 = (THnSparseF*)histList_4->FindObject("hsDphiz");
    sparse_4->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD_4=(TH3D*)sparse_4->Projection(3,1,2);
    hInvMassptD_4->SetName("hInvMassptD_4");
    hInvMassptD->Add(hInvMassptD_4);
*/
    //if(!hInvMassptD || !hInvMassptD_1 || !hInvMassptD_2 || !hInvMassptD_3 || !hInvMassptD_4) return kFALSE;
    if(!hInvMassptD || !hInvMassptD_1 || !hInvMassptD_2 || !hInvMassptD_3) return kFALSE;

    TH1F* hmassjet[fnDbins];
    TH1F* hmassjet_scale[fnDbins];
    TH1F *hmass = 0;

    for(int j=0; j<fnDbins; j++){

       hmassjet[j]=(TH1F*)hInvMassptD->ProjectionX(Form("hmassjet%d",j),hInvMassptD->GetYaxis()->FindBin(fpTmin),hInvMassptD->GetYaxis()->FindBin(fpTmax)-1,hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[j]),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[j+1])-1);

       hmassjet_scale[j] = (TH1F*)hmassjet[j]->Clone(Form("hmassjet%d_scale",j));
       hmassjet_scale[j]->Scale(1./fDEffValues[j]);

       if(!j) hmass = (TH1F*)hmassjet_scale[j]->Clone("hmass");
       else hmass->Add(hmassjet_scale[j]);

    }

    fMassPlot = (TH1D*)hmass->Clone("inputSpectrum");

    if(!fMassPlot) {
      std::cout << "Error in extracting the mass plot! Exiting..." << std::endl;
      return kFALSE;
    }

    return kTRUE;

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::ExtractInputMassPlotDstarSideband() {

    std::cout << "Extracting input mass plot: " << fpTmin << " to " << fpTmax << std::endl;

    double jetmin = 0, jetmax = 50;

    TDirectoryFile* dir=(TDirectoryFile*)fFileInput->Get(fDirName.Data());

    TList *histList =  (TList*)dir->Get(Form("%s0",fListName.Data()));
    THnSparseF *sparse = (THnSparseF*)histList->FindObject(fObjectName.Data());
    sparse->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD=(TH3D*)sparse->Projection(3,1,2);
    hInvMassptD->SetName("hInvMassptD");

    TList *histList_1 =  (TList*)dir->Get(Form("%s1",fListName.Data()));
    THnSparseF *sparse_1 = (THnSparseF*)histList_1->FindObject("hsDphiz");
    sparse_1->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD_1=(TH3D*)sparse_1->Projection(3,1,2);
    hInvMassptD_1->SetName("hInvMassptD_1");
    hInvMassptD->Add(hInvMassptD_1);

    TList *histList_2 =  (TList*)dir->Get(Form("%s2",fListName.Data()));
    THnSparseF *sparse_2 = (THnSparseF*)histList_2->FindObject("hsDphiz");
    sparse_2->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD_2=(TH3D*)sparse_2->Projection(3,1,2);
    hInvMassptD_2->SetName("hInvMassptD_2");
    hInvMassptD->Add(hInvMassptD_2);

    TList *histList_3 =  (TList*)dir->Get(Form("%s3",fListName.Data()));
    THnSparseF *sparse_3 = (THnSparseF*)histList_3->FindObject("hsDphiz");
    sparse_3->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD_3=(TH3D*)sparse_3->Projection(3,1,2);
    hInvMassptD_3->SetName("hInvMassptD_3");
    hInvMassptD->Add(hInvMassptD_3);

/*    TList *histList_4 =  (TList*)dir->Get(Form("%s4",fListName.Data()));
    THnSparseF *sparse_4 = (THnSparseF*)histList_4->FindObject("hsDphiz");
    sparse_4->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD_4=(TH3D*)sparse_4->Projection(3,1,2);
    hInvMassptD_4->SetName("hInvMassptD_4");
    hInvMassptD->Add(hInvMassptD_4);
*/
    //if(!hInvMassptD || !hInvMassptD_1 || !hInvMassptD_2 || !hInvMassptD_3 || !hInvMassptD_4) return kFALSE;
    if(!hInvMassptD || !hInvMassptD_1 || !hInvMassptD_2 || !hInvMassptD_3) return kFALSE;

	 TCanvas *c = new TCanvas;
  hInvMassptD_1->Draw();

    fMassPlot = (TH1D*)(hInvMassptD->ProjectionX("projX",hInvMassptD->GetYaxis()->FindBin(jetmin),hInvMassptD->GetYaxis()->FindBin(jetmax)-1,hInvMassptD->GetZaxis()->FindBin(fpTmin), hInvMassptD->GetZaxis()->FindBin(fpTmax)-1))->Clone("inputSpectrum");

    if(!fMassPlot) {
      std::cout << "Error in extracting the mass plot! Exiting..." << std::endl;
      return kFALSE;
    }

    return kTRUE;

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::RunMultiTrial(){

  std::cout << "Running MultiTrial on pT bin" << fpTmin << " to " << fpTmax << std::endl;

  Double_t massD;
  if(fDmesonSpecie==kD0toKpi) massD = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  if(fDmesonSpecie==kDStarD0pi) massD = TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass();
  if(fDebug) std::cout << "D-meson mass: " << massD << std::endl;

  TString outfilnam=Form("RawYieldVariations_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fpTmin,fpTmax);
  AliHFMultiTrials* mt = new AliHFMultiTrials();
  mt->SetSuffixForHistoNames("");
  mt->SetMass(massD);
  mt->SetSigmaGaussMC(fSigmaToFix);
  mt->SetUseFixSigFreeMean(fMeanSigmaVar[0]);
  mt->SetUseFixSigUpFreeMean(fMeanSigmaVar[1]);
  mt->SetUseFixSigDownFreeMean(fMeanSigmaVar[2]);
  mt->SetUseFreeS(fMeanSigmaVar[3]);
  mt->SetUseFixedMeanFreeS(fMeanSigmaVar[4]);
  mt->SetUseFixSigFixMean(fMeanSigmaVar[5]);
  mt->SetUseExpoBackground(fBkgVar[0]);
  mt->SetUseLinBackground(fBkgVar[1]);
  mt->SetUsePol2Background(fBkgVar[2]);
  mt->SetUsePol3Background(fBkgVar[3]);
  mt->SetUsePol4Background(fBkgVar[4]);
  mt->SetUsePol5Background(fBkgVar[5]);
  mt->SetUsePowerLawBackground(fBkgVar[6]);
  mt->SetUsePowerLawTimesExpoBackground(fBkgVar[7]);
  mt->ConfigureRebinSteps(fnRebinSteps,fRebinSteps);
  mt->ConfigureLowLimFitSteps(fnMinMassSteps,fMinMassSteps);
  mt->ConfigureUpLimFitSteps(fnMaxMassSteps,fMaxMassSteps);
printf("sigmaBC=%d\n",fnSigmaBC);
  if (fnSigmaBC) mt->ConfigurenSigmaBinCSteps(fnSigmaBC,fSigmaBC);
  mt->SetSaveBkgValue(kTRUE,fnSigmaSignReg);
  mt->SetDrawIndividualFits(kFALSE);
  if(fDebug>=2) mt->SetDrawIndividualFits(kTRUE);

  TCanvas* c0 = new TCanvas("c0","MassFit");

  if(fFitRefl) { //reflection treatment

    cout<< " Reflection template file: " << fReflFilenameInput.Data() << " - Reflection histo name:  " << fReflHistoName.Data() << endl;
    TFile *fRefl=TFile::Open(fReflFilenameInput.Data(),"read");
    if(!fRefl) {cout << " Reflection file not found! Exiting... " << endl; return kFALSE;}
    TH1F *hTemplRefl=(TH1F*)fRefl->Get(fReflHistoName.Data());
    if(!hTemplRefl) {cout << " Reflection histo not found! Exiting... " << endl; return kFALSE;}
    hTemplRefl->Draw();

        if(fFixRiflOverS<0){
          cout<< " MC signal file: " << fSigMCFilenameInput.Data() << " - MC signal histo name:  " << fReflHistoName.Data() << endl;
          TFile *fSigMC=TFile::Open(fSigMCFilenameInput.Data(),"read"); //needed if refl/sigMC ratio is not fixed
          if(!fSigMC) {cout << " MC signal file not found! Exiting... " << endl; return kFALSE;}
          TH1F *hSignMC=(TH1F*)fSigMC->Get(fSigMCHistoName.Data());
	  if(!hSignMC) {cout << " MC signal histo not found! Exiting... " << endl; return kFALSE;}
          fFixRiflOverS=hTemplRefl->Integral(hTemplRefl->FindBin(fReflRangeL+0.00001),hTemplRefl->FindBin(fReflRangeR-0.00001))/hSignMC->Integral(hSignMC->FindBin(fReflRangeL+0.00001),hSignMC->FindBin(fReflRangeR-0.00001));
          printf("Sign over Refl bin %1.1f-%1.1f = %f",fpTmin,fpTmax,fFixRiflOverS);
        }

    mt->SetTemplateRefl(hTemplRefl);
    mt->SetFixRefoS(fFixRiflOverS);

  }

  Bool_t isOK = mt->DoMultiTrials(fMassPlot,c0);

  TCanvas* cOut =new TCanvas("cry","All Trials");
  std::cout << "Has MultiTrial suceeded? " << isOK << std::endl;

  if(isOK){
    mt->DrawHistos(cOut);
    mt->SaveToRoot(outfilnam.Data(),"recreate");
  }
  else return kFALSE;

  CombineMultiTrialOutcomes();

  return isOK;

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::CombineMultiTrialOutcomes(){

  TString infilnam=Form("RawYieldVariations_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fpTmin,fpTmax);
  TFile* fil=new TFile(infilnam.Data());

  TString confCase[6]={"","","","","",""};
  TString bkgFunc[8]={"","","","","","","",""};

  Int_t nConfigCases=0;
  for(int i=0;i<6;i++) {
    if(i==0) confCase[nConfigCases] = "FixedS";
    if(i==1) confCase[nConfigCases] = "FixedSp20";
    if(i==2) confCase[nConfigCases] = "FixedSm20";
    if(i==3) confCase[nConfigCases] = "FreeS";
    if(i==4) confCase[nConfigCases] = "FixedMeanFreeS";
    if(i==5) confCase[nConfigCases] = "FixedMeanFixedS";
    if(fMeanSigmaVar[i]) nConfigCases++;
  }

  Int_t nBackFuncCases=0;
  for(int i=0;i<8;i++) {
    if(i==0) bkgFunc[nBackFuncCases] = "Expo";
    if(i==1) bkgFunc[nBackFuncCases] = "Lin";
    if(i==2) bkgFunc[nBackFuncCases] = "Pol2";
    if(i==3) bkgFunc[nBackFuncCases] = "Pol3";
    if(i==4) bkgFunc[nBackFuncCases] = "Pol4";
    if(i==5) bkgFunc[nBackFuncCases] = "Pol5";
    if(i==6) bkgFunc[nBackFuncCases] = "PowLaw";
    if(i==7) bkgFunc[nBackFuncCases] = "PowLawExpo";
    if(fBkgVar[i]) nBackFuncCases++;
  }

  Int_t totCases=nBackFuncCases*nConfigCases;
  if(totCases!=fnMask) {
    std::cout << "Error in the configuration of the mask! Mismatch with the number of active sigma/mean and bkg types settings!" << std::endl;
    return kFALSE;
  }

  TH1F* histo[totCases];
  std::cout << " Total cases (sigma/mean * bkg configs): " << totCases << std::endl;

  Int_t jh=0;
  for(Int_t iConf=0; iConf<nConfigCases; iConf++){
    for(Int_t iType=0; iType<nBackFuncCases; iType++){
      histo[jh++]=(TH1F*)fil->Get(Form("hRawYieldTrial%s%s",bkgFunc[iType].Data(),confCase[iConf].Data()));
      if(fDebug) std::cout << "Loading histo: " << Form("hRawYieldTrial%s%s",bkgFunc[iType].Data(),confCase[iConf].Data()) << std::endl;
    }
  }

  Int_t totTrials=0;
  Int_t totHistos=0;
  Int_t first[totCases];
  for(Int_t j=0; j<totCases; j++){
    if(fMask[j]){
      if(fDebug) std::cout << " case (with fMask active): " << j << std::endl;
      if(histo[j]){
	first[j]=totTrials;
	totTrials+=histo[j]->GetNbinsX();
	++totHistos;
      }else{
	fMask[j]=0;
      }
    }
  }

  TLine **vlines=new TLine*[totCases];

  if(fDebug) printf("Histos merged = %d, totTrials = %d\n",totHistos,totTrials);
  for(Int_t j=0; j<totCases; j++){
    if(fMask[j]){
      if(fDebug) printf("  %d) %s  -- %d \n",j,histo[j]->GetName(),first[j]);
      vlines[j]=new TLine(first[j],0.,first[j],50000.);
      vlines[j]->SetLineColor(kMagenta+2);
      vlines[j]->SetLineStyle(2);
    }
  }

  TH1F* hRawYieldAll=new TH1F("hRawYieldAll"," ; Trial # ; Raw Yield",totTrials,-0.5,totTrials-0.5);
  TH1F* hMeanAll=new TH1F("hMeanAll"," ; Trial # ; Gaussian mean",totTrials,-0.5,totTrials-0.5);
  TH1F* hSigmaAll=new TH1F("hSigmaAll"," ; Trial # ; Gaussian #sigma",totTrials,-0.5,totTrials-0.5);
  TH1F* hChi2All=new TH1F("hChi2All"," ; Trial # ; #chi^{2}",totTrials,-0.5,totTrials-0.5);
  TH1F* hBkgAll=new TH1F("hBkgAll"," ; Trial # ; Background under peak",totTrials,-0.5,totTrials-0.5);
  TH1F* hRawYieldDistAll=new TH1F("hRawYieldDistAll","  ; Raw Yield",25000,0.,25000.);
  TH1F* hStatErrDistAll=new TH1F("hStatErrDistAll","  ; Stat Unc on Yield",1000,0.,20000.);
  TH1F* hRelStatErrDistAll=new TH1F("hRelStatErrDistAll","  ; Rel Stat Unc on Yield",100,0.,1.);

  Double_t minYield=999999.;
  Double_t maxYield=0.;
  Double_t sumy[4]={0.,0.,0.,0.};
  Double_t sumwei[4]={0.,0.,0.,0.};
  Double_t sumerr[4]={0.,0.,0.,0.};
  Double_t counts=0.;
  Double_t wei[4];

  if(fDebug) printf("Building overall variation plot\n");
  for(Int_t j=0; j<totCases; j++){
    if(fMask[j]){
      if(fDebug) printf("-> Case %d\n",j);
      TString hmeanname=histo[j]->GetName();
      hmeanname.ReplaceAll("RawYield","Mean");
      TH1F* hmeant=(TH1F*)fil->Get(hmeanname.Data());

      TString hsigmaname=histo[j]->GetName();
      hsigmaname.ReplaceAll("RawYield","Sigma");
      TH1F* hsigmat=(TH1F*)fil->Get(hsigmaname.Data());

      TString hchi2name=histo[j]->GetName();
      hchi2name.ReplaceAll("RawYield","Chi2");
      TH1F* hchi2t=(TH1F*)fil->Get(hchi2name.Data());

      TString hbkgname=histo[j]->GetName();
      if(fDmesonSpecie==kD0toKpi) {
        hbkgname.ReplaceAll("RawYield","Bkg");
	if(j==0) printf("*** Using bkg values in real nSigma region, not in bin edges! (Salvatore) ***");
	if(j==0 && fDebug) getchar();
      } else {
        hbkgname.ReplaceAll("RawYield","BkgInBinEdges");
	if(j==0) printf("*** Using bkg values in bin edges! (Barbara/Antonio) ***");
	if(j==0 && fDebug) getchar();
      }
      TH1F* hbkg=(TH1F*)fil->Get(hbkgname.Data());

      //TCanvas *c=new TCanvas(Form("ciao%d",j));
      //histo[j]->Draw();
      for(Int_t ib=1; ib<=histo[j]->GetNbinsX(); ib++){
	Double_t ry=histo[j]->GetBinContent(ib);//rawyield
	Double_t ery=histo[j]->GetBinError(ib);

	Double_t pos=hmeant->GetBinContent(ib);
	Double_t epos=hmeant->GetBinError(ib);

	Double_t sig=hsigmat->GetBinContent(ib);
	Double_t esig=hsigmat->GetBinError(ib);

	Double_t chi2=hchi2t->GetBinContent(ib);
	Double_t bkg=hbkg->GetBinContent(ib);
	Double_t ebkg=hbkg->GetBinError(ib);

 	if(fDebug) std::cout<< " ry " << ry << " ery " << ery << " chi2 " << chi2 << std::endl;
	if(ry>0.001 && ery>(0.01*ry) && ery<(0.5*ry) && chi2<fChi2Cut){
	  hRawYieldDistAll->Fill(ry);
	  hStatErrDistAll->Fill(ery);
	  hRelStatErrDistAll->Fill(ery/ry);
	  hRawYieldAll->SetBinContent(first[j]+ib,ry);
	  hRawYieldAll->SetBinError(first[j]+ib,ery);
	  if(ry<minYield) minYield=ry;
	  if(ry>maxYield) maxYield=ry;
	  wei[0]=1.;
	  wei[1]=1./(ery*ery);
	  wei[2]=1./(ery*ery/(ry*ry));
	  wei[3]=1./(ery*ery/ry);
	  for(Int_t kw=0; kw<4; kw++){
	    sumy[kw]+=wei[kw]*ry;
	    sumerr[kw]+=wei[kw]*wei[kw]*ery*ery;
	    sumwei[kw]+=wei[kw];
	  }
	  counts+=1.;
	  hSigmaAll->SetBinContent(first[j]+ib,sig);
	  hSigmaAll->SetBinError(first[j]+ib,esig);
	  hMeanAll->SetBinContent(first[j]+ib,pos);
	  hMeanAll->SetBinError(first[j]+ib,epos);
	  hChi2All->SetBinContent(first[j]+ib,chi2);
	  hChi2All->SetBinError(first[j]+ib,0.000001);
	  hBkgAll->SetBinContent(first[j]+ib,bkg);
	  hBkgAll->SetBinError(first[j]+ib,ebkg);
	}
      }
    }
  }

  Double_t weiav[4]={0.,0.,0.,0.};
  Double_t eweiav[4]={0.,0.,0.,0.};
  for(Int_t kw=0; kw<4; kw++){
    if(sumwei[kw]>0.){
      weiav[kw]=sumy[kw]/sumwei[kw];
      eweiav[kw]=TMath::Sqrt(sumerr[kw])/sumwei[kw];
    }
  }

  hRawYieldAll->SetStats(0);
  hMeanAll->SetStats(0);
  hSigmaAll->SetStats(0);
  hChi2All->SetStats(0);
  hChi2All->SetMarkerStyle(7);
  hBkgAll->SetStats(0);
  //hMeanAll->SetMinimum(1.85);
  //hMeanAll->SetMaximum(1.88);
  hSigmaAll->SetMinimum(0.);
  //if(hSigmaAll->GetMaximum()<0.018) hSigmaAll->SetMaximum(0.018);

  TCanvas* call=new TCanvas("call","All",1400,800);
  call->Divide(4,2);
  call->cd(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.06);
  hSigmaAll->GetYaxis()->SetTitleOffset(1.7);
  hSigmaAll->Draw();
  call->cd(2);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.06);
  hMeanAll->GetYaxis()->SetTitleOffset(1.7);
  hMeanAll->Draw();
  call->cd(3);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.06);
  hChi2All->GetYaxis()->SetTitleOffset(1.7);
  hChi2All->Draw();
  call->cd(4);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.06);
  hBkgAll->GetYaxis()->SetTitleOffset(1.7);
  hBkgAll->Draw();
  call->cd(5);
  hRawYieldAll->SetTitle(Form("%s",infilnam.Data()));
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.06);
  Double_t newmax=1.25*(hRawYieldAll->GetMaximum()+hRawYieldAll->GetBinError(1));
  hRawYieldAll->GetYaxis()->SetTitleOffset(1.7);
  //hRawYieldAll->SetMaximum(newmax);
  hRawYieldAll->Draw();
  TLatex* tweimean[4];
  for(Int_t kw=0; kw<4; kw++){
    tweimean[kw]=new TLatex(0.16,0.84-0.06*kw,Form("<Yield>_{wei%d} = %.1f #pm %.1f\n",kw,weiav[kw],eweiav[kw]*sqrt(counts)));
    tweimean[kw]->SetNDC();
    tweimean[kw]->SetTextColor(4);
    //    tweimean[kw]->Draw();
  }

  for(Int_t j=1; j<totCases; j++){
    if(fMask[j]){
      vlines[j]->SetY2(newmax);
      vlines[j]->Draw("same");
    }
  }

  call->cd(6);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  hRawYieldDistAll->SetTitle(Form("%s",infilnam.Data()));
  hRawYieldDistAll->Draw();
  hRawYieldDistAll->GetXaxis()->SetRangeUser(minYield*0.8,maxYield*1.2);
  Double_t perc[3]={0.15,0.5,0.85}; // quantiles for +-1 sigma
  Double_t lim70[3];
  hRawYieldDistAll->GetQuantiles(3,lim70,perc);
  call->cd(7);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  Double_t aver=hRawYieldDistAll->GetMean();
  TLatex* tmean=new TLatex(0.15,0.94,Form("mean=%.1f",aver));
  tmean->SetNDC();
  tmean->Draw();
  TLatex* tmedian=new TLatex(0.15,0.86,Form("median=%.1f",lim70[1]));
  tmedian->SetNDC();
  tmedian->Draw();
  Double_t val=hRawYieldDistAll->GetRMS();
  TLatex* thrms=new TLatex(0.15,0.78,Form("rms=%.1f  (%.1f%%)",val,val/aver*100.));
  Double_t rms=val/aver*100.;
  thrms->SetNDC();
  thrms->Draw();
  TLatex* tmin=new TLatex(0.15,0.64,Form("min=%.1f      max=%.1f",minYield,maxYield));
  tmin->SetNDC();
  tmin->Draw();
  val=(maxYield-minYield)/sqrt(12);
  TLatex* trms=new TLatex(0.15,0.56,Form("(max-min)/sqrt(12)=%.1f  (%.1f%%)",val,val/aver*100.));
  trms->SetNDC();
  trms->Draw();
  val=(maxYield-aver)/sqrt(3);
  TLatex* tup=new TLatex(0.15,0.48,Form("(max-mean)/sqrt(3)=%.1f  (%.1f%%)",val,val/aver*100.));
  tup->SetNDC();
  tup->Draw();
  val=(aver-minYield)/sqrt(3);
  TLatex* tdw=new TLatex(0.15,0.40,Form("(mean-min)/sqrt(3)=%.1f  (%.1f%%)",val,val/aver*100.));
  tdw->SetNDC();
  tdw->Draw();
  TLatex* tl15=new TLatex(0.15,0.26,Form("15 percentile=%.1f",lim70[0]));
  tl15->SetNDC();
  tl15->Draw();
  TLatex* tl85=new TLatex(0.15,0.18,Form("85 percentile=%.1f",lim70[2]));
  tl85->SetNDC();
  tl85->Draw();
  val=(lim70[2]-lim70[0])/2.;
  TLatex* t1s=new TLatex(0.15,0.1,Form("70%% range =%.1f  (%.1f%%)",val,val/aver*100.));
  t1s->SetNDC();
  t1s->Draw();

  for(Int_t kw=0; kw<4; kw++){
    if(fDebug) printf("Weight %d: %.1f +- %.1f(stat) +- %.1f (syst)\n",kw,weiav[kw],eweiav[kw]*sqrt(counts),(maxYield-minYield)/sqrt(12));
  }

  TString outfilnam=Form("RawYieldSyst_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fpTmin,fpTmax);
  call->SaveAs(Form("%s",outfilnam.Data()));

  TH1D *hUnc = new TH1D(Form("hUnc_%1.1fto%1.1f",fpTmin,fpTmax),"hUnc",1,0.,1.);
  hUnc->SetBinContent(1,aver);
  hUnc->SetBinError(1,rms*aver/100.);
  hUnc->SaveAs(Form("Hist_%s",outfilnam.Data()));

  return kTRUE;

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::EvaluateUncertainty(){

  Bool_t success = kTRUE;
  if(fDmesonSpecie==kD0toKpi && fYieldApproach==kEffScale) success = EvaluateUncertaintyDzeroEffScale();
  if(fDmesonSpecie==kD0toKpi && fYieldApproach==kSideband) success = EvaluateUncertaintyDzeroSideband();
  if(fDmesonSpecie==kDStarD0pi && fYieldApproach==kEffScale) success = EvaluateUncertaintyDstarEffScale();
  if(fDmesonSpecie==kDStarD0pi && fYieldApproach==kSideband) success = EvaluateUncertaintyDstarSideband();

  if(success) std::cout << "Evaluated raw yield uncertainty" << std::endl;

  return success;

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::EvaluateUncertainty_CoherentTrialChoice(){

  Bool_t success = kTRUE;
  if(fDmesonSpecie==kD0toKpi) success = EvaluateUncertaintyDzeroSideband_CoherentTrialChoice();
  if(fDmesonSpecie==kDStarD0pi) success = EvaluateUncertaintyDstarSideband_CoherentTrialChoice();

  if(success) std::cout << "Evaluated raw yield uncertainty (CoherentTrialChoice)" << std::endl;

  return success;

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::EvaluateUncertaintyDzeroEffScale() {

  std::cout << "Jet spectrum pT bin edges: ";
  for(int i=0; i<fnJetbins; i++) std::cout << fJetbinpTedges[i] << " - ";
  std::cout << fJetbinpTedges[fnJetbins] << std::endl;

  fJetYieldUnc = new TH1D("JetRawYieldUncert","Raw yield uncertainty on jet spectrum - Dzero - Eff.scaling",fnJetbins,fJetbinpTedges);
  fJetYieldCentral = new TH1D("JetRawYieldCentral","Jet spectrum central values + syst yield uncertainty - Dzero - Eff.scaling",fnJetbins,fJetbinpTedges);

  TH1D *hUnc;

  //loop over the jet pT bins already extracted
  for(int ibin=0; ibin<fnJetbins; ibin++) {
    if(fDebug) std::cout << Form("Hist_RawYieldSyst_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fJetbinpTedges[ibin],fJetbinpTedges[ibin+1]) << std::endl;
    TFile *f = TFile::Open(Form("Hist_RawYieldSyst_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fJetbinpTedges[ibin],fJetbinpTedges[ibin+1]),"read");
    if(!f){
      std::cout << "Uncertainty file for bin " << fJetbinpTedges[ibin] <<" - "<< fJetbinpTedges[ibin+1] << " cannot be opened! Did you already evaluate it?" << std::endl;
      return kFALSE;
    }
    else {
      if(fDebug) std::cout << "Building uncertainty for bin " << fJetbinpTedges[ibin] <<" - "<< fJetbinpTedges[ibin+1] << std::endl;
    }

    hUnc = (TH1D*)f->Get(Form("hUnc_%1.1fto%1.1f",fJetbinpTedges[ibin],fJetbinpTedges[ibin+1]));
    if(!hUnc){
      std::cout << "Histogram with uncertainty from mass plot not found! Returning..." << std::endl;
      return kFALSE;
    }

    Double_t centrYield = hUnc->GetBinContent(1);
    Double_t rmsPct = hUnc->GetBinError(1);

    fJetYieldUnc->SetBinContent(ibin+1,rmsPct);
    fJetYieldCentral->SetBinContent(ibin+1,centrYield);
    fJetYieldCentral->SetBinError(ibin+1,rmsPct);
  }

  fJetYieldUnc->SetStats(kFALSE);
  fJetYieldUnc->Draw();
  fJetYieldUnc->SaveAs(Form("FinalRawYieldUncertainty_%s.root",fDmesonLabel.Data()));
  fJetYieldCentral->SetStats(kFALSE);
  fJetYieldCentral->Draw();
  fJetYieldCentral->SaveAs(Form("FinalRawYieldCentralPlusSystUncertainty_%s.root",fDmesonLabel.Data()));

  //print distribution of yields for each variation
  for(int ibin=0; ibin<fnJetbins; ibin++) {
    TFile *f2 = TFile::Open(Form("RawYieldSyst_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fJetbinpTedges[ibin],fJetbinpTedges[ibin+1]),"read");
    TCanvas *c = (TCanvas*)f2->Get("call");
    TH1F *hDist = (TH1F*)c->FindObject("hRawYieldDistAll");
    hDist->SetStats(kTRUE);
    hDist->SaveAs(Form("YieldDistribution_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fJetbinpTedges[ibin],fJetbinpTedges[ibin+1]));
  }

  return kTRUE;

}


//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::EvaluateUncertaintyDzeroSideband() {


	std::cout << "Jet spectrum pT bin edges: ";
	for (int i = 0; i<fnJetbins; i++) std::cout << fJetbinpTedges[i] << " - ";
	std::cout << fJetbinpTedges[fnJetbins] << std::endl;

	//Reload input object to extract projections in jet pT for each D-mass pTrange
	fFileInput = TFile::Open(fFileNameInput.Data(), "read");
	if (!fFileInput) {
		std::cout << "File " << fFileInput << " cannot be opened! check your file path!" << std::endl; return kFALSE;
	}

  double jetmin = 0, jetmax = 50;
  double fzmin = -2, fzmax = 2;
  const int ND = 4;
  TDirectoryFile* dir=(TDirectoryFile*)fFileInput->Get(fDirName.Data());
  TList *histList;
  THnSparseF *sparse;
  TH3D* hInvMassptD;
  for(int i=0;i<ND; i++){
      histList =  (TList*)dir->Get(Form("%s%d",fListName.Data(),i));
      sparse = (THnSparseF*)histList->FindObject(fObjectName.Data());
      sparse->GetAxis(0)->SetRangeUser(fzmin,fzmax);
      //sparse->GetAxis(1)->SetRangeUser(jetmin,jetmax);
      if(i==0) hInvMassptD=(TH3D*)sparse->Projection(3,1,2);
      else hInvMassptD->Add((TH3D*)sparse->Projection(3,1,2));
  }
  if(!hInvMassptD) return kFALSE;


	TRandom3 *gen = new TRandom3();
	gen->SetSeed(0);

	//define list of histograms, one per variation
	fJetSpectrSBVars = new TH1F*[fnMaxTrials];

	//for debug and thorough studies
	Double_t arrYldBinPerBin[fnDbins][fnJetbins][fnMaxTrials];
	for (Int_t i = 0; i<fnDbins; i++) for (Int_t j = 0; j<fnJetbins; j++) for (Int_t k = 0; k<fnMaxTrials; k++) arrYldBinPerBin[i][j][k] = 0;

	for (int iDbin = 0; iDbin<fnDbins; iDbin++) {

		//open file with summary of variations from MultiTrial - get histos of variations
		TFile *fileMult = TFile::Open(Form("RawYieldSyst_%s_%1.1fto%1.1f.root", fDmesonLabel.Data(), fDbinpTedges[iDbin], fDbinpTedges[iDbin + 1]), "read");
		if (!fileMult) {
			std::cout << "Uncertainty file for bin " << fDbinpTedges[iDbin] << " - " << fDbinpTedges[iDbin + 1] << " cannot be opened! Did you already evaluate it?" << std::endl;
			return kFALSE;
		}
		TCanvas *c = (TCanvas*)fileMult->Get("call");
		TH1F *hMean = (TH1F*)c->FindObject("hMeanAll");
		TH1F *hSigma = (TH1F*)c->FindObject("hSigmaAll");
		TH1F *hBkg = (TH1F*)c->FindObject("hBkgAll");

		Int_t extracted[fnMaxTrials];
		for (int iTrial = 0; iTrial<fnMaxTrials; iTrial++) extracted[iTrial] = -1;

		TH1F* hjetpt = new TH1F(Form("hjetpt%d", iDbin), "hJetPt_signReg", fnJetbins, fJetbinpTedges);
		TH1F* hjetpt_s1 = new TH1F(Form("hjetpt_s1%d", iDbin), "hJetPt_sb1", fnJetbins, fJetbinpTedges);
		TH1F* hjetpt_s2 = new TH1F(Form("hjetpt_s2%d", iDbin), "hJetPt_sb2", fnJetbins, fJetbinpTedges);
		TH1F* hjetpt_s = new TH1F(Form("hjetpt_s%d", iDbin), "hJetPt_sb2", fnJetbins, fJetbinpTedges);
		hjetpt->Sumw2();
		hjetpt_s1->Sumw2();
		hjetpt_s2->Sumw2();

		if (fDebug) std::cout << "Running bin pT(D) " << iDbin << std::endl;

		for (int iTrial = 0; iTrial<fnMaxTrials; iTrial++) {

			Bool_t extractOk = kFALSE;
			Int_t rnd = -1;

			do {  //just one time if fAllowRepetitions==kTRUE, repeat extraction till new number is obtained if fAllowRepetitions==kFALSE

				rnd = gen->Integer(hMean->GetNbinsX()) + 1;
				if (hSigma->GetBinContent(rnd)>0) extractOk = kTRUE;  //avoid 'empty' cases

				if (!fAllowRepetitions) { //check if already extracted for this pT(D) bin

					if (fnMaxTrials>hMean->GetNbinsX()) {
						std::cout << "Error! you set more set spectrum total variations than those done for pT(D) bin" << fDbinpTedges[iDbin] << " - " << fDbinpTedges[iDbin + 1] << "! ";
						std::cout << "Impossible to do without allowing repetitions! Exiting..." << std::endl;
						return kFALSE;
					}

					for (int j = 0; j<iTrial; j++) {
						if (rnd == extracted[j]) extractOk = kFALSE;
					}

				} //end of if(!fAllowRepetitions)

				extracted[iTrial] = rnd;

			} while (extractOk == kFALSE);

			Double_t mean = hMean->GetBinContent(rnd);
			Double_t sigma = hSigma->GetBinContent(rnd);
			Double_t bkg = hBkg->GetBinContent(rnd);

      TH1F *hmass = (TH1F*)(hInvMassptD->ProjectionX("projX",hInvMassptD->GetYaxis()->FindBin(jetmin),hInvMassptD->GetYaxis()->FindBin(jetmax)-1,hInvMassptD->GetZaxis()->FindBin(fpTmin), hInvMassptD->GetZaxis()->FindBin(fpTmax)-1));
      Float_t hmin = hmass->GetBinLowEdge(2);
      Float_t hmax = hmass->GetBinLowEdge(hmass->GetNbinsX());
			Float_t signal_l_min = mean - 9 * sigma;
			Float_t signal_l_max = mean - 4 * sigma;
			Float_t signal_u_min = mean + 4 * sigma;
			Float_t signal_u_max = mean + 9 * sigma;
			Float_t signal_c_min = mean - fnSigmaSignReg*sigma;
			Float_t signal_c_max = mean + fnSigmaSignReg*sigma;
      if(signal_l_min<hmin) signal_l_min = hmin;
      if(signal_u_max>hmax) signal_u_max = hmax;

      Double_t scaling;

      if(!fRebinDstarSB) { //do not rebin, just use projections

        //extract signal and sideband region spectra
        hjetpt=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_c_min), hInvMassptD->GetXaxis()->FindBin(signal_c_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
        hjetpt_s1=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_s1%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_l_min), hInvMassptD->GetXaxis()->FindBin(signal_l_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
        hjetpt_s2=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_s2%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_u_min), hInvMassptD->GetXaxis()->FindBin(signal_u_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
        hjetpt_s = (TH1F*)hjetpt_s1->Clone(Form("hjetpt_s%d",iDbin));
        hjetpt_s->Add(hjetpt_s2);

        // scale background from side bands to the background under the peak
      if(hjetpt_s->Integral()==0) {
            std::cout << "Error! At least one variation with no entries! Exiting..." << std::endl;
        return kFALSE;
            }
            scaling = bkg / hjetpt_s->Integral(hjetpt_s->FindBin(jetmin),hjetpt_s->FindBin(jetmax));

      }
      else { //rebin spectrum

        hjetpt = new TH1F(Form("hjetpt%d",iDbin),"hJetPt_signReg_Rebinned",fnJetbins,fJetbinpTedges);
        hjetpt_s1 = new TH1F(Form("hjetpt_s1%d",iDbin),"hJetPt_sb1_Rebinned",fnJetbins,fJetbinpTedges);
        hjetpt_s2 = new TH1F(Form("hjetpt_s2%d",iDbin),"hJetPt_sb2_Rebinned",fnJetbins,fJetbinpTedges);
        hjetpt_s = new TH1F(Form("hjetpt_s%d",iDbin),"hJetPt_sb2_Rebinned",fnJetbins,fJetbinpTedges);
        hjetpt->Sumw2();
        hjetpt_s1->Sumw2();
        hjetpt_s2->Sumw2();

        //extract signal and sideband region spectra
        TH1F* tmphjetpt=(TH1F*)hInvMassptD->ProjectionY(Form("tmphjetpt%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_c_min), hInvMassptD->GetXaxis()->FindBin(signal_c_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
      //  TH1F* tmphjetpt=(TH1F*)hInvMassptD->ProjectionY(Form("tmphjetpt%d",iDbin),hInvMassptD->GetXaxis()->FindBin(1.83), hInvMassptD->GetXaxis()->FindBin(1.93),hInvMassptD->GetZaxis()->FindBin(3), hInvMassptD->GetZaxis()->FindBin(4)-1);
        TH1F* tmphjetpt_s1=(TH1F*)hInvMassptD->ProjectionY(Form("tmphjetpt_s1%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_l_min), hInvMassptD->GetXaxis()->FindBin(signal_l_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
        TH1F* tmphjetpt_s2=(TH1F*)hInvMassptD->ProjectionY(Form("tmphjetpt_s2%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_u_min), hInvMassptD->GetXaxis()->FindBin(signal_u_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
        TH1F* tmphjetpt_s = (TH1F*)tmphjetpt_s1->Clone(Form("tmphjetpt_s%d",iDbin));
        tmphjetpt_s->Add(tmphjetpt_s2);

        // scale background from side bands to the background under the peak
        if(tmphjetpt_s->Integral()==0) {
          std::cout << "Error! At least one variation with no entries! Exiting reb..." << iDbin << std::endl;
          return kFALSE;
        }
        scaling = bkg / tmphjetpt_s->Integral(tmphjetpt_s->FindBin(jetmin),tmphjetpt_s->FindBin(jetmax)); //integral btw 0 and 50 pTJet (where you get the bkg from the mass plot)

        for(int j=1; j<=tmphjetpt->GetNbinsX(); j++) {
          Double_t centerbin = tmphjetpt->GetBinCenter(j); //bin of hjetpt corresponding to j-th bin of THnSparse projection
          hjetpt->Fill(centerbin,tmphjetpt->GetBinContent(j));
          hjetpt_s1->Fill(centerbin,tmphjetpt_s1->GetBinContent(j));
          hjetpt_s2->Fill(centerbin,tmphjetpt_s2->GetBinContent(j));
          hjetpt_s->Fill(centerbin,tmphjetpt_s->GetBinContent(j));
        }
        for(int j=1; j<=hjetpt->GetNbinsX(); j++) {
          hjetpt->SetBinError(j,TMath::Sqrt(hjetpt->GetBinContent(j)));
          hjetpt_s1->SetBinError(j,TMath::Sqrt(hjetpt_s1->GetBinContent(j)));
          hjetpt_s2->SetBinError(j,TMath::Sqrt(hjetpt_s2->GetBinContent(j)));
          hjetpt_s->SetBinError(j,TMath::Sqrt(hjetpt_s->GetBinContent(j)));
        }

      } //end of rebinning of spectrum

			hjetpt_s->Scale(scaling);

			// subtract background from signal jet
			hjetpt->Add(hjetpt_s, -1);
			hjetpt->SetMarkerColor(kGreen + 3);
			hjetpt->SetLineColor(kGreen + 3);

			// correct for D* efficiency
			hjetpt->Scale(1. / fDEffValues[iDbin]); // D efficiency
			hjetpt->SetMarkerColor(kBlue + 3);
			hjetpt->SetLineColor(kBlue + 3);

			// for every trial of every pT(D) bin, save the value of the yield, after eff scaling, in each pT(jet) bin (to study pT(D)->pT(jet) yield correlations)
			for (int l = 0; l<hjetpt->GetNbinsX(); l++) arrYldBinPerBin[iDbin][l][iTrial] = hjetpt->GetBinContent(l + 1);

			// add 'iDbin' pT(D) bin to total spectrum for variation 'iTrial'
			if (!iDbin) fJetSpectrSBVars[iTrial] = (TH1F*)hjetpt->Clone(Form("JetRawYieldUncert_%d", iTrial));
			else fJetSpectrSBVars[iTrial]->Add(hjetpt);

			hjetpt->Reset();
			hjetpt_s1->Reset();
			hjetpt_s2->Reset();
			hjetpt_s->Reset();

		} //end loop on trials for sideband approach

	} //end loop on pT(D) bins

	  //Now evaluate central value + rms in each pT(jet) bin to build the uncertainty
	Int_t nJetBins = fJetSpectrSBVars[0]->GetNbinsX();
	Double_t arrYld[nJetBins][fnMaxTrials];
	for (Int_t i = 0; i<nJetBins; i++) for (Int_t j = 0; j<fnMaxTrials; j++) arrYld[i][j] = 0;

	fJetPtBinYieldDistribution = new TH1F*[nJetBins];

	fJetYieldUnc = (TH1D*)fJetSpectrSBVars[0]->Clone("JetRawYieldUncert");
	fJetYieldUnc->Reset();
	fJetYieldUnc->SetTitle("Raw yield uncertainty on jet spectrum - Dzero - Sideband subtraction");
	fJetYieldCentral = (TH1D*)fJetSpectrSBVars[0]->Clone("JetRawYieldCentral");
	fJetYieldCentral->Reset();
	fJetYieldCentral->SetTitle("Jet spectrum central values + syst yield uncertainty - Dzero - Sideband subtraction");

	for (Int_t iJetbin = 0; iJetbin<nJetBins; iJetbin++) { //loop on jet spectrum pT bins

		fJetPtBinYieldDistribution[iJetbin] = new TH1F(Form("fJetPtBinYieldDistribution_Bin%d", iJetbin), "  ; Yield distribution", 50000, 0., 50000.);

		for (Int_t iTrial = 0; iTrial<fnMaxTrials; iTrial++) { //loop on trials and build array of variations for a given pT(jet) bin
			arrYld[iJetbin][iTrial] = fJetSpectrSBVars[iTrial]->GetBinContent(iJetbin + 1);
			fJetPtBinYieldDistribution[iJetbin]->Fill(arrYld[iJetbin][iTrial]);
		}

		Double_t mean = TMath::Mean(fnMaxTrials, arrYld[iJetbin]);
		Double_t rms = TMath::RMS(fnMaxTrials, arrYld[iJetbin]);
		if (fDebug) {
			std::cout << "Jet bin " << iJetbin << " (" << fJetSpectrSBVars[0]->GetXaxis()->GetBinLowEdge(iJetbin + 1) << "-" << fJetSpectrSBVars[0]->GetXaxis()->GetBinUpEdge(iJetbin + 1) << ")";
			std::cout << ": Mean = " << mean << ", RMS = " << rms << std::endl;
		}

		fJetYieldUnc->SetBinContent(iJetbin + 1, rms);
		fJetYieldCentral->SetBinContent(iJetbin + 1, mean);
		fJetYieldCentral->SetBinError(iJetbin + 1, rms);

		if (fRebinDstarSB) fJetPtBinYieldDistribution[iJetbin]->SaveAs(Form("YieldDistribution_%s_%1.1fto%1.1f.root", fDmesonLabel.Data(), fJetbinpTedges[iJetbin], fJetbinpTedges[iJetbin + 1]));

	}

	//fOut->Close();

	fJetYieldUnc->SetStats(kFALSE);
	fJetYieldUnc->Draw();
	fJetYieldUnc->SaveAs(Form("FinalRawYieldUncertainty_%s.root", fDmesonLabel.Data()));
	fJetYieldCentral->SetStats(kFALSE);
	fJetYieldCentral->Draw();
	fJetYieldCentral->SaveAs(Form("FinalRawYieldCentralPlusSystUncertainty_%s.root", fDmesonLabel.Data()));

	if (fDebug) {
		//ADVANCED - save distribution of final jet yields (summing all pT(D) bins) in a single plot
		TCanvas *cDistr = new TCanvas("cDistr", "cDistr", 900, 600);
		for (Int_t iTrial = 0; iTrial<fnMaxTrials; iTrial++) {
			for (int l = 0; l<fJetSpectrSBVars[iTrial]->GetNbinsX(); l++) fJetSpectrSBVars[iTrial]->SetBinError(l + 1, 0.0001);
			fJetSpectrSBVars[iTrial]->SetMarkerColor(iTrial + 1);
			fJetSpectrSBVars[iTrial]->SetLineColor(iTrial + 1);
			if (!iTrial) fJetSpectrSBVars[iTrial]->Draw();
			else fJetSpectrSBVars[iTrial]->Draw("same");
		}
		cDistr->SaveAs(Form("DistributionOfFinalYields_SBApproach_%s_AfterDbinSum.root", fDmesonLabel.Data()));

		//ADVANCED - save distribution of final jet yields from each single pT(D) bin in a single plot (one per each pT(D) bin)
		for (int iDbin = 0; iDbin<fnDbins; iDbin++) {
			TCanvas *cDistr1 = new TCanvas(Form("cDistr%d", iDbin), Form("cDistr%d", iDbin), 900, 600);
			TH1F** hJetSpectrFromSingleDbin = new TH1F*[fnMaxTrials];
			for (Int_t iTrial = 0; iTrial<fnMaxTrials; iTrial++) {
				hJetSpectrFromSingleDbin[iTrial] = (TH1F*)fJetSpectrSBVars[0]->Clone(Form("JetRawYieldDistr_Dbin%d", iDbin));
				for (int l = 0; l<hJetSpectrFromSingleDbin[iTrial]->GetNbinsX(); l++) {
					hJetSpectrFromSingleDbin[iTrial]->SetBinContent(l + 1, arrYldBinPerBin[iDbin][l][iTrial]);
					hJetSpectrFromSingleDbin[iTrial]->SetBinError(l + 1, 0.0001);
				}
				hJetSpectrFromSingleDbin[iTrial]->SetMarkerColor(iTrial + 1);
				hJetSpectrFromSingleDbin[iTrial]->SetLineColor(iTrial + 1);
				if (!iTrial) hJetSpectrFromSingleDbin[iTrial]->Draw();
				else hJetSpectrFromSingleDbin[iTrial]->Draw("same");
			}
			cDistr1->SaveAs(Form("DistributionOfFinalYields_SBApproach_%s_Bin%d.root", fDmesonLabel.Data(), iDbin));
		}

		//ADVANCED - save averages of final jet yields from each single pT(D) bin, with their RMS, without summing them, in a single plot
		TCanvas *cDistr2 = new TCanvas("cDistrAllAvgs", "cDistrAllAvgs", 900, 600);

		fJetYieldCentral->SetLineWidth(3);
		fJetYieldCentral->Draw();

		TH1F** hJetSpectrFromSingleDbin_Avg = new TH1F*[fnDbins];

		for (Int_t iDbin = 0; iDbin<fnDbins; iDbin++) {
			hJetSpectrFromSingleDbin_Avg[iDbin] = (TH1F*)fJetSpectrSBVars[0]->Clone(Form("JetRawYieldAvgDistr_%d", iDbin));

			for (Int_t iJetbin = 0; iJetbin<nJetBins; iJetbin++) { //loop on jet spectrum pT bins

				Double_t mean = TMath::Mean(fnMaxTrials, arrYldBinPerBin[iDbin][iJetbin]);
				Double_t rms = TMath::RMS(fnMaxTrials, arrYldBinPerBin[iDbin][iJetbin]);

				fJetYieldUnc->SetBinContent(iJetbin + 1, rms);
				hJetSpectrFromSingleDbin_Avg[iDbin]->SetBinContent(iJetbin + 1, mean);
				hJetSpectrFromSingleDbin_Avg[iDbin]->SetBinError(iJetbin + 1, rms);
			}

			hJetSpectrFromSingleDbin_Avg[iDbin]->SetMarkerColor(iDbin + 1);
			hJetSpectrFromSingleDbin_Avg[iDbin]->SetLineColor(iDbin + 1);
			hJetSpectrFromSingleDbin_Avg[iDbin]->Draw("same");
		}

		TLegend* leg = new TLegend(0.1, 0.7, 0.48, 0.9);
		leg->AddEntry(fJetYieldCentral, "Average after pT(D) bin sum", "pl");
		for (Int_t iDbin = 0; iDbin<fnDbins; iDbin++) leg->AddEntry(hJetSpectrFromSingleDbin_Avg[iDbin], Form("pt(D) %1.1f - %1.1f", fDbinpTedges[iDbin], fDbinpTedges[iDbin + 1]), "pl");
		leg->Draw();

		cDistr2->SaveAs(Form("AverageOfFinalYields_SBApproach_%s_AllDBins.root", fDmesonLabel.Data()));

	} //end of advanced plots

	return kTRUE;

}


//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::EvaluateUncertaintyDzeroSideband_CoherentTrialChoice() {

  std::cout << "Jet spectrum pT bin edges: ";
  for(int i=0; i<fnJetbins; i++) std::cout << fJetbinpTedges[i] << " - ";
  std::cout << fJetbinpTedges[fnJetbins] << std::endl;

  //Reload input object to extract projections in jet pT for each D-mass pTrange
  fFileInput = TFile::Open(fFileNameInput.Data(),"read");
  if(!fFileInput){
    std::cout << "File " << fFileInput << " cannot be opened! check your file path!" << std::endl; return kFALSE;
  }

  TTree *tree = (TTree*)fFileInput->Get(fTreeName);
  AliAnalysisTaskDmesonJets::AliD0InfoSummary *brD = 0;
  AliAnalysisTaskDmesonJets::AliJetInfoSummary *brJet = 0;
  tree->SetBranchAddress(fDBranchName,&brD);
  tree->SetBranchAddress(fJetBranchName,&brJet);

  if(!tree || !brD || !brJet) {
    std::cout << "Error in setting the tree/branch names! Exiting..." << std::endl;
    return kFALSE;
  }

  TRandom3 *gen = new TRandom3();
  gen->SetSeed(0);

  //define list of histograms, one per variation
  fJetSpectrSBVars = new TH1F*[fnMaxTrials];

  //for debug and thorough studies
  Double_t arrYldBinPerBin[fnDbins][fnJetbins][fnMaxTrials];
  for (Int_t i = 0; i<fnDbins; i++) for (Int_t j = 0; j<fnJetbins; j++) for (Int_t k = 0; k<fnMaxTrials; k++) arrYldBinPerBin[i][j][k] = 0;

  for(int iDbin=0; iDbin<fnDbins; iDbin++) {

    //open file with summary of variations from MultiTrial - get histos of variations
    TFile *fileMult = TFile::Open(Form("RawYieldSyst_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fDbinpTedges[iDbin],fDbinpTedges[iDbin+1]),"read");
    if(!fileMult) {
      std::cout << "Uncertainty file for bin " << fDbinpTedges[iDbin] <<" - "<< fDbinpTedges[iDbin+1] << " cannot be opened! Did you already evaluate it?" << std::endl;
      return kFALSE;
    }
    TCanvas *c = (TCanvas*)fileMult->Get("call");
    TH1F *hMean = (TH1F*)c->FindObject("hMeanAll");
    TH1F *hSigma = (TH1F*)c->FindObject("hSigmaAll");
    TH1F *hBkg = (TH1F*)c->FindObject("hBkgAll");

    Int_t extracted[fnMaxTrials];
    for(int iTrial=0; iTrial<fnMaxTrials; iTrial++) extracted[iTrial] = -1;

    TH1F* hjetpt = new TH1F(Form("hjetpt%d",iDbin),"hJetPt_signReg",fnJetbins,fJetbinpTedges);
    TH1F* hjetpt_s1 = new TH1F(Form("hjetpt_s1%d",iDbin),"hJetPt_sb1",fnJetbins,fJetbinpTedges);
    TH1F* hjetpt_s2 = new TH1F(Form("hjetpt_s2%d",iDbin),"hJetPt_sb2",fnJetbins,fJetbinpTedges);
    TH1F* hjetpt_s = new TH1F(Form("hjetpt_s%d",iDbin),"hJetPt_sb2",fnJetbins,fJetbinpTedges);
    hjetpt->Sumw2();
    hjetpt_s1->Sumw2();
    hjetpt_s2->Sumw2();

    if(fDebug) std::cout << "Running bin pT(D) " << iDbin << std::endl;

    for(int iTrial=0; iTrial<fnMaxTrials; iTrial++) {

        Double_t mean = hMean->GetBinContent(iTrial+1);
        Double_t sigma = hSigma->GetBinContent(iTrial+1);
        Double_t bkg = hBkg->GetBinContent(iTrial+1);

  	Float_t signal_l_min = mean-8*sigma;
  	Float_t signal_l_max = mean-4*sigma;
  	Float_t signal_u_min = mean+4*sigma;
  	Float_t signal_u_max = mean+8*sigma;
  	Float_t signal_c_min = mean-fnSigmaSignReg*sigma;
  	Float_t signal_c_max = mean+fnSigmaSignReg*sigma;

        for(int k=0; k<tree->GetEntries(); k++) {
          tree->GetEntry(k);
	  if(TMath::Abs(brJet->fEta)>0.5) continue;
          if(brD->fInvMass > signal_c_min && brD->fInvMass < signal_c_max && brD->fPt > fDbinpTedges[iDbin] && brD->fPt <= fDbinpTedges[iDbin+1]) hjetpt->Fill(brJet->fPt);
  	  if(brD->fInvMass > signal_l_min && brD->fInvMass < signal_l_max && brD->fPt > fDbinpTedges[iDbin] && brD->fPt <= fDbinpTedges[iDbin+1]) hjetpt_s1->Fill(brJet->fPt);
	  if(brD->fInvMass > signal_u_min && brD->fInvMass < signal_u_max && brD->fPt > fDbinpTedges[iDbin] && brD->fPt <= fDbinpTedges[iDbin+1]) hjetpt_s2->Fill(brJet->fPt);
        }

        hjetpt_s = (TH1F*)hjetpt_s1->Clone(Form("hjetpt_s%d",iDbin));
        hjetpt_s->Add(hjetpt_s2);

       // scale background from side bands to the background under the peak
	if(hjetpt_s->Integral(0,hjetpt_s->GetNbinsX()+1)==0) {
     	  std::cout << "Error! At least one variation with integral of signal = 0! Exiting..." << std::endl;
	  return kFALSE;
        }

        Double_t scaling = bkg / hjetpt_s->Integral(0,hjetpt_s->GetNbinsX()+1);
        hjetpt_s->Scale(scaling);

        // subtract background from signal jet
        hjetpt->Add(hjetpt_s,-1);
        hjetpt->SetMarkerColor(kGreen+3);
        hjetpt->SetLineColor(kGreen+3);

        // correct for D* efficiency
        hjetpt->Scale(1./fDEffValues[iDbin]); // D efficiency
        hjetpt->SetMarkerColor(kBlue+3);
        hjetpt->SetLineColor(kBlue+3);

	// for every trial of every pT(D) bin, save the value of the yield, after eff scaling, in each pT(jet) bin (to study pT(D)->pT(jet) yield correlations)
	for (int l = 0; l<hjetpt->GetNbinsX(); l++) arrYldBinPerBin[iDbin][l][iTrial] = hjetpt->GetBinContent(l + 1);

	// add 'iDbin' pT(D) bin to total spectrum for variation 'iTrial'
        if(!iDbin) fJetSpectrSBVars[iTrial] = (TH1F*)hjetpt->Clone(Form("JetRawYieldUncert_%d",iTrial));
        else fJetSpectrSBVars[iTrial]->Add(hjetpt);

	hjetpt->Reset();
	hjetpt_s1->Reset();
	hjetpt_s2->Reset();
	hjetpt_s->Reset();

    } //end loop on trials for sideband approach

  } //end loop on pT(D) bins

  //Now evaluate central value + rms in each pT(jet) bin to build the uncertainty
  Int_t nJetBins = fJetSpectrSBVars[0]->GetNbinsX();
  Double_t arrYld[nJetBins][fnMaxTrials];
  for(Int_t i=0; i<nJetBins; i++) for(Int_t j=0; j<fnMaxTrials; j++) arrYld[i][j] = 0;

  fJetPtBinYieldDistribution = new TH1F*[nJetBins];

  fJetYieldUnc = (TH1D*)fJetSpectrSBVars[0]->Clone("JetRawYieldUncert");
  fJetYieldUnc->Reset();
  fJetYieldUnc->SetTitle("Raw yield uncertainty on jet spectrum - Dzero - Sideband subtraction");
  fJetYieldCentral = (TH1D*)fJetSpectrSBVars[0]->Clone("JetRawYieldCentral");
  fJetYieldCentral->Reset();
  fJetYieldCentral->SetTitle("Jet spectrum central values + syst yield uncertainty - Dzero - Sideband subtraction");

  for(Int_t iJetbin=0; iJetbin<nJetBins; iJetbin++) { //loop on jet spectrum pT bins

    fJetPtBinYieldDistribution[iJetbin] = new TH1F(Form("fJetPtBinYieldDistribution_Bin%d",iJetbin),"  ; Yield distribution",50000,0.,50000.);

    for(Int_t iTrial=0; iTrial<fnMaxTrials; iTrial++) { //loop on trials and build array of variations for a given pT(jet) bin
      arrYld[iJetbin][iTrial] = fJetSpectrSBVars[iTrial]->GetBinContent(iJetbin+1);
      fJetPtBinYieldDistribution[iJetbin]->Fill(arrYld[iJetbin][iTrial]);
    }

    Double_t mean = TMath::Mean(fnMaxTrials,arrYld[iJetbin]);
    Double_t rms = TMath::RMS(fnMaxTrials,arrYld[iJetbin]);
    if(fDebug) {
      std::cout << "Jet bin " << iJetbin << " (" << fJetSpectrSBVars[0]->GetXaxis()->GetBinLowEdge(iJetbin+1) << "-" << fJetSpectrSBVars[0]->GetXaxis()->GetBinUpEdge(iJetbin+1) << ")";
      std::cout << ": Mean = " << mean << ", RMS = " << rms << std::endl;
    }

    fJetYieldUnc->SetBinContent(iJetbin+1,rms);
    fJetYieldCentral->SetBinContent(iJetbin+1,mean);
    fJetYieldCentral->SetBinError(iJetbin+1,rms);

    if(fRebinDstarSB) fJetPtBinYieldDistribution[iJetbin]->SaveAs(Form("YieldDistribution_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fJetbinpTedges[iJetbin],fJetbinpTedges[iJetbin+1]));

  }

  //fOut->Close();

  fJetYieldUnc->SetStats(kFALSE);
  fJetYieldUnc->Draw();
  fJetYieldUnc->SaveAs(Form("FinalRawYieldUncertainty_%s.root",fDmesonLabel.Data()));
  fJetYieldCentral->SetStats(kFALSE);
  fJetYieldCentral->Draw();
  fJetYieldCentral->SaveAs(Form("FinalRawYieldCentralPlusSystUncertainty_%s.root",fDmesonLabel.Data()));

  if (fDebug) {
	  //ADVANCED - save distribution of final jet yields (summing all pT(D) bins) in a single plot
	  TCanvas *cDistr = new TCanvas("cDistr", "cDistr", 900, 600);
	  for (Int_t iTrial = 0; iTrial<fnMaxTrials; iTrial++) {
		  for (int l = 0; l<fJetSpectrSBVars[iTrial]->GetNbinsX(); l++) fJetSpectrSBVars[iTrial]->SetBinError(l + 1, 0.0001);
		  fJetSpectrSBVars[iTrial]->SetMarkerColor(iTrial + 1);
		  fJetSpectrSBVars[iTrial]->SetLineColor(iTrial + 1);
		  if (!iTrial) fJetSpectrSBVars[iTrial]->Draw();
		  else fJetSpectrSBVars[iTrial]->Draw("same");
	  }
	  cDistr->SaveAs(Form("DistributionOfFinalYields_SBApproach_%s_AfterDbinSum.root", fDmesonLabel.Data()));

	  //ADVANCED - save distribution of final jet yields from each single pT(D) bin in a single plot (one per each pT(D) bin)
	  for (int iDbin = 0; iDbin<fnDbins; iDbin++) {
		  TCanvas *cDistr1 = new TCanvas(Form("cDistr%d", iDbin), Form("cDistr%d", iDbin), 900, 600);
		  TH1F** hJetSpectrFromSingleDbin = new TH1F*[fnMaxTrials];
		  for (Int_t iTrial = 0; iTrial<fnMaxTrials; iTrial++) {
			  hJetSpectrFromSingleDbin[iTrial] = (TH1F*)fJetSpectrSBVars[0]->Clone(Form("JetRawYieldDistr_Dbin%d", iDbin));
			  for (int l = 0; l<hJetSpectrFromSingleDbin[iTrial]->GetNbinsX(); l++) {
				  hJetSpectrFromSingleDbin[iTrial]->SetBinContent(l + 1, arrYldBinPerBin[iDbin][l][iTrial]);
				  hJetSpectrFromSingleDbin[iTrial]->SetBinError(l + 1, 0.0001);
			  }
			  hJetSpectrFromSingleDbin[iTrial]->SetMarkerColor(iTrial + 1);
			  hJetSpectrFromSingleDbin[iTrial]->SetLineColor(iTrial + 1);
			  if (!iTrial) hJetSpectrFromSingleDbin[iTrial]->Draw();
			  else hJetSpectrFromSingleDbin[iTrial]->Draw("same");
		  }
		  cDistr1->SaveAs(Form("DistributionOfFinalYields_SBApproach_%s_Bin%d.root", fDmesonLabel.Data(), iDbin));
	  }

	  //ADVANCED - save averages of final jet yields from each single pT(D) bin, with their RMS, without summing them, in a single plot
	  TCanvas *cDistr2 = new TCanvas("cDistrAllAvgs", "cDistrAllAvgs", 900, 600);

	  fJetYieldCentral->SetLineWidth(3);
	  fJetYieldCentral->Draw();

	  TH1F** hJetSpectrFromSingleDbin_Avg = new TH1F*[fnDbins];

	  for (Int_t iDbin = 0; iDbin<fnDbins; iDbin++) {
		  hJetSpectrFromSingleDbin_Avg[iDbin] = (TH1F*)fJetSpectrSBVars[0]->Clone(Form("JetRawYieldAvgDistr_%d", iDbin));

		  for (Int_t iJetbin = 0; iJetbin<nJetBins; iJetbin++) { //loop on jet spectrum pT bins

			  Double_t mean = TMath::Mean(fnMaxTrials, arrYldBinPerBin[iDbin][iJetbin]);
			  Double_t rms = TMath::RMS(fnMaxTrials, arrYldBinPerBin[iDbin][iJetbin]);

			  fJetYieldUnc->SetBinContent(iJetbin + 1, rms);
			  hJetSpectrFromSingleDbin_Avg[iDbin]->SetBinContent(iJetbin + 1, mean);
			  hJetSpectrFromSingleDbin_Avg[iDbin]->SetBinError(iJetbin + 1, rms);
		  }

		  hJetSpectrFromSingleDbin_Avg[iDbin]->SetMarkerColor(iDbin + 1);
		  hJetSpectrFromSingleDbin_Avg[iDbin]->SetLineColor(iDbin + 1);
		  hJetSpectrFromSingleDbin_Avg[iDbin]->Draw("same");
	  }

	  TLegend* leg = new TLegend(0.1, 0.7, 0.48, 0.9);
	  leg->AddEntry(fJetYieldCentral, "Average after pT(D) bin sum", "pl");
	  for (Int_t iDbin = 0; iDbin<fnDbins; iDbin++) leg->AddEntry(hJetSpectrFromSingleDbin_Avg[iDbin], Form("pt(D) %1.1f - %1.1f", fDbinpTedges[iDbin], fDbinpTedges[iDbin + 1]), "pl");
	  leg->Draw();

	  cDistr2->SaveAs(Form("AverageOfFinalYields_SBApproach_%s_AllDBins.root", fDmesonLabel.Data()));

  } //end of advanced plots

  return kTRUE;

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::EvaluateUncertaintyDstarEffScale() {

  std::cout << "Jet spectrum pT bin edges: ";
  for(int i=0; i<fnJetbins; i++) std::cout << fJetbinpTedges[i] << " - ";
  std::cout << fJetbinpTedges[fnJetbins] << std::endl;

  fJetYieldUnc = new TH1D("JetRawYieldUncert","Raw yield uncertainty on jet spectrum - Dstar - Eff.scaling",fnJetbins,fJetbinpTedges);
  fJetYieldCentral = new TH1D("JetRawYieldCentral","Jet spectrum central values + syst yield uncertainty - Dstar - Eff.scaling",fnJetbins,fJetbinpTedges);

  TH1D *hUnc;

  //loop over the jet pT bins already extracted
  for(int ibin=0; ibin<fnJetbins; ibin++) {
    if(fDebug) std::cout << Form("Hist_RawYieldSyst_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fJetbinpTedges[ibin],fJetbinpTedges[ibin+1]) << std::endl;
    TFile *f = TFile::Open(Form("Hist_RawYieldSyst_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fJetbinpTedges[ibin],fJetbinpTedges[ibin+1]),"read");
    if(!f){
      std::cout << "Uncertainty file for bin " << fJetbinpTedges[ibin] <<" - "<< fJetbinpTedges[ibin+1] << " cannot be opened! Did you already evaluate it?" << std::endl;
      return kFALSE;
    }
    else {
      if(fDebug) std::cout << "Building uncertainty for bin " << fJetbinpTedges[ibin] <<" - "<< fJetbinpTedges[ibin+1] << std::endl;
    }

    hUnc = (TH1D*)f->Get(Form("hUnc_%1.1fto%1.1f",fJetbinpTedges[ibin],fJetbinpTedges[ibin+1]));
    if(!hUnc){
      std::cout << "Histogram with uncertainty from mass plot not found! Returning..." << std::endl;
      return kFALSE;
    }

    Double_t centrYield = hUnc->GetBinContent(1);
    Double_t rmsPct = hUnc->GetBinError(1);

    fJetYieldUnc->SetBinContent(ibin+1,rmsPct);
    fJetYieldCentral->SetBinContent(ibin+1,centrYield);
    fJetYieldCentral->SetBinError(ibin+1,rmsPct);
  }

  fJetYieldUnc->SetStats(kFALSE);
  fJetYieldUnc->Draw();
  fJetYieldUnc->SaveAs(Form("FinalRawYieldUncertainty_%s.root",fDmesonLabel.Data()));
  fJetYieldCentral->SetStats(kFALSE);
  fJetYieldCentral->Draw();
  fJetYieldCentral->SaveAs(Form("FinalRawYieldCentralPlusSystUncertainty_%s.root",fDmesonLabel.Data()));

  //print distribution of yields for each variation
  for(int ibin=0; ibin<fnJetbins; ibin++) {
    TFile *f2 = TFile::Open(Form("RawYieldSyst_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fJetbinpTedges[ibin],fJetbinpTedges[ibin+1]),"read");
    TCanvas *c = (TCanvas*)f2->Get("call");
    TH1F *hDist = (TH1F*)c->FindObject("hRawYieldDistAll");
    hDist->SetStats(kTRUE);
    hDist->SaveAs(Form("YieldDistribution_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fJetbinpTedges[ibin],fJetbinpTedges[ibin+1]));
  }

  return kTRUE;

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::EvaluateUncertaintyDstarSideband() {

  double jetmin = 0, jetmax = 50;
  fzmin = -2;
  fzmax = 2;

  //Reload input object to extract projections in jet pT for each D-mass pTrange
  fFileInput = TFile::Open(fFileNameInput.Data(),"read");
  if(!fFileInput){
    std::cout << "File " << fFileInput << " cannot be opened! check your file path!" << std::endl; return kFALSE;
  }

  TDirectoryFile* dir=(TDirectoryFile*)fFileInput->Get(fDirName.Data());

    TList *histList =  (TList*)dir->Get(Form("%s0",fListName.Data()));
    THnSparseF *sparse = (THnSparseF*)histList->FindObject(fObjectName.Data());
    sparse->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD=(TH3D*)sparse->Projection(3,1,2);
    hInvMassptD->SetName("hInvMassptD");

   TCanvas *c = new TCanvas;
  hInvMassptD->Draw("colz");


  std::cout << "file: " << fFileNameInput.Data() << std::endl;
  std::cout << "dir: " << fDirName.Data() << std::endl;
  std::cout << "list: " << fListName.Data() << std::endl;
  std::cout << "object: " << fObjectName.Data() << std::endl;


    TList *histList_1 =  (TList*)dir->Get(Form("%s1",fListName.Data()));
    THnSparseF *sparse_1 = (THnSparseF*)histList_1->FindObject("hsDphiz");
    sparse_1->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD_1=(TH3D*)sparse_1->Projection(3,1,2);
    hInvMassptD_1->SetName("hInvMassptD_1");
    hInvMassptD->Add(hInvMassptD_1);

    TList *histList_2 =  (TList*)dir->Get(Form("%s2",fListName.Data()));
    THnSparseF *sparse_2 = (THnSparseF*)histList_2->FindObject("hsDphiz");
    sparse_2->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD_2=(TH3D*)sparse_2->Projection(3,1,2);
    hInvMassptD_2->SetName("hInvMassptD_2");
    hInvMassptD->Add(hInvMassptD_2);

    TList *histList_3 =  (TList*)dir->Get(Form("%s3",fListName.Data()));
    THnSparseF *sparse_3 = (THnSparseF*)histList_3->FindObject("hsDphiz");
    sparse_3->GetAxis(0)->SetRangeUser(fzmin,fzmax);
    TH3D* hInvMassptD_3=(TH3D*)sparse_3->Projection(3,1,2);
    hInvMassptD_3->SetName("hInvMassptD_3");
    hInvMassptD->Add(hInvMassptD_3);

 /* TList *histList_4 =  (TList*)dir->Get(Form("%s4",fListName.Data()));
  THnSparseF *sparse_4 = (THnSparseF*)histList_4->FindObject("hsDphiz");
  sparse_4->GetAxis(0)->SetRangeUser(fzmin,fzmax);
  TH3D* hInvMassptD_4=(TH3D*)sparse_4->Projection(3,1,2);
  hInvMassptD_4->SetName("hInvMassptD_4");
  hInvMassptD->Add(hInvMassptD_4);
*/
  //if(!hInvMassptD || !hInvMassptD_1 || !hInvMassptD_2 || !hInvMassptD_3 || !hInvMassptD_4) return kFALSE;
  if(!hInvMassptD || !hInvMassptD_1 || !hInvMassptD_2 || !hInvMassptD_3) return kFALSE;


  TRandom3 *gen = new TRandom3();
  gen->SetSeed(0);

  //define list of histograms, one per variation
  fJetSpectrSBVars = new TH1F*[fnMaxTrials];

  //for debug and thorough studies
  Double_t arrYldBinPerBin[fnDbins][fnJetbins][fnMaxTrials];
  for(Int_t i=0; i<fnDbins; i++) for(Int_t j=0; j<fnJetbins; j++) for(Int_t k=0; k<fnMaxTrials; k++) arrYldBinPerBin[i][j][k] = 0;

  for(int iDbin=0; iDbin<fnDbins; iDbin++) {

    //open file with summary of variations from MultiTrial - get histos of variations
    TFile *fileMult = TFile::Open(Form("RawYieldSyst_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fDbinpTedges[iDbin],fDbinpTedges[iDbin+1]),"read");
    if(!fileMult) {
      std::cout << "Uncertainty file for bin " << fDbinpTedges[iDbin] <<" - "<< fDbinpTedges[iDbin+1] << " cannot be opened! Did you already evaluate it?" << std::endl;
      return kFALSE;
    }
    TCanvas *c = (TCanvas*)fileMult->Get("call");
    TH1F *hMean = (TH1F*)c->FindObject("hMeanAll");
    TH1F *hSigma = (TH1F*)c->FindObject("hSigmaAll");
    TH1F *hBkg = (TH1F*)c->FindObject("hBkgAll");

    Int_t extracted[fnMaxTrials];
    for(int iTrial=0; iTrial<fnMaxTrials; iTrial++) extracted[iTrial] = -1;

    TH1F* hjetpt;
    TH1F* hjetpt_s1;
    TH1F* hjetpt_s2;
    TH1F* hjetpt_s;

    if(fDebug) std::cout << "Running bin pT(D) " << iDbin << std::endl;

    for(int iTrial=0; iTrial<fnMaxTrials; iTrial++) {

        Bool_t extractOk = kFALSE;
        Int_t rnd = -1;

	       do {  //just one time if fAllowRepetitions==kTRUE, repeat extraction till new number is obtained if fAllowRepetitions==kFALSE

          rnd = gen->Integer(hMean->GetNbinsX())+1;
	        if(hSigma->GetBinContent(rnd)>0) extractOk = kTRUE;  //avoid 'empty' cases

          if(!fAllowRepetitions) { //check if already extracted for this pT(D) bin

            if(fnMaxTrials>hMean->GetNbinsX()) {
       	      std::cout << "Error! you set more set spectrum total variations than those done for pT(D) bin" << fDbinpTedges[iDbin] <<" - "<< fDbinpTedges[iDbin+1] << "! ";
       	      std::cout << "Impossible to do without allowing repetitions! Exiting..." << std::endl;
  	          return kFALSE;
            }

            for(int j=0; j<iTrial; j++) {
	             if(rnd==extracted[j]) extractOk=kFALSE;
	           }

          } //end of if(!fAllowRepetitions)

          extracted[iTrial] = rnd;

        } while (extractOk==kFALSE);

        Double_t mean = hMean->GetBinContent(rnd);
        Double_t sigma = hSigma->GetBinContent(rnd);
        Double_t bkg = hBkg->GetBinContent(rnd);

      	Float_t signal_l_min = mean-8*sigma;
      	Float_t signal_l_max = mean-5*sigma;
      	Float_t signal_u_min = mean+5*sigma;
      	Float_t signal_u_max = mean+13*sigma;
      	Float_t signal_c_min = mean-fnSigmaSignReg*sigma;
      	Float_t signal_c_max = mean+fnSigmaSignReg*sigma;

        Double_t scaling;

        if(!fRebinDstarSB) { //do not rebin, just use projections

          //extract signal and sideband region spectra
          hjetpt=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_c_min), hInvMassptD->GetXaxis()->FindBin(signal_c_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
          hjetpt_s1=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_s1%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_l_min), hInvMassptD->GetXaxis()->FindBin(signal_l_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
          hjetpt_s2=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_s2%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_u_min), hInvMassptD->GetXaxis()->FindBin(signal_u_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
          hjetpt_s = (TH1F*)hjetpt_s1->Clone(Form("hjetpt_s%d",iDbin));
          hjetpt_s->Add(hjetpt_s2);

          // scale background from side bands to the background under the peak
	         if(hjetpt_s->Integral()==0) {
     	        std::cout << "Error! At least one variation with no entries! Exiting..." << std::endl;
	             return kFALSE;
          }
	         scaling = bkg / hjetpt_s->Integral(hjetpt_s->FindBin(jetmin),hjetpt_s->FindBin(jetmax));

        }
        else { //rebin spectrum

 	        hjetpt = new TH1F(Form("hjetpt%d",iDbin),"hJetPt_signReg_Rebinned",fnJetbins,fJetbinpTedges);
   	      hjetpt_s1 = new TH1F(Form("hjetpt_s1%d",iDbin),"hJetPt_sb1_Rebinned",fnJetbins,fJetbinpTedges);
      	  hjetpt_s2 = new TH1F(Form("hjetpt_s2%d",iDbin),"hJetPt_sb2_Rebinned",fnJetbins,fJetbinpTedges);
    	    hjetpt_s = new TH1F(Form("hjetpt_s%d",iDbin),"hJetPt_sb2_Rebinned",fnJetbins,fJetbinpTedges);
    	    hjetpt->Sumw2();
	        hjetpt_s1->Sumw2();
	        hjetpt_s2->Sumw2();

          //extract signal and sideband region spectra
          TH1F* tmphjetpt=(TH1F*)hInvMassptD->ProjectionY(Form("tmphjetpt%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_c_min), hInvMassptD->GetXaxis()->FindBin(signal_c_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
          TH1F* tmphjetpt_s1=(TH1F*)hInvMassptD->ProjectionY(Form("tmphjetpt_s1%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_l_min), hInvMassptD->GetXaxis()->FindBin(signal_l_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
          TH1F* tmphjetpt_s2=(TH1F*)hInvMassptD->ProjectionY(Form("tmphjetpt_s2%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_u_min), hInvMassptD->GetXaxis()->FindBin(signal_u_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
          TH1F* tmphjetpt_s = (TH1F*)tmphjetpt_s1->Clone(Form("tmphjetpt_s%d",iDbin));
          tmphjetpt_s->Add(tmphjetpt_s2);

		//tmphjetpt->Draw();

          // scale background from side bands to the background under the peak
	        if(tmphjetpt_s->Integral()==0) {
     	      std::cout << "Error! At least one variation with no entries! Exiting reb..." << std::endl;
	          return kFALSE;
          }
	        scaling = bkg / tmphjetpt_s->Integral(tmphjetpt_s->FindBin(jetmin),tmphjetpt_s->FindBin(jetmax)); //integral btw 0 and 50 pTJet (where you get the bkg from the mass plot)

      	  for(int j=1; j<=tmphjetpt->GetNbinsX(); j++) {
      	    Double_t centerbin = tmphjetpt->GetBinCenter(j); //bin of hjetpt corresponding to j-th bin of THnSparse projection
      	    hjetpt->Fill(centerbin,tmphjetpt->GetBinContent(j));
      	    hjetpt_s1->Fill(centerbin,tmphjetpt_s1->GetBinContent(j));
      	    hjetpt_s2->Fill(centerbin,tmphjetpt_s2->GetBinContent(j));
      	    hjetpt_s->Fill(centerbin,tmphjetpt_s->GetBinContent(j));
          }
          for(int j=1; j<=hjetpt->GetNbinsX(); j++) {
      	    hjetpt->SetBinError(j,TMath::Sqrt(hjetpt->GetBinContent(j)));
      	    hjetpt_s1->SetBinError(j,TMath::Sqrt(hjetpt_s1->GetBinContent(j)));
      	    hjetpt_s2->SetBinError(j,TMath::Sqrt(hjetpt_s2->GetBinContent(j)));
      	    hjetpt_s->SetBinError(j,TMath::Sqrt(hjetpt_s->GetBinContent(j)));
    	    }

        } //end of rebinning of spectrum

        hjetpt_s->Scale(scaling);

        // subtract background from signal jet
        hjetpt->Add(hjetpt_s,-1);
        hjetpt->SetMarkerColor(kGreen+3);
        hjetpt->SetLineColor(kGreen+3);

        // correct for D* efficiency
        hjetpt->Scale(1./fDEffValues[iDbin]); // D efficiency
        hjetpt->SetMarkerColor(kBlue+3);
        hjetpt->SetLineColor(kBlue+3);

	// for every trial of every pT(D) bin, save the value of the yield, after eff scaling, in each pT(jet) bin (to study pT(D)->pT(jet) yield correlations)
        for(int l=0;l<hjetpt->GetNbinsX();l++) arrYldBinPerBin[iDbin][l][iTrial] = hjetpt->GetBinContent(l+1);

	// add 'iDbin' pT(D) bin to total spectrum for variation 'iTrial'
        if(!iDbin) fJetSpectrSBVars[iTrial] = (TH1F*)hjetpt->Clone(Form("JetRawYieldUncert_%d",iTrial));
        else fJetSpectrSBVars[iTrial]->Add(hjetpt);

	hjetpt->Reset();
	hjetpt_s1->Reset();
	hjetpt_s2->Reset();
	hjetpt_s->Reset();

    } //end loop on trials for sideband approach

  } //end loop on pT(D) bins

  //Now evaluate central value + rms in each pT(jet) bin to build the uncertainty
  Int_t nJetBins = fJetSpectrSBVars[0]->GetNbinsX();
  Double_t arrYld[nJetBins][fnMaxTrials];
  for(Int_t i=0; i<nJetBins; i++) for(Int_t j=0; j<fnMaxTrials; j++) arrYld[i][j] = 0;

  fJetPtBinYieldDistribution = new TH1F*[nJetBins];

  fJetYieldUnc = (TH1D*)fJetSpectrSBVars[0]->Clone("JetRawYieldUncert");
  fJetYieldUnc->Reset();
  fJetYieldUnc->SetTitle("Raw yield uncertainty on jet spectrum - Dstar - Sideband subtraction");
  fJetYieldCentral = (TH1D*)fJetSpectrSBVars[0]->Clone("JetRawYieldCentral");
  fJetYieldCentral->Reset();
  fJetYieldCentral->SetTitle("Jet spectrum central values + syst yield uncertainty - Dstar - Sideband subtraction");

  for(Int_t iJetbin=0; iJetbin<nJetBins; iJetbin++) { //loop on jet spectrum pT bins

    fJetPtBinYieldDistribution[iJetbin] = new TH1F(Form("fJetPtBinYieldDistribution_Bin%d",iJetbin),"  ; Yield distribution",50000,0.,50000.);

    for(Int_t iTrial=0; iTrial<fnMaxTrials; iTrial++) { //loop on trials and build array of variations for a given pT(jet) bin
      arrYld[iJetbin][iTrial] = fJetSpectrSBVars[iTrial]->GetBinContent(iJetbin+1);
      fJetPtBinYieldDistribution[iJetbin]->Fill(arrYld[iJetbin][iTrial]);
    }

    Double_t mean = TMath::Mean(fnMaxTrials,arrYld[iJetbin]);
    Double_t rms = TMath::RMS(fnMaxTrials,arrYld[iJetbin]);
    if(fDebug) {std::cout << "Jet bin " << iJetbin << " (" << fJetSpectrSBVars[0]->GetXaxis()->GetBinLowEdge(iJetbin+1) << "-" << fJetSpectrSBVars[0]->GetXaxis()->GetBinUpEdge(iJetbin+1) << ")";
      std::cout << ": Mean = " << mean << ", RMS = " << rms << std::endl;
    }

    fJetYieldUnc->SetBinContent(iJetbin+1,rms);
    fJetYieldCentral->SetBinContent(iJetbin+1,mean);
    fJetYieldCentral->SetBinError(iJetbin+1,rms);

    fJetPtBinYieldDistribution[iJetbin]->SaveAs(Form("YieldDistribution_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fJetbinpTedges[iJetbin],fJetbinpTedges[iJetbin+1]));

  }

  //fOut->Close();

  fJetYieldUnc->SetStats(kFALSE);
  fJetYieldUnc->Draw();
  fJetYieldUnc->SaveAs(Form("FinalRawYieldUncertainty_%s.root",fDmesonLabel.Data()));
  fJetYieldCentral->SetStats(kFALSE);
  fJetYieldCentral->Draw();
  fJetYieldCentral->SaveAs(Form("FinalRawYieldCentralPlusSystUncertainty_%s.root",fDmesonLabel.Data()));

  if(fDebug) {
    //ADVANCED - save distribution of final jet yields (summing all pT(D) bins) in a single plot
    TCanvas *cDistr = new TCanvas("cDistr","cDistr",900,600);
    for(Int_t iTrial=0; iTrial<fnMaxTrials; iTrial++) {
      for(int l=0; l<fJetSpectrSBVars[iTrial]->GetNbinsX();l++) fJetSpectrSBVars[iTrial]->SetBinError(l+1,0.0001);
      fJetSpectrSBVars[iTrial]->SetMarkerColor(iTrial+1);
      fJetSpectrSBVars[iTrial]->SetLineColor(iTrial+1);
      if(!iTrial) fJetSpectrSBVars[iTrial]->Draw();
      else fJetSpectrSBVars[iTrial]->Draw("same");
    }
    cDistr->SaveAs(Form("DistributionOfFinalYields_SBApproach_%s_AfterDbinSum.root",fDmesonLabel.Data()));

    //ADVANCED - save distribution of final jet yields from each single pT(D) bin in a single plot (one per each pT(D) bin)
    for(int iDbin=0; iDbin<fnDbins; iDbin++) {
      TCanvas *cDistr1 = new TCanvas(Form("cDistr%d",iDbin),Form("cDistr%d",iDbin),900,600);
      TH1F** hJetSpectrFromSingleDbin = new TH1F*[fnMaxTrials];
      for(Int_t iTrial=0; iTrial<fnMaxTrials; iTrial++) {
        hJetSpectrFromSingleDbin[iTrial] = (TH1F*)fJetSpectrSBVars[0]->Clone(Form("JetRawYieldDistr_Dbin%d",iDbin));
        for(int l=0; l<hJetSpectrFromSingleDbin[iTrial]->GetNbinsX();l++) {
          hJetSpectrFromSingleDbin[iTrial]->SetBinContent(l+1,arrYldBinPerBin[iDbin][l][iTrial]);
          hJetSpectrFromSingleDbin[iTrial]->SetBinError(l+1,0.0001);
        }
        hJetSpectrFromSingleDbin[iTrial]->SetMarkerColor(iTrial+1);
        hJetSpectrFromSingleDbin[iTrial]->SetLineColor(iTrial+1);
        if(!iTrial) hJetSpectrFromSingleDbin[iTrial]->Draw();
        else hJetSpectrFromSingleDbin[iTrial]->Draw("same");
      }
      cDistr1->SaveAs(Form("DistributionOfFinalYields_SBApproach_%s_Bin%d.root",fDmesonLabel.Data(),iDbin));
    }

    //ADVANCED - save averages of final jet yields from each single pT(D) bin, with their RMS, without summing them, in a single plot
    TCanvas *cDistr2 = new TCanvas("cDistrAllAvgs","cDistrAllAvgs",900,600);

    fJetYieldCentral->SetLineWidth(3);
    fJetYieldCentral->Draw();

    TH1F** hJetSpectrFromSingleDbin_Avg = new TH1F*[fnDbins];

    for(Int_t iDbin=0; iDbin<fnDbins; iDbin++) {
      hJetSpectrFromSingleDbin_Avg[iDbin] = (TH1F*)fJetSpectrSBVars[0]->Clone(Form("JetRawYieldAvgDistr_%d",iDbin));

        for(Int_t iJetbin=0; iJetbin<nJetBins; iJetbin++) { //loop on jet spectrum pT bins

          Double_t mean = TMath::Mean(fnMaxTrials,arrYldBinPerBin[iDbin][iJetbin]);
          Double_t rms = TMath::RMS(fnMaxTrials,arrYldBinPerBin[iDbin][iJetbin]);

          fJetYieldUnc->SetBinContent(iJetbin+1,rms);
          hJetSpectrFromSingleDbin_Avg[iDbin]->SetBinContent(iJetbin+1,mean);
          hJetSpectrFromSingleDbin_Avg[iDbin]->SetBinError(iJetbin+1,rms);
        }

      hJetSpectrFromSingleDbin_Avg[iDbin]->SetMarkerColor(iDbin+1);
      hJetSpectrFromSingleDbin_Avg[iDbin]->SetLineColor(iDbin+1);
      hJetSpectrFromSingleDbin_Avg[iDbin]->Draw("same");
    }

    TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->AddEntry(fJetYieldCentral,"Average after pT(D) bin sum","pl");
    for(Int_t iDbin=0; iDbin<fnDbins; iDbin++) leg->AddEntry(hJetSpectrFromSingleDbin_Avg[iDbin],Form("pt(D) %1.1f - %1.1f",fDbinpTedges[iDbin],fDbinpTedges[iDbin+1]),"pl");
    leg->Draw();

    cDistr2->SaveAs(Form("AverageOfFinalYields_SBApproach_%s_AllDBins.root",fDmesonLabel.Data()));

  } //end of advanced plots

  return kTRUE;

}

//___________________________________________________________________________________________
Bool_t AliDJetRawYieldUncertaintyLocal::EvaluateUncertaintyDstarSideband_CoherentTrialChoice() {

  double jetmin = 0, jetmax = 50;

  //Reload input object to extract projections in jet pT for each D-mass pTrange
  fFileInput = TFile::Open(fFileNameInput.Data(),"read");
  if(!fFileInput){
    std::cout << "File " << fFileInput << " cannot be opened! check your file path!" << std::endl; return kFALSE;
  }

  TDirectoryFile* dir=(TDirectoryFile*)fFileInput->Get(fDirName.Data());

  TList *histList =  (TList*)dir->Get(Form("%s0",fListName.Data()));
  THnSparseF *sparse = (THnSparseF*)histList->FindObject(fObjectName.Data());
  sparse->GetAxis(0)->SetRangeUser(fzmin,fzmax);
  TH3D* hInvMassptD=(TH3D*)sparse->Projection(3,1,2);
  hInvMassptD->SetName("hInvMassptD");

  TList *histList_1 =  (TList*)dir->Get(Form("%s1",fListName.Data()));
  THnSparseF *sparse_1 = (THnSparseF*)histList_1->FindObject("hsDphiz");
  sparse_1->GetAxis(0)->SetRangeUser(fzmin,fzmax);
  TH3D* hInvMassptD_1=(TH3D*)sparse_1->Projection(3,1,2);
  hInvMassptD_1->SetName("hInvMassptD_1");
  hInvMassptD->Add(hInvMassptD_1);

  TList *histList_2 =  (TList*)dir->Get(Form("%s2",fListName.Data()));
  THnSparseF *sparse_2 = (THnSparseF*)histList_2->FindObject("hsDphiz");
  sparse_2->GetAxis(0)->SetRangeUser(fzmin,fzmax);
  TH3D* hInvMassptD_2=(TH3D*)sparse_2->Projection(3,1,2);
  hInvMassptD_2->SetName("hInvMassptD_2");
  hInvMassptD->Add(hInvMassptD_2);

  TList *histList_3 =  (TList*)dir->Get(Form("%s3",fListName.Data()));
  THnSparseF *sparse_3 = (THnSparseF*)histList_3->FindObject("hsDphiz");
  sparse_3->GetAxis(0)->SetRangeUser(fzmin,fzmax);
  TH3D* hInvMassptD_3=(TH3D*)sparse_3->Projection(3,1,2);
  hInvMassptD_3->SetName("hInvMassptD_3");
  hInvMassptD->Add(hInvMassptD_3);

 /* TList *histList_4 =  (TList*)dir->Get(Form("%s4",fListName.Data()));
  THnSparseF *sparse_4 = (THnSparseF*)histList_4->FindObject("hsDphiz");
  sparse_4->GetAxis(0)->SetRangeUser(fzmin,fzmax);
  TH3D* hInvMassptD_4=(TH3D*)sparse_4->Projection(3,1,2);
  hInvMassptD_4->SetName("hInvMassptD_4");
  hInvMassptD->Add(hInvMassptD_4);
*/
  //if(!hInvMassptD || !hInvMassptD_1 || !hInvMassptD_2 || !hInvMassptD_3 || !hInvMassptD_4) return kFALSE;
  if(!hInvMassptD || !hInvMassptD_1 || !hInvMassptD_2 || !hInvMassptD_3 ) return kFALSE;

  TRandom3 *gen = new TRandom3();
  gen->SetSeed(0);

  //define list of histograms, one per variation
  fJetSpectrSBVars = new TH1F*[fnMaxTrials];

  //for debug and thorough studies
  Double_t arrYldBinPerBin[fnDbins][fnJetbins][fnMaxTrials];
  for(Int_t i=0; i<fnDbins; i++) for(Int_t j=0; j<fnJetbins; j++) for(Int_t k=0; k<fnMaxTrials; k++) arrYldBinPerBin[i][j][k] = 0;

  for(int iDbin=0; iDbin<fnDbins; iDbin++) {

    //open file with summary of variations from MultiTrial - get histos of variations
    TFile *fileMult = TFile::Open(Form("RawYieldSyst_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fDbinpTedges[iDbin],fDbinpTedges[iDbin+1]),"read");
    if(!fileMult) {
      std::cout << "Uncertainty file for bin " << fDbinpTedges[iDbin] <<" - "<< fDbinpTedges[iDbin+1] << " cannot be opened! Did you already evaluate it?" << std::endl;
      return kFALSE;
    }
    TCanvas *c = (TCanvas*)fileMult->Get("call");
    TH1F *hMean = (TH1F*)c->FindObject("hMeanAll");
    TH1F *hSigma = (TH1F*)c->FindObject("hSigmaAll");
    TH1F *hBkg = (TH1F*)c->FindObject("hBkgAll");

    Int_t extracted[fnMaxTrials];
    for(int iTrial=0; iTrial<fnMaxTrials; iTrial++) extracted[iTrial] = -1;

    TH1F* hjetpt;
    TH1F* hjetpt_s1;
    TH1F* hjetpt_s2;
    TH1F* hjetpt_s;

    if(fDebug) std::cout << "Running bin pT(D) " << iDbin << std::endl;

    for(int iTrial=0; iTrial<fnMaxTrials; iTrial++) {

        Double_t mean = hMean->GetBinContent(iTrial+1);
        Double_t sigma = hSigma->GetBinContent(iTrial+1);
        Double_t bkg = hBkg->GetBinContent(iTrial+1);

  	Float_t signal_l_min = mean-8*sigma;
  	Float_t signal_l_max = mean-5*sigma;
  	Float_t signal_u_min = mean+5*sigma;
  	Float_t signal_u_max = mean+13*sigma;
  	Float_t signal_c_min = mean-fnSigmaSignReg*sigma;
  	Float_t signal_c_max = mean+fnSigmaSignReg*sigma;

        Double_t scaling;

        if(!fRebinDstarSB) { //do not rebin, just use projections

          //extract signal and sideband region spectra
          hjetpt=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_c_min), hInvMassptD->GetXaxis()->FindBin(signal_c_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
          hjetpt_s1=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_s1%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_l_min), hInvMassptD->GetXaxis()->FindBin(signal_l_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
          hjetpt_s2=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_s2%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_u_min), hInvMassptD->GetXaxis()->FindBin(signal_u_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
          hjetpt_s = (TH1F*)hjetpt_s1->Clone(Form("hjetpt_s%d",iDbin));
          hjetpt_s->Add(hjetpt_s2);

          // scale background from side bands to the background under the peak
	  if(hjetpt_s->Integral()==0) {
     	    std::cout << "Error! At least one variation with no entries! Exiting..." << std::endl;
	    return kFALSE;
          }
	  scaling = bkg / hjetpt_s->Integral(hjetpt_s->FindBin(jetmin),hjetpt_s->FindBin(jetmax));

        } else { //rebin spectrum

 	  hjetpt = new TH1F(Form("hjetpt%d",iDbin),"hJetPt_signReg_Rebinned",fnJetbins,fJetbinpTedges);
   	  hjetpt_s1 = new TH1F(Form("hjetpt_s1%d",iDbin),"hJetPt_sb1_Rebinned",fnJetbins,fJetbinpTedges);
      	  hjetpt_s2 = new TH1F(Form("hjetpt_s2%d",iDbin),"hJetPt_sb2_Rebinned",fnJetbins,fJetbinpTedges);
    	  hjetpt_s = new TH1F(Form("hjetpt_s%d",iDbin),"hJetPt_sb2_Rebinned",fnJetbins,fJetbinpTedges);
    	  hjetpt->Sumw2();
	  hjetpt_s1->Sumw2();
	  hjetpt_s2->Sumw2();

          //extract signal and sideband region spectra
          TH1F* tmphjetpt=(TH1F*)hInvMassptD->ProjectionY(Form("tmphjetpt%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_c_min), hInvMassptD->GetXaxis()->FindBin(signal_c_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
          TH1F* tmphjetpt_s1=(TH1F*)hInvMassptD->ProjectionY(Form("tmphjetpt_s1%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_l_min), hInvMassptD->GetXaxis()->FindBin(signal_l_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
          TH1F* tmphjetpt_s2=(TH1F*)hInvMassptD->ProjectionY(Form("tmphjetpt_s2%d",iDbin),hInvMassptD->GetXaxis()->FindBin(signal_u_min), hInvMassptD->GetXaxis()->FindBin(signal_u_max),hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[iDbin+1])-1);
          TH1F* tmphjetpt_s = (TH1F*)tmphjetpt_s1->Clone(Form("tmphjetpt_s%d",iDbin));
          tmphjetpt_s->Add(tmphjetpt_s2);

          // scale background from side bands to the background under the peak
	  if(tmphjetpt_s->Integral()==0) {
     	    std::cout << "Error! At least one variation with no entries! Exiting..." << std::endl;
	    return kFALSE;
          }
	  scaling = bkg / tmphjetpt_s->Integral(tmphjetpt_s->FindBin(jetmin),tmphjetpt_s->FindBin(jetmax)); //integral btw 0 and 50 pTJet (where you get the bkg from the mass plot)

	  for(int j=1; j<=tmphjetpt->GetNbinsX(); j++) {
	    Double_t centerbin = tmphjetpt->GetBinCenter(j); //bin of hjetpt corresponding to j-th bin of THnSparse projection
	    hjetpt->Fill(centerbin,tmphjetpt->GetBinContent(j));
	    hjetpt_s1->Fill(centerbin,tmphjetpt_s1->GetBinContent(j));
	    hjetpt_s2->Fill(centerbin,tmphjetpt_s2->GetBinContent(j));
	    hjetpt_s->Fill(centerbin,tmphjetpt_s->GetBinContent(j));
          }
          for(int j=1; j<=hjetpt->GetNbinsX(); j++) {
	    hjetpt->SetBinError(j,TMath::Sqrt(hjetpt->GetBinContent(j)));
	    hjetpt_s1->SetBinError(j,TMath::Sqrt(hjetpt_s1->GetBinContent(j)));
	    hjetpt_s2->SetBinError(j,TMath::Sqrt(hjetpt_s2->GetBinContent(j)));
	    hjetpt_s->SetBinError(j,TMath::Sqrt(hjetpt_s->GetBinContent(j)));
	  }

        } //end of rebinning of spectrum

        hjetpt_s->Scale(scaling);

        // subtract background from signal jet
        hjetpt->Add(hjetpt_s,-1);
        hjetpt->SetMarkerColor(kGreen+3);
        hjetpt->SetLineColor(kGreen+3);

        // correct for D* efficiency
        hjetpt->Scale(1./fDEffValues[iDbin]); // D efficiency
        hjetpt->SetMarkerColor(kBlue+3);
        hjetpt->SetLineColor(kBlue+3);

	// for every trial of every pT(D) bin, save the value of the yield, after eff scaling, in each pT(jet) bin (to study pT(D)->pT(jet) yield correlations)
        for(int l=0;l<hjetpt->GetNbinsX();l++) arrYldBinPerBin[iDbin][l][iTrial] = hjetpt->GetBinContent(l+1);

	// add 'iDbin' pT(D) bin to total spectrum for variation 'iTrial'
        if(!iDbin) fJetSpectrSBVars[iTrial] = (TH1F*)hjetpt->Clone(Form("JetRawYieldUncert_%d",iTrial));
        else fJetSpectrSBVars[iTrial]->Add(hjetpt);

	hjetpt->Reset();
	hjetpt_s1->Reset();
	hjetpt_s2->Reset();
	hjetpt_s->Reset();

    } //end loop on trials for sideband approach

  } //end loop on pT(D) bins

  //Now evaluate central value + rms in each pT(jet) bin to build the uncertainty
  Int_t nJetBins = fJetSpectrSBVars[0]->GetNbinsX();
  Double_t arrYld[nJetBins][fnMaxTrials];
  for(Int_t i=0; i<nJetBins; i++) for(Int_t j=0; j<fnMaxTrials; j++) arrYld[i][j] = 0;

  fJetPtBinYieldDistribution = new TH1F*[nJetBins];

  fJetYieldUnc = (TH1D*)fJetSpectrSBVars[0]->Clone("JetRawYieldUncert");
  fJetYieldUnc->Reset();
  fJetYieldUnc->SetTitle("Raw yield uncertainty on jet spectrum - Dstar - Sideband subtraction");
  fJetYieldCentral = (TH1D*)fJetSpectrSBVars[0]->Clone("JetRawYieldCentral");
  fJetYieldCentral->Reset();
  fJetYieldCentral->SetTitle("Jet spectrum central values + syst yield uncertainty - Dstar - Sideband subtraction");

  for(Int_t iJetbin=0; iJetbin<nJetBins; iJetbin++) { //loop on jet spectrum pT bins

    fJetPtBinYieldDistribution[iJetbin] = new TH1F(Form("fJetPtBinYieldDistribution_Bin%d",iJetbin),"  ; Yield distribution",50000,0.,50000.);

    for(Int_t iTrial=0; iTrial<fnMaxTrials; iTrial++) { //loop on trials and build array of variations for a given pT(jet) bin
      arrYld[iJetbin][iTrial] = fJetSpectrSBVars[iTrial]->GetBinContent(iJetbin+1);
      fJetPtBinYieldDistribution[iJetbin]->Fill(arrYld[iJetbin][iTrial]);
    }

    Double_t mean = TMath::Mean(fnMaxTrials,arrYld[iJetbin]);
    Double_t rms = TMath::RMS(fnMaxTrials,arrYld[iJetbin]);
    if(fDebug) {std::cout << "Jet bin " << iJetbin << " (" << fJetSpectrSBVars[0]->GetXaxis()->GetBinLowEdge(iJetbin+1) << "-" << fJetSpectrSBVars[0]->GetXaxis()->GetBinUpEdge(iJetbin+1) << ")";
      std::cout << ": Mean = " << mean << ", RMS = " << rms << std::endl;
    }

    fJetYieldUnc->SetBinContent(iJetbin+1,rms);
    fJetYieldCentral->SetBinContent(iJetbin+1,mean);
    fJetYieldCentral->SetBinError(iJetbin+1,rms);

    fJetPtBinYieldDistribution[iJetbin]->SaveAs(Form("YieldDistribution_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fJetbinpTedges[iJetbin],fJetbinpTedges[iJetbin+1]));

  }

  //fOut->Close();

  fJetYieldUnc->SetStats(kFALSE);
  fJetYieldUnc->Draw();
  fJetYieldUnc->SaveAs(Form("FinalRawYieldUncertainty_%s.root",fDmesonLabel.Data()));
  fJetYieldCentral->SetStats(kFALSE);
  fJetYieldCentral->Draw();
  fJetYieldCentral->SaveAs(Form("FinalRawYieldCentralPlusSystUncertainty_%s.root",fDmesonLabel.Data()));

  if(fDebug) {
    //ADVANCED - save distribution of final jet yields (summing all pT(D) bins) in a single plot
    TCanvas *cDistr = new TCanvas("cDistr","cDistr",900,600);
    for(Int_t iTrial=0; iTrial<fnMaxTrials; iTrial++) {
      for(int l=0; l<fJetSpectrSBVars[iTrial]->GetNbinsX();l++) fJetSpectrSBVars[iTrial]->SetBinError(l+1,0.0001);
      fJetSpectrSBVars[iTrial]->SetMarkerColor(iTrial+1);
      fJetSpectrSBVars[iTrial]->SetLineColor(iTrial+1);
      if(!iTrial) fJetSpectrSBVars[iTrial]->Draw();
      else fJetSpectrSBVars[iTrial]->Draw("same");
    }
    cDistr->SaveAs(Form("DistributionOfFinalYields_SBApproach_%s_AfterDbinSum.root",fDmesonLabel.Data()));

    //ADVANCED - save distribution of final jet yields from each single pT(D) bin in a single plot (one per each pT(D) bin)
    for(int iDbin=0; iDbin<fnDbins; iDbin++) {
      TCanvas *cDistr1 = new TCanvas(Form("cDistr%d",iDbin),Form("cDistr%d",iDbin),900,600);
      TH1F** hJetSpectrFromSingleDbin = new TH1F*[fnMaxTrials];
      for(Int_t iTrial=0; iTrial<fnMaxTrials; iTrial++) {
        hJetSpectrFromSingleDbin[iTrial] = (TH1F*)fJetSpectrSBVars[0]->Clone(Form("JetRawYieldDistr_Dbin%d",iDbin));
        for(int l=0; l<hJetSpectrFromSingleDbin[iTrial]->GetNbinsX();l++) {
          hJetSpectrFromSingleDbin[iTrial]->SetBinContent(l+1,arrYldBinPerBin[iDbin][l][iTrial]);
          hJetSpectrFromSingleDbin[iTrial]->SetBinError(l+1,0.0001);
        }
        hJetSpectrFromSingleDbin[iTrial]->SetMarkerColor(iTrial+1);
        hJetSpectrFromSingleDbin[iTrial]->SetLineColor(iTrial+1);
        if(!iTrial) hJetSpectrFromSingleDbin[iTrial]->Draw();
        else hJetSpectrFromSingleDbin[iTrial]->Draw("same");
      }
      cDistr1->SaveAs(Form("DistributionOfFinalYields_SBApproach_%s_Bin%d.root",fDmesonLabel.Data(),iDbin));
    }

    //ADVANCED - save averages of final jet yields from each single pT(D) bin, with their RMS, without summing them, in a single plot
    TCanvas *cDistr2 = new TCanvas("cDistrAllAvgs","cDistrAllAvgs",900,600);

    fJetYieldCentral->SetLineWidth(3);
    fJetYieldCentral->Draw();

    TH1F** hJetSpectrFromSingleDbin_Avg = new TH1F*[fnDbins];

    for(Int_t iDbin=0; iDbin<fnDbins; iDbin++) {
      hJetSpectrFromSingleDbin_Avg[iDbin] = (TH1F*)fJetSpectrSBVars[0]->Clone(Form("JetRawYieldAvgDistr_%d",iDbin));

        for(Int_t iJetbin=0; iJetbin<nJetBins; iJetbin++) { //loop on jet spectrum pT bins

          Double_t mean = TMath::Mean(fnMaxTrials,arrYldBinPerBin[iDbin][iJetbin]);
          Double_t rms = TMath::RMS(fnMaxTrials,arrYldBinPerBin[iDbin][iJetbin]);

          fJetYieldUnc->SetBinContent(iJetbin+1,rms);
          hJetSpectrFromSingleDbin_Avg[iDbin]->SetBinContent(iJetbin+1,mean);
          hJetSpectrFromSingleDbin_Avg[iDbin]->SetBinError(iJetbin+1,rms);
        }

      hJetSpectrFromSingleDbin_Avg[iDbin]->SetMarkerColor(iDbin+1);
      hJetSpectrFromSingleDbin_Avg[iDbin]->SetLineColor(iDbin+1);
      hJetSpectrFromSingleDbin_Avg[iDbin]->Draw("same");
    }

    TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->AddEntry(fJetYieldCentral,"Average after pT(D) bin sum","pl");
    for(Int_t iDbin=0; iDbin<fnDbins; iDbin++) leg->AddEntry(hJetSpectrFromSingleDbin_Avg[iDbin],Form("pt(D) %1.1f - %1.1f",fDbinpTedges[iDbin],fDbinpTedges[iDbin+1]),"pl");
    leg->Draw();

    cDistr2->SaveAs(Form("AverageOfFinalYields_SBApproach_%s_AllDBins.root",fDmesonLabel.Data()));

  } //end of advanced plots

  return kTRUE;

}


//___________________________________________________________________________________________
void AliDJetRawYieldUncertaintyLocal::FitReflDistr(Int_t nPtBins, TString inputfile, TString fitType) {

  TFile *fReflections = new TFile(inputfile.Data());
  TFile *fFitReflection = new TFile(Form("reflections_fitted_%s.root",fitType.Data()), "recreate");

  for (Int_t i=0; i<nPtBins; i++){
    TH1F *hSignMC= (TH1F*)fReflections->Get(Form("histSgn_%d",i));
    if(hSignMC) {
      fFitReflection->cd();
      hSignMC->Write(Form("histSgn_%d",i));
    }
  }

  TCanvas *cy = new TCanvas("fitCanv", "fitCanv");
  cy->Divide(4,3);
  TCanvas *cy2 = new TCanvas("fitCanv2", "fitCanv2");
  cy2->Divide(4,3);
  TCanvas *cyRatio = new TCanvas("fitCanvRatio", "fitCanvRatio");
  cyRatio->Divide(4,3);


  for (Int_t iBin=0; iBin<nPtBins; iBin++){
    TH1F *hfitRefl= (TH1F*)fReflections->Get(Form("histRfl_%d", iBin));
    hfitRefl->SetName(Form("histoRfl_%d",iBin));
    hfitRefl->SetMarkerStyle(1);
    hfitRefl->SetLineStyle(1);
    hfitRefl->Rebin(2);
    cy->cd(iBin+1);
    hfitRefl->Draw();
    hfitRefl->Sumw2();

    TF1 *finput=0x0;
    if (fitType == "DoubleGaus"){//DoubleGaus
      finput= new TF1 ("finput","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]/( TMath::Sqrt(2.*TMath::Pi())*[5])*TMath::Exp(-(x-[4])*(x-[4])/(2.*[5]*[5]))");

      finput->SetParameter(0, 1);
      finput->SetParameter(1, 1);
      finput->SetParameter(2, 1);
      finput->SetParameter(3, 1);
      finput->SetParameter(4, 1);
      finput->SetParameter(5, 1);

    }
    else if (fitType == "pol3"){
      finput= new TF1 ("finput", "pol3");
    }
    else if (fitType == "pol6"){
      finput= new TF1 ("finput", "pol6");
    }
    else if (fitType == "gaus"){
      finput= new TF1 ("finput", "gaus");
    }

    hfitRefl->Fit("finput", "MLFI");
    TF1 *fFitRefl = hfitRefl->GetFunction("finput");

    TH1F *hFitReflNewTemp = (TH1F*)hfitRefl->Clone(Form("histRflFitted%s_ptBin%d", fitType.Data(), iBin));
    TH1F *ratio = (TH1F*)hfitRefl->Clone(Form("ratioRelDistr_%s_bin%d", fitType.Data(), iBin));
    for(Int_t iBin2=1; iBin2<=hfitRefl->GetNbinsX(); iBin2++){
      hFitReflNewTemp->SetBinContent(iBin2, 0.);
      ratio->SetBinContent(iBin2, 0.);

      hFitReflNewTemp->SetBinContent(iBin2, fFitRefl->Eval(hfitRefl->GetBinCenter(iBin2)));
      if(fFitRefl->Eval(hfitRefl->GetBinCenter(iBin2))>=0.) hFitReflNewTemp->SetBinError(iBin2, TMath::Sqrt(fFitRefl->Eval(hfitRefl->GetBinCenter(iBin2))));
      else hFitReflNewTemp->SetBinError(iBin2, 0.);
      ratio->SetBinContent(iBin2, (hfitRefl->GetBinContent(iBin2) / fFitRefl->Eval(hfitRefl->GetBinCenter(iBin2))));
    }

    cy2->cd(iBin+1);
    hFitReflNewTemp->Draw();

    //gStyle->SetOptFit(111111);
    cyRatio->cd(iBin+1);
    ratio->GetYaxis()->SetRangeUser(-1.5, 3.);
    ratio->SetMarkerStyle(20);
    ratio->Fit("pol0", "FM");
    ratio->Draw("p");
    gPad->Update();

    fFitReflection->cd();
    hFitReflNewTemp->Write();
    ratio->Write();

  }

  return;

}



//___________________________________________________________________________________________
void AliDJetRawYieldUncertaintyLocal::SetDmesonPtBins(Int_t nbins, Double_t* ptedges) {

  if(!nbins) return;
  fnDbins=nbins;
  fDbinpTedges = new Double_t[fnDbins+1];
  for(int i=0;i<fnDbins+1;i++) {
    fDbinpTedges[i]=ptedges[i];
  }

  return;
}

//___________________________________________________________________________________________
void AliDJetRawYieldUncertaintyLocal::SetJetPtBins(Int_t nbins, Double_t* ptedges) {

  if(!fnJetbins) return;
  fnJetbins=nbins;
  fJetbinpTedges = new Double_t[fnJetbins+1];
  for(int i=0;i<fnJetbins+1;i++) {
    fJetbinpTedges[i]=ptedges[i];
  }

  return;
}

//___________________________________________________________________________________________
void AliDJetRawYieldUncertaintyLocal::SetDmesonEfficiency(Double_t* effvalues) {

  if(!fnDbins) return;
  fDEffValues = new Double_t[fnDbins];
  for(int i=0;i<fnDbins;i++) {
    fDEffValues[i]=effvalues[i];
  }

  return;
}

//___________________________________________________________________________________________
void AliDJetRawYieldUncertaintyLocal::SetMeanSigmaVariations(Bool_t* cases) {

  fMeanSigmaVar = new Bool_t[6];
  for(int i=0;i<6;i++) {
    fMeanSigmaVar[i]=cases[i];
  }

  return;
}

//___________________________________________________________________________________________
void AliDJetRawYieldUncertaintyLocal::SetBkgVariations(Bool_t* cases) {

  fBkgVar = new Bool_t[8];
  for(int i=0;i<8;i++) {
    fBkgVar[i]=cases[i];
  }

  return;
}

//___________________________________________________________________________________________
void AliDJetRawYieldUncertaintyLocal::SetRebinSteps(Int_t nsteps, Int_t* cases) {

  fnRebinSteps = nsteps;
  fRebinSteps = new Int_t[fnRebinSteps];
  for(int i=0;i<fnRebinSteps;i++) {
    fRebinSteps[i]=cases[i];
  }

  return;
}

//___________________________________________________________________________________________
void AliDJetRawYieldUncertaintyLocal::SetMinMassSteps(Int_t nsteps, Double_t* cases) {

  fnMinMassSteps = nsteps;
  fMinMassSteps = new Double_t[fnMinMassSteps];
  for(int i=0;i<fnMinMassSteps;i++) {
    fMinMassSteps[i]=cases[i];
  }

  return;
}

//___________________________________________________________________________________________
void AliDJetRawYieldUncertaintyLocal::SetMaxMassSteps(Int_t nsteps, Double_t* cases) {

  fnMaxMassSteps = nsteps;
  fMaxMassSteps = new Double_t[fnMaxMassSteps];
  for(int i=0;i<fnMaxMassSteps;i++) {
    fMaxMassSteps[i]=cases[i];
  }

  return;
}

//___________________________________________________________________________________________
void AliDJetRawYieldUncertaintyLocal::SetSigmaBinCounting(Int_t nsteps, Double_t* cases) {

  fnSigmaBC = nsteps;
  if(!fnSigmaBC) return;
  fSigmaBC = new Double_t[fnSigmaBC];
  for(int i=0;i<fnSigmaBC;i++) {
    fSigmaBC[i]=cases[i];
  }

  return;
}

//___________________________________________________________________________________________
void AliDJetRawYieldUncertaintyLocal::SetMaskOfVariations(Int_t ncases, Bool_t* cases) {

  fnMask = ncases;
  fMask = new Bool_t[fnMask];
  for(int i=0;i<fnMask;i++) {
    fMask[i]=cases[i];
  }

  return;
}

//__________________________________________
void AliDJetRawYieldUncertaintyLocal::ClearObjects() {

  if(fDbinpTedges) delete[] fDbinpTedges;
  if(fJetbinpTedges) delete[] fJetbinpTedges;
  if(fDEffValues) delete[] fDEffValues;
  if(fMeanSigmaVar) delete[] fMeanSigmaVar;
  if(fBkgVar) delete[] fBkgVar;
  if(fRebinSteps) delete[] fRebinSteps;
  if(fMinMassSteps) delete[] fMinMassSteps;
  if(fMaxMassSteps) delete[] fMaxMassSteps;
  if(fMask) delete[] fMask;

  return;
}
