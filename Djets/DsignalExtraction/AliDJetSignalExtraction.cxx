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

#include "AliDJetSignalExtraction.h"

//___________________________________________________________________________________________
AliDJetSignalExtraction::AliDJetSignalExtraction():
fFileInput(0x0),
fFileRawOutput(0x0),
fDmesonSpecie(kDStarD0pi),
fDmesonLabel("Dstar"),
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
AliDJetSignalExtraction::AliDJetSignalExtraction(const AliDJetRawYieldUncertaintyLocal &source):
// copy constructor
fFileInput(source.fFileInput),
fFileRawOutput(source.fFileRawOutput),
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
AliDJetSignalExtraction::~AliDJetSignalExtraction() {
//destructor

}

//___________________________________________________________________________________________
Bool_t AliDJetSignalExtraction::SetDmesonSpecie(DMesonSpecies k){

  if(k<0 || k>1) {
    printf("Error! D meson specie not correctly set!\n");
    return kFALSE;
  } else if(k==0) fDmesonLabel="Dzero";
  else fDmesonLabel="Dstar";

  fDmesonSpecie=k;
  return kTRUE;

}

//___________________________________________________________________________________________
Bool_t AliDJetSignalExtraction::ExtractInputMassPlot(){

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

  if(success) std::cout << "Extracted mass spectrum" << std::endl;

 
  //fMassPlot->SaveAs("spettrotemp.root"); //DEBUG TEMP

  //fFileInput->Close(); 

  return success;

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
    TH1F *hmass;

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
void AliDJetSignalExtraction::SetDmesonPtBins(Int_t nbins, Double_t* ptedges) {

  if(!nbins) return;
  fnDbins=nbins;
  fDbinpTedges = new Double_t[fnDbins+1];
  for(int i=0;i<fnDbins+1;i++) {
    fDbinpTedges[i]=ptedges[i];
  }

  return;
}

//___________________________________________________________________________________________
void AliDJetSignalExtraction::SetJetPtBins(Int_t nbins, Double_t* ptedges) {

  if(!fnJetbins) return;
  fnJetbins=nbins;
  fJetbinpTedges = new Double_t[fnJetbins+1];
  for(int i=0;i<fnJetbins+1;i++) {
    fJetbinpTedges[i]=ptedges[i];
  }

  return;
}

//___________________________________________________________________________________________
void AliDJetSignalExtraction::SetDmesonEfficiency(Double_t* effvalues) {

  if(!fnDbins) return;
  fDEffValues = new Double_t[fnDbins];
  for(int i=0;i<fnDbins;i++) {
    fDEffValues[i]=effvalues[i];
  }

  return;
}


//__________________________________________
void AliDJetSignalExtraction::ClearObjects() {
   
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

