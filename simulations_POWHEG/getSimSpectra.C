//
// Macro to extract jet pt spectra from simulation, prompt or non-prompt D
// 
//
// Author: B.Trzeciak (barbara.antonina.trzeciak@cern.ch)
//
/*
#include <cstring>
#include <string>
#include <iostream>
#include <vector>

#include <TMath.h>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <THnSparse.h>
#include <TDatabasePDG.h>
#include <TChain.h>
#include <TNtuple.h>
#include <Riostream.h>
*/
//#include <AliAnalysisTaskDmesonJets.h>

double jetptmin = 5, jetptmax = 30; // for D pT spectra
double jetEta = 0.5;

const Int_t nDbins = 10;
double ptDbins[nDbins+1] = {3,4,5,6,7,8,10,12,16,24,36}; // for jet pt spectra, or for efficiency scaling


 TString runB[] = { "1496995621", "1497227464", "1497228262", "1497207414", "1497173624", "1497172782", "1497132057", "1497130041", "1497121734" };
 TString descB[] = {"central", "m_{b}=4.5", "m_{b}=5", "muR=2,muF=2", "muR=1,muF=2", "muR=2,muF=1", "muR=0.5,muF=0.5", "muR=1,muF=0.5" ,"muR=0.5,muF=1" };
 
TString runC[] = { "1496999831", "1497979153", "1497983041", "1497891573", "1497817179", "1497857391", "1497726041", "1497699194", "1497713588" };
//Tstring descC[] = {"central", "m_{c}=1.3", "m_{c}=1.7", "muR=2,muF=2", "muR=1,muF=2", "muR=2,muF=1", "muR=0.5,muF=0.5", "muR=1,muF=0.5" ,"muR=0.5,muF=1"};

 
// quark: 1 = beauty, 0 = charm 

void getSimSpectra(int simNr = 0, int quark = 1, bool jet = 1, bool isjetptcut = 0, bool isDptcut = 0, bool isEff = 0, TString effFilePrompt = "/home/basia/Work/alice/analysis/pPb_run2/Efficiencies/out_806Base/DjetEff_prompt_jetpt2_50.root", TString effFileNonPrompt = "/home/basia/Work/alice/analysis/pPb_run2/Efficiencies/out_806Base_FD/DjetEff_nonPrompt_jetpt2_50.root", TString outHistName = "out/ptSpectrumSim_" );


double GetEfficiency(TH1 *hh, double Dpt);

void ScaleHist(TH1 *hh, int full = 0);
void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, Width_t width, string name);
void SaveCanvas(TCanvas *c, string name = "tmp");

// quark: 1 = beauty, 0 = charm 
void getSimSpectra(int simNr, int quark, bool jet, bool isjetptcut, bool isDptcut, bool isEff, TString effFilePrompt, TString effFileNonPrompt, TString outHistName )
{

    TString simFile = "/home/basia/Work/alice/analysis/fastSim_pPb5TeV/files/AnalysisResults_FastSim_powheg_";
    if(quark == 1) simFile += "beauty";
    else simFile += "charm";
    simFile += "_";
    if(quark == 1) simFile += runB[simNr];
    else if(quark == 0) simFile += runC[simNr];
    simFile += ".root";


    TH1D *hPt;
    if(jet) hPt = (TH1D*) GetInputSimHistJet(simFile.Data(),hPt, isEff, effFilePrompt.Data(), effFileNonPrompt.Data(),isDptcut);
    else hPt = (TH1D*) GetInputSimHistD(simFile.Data(),hPt,isjetptcut);
    
    //if(quark == 1){ hPt->Scale(1./eventsB[simNr]); hPt->Scale(sigma_B[simNr]); }
    //else if(quark == 0) { hPt->Scale(1./eventsC[simNr]); hPt->Scale(sigma_C[simNr]); }
 
    TString out = outHistName;
    if(jet) out += "JetPt_";
    else out += "DPt_";
    if(quark == 1) out += "beauty";
    else out += "charm";
    out += "_";
    if(quark == 1) out += runB[simNr];
    else if(quark == 0) out += runC[simNr];
    if(jet) { 
        if(isDptcut) { out += "_Dpt"; out +=  ptDbins[0]; out += "_"; out += ptDbins[nDbins];  }
    }
    else{ 
        if(isjetptcut){ out += "_Jetpt"; out +=  jetptmin; out += "_"; out += jetptmax; }
    }
    if(jet && isEff) out += "_effScaled";
    out += ".root";
    
    
    TFile *ofile = new TFile( out.Data() ,"RECREATE");
    hPt->Write();
    ofile->Close();
    
 
    return;
}


TH1* GetInputHist(TString inFile, string histName,TH1 *hh){

	TFile *jetPtFile = new TFile(inFile,"read");  
    hh = (TH1F*)jetPtFile->Get(histName.c_str());
 
    //hh = (TH1F*)hh_tmp->Rebin(ptbinsJetN,"hh",ptbinsJet);
  
  return hh;
  
}


TH1* GetInputSimHistJet(TString inFile, TH1 *hPt, bool isEff, TString effFilePrompt, TString effFileNonPrompt,bool isDptcut){

  TFile *fileInput = new TFile(inFile,"read");
  if(!fileInput){
    std::cout << "File " << fileInput << " cannot be opened! check your file path!" << std::endl; return kFALSE;
  }
  
 
 
  TList* dir=(TList*)fileInput->Get("AliAnalysisTaskDmesonJets_histos");
  if(!dir) {
    std::cout << "Error in getting dir! Exiting..." << std::endl;
    return kFALSE;
  }
  
  
  
  TH1D *hxsection = (TH1D*)dir->FindObject("fHistXsectionVsPtHard");
   if(!hxsection) {
    std::cout << "Error in getting x-section hist! Exiting..." << std::endl;
    return kFALSE;
  }
  double xsection = hxsection->GetMean(2);
  double events = (double)hxsection->GetEntries();
  double scaling = xsection/events;
  
 
  
  TTree *tree = (TTree*)fileInput->Get("AliAnalysisTaskDmesonJets_DStar_MCTruth");
  //TTree *tree = (TTree*)fileInput->Get("AliAnalysisTaskDmesonJets_D0_MCTruth");
  AliAnalysisTaskDmesonJets::AliDmesonInfoSummary *brD = 0;
  AliAnalysisTaskDmesonJets::AliJetInfoSummary *brJet = 0;
  tree->SetBranchAddress("DmesonJet",&brD);
  tree->SetBranchAddress("Jet_AKTChargedR040_pt_scheme",&brJet);

  if(!tree || !brD || !brJet) {
    std::cout << "Error in setting the tree/branch names! Exiting..." << std::endl;
    return kFALSE;
}
    
    TH1D *hPromptEff;
    if(isEff) hPromptEff = (TH1D*)GetInputHist(effFilePrompt,"hEff_reb",hPromptEff);
    TH1D *hNonPromptEff;
    if(isEff) hNonPromptEff = (TH1D*)GetInputHist(effFileNonPrompt,"hEff_reb",hNonPromptEff);
    
    
    
    TH1D *hjetpt[nDbins];
    for (int j=0; j<nDbins; j++) {
        hjetpt[j] = new TH1D(Form("hjetpt_%d",j),"hjetpt",100,0,100);
        hjetpt[j]->Sumw2();
    }
    double effC, effB, eff;
    
  
     
    for (int k=0; k<tree->GetEntries(); k++) {
    tree->GetEntry(k);
    if (brJet->fEta < -0.5 || brJet->fEta >= 0.5) continue;
    for (int j=0; j<nDbins; j++) {
    if(isDptcut){  if (brD->fPt < ptDbins[j] || brD->fPt >= ptDbins[j+1]) continue; }
      hjetpt[j]->Fill(brJet->fPt);
    }//end of D-meson pT bin loop
    }

  for (int j=0; j<nDbins; j++) {
      
    double pt = (ptDbins[j]+ptDbins[j+1])/2.;
    if(isEff)  {
        effC = GetEfficiency(hPromptEff,pt);
        effB = GetEfficiency(hNonPromptEff,pt);
        eff = effB / effC;
    }
    else eff = 1;  
     
    if (!j) { 
        hPt = (TH1D*)hjetpt[j]->Clone("hPt");
        hPt->Scale(eff);
    }
    else hPt->Add(hjetpt[j],eff);
    //else hJetPt_B->Add(hjetpt[j]);
  
  }
 // return hPt;
 hPt->Scale(scaling);

  if(!hPt) {
    std::cout << "Error in extracting the mass plot! Exiting..." << std::endl;
    return kFALSE;
  }
  
    
return hPt;

}


TH1* GetInputSimHistD(TString inFile, TH1 *hPt, bool isjetptcut){

  TFile *fileInput = new TFile(inFile,"read");
  if(!fileInput){
    std::cout << "File " << fileInput << " cannot be opened! check your file path!" << std::endl; return kFALSE;
  }
  
  TList* dir=(TList*)fileInput->Get("AliAnalysisTaskDmesonJets_histos");
  if(!dir) {
    std::cout << "Error in getting dir! Exiting..." << std::endl;
    return kFALSE;
  }
  
  TH1D *hxsection = (TH1D*)dir->FindObject("fHistXsectionVsPtHard");
   if(!hxsection) {
    std::cout << "Error in getting x-section hist! Exiting..." << std::endl;
    return kFALSE;
  }
  double xsection = hxsection->GetMean(2);
  double events = (double)hxsection->GetEntries();
  double scaling = xsection/events;
  
  
  TTree *tree = (TTree*)fileInput->Get("AliAnalysisTaskDmesonJets_DStar_MCTruth");
  //TTree *tree = (TTree*)fileInput->Get("AliAnalysisTaskDmesonJets_D0_MCTruth");
  AliAnalysisTaskDmesonJets::AliDmesonInfoSummary *brD = 0;
  AliAnalysisTaskDmesonJets::AliJetInfoSummary *brJet = 0;
  tree->SetBranchAddress("DmesonJet",&brD);
  tree->SetBranchAddress("Jet_AKTChargedR040_pt_scheme",&brJet);

  if(!tree || !brD || !brJet) {
    std::cout << "Error in setting the tree/branch names! Exiting..." << std::endl;
    return kFALSE;
}
    
   hPt = new TH1D("hPt","hDpt",100,0,100);
   // hPt = new TH1D("hPt","hDpt",20,-2,2);
    for (int k=0; k<tree->GetEntries(); k++) {
    tree->GetEntry(k);
    if (brJet->fEta < -0.5 || brJet->fEta >= 0.5) continue;
   if(isjetptcut) { if (brJet->fPt < jetptmin || brJet->fPt >= jetptmax) continue; }
      hPt->Fill(brD->fPt);
      //hPt->Fill(brD->fEta);
    }
    
    hPt->Scale(scaling);

  if(!hPt) {
    std::cout << "Error in extracting the mass plot! Exiting..." << std::endl;
    return kFALSE;
  }
  
return hPt;

}

double GetEfficiency(TH1 *hh, double Dpt){
    
    return hh->GetBinContent(hh->GetXaxis()->FindBin(Dpt));
    
}


void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, Width_t width, string name){

    hh->SetMarkerColor(color);
    hh->SetMarkerStyle(Mstyle);;
    hh->SetLineColor(color);
    hh->SetLineWidth(2);
    hh->SetMarkerSize(1.1);
    hh->SetName(name.c_str());
    
    hh->SetTitle();
    hh->GetXaxis()->SetTitle("p_{T}^{ch,jet} (GeV/c)");
    
}

void SaveCanvas(TCanvas *c, string name){
    
    c->SaveAs(Form("%s.png",name.c_str()));
    c->SaveAs(Form("%s.pdf",name.c_str()));
   
}
