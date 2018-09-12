//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------
//
// Macro to extract jet pt spectra from simulation, prompt or non-prompt D
//

#include "config.h"

//fDmesonSpecie = 1;
//fRpar = 0.4;
//Rpar = 4;

double jetptmin = 5, jetptmax = 30; // for D pT spectra
double jetEta = 0.9-fRpar;
int BFDsim;



//quark: 1 = beauty, 0 = charm
void getSimSpectra(
TString simFile = "./",
int simNr = 0, int quark = 1, bool jet = 1,  bool isDptcut = 1, bool isEff = 0,
TString effFilePrompt = "$HOME/file.root",
TString effFileNonPrompt = "$HOME/file.root",
TString outFileDir = "outR03Test/",
bool isjetptcut = 0, double jetmin = 5, double jetmax = 30 );

double GetEfficiency(TH1 *hh, double Dpt);
void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, Width_t width, string name);
void SaveCanvas(TCanvas *c, string name = "tmp");

//quark: 1 = beauty, 0 = charm
void getSimSpectra(TString simFile, int simNr,
  int quark, bool jet, bool isDptcut,
  bool isEff, TString effFilePrompt, TString effFileNonPrompt,
  TString outFileDir,
  bool isjetptcut, double jetmin, double jetmax )
{

    gSystem->Exec(Form("mkdir -p %s",outFileDir.Data()));

    jetptmin = jetmin;
    jetptmax = jetmax;
    BFDsim = quark;

    simFile += "/";
    if(quark == 1) simFile += fRunB[simNr];
    else if(quark == 0) simFile += fRunC[simNr];
  /*  if(quark == 1) simFile += "beauty";
    else simFile += "charm";
    simFile += "_";
    if(quark == 1) simFile += fRunB[simNr];
    else if(quark == 0) simFile += fRunC[simNr];*/
    simFile += ".root";

    TH1D *hPt;
    if(jet) hPt = (TH1D*) GetInputSimHistJet(simFile.Data(),hPt, isEff, effFilePrompt.Data(), effFileNonPrompt.Data(),isDptcut,quark,simNr);
    else hPt = (TH1D*) GetInputSimHistD(simFile.Data(),hPt,isjetptcut);

    TString out = outFileDir;
    if(jet) out += "/JetPt_";
    else out += "/DPt_";
  //  if(quark == 1) out += "beauty";
  //  else out += "charm";
  //  out += "_";
    if(quark == 1) out += fRunB[simNr];
    else if(quark == 0) out += fRunC[simNr];

    if(jet){
        if(isDptcut) { out += "_Dpt"; out +=  fptbinsDA[0]; out += "_"; out += fptbinsDA[fptbinsDN];  }
    }
    else{
        if(isjetptcut){ out += "_Jetpt"; out +=  jetptmin; out += "_"; out += jetptmax; }
    }

    if(jet && isEff) out += "_effScaled";
    if(fDmesonSpecie) out += "_Dstar";
    else out += "_Dzero";
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

TH1* GetInputSimHistJet(TString inFile, TH1 *hPt, bool isEff, TString effFilePrompt, TString effFileNonPrompt,bool isDptcut, int isNPrompt, int SimNr){

    TFile *fileInput = new TFile(inFile,"read");
    if(!fileInput){
      std::cout << "File " << fileInput << " cannot be opened! check your file path!" << std::endl; return kFALSE;
    }

    TList* dir=(TList*)fileInput->Get("AliAnalysisTaskDmesonJets_histos");
    if(!dir) {
      std::cout << "Error in getting dir! Exiting..." << std::endl;
      return kFALSE;
    }

    TH1D *hxsection;
    if(!isNPrompt)	hxsection = (TH1D*)dir->FindObject("fHistXsection");
    else { 
          if(SimNr == 11)  hxsection = (TH1D*)dir->FindObject("fHistXsection");
          else hxsection = (TH1D*)dir->FindObject("fHistXsectionVsPtHard");
    }   
 
    if(!hxsection) {
      std::cout << "Error in getting x-section hist! Exiting..." << std::endl;
      return kFALSE;
    }
    double xsection = hxsection->GetMean(2);
    double events = (double)hxsection->GetEntries();
    double scaling = xsection/events;

    TTree *tree;
    if(!fDmesonSpecie) tree = (TTree*)fileInput->Get("AliAnalysisTaskDmesonJets_D0_MCTruth");
    else tree = (TTree*)fileInput->Get("AliAnalysisTaskDmesonJets_DStar_MCTruth");
    AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary *brD = 0;

    AliAnalysisTaskDmesonJets::AliJetInfoSummary *brJet = 0;
    tree->SetBranchAddress("DmesonJet",&brD);
    tree->SetBranchAddress(Form("Jet_AKTChargedR0%d0_pt_scheme",Rpar),&brJet);

    if(!tree || !brD || !brJet) {
      std::cout << "Error in setting the tree/branch names! Exiting..." << std::endl;
      return NULL;
    }

    TH1D *hPromptEff = NULL;
    if(isEff) hPromptEff = (TH1D*)GetInputHist(effFilePrompt,"hEff_reb",hPromptEff);
    TH1D *hNonPromptEff = NULL;
    if(isEff) hNonPromptEff = (TH1D*)GetInputHist(effFileNonPrompt,"hEff_reb",hNonPromptEff);


    TH1D *hjetpt[fptbinsDN];
    for (int j=0; j<fptbinsDN; j++) {
        hjetpt[j] = new TH1D(Form("hjetpt_%d",j),"hjetpt",100,0,100);
        hjetpt[j]->Sumw2();
    }
    double effC, effB, eff;
    hPt = new TH1D("hPt","hjetpt",100,0,100);

    for (int k=0; k<tree->GetEntries(); k++) {
    tree->GetEntry(k);
    if (brJet->fEta < -jetEta || brJet->fEta > jetEta) continue;


    if(BFDsim){
      if(brD->fPartonType != 5) continue;
    }
    else if(brD->fPartonType != 4) continue;
    if(brD->fAncestorPDG == 2212) continue; // check if not coming from proton

    if(isDptcut){
      for (int j=0; j<fptbinsDN; j++) {
        if (brD->fPt < fptbinsDA[j] || brD->fPt >= fptbinsDA[j+1]) continue;
        hjetpt[j]->Fill(brJet->fPt);
      }//end of D-meson pT bin loop
    }
    else hPt->Fill(brJet->fPt);
    }

if(isDptcut){
  for (int j=0; j<fptbinsDN; j++){
    double pt = (fptbinsDA[j]+fptbinsDA[j+1])/2.;
    if(isEff)  {
        effC = GetEfficiency(hPromptEff,pt);
        effB = GetEfficiency(hNonPromptEff,pt);
        eff = effB / effC;
    }
    else eff = 1;
    if (!j){
        hPt = (TH1D*)hjetpt[j]->Clone("hPt");
        hPt->Scale(eff);
    }
    else hPt->Add(hjetpt[j],eff);
  }
}

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
    if(!dir){
      std::cout << "Error in getting dir! Exiting..." << std::endl;
      return NULL;
    }

    TH1D *hxsection = (TH1D*)dir->FindObject("fHistXsectionVsPtHard");
    if(!hxsection){
      std::cout << "Error in getting x-section hist! Exiting..." << std::endl;
      return NULL;
    }
    double xsection = hxsection->GetMean(2);
    double events = (double)hxsection->GetEntries();
    double scaling = xsection/events;



    TTree *tree;
    if(!fDmesonSpecie) tree = (TTree*)fileInput->Get("AliAnalysisTaskDmesonJets_D0_MCTruth");
    else tree = (TTree*)fileInput->Get("AliAnalysisTaskDmesonJets_DStar_MCTruth");
    AliAnalysisTaskDmesonJets::AliDmesonInfoSummary *brD = 0;

    AliAnalysisTaskDmesonJets::AliJetInfoSummary *brJet = 0;
    tree->SetBranchAddress("DmesonJet",&brD);

    tree->SetBranchAddress(Form("Jet_AKTChargedR0%d0_pt_scheme",Rpar),&brJet);

    if(!tree || !brD || !brJet) {
      std::cout << "Error in setting the tree/branch names! Exiting..." << std::endl;
      return NULL;
    }

    hPt = new TH1D("hPt","hDpt",100,0,100);
    // hPt = new TH1D("hPt","hDpt",20,-2,2);
    for (int k=0; k<tree->GetEntries(); k++) {
    tree->GetEntry(k);
    //if(isjetptcut) { if (brJet->fEta < -jetEta || brJet->fEta > jetEta) continue; }
    //if (brJet->fEta < -jetEta || brJet->fEta > jetEta) continue;
    if(isjetptcut) { if (brJet->fPt < jetptmin || brJet->fPt >= jetptmax) continue; }
    if (brD->fEta < -1 || brD->fEta > 1) continue;
      hPt->Fill(brD->fPt);
      //hPt->Fill(brD->fEta);
    }

    hPt->Scale(scaling);

    if(!hPt) {
      std::cout << "Error in extracting the mass plot! Exiting..." << std::endl;
      return NULL;
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
