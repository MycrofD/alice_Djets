//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------
//
//-----------------------------------------------------------------------
//  [Modified] A.Mohanty
//  Utrecht University
//  auro.mohanty@cern.ch
//-----------------------------------------------------------------------
//
// Macro to extract z_||,ch spectra from simulation, prompt or non-prompt D
//

#include "../DsignalExtraction/configDzero_ppz.h"

//fDmesonSpecie = 1;
//fRpar = 0.4;
//Rpar = 4;

double jetEta = 0.9-fRpar;
int BFDsim;

//quark: 1 = beauty, 0 = charm
void getSimSpectra_z(
TString simFile = "/home/jackbauer/ALICE_HeavyFlavour/work/Djets/out/outMC/beauty/",//$2
int simNr = 1, //$count
int quark = 0, //$3
bool zfrac = 1,  //$4
bool isDptcut = 1, //$5
bool isEff = 0,//$6
//TString effFilePrompt = "/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results/DzeroR03_pPbCuts/Default/efficiency/DjetEff_prompt_jetpt5_50.root",//$7
//TString effFileNonPrompt = "/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results/DzeroR03_pPbCuts/Default/efficiency/DjetEff_nonPrompt_jetpt5_50.root",//$8
TString effFilePrompt = "",//$7
TString effFileNonPrompt = "",//$8
TString outFileDir = "/media/jackbauer/data/z_out/SimFiles/"//"outR03Test/"//$9
);

double GetEfficiency(TH1 *hh, double Dpt);
void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, Width_t width, string name);
void SaveCanvas(TCanvas *c, string name = "tmp");
TH1* GetInputSimHistJet(TString inFile, TH1 *hPt, bool isEff, TString effFilePrompt, TString effFileNonPrompt, bool isDptcut, int isNPrompt, int SimNr);

//quark: 1 = beauty, 0 = charm
void getSimSpectra_z(
TString simFile, 
int simNr,
int quark, 
bool zfrac, 
bool isDptcut,
bool isEff, 
TString effFilePrompt, 
TString effFileNonPrompt,
TString outFileDir
)
{
effFilePrompt+=Form("%d_%d.root",(int)fptbinsJetA[(int)zjetbin-1], (int)fptbinsJetA[(int)zjetbin]);
effFileNonPrompt+=Form("%d_%d.root",(int)fptbinsJetA[(int)zjetbin-1], (int)fptbinsJetA[(int)zjetbin]);

    TString subdir = Form("/Jetpt%d_%d",(int)fptbinsJetA[(int)zjetbin-1],(int)fptbinsJetA[(int)zjetbin]);
//    gSystem->Exec(Form("mkdir -p %s",outFileDir.Data()));
    outFileDir+=subdir;
    gSystem->Exec(Form("mkdir -p %s",outFileDir.Data()));

    BFDsim = quark;

    simFile += "/";
    if(quark == 1) simFile += fRunB[simNr];
    else if(quark == 0) simFile += fRunC[simNr];
    simFile += ".root";

    TString out = outFileDir;
    if(zfrac) out += "/Z_";
    else out += "/DPt_";
    if(quark == 1) out += fRunB[simNr];
    else if(quark == 0) out += fRunC[simNr];

    if(zfrac){
        if(isDptcut) { out += "_Jetpt"; out +=  (int)fptbinsJetA[(int)zjetbin-1]; out += "_"; out += (int)fptbinsJetA[(int)zjetbin];
			out += "_Dpt"; out +=  fptbinsDA[0]; out += "_"; out += fptbinsDA[fptbinsDN];  
	}
    }
    else{return NULL;
    }

    if(isEff) out += "_effScaled";
    if(fDmesonSpecie) out += "_Dstar";
    else out += "_Dzero";
    out += ".root";

// sanity-check if file exists?
if(!gSystem->AccessPathName(Form("%s",out.Data()))){return;}

    TH1D *hPt;
    if(zfrac) hPt = (TH1D*) GetInputSimHistJet(simFile.Data(),hPt, isEff, effFilePrompt.Data(), effFileNonPrompt.Data(),isDptcut,quark,simNr);
    else {cout<<"---------- NOT DOING D analysis here. This is Z analysis. Check older macro--------------"<<endl; return NULL;}//hPt = (TH1D*) GetInputSimHistD(simFile.Data(),hPt,isjetptcut);

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
      std::cout << "File " << fileInput << " cannot be opened! check your file path!" << std::endl; return NULL;
    }

    TList* dir=(TList*)fileInput->Get("AliAnalysisTaskDmesonJets_histos");
    if(!dir) {
      std::cout << "Error in getting dir! Exiting..." << std::endl;
      return NULL;
    }

    TH1D *hxsection;
    if(!isNPrompt)	hxsection = (TH1D*)dir->FindObject("fHistXsection");
    else { 
          //if(SimNr == 0 || SimNr > 11)  
		hxsection = (TH1D*)dir->FindObject("fHistXsection");//for eventgen
          //else hxsection = (TH1D*)dir->FindObject("fHistXsectionVsPtHard");//for normal
          ////else hxsection = (TH1D*)dir->FindObject("fHistXsection");
          //hxsection = (TH1D*)dir->FindObject("fHistXsectionVsPtHard");
    }   
 
    if(!hxsection) {
      std::cout << "Error in getting x-section hist! Exiting..." << std::endl;
      return NULL;
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
        //hjetpt[j] = new TH1D(Form("hjetpt_%d",j),"hjetpt",10,0,1.0);
        hjetpt[j] = new TH1D(Form("hjetpt_%d",j),"hjetpt",fptbinsZMeasN,fptbinsZMeasA);
        hjetpt[j]->Sumw2();
    }
    double effC, effB, eff;
    //hPt = new TH1D("hPt","hjetpt",10,0.0,1.0);
    hPt = new TH1D("hPt","hjetpt",fptbinsZMeasN,fptbinsZMeasA);

    for (int k=0; k<tree->GetEntries(); k++) {
    tree->GetEntry(k);
    if (brJet->fEta < -jetEta || brJet->fEta > jetEta) continue;


    if(BFDsim){
      if(brD->fPartonType != 5) continue;
    }
    else if(brD->fPartonType != 4) continue;
    if(brD->fAncestorPDG == 2212) continue; // check if not coming from proton

    if (brJet->fPt < (int)fptbinsJetA[(int)zjetbin-1] || brJet->fPt >= (int)fptbinsJetA[(int)zjetbin]) continue;
    if(isDptcut){
      for (int j=0; j<fptbinsDN; j++) {
        if  (brD->fPt < fptbinsDA[j] || brD->fPt >= fptbinsDA[j+1])
	    continue;
        hjetpt[j]->Fill(brJet->fZ);
      }//end of D-meson pT bin loop
    }
    else hPt->Fill(brJet->fZ);
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
    hh->SetTitle("");
    hh->GetXaxis()->SetTitle("p_{T}^{ch,jet} (GeV/c)");
}

void SaveCanvas(TCanvas *c, string name){
    c->SaveAs(Form("%s.png",name.c_str()));
    c->SaveAs(Form("%s.pdf",name.c_str()));
}
