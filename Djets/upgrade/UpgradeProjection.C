//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------

#include "config.h"

using namespace std;

void Init() {

  fIsHBins = false;

  fJetptmin = 0;
  fJetptmax = 100;

  fRaaC = new double[fPtbinsDN];
  fRaaB = new double[fPtbinsDN];
  for(int i=0;i<fPtbinsDN;i++){
    fRaaC[i] = 0.2;
    fRaaB[i] = 0.4;
  }

  double weightC[6] = {1,1,1,1,1,1};
  double weightB[6] = {1,1,1,1,1,1};
  for(int i=0; i<6; i++){
    fHBWeightsC[i] = weightC[i];
    fHBWeightsB[i] = weightB[i];
  }


}

void UpgradeProjection() {

  gStyle->SetOptStat(0000); //Mean and RMS shown
	//gSystem->Exec(Form("mkdir -p %s",OUTDIRECTORY.Data()));

  Init();

  TH1D *hJetMCPt;
  if(fIsHBins) hJetMcPt = (TH1D*)getMCJetPt(); // x-check for hard-bin merging

  // get prompt efficiency
  fHEff_prompt = (TH1D*)getEfficiency(fFileEff_prompt);
  if(!fHEff_prompt) { cout << "\n!!! Prompt efficincy not extracted !!! " << endl; return; }
  setHistoDetails(fHEff_prompt,2,20,1.2);

  // get non-prompt getEfficiency
  fHEff_nprompt = (TH1D*)getEfficiency(fFileEff_nprompt);
  if(!fHEff_nprompt) { cout << "\n!!! Non-Prompt efficincy not extracted !!! " << endl; return; }
  setHistoDetails(fHEff_nprompt,4,25,1.2);

  // get invariant mass histogram for the prompt D-jet signal
  TH3D *hInvMassSignal = (TH3D*)getInvMass(fFileSignal_prompt+fSMBBin);
  if(!hInvMassSignal) { cout << "\n!!! Prompt signal not extracted !!! " << endl; return; }
  TH3D *hInvMassBkg = (TH3D*)getInvMass(fFileBkg,1);
  if(!hInvMassBkg) { cout << "\n!!! Bkg not extracted !!! " << endl; return; }

  // get D x-section from POWHEG
  fHPowhegDPrompt = (TH1D*)getPowhegDSpectra(fFilePowhegPrompt);
  if(!fHPowhegDPrompt) { cout << "\n!!!Prompt D POWHEG spectrum not extracted !!!" << endl; return; }
  setHistoDetails(fHPowhegDPrompt,2,20,1.2);
  fHPowhegDNonPrompt = (TH1D*)getPowhegDSpectra(fFilePowhegNonPrompt,0);
  if(!fHPowhegDNonPrompt) { cout << "\n!!! Non-Prompt D POWHEG spectrum not extracted !!!" << endl; return; }
  setHistoDetails(fHPowhegDNonPrompt,4,25,1.2);

  extractSignal(hInvMassSignal, hInvMassBkg);

//----------- x-check drawing //

  TCanvas* cEffReb = new TCanvas();
  cEffReb->SetLogy();
  fHEff_prompt->Draw("ep");
  fHEff_nprompt->Draw("epsame");

  TH1D *hPowhegDInclusive = (TH1D*)fHPowhegDPrompt->Clone("hPowhegDInclusive");
  hPowhegDInclusive->Add(fHPowhegDNonPrompt);
  setHistoDetails(hPowhegDInclusive,8,25,1.2);
  fHPowhegDPrompt->SetMinimum(1e-6);
  TCanvas *cp = new TCanvas();
  cp->SetLogy();
  fHPowhegDPrompt->Draw();
  fHPowhegDNonPrompt->Draw("same");
  hPowhegDInclusive->Draw("same");

  TH1D* hPowhegRatio = (TH1D*)fHPowhegDNonPrompt->Clone("hPowhegRatio");
  hPowhegRatio->Divide(hPowhegDInclusive);
  hPowhegRatio->GetYaxis()->SetRangeUser(0,0.5);
  TCanvas *cpr = new TCanvas();
  hPowhegRatio->Draw();

//---------------------- //

return;

  // get POWHEG D-jet x-sections
  fHPowhegPrompt = getPowhegJetSpectra(fFilePowhegPrompt);
  fHPowhegNonPrompt = getPowhegJetSpectra(fFilePowhegNonPrompt,0);
  TH1D *hPowhegNonPromptScaled = getPowhegJetSpectra(fFilePowhegNonPrompt,0,1);



return;
}

void extractSignal(TH3D *hsignal, TH3D *hbkg) {

  double bkgscale = getBkgScaling();

  TCanvas *cmass = new TCanvas("cmass","cmass",1200,1200);
  cmass->Divide(3,4);

  TCanvas *cpt = new TCanvas("cpt","cpt",1200,1200);
  cpt->Divide(3,4);

  TH1D *hmean = new TH1D("hmean","hmean",fPtbinsDN,fPtbinsDA);
  TH1D *hsigma = new TH1D("hsigma","hsigma",fPtbinsDN,fPtbinsDA);

  for(int i=0; i<fPtbinsDN; i++){

      double pt = (fPtbinsDA[i]+fPtbinsDA[i+1])/2.;
      double effP = GetContent(fHEff_prompt,pt);
      double effNP = GetContent(fHEff_nprompt,pt);
      double sigP = GetContent(fHPowhegDPrompt,pt) * fSimScalingC * fRaaC[i] * effP * fDataEv;
      double sigNP = GetContent(fHPowhegDNonPrompt,pt) * fSimScalingB *fRaaB[i] * effNP * fDataEv;

      cout << "\n\npT: " << pt << "\ttotalC: " << sigP << "\ttotalB: " << sigNP << "\ttotal: " << sigP+sigNP << endl << endl;

      TH1F *hhsignal = (TH1F*)hsignal->ProjectionX(Form("hh_%d",i),hsignal->GetYaxis()->FindBin(fJetptmin), hsignal->GetYaxis()->FindBin(fJetptmax)-1,hsignal->GetZaxis()->FindBin(fPtbinsDA[i]), hsignal->GetZaxis()->FindBin(fPtbinsDA[i+1])-1);
      //hhsignal->Rebin(fRebinMass);
      hhsignal->GetXaxis()->SetRangeUser(minf,maxf);
      hhsignal->SetTitle(Form("%.1lf < pt^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));
      hhsignal->Scale(1./hhsignal->Integral());
      //hhsignal->Scale(sigP);
      hhsignal->Scale(sigP+sigNP);
      for(int k=0; k<hhsignal->GetNbinsX();k++){
        hhsignal->SetBinError(k+1,TMath::Sqrt(hhsignal->GetBinContent(k+1)));
      }

      /*
      TH1F *hmassfit = (TH1F*)hhsignal->Clone("hmassfit");
      float hmin = TMath::Max(minf,hmassfit->GetBinLowEdge(2));
      float hmax = TMath::Min(maxf,hmassfit->GetBinLowEdge(hmassfit->GetNbinsX()));
      */

      TF1 *fitG = new TF1("fitG","gaus",1.78,1.92);
      hhsignal->Fit("fitG","0RM");
      fitG->SetLineColor(2);

      double sigma = fitG->GetParameter(2);
      double mean = fitG->GetParameter(1);
      hsigma->SetBinContent(i+1,fitG->GetParameter(2));
      hsigma->SetBinError(i+1,fitG->GetParError(2));
      hmean->SetBinContent(i+1,fitG->GetParameter(1));
      hmean->SetBinError(i+1,fitG->GetParError(1));

      double massmin = 1.78; // mean - 3*sigma;
      double massmax = 1.92; // mean + 3*sigma;
      double sbmin1 = 1.71, sbmax1 = 1.77;
      double sbmin2 = 1.93, sbmax2 = 2.1;

      TH1F *hhbkg = (TH1F*)hbkg->ProjectionX(Form("hhbkg_%d",i),hbkg->GetYaxis()->FindBin(fJetptmin), hbkg->GetYaxis()->FindBin(fJetptmax)-1,hbkg->GetZaxis()->FindBin(fPtbinsDA[i]), hbkg->GetZaxis()->FindBin(fPtbinsDA[i+1])-1);
      //hhsignal->Rebin(fRebinMass);
      hhbkg->GetXaxis()->SetRangeUser(minf,maxf);
      hhbkg->SetTitle(Form("%.1lf < pt^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));
      //hhbkg->Scale(1./hhbkg->Integral());
      hhbkg->Scale(bkgscale);
      for(int k=0; k<hhbkg->GetNbinsX();k++){
        hhbkg->SetBinError(k+1,TMath::Sqrt(hhbkg->GetBinContent(k+1)));
      }


      TF1 *fitB = new TF1("fitB","[0]*[1]/(TMath::Exp([1]*2.1)-TMath::Exp([1]*1.71))*TMath::Exp([1]*x)",1.78,1.92);
      hhbkg->Fit("fitBkg","0RM");
      fitB->SetLineColor(2);

    //  TF1* funcbkg =  new TF1(fname.Data(),this,&FitFunction4Bkg,fMinMass,fMaxMass,fNParsBkg,"AliHFInvMassFitter","FitFunction4Bkg");
    //  funcbkg->SetParNames("BkgInt","Slope");
    //  funcbkg->SetParameters(integral,-2.);

      cmass->cd(i+1);
      TVirtualPad *pad = (TVirtualPad*)cmass->GetPad(i+1);
      //hhsignal->Draw();
      hhbkg->Draw();
      fitG->Draw("same");

      TH1F *hhsignalPt = (TH1F*)hsignal->ProjectionY(Form("hhsignalPt_%d",i),hsignal->GetXaxis()->FindBin(massmin), hsignal->GetXaxis()->FindBin(massmax)-1,hsignal->GetZaxis()->FindBin(fPtbinsDA[i]), hsignal->GetZaxis()->FindBin(fPtbinsDA[i+1])-1);
      //hhsignal->Rebin(fRebinMass);
      hhsignalPt->GetXaxis()->SetRangeUser(0,50);
      hhsignalPt->SetTitle(Form("%.1lf < pt^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));
      hhsignalPt->Scale(1./hhsignalPt->Integral());
      //hhsignal->Scale(sigP);
      hhsignalPt->Scale(sigP+sigNP);

      for(int k=0; k<hhsignalPt->GetNbinsX();k++){
        hhsignalPt->SetBinError(k+1,TMath::Sqrt(hhsignalPt->GetBinContent(k+1)));
      }

      cpt->cd(i+1);
      TVirtualPad *pad = (TVirtualPad*)cpt->GetPad(i+1);
      hhsignalPt->Draw();

      //cout << "\n mass int: " << hhsignal->Integral() << "\tpt int: " << hhsignalPt->Integral() << endl << endl;

  }


  return;
}

Double_t FitFunction4Bkg(Double_t *x, Double_t *par){

  Double_t maxDeltaM = fNSigma4SideBands*fSigmaSgn;
  if(fOnlySideBands && TMath::Abs(x[0]-fMass) < maxDeltaM) {
    TF1::RejectPoint();
    return 0;
  }
  if(fOnlySideBands && fSecondPeak && TMath::Abs(x[0]-fSecMass) < (fNSigma4SideBands*fSecWidth)){
    TF1::RejectPoint();
    return 0;
  }
  Double_t total=0;
  //exponential
   //exponential = A*exp(B*x) -> integral(exponential)=A/B*exp(B*x)](min,max)
   //-> A = B*integral/(exp(B*max)-exp(B*min)) where integral can be written
   //as integralTot- integralGaus (=par [2])
   //Par:
   // * [0] = integralBkg;
   // * [1] = B;
   //exponential = [1]*[0]/(exp([1]*max)-exp([1]*min))*exp([1]*x)
   total = par[0]*par[1]/(TMath::Exp(par[1]*fMaxMass)-TMath::Exp(par[1]*fMinMass))*TMath::Exp(par[1]*x[0]);
   //    AliInfo("Background function set to: exponential");

}

double getBkgScaling(){

  TFile *File = new TFile(fFileBkg,"read");;
  if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return ;}
  TDirectoryFile *dir = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
  TList *histList = (TList*)dir->Get("histosD0UpgradeN0MCrec");
  TH1F* hEvents = (TH1F*)histList->FindObject("hstat");
  //double nEvSel = hEvents->GetBinContent(2);
  double nEvAn = hEvents->GetBinContent(1); // number of analyzed events
  double scale = fDataEv / nEvAn;

  return scale;

}

TH1* getPowhegJetSpectra(TString data, Bool_t isPrompt=1, Bool_t isEffScale=0){

  TFile *fileInput = new TFile(data,"read");
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

  TTree *tree;
  if(!fDmesonSpecie) tree = (TTree*)fileInput->Get("AliAnalysisTaskDmesonJets_D0_MCTruth");
  else tree = (TTree*)fileInput->Get("AliAnalysisTaskDmesonJets_DStar_MCTruth");
  AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary *brD = 0;
  AliAnalysisTaskDmesonJets::AliJetInfoSummary *brJet = 0;
  tree->SetBranchAddress("DmesonJet",&brD);
  //tree->SetBranchAddress(Form("Jet_AKTChargedR0%d0_pt_scheme",Rpar),&brJet);
  tree->SetBranchAddress(Form("Jet_AKTChargedR0%d0_pt_scheme",3),&brJet);
  if(!tree || !brD || !brJet) {
    std::cout << "Error in setting the tree/branch names! Exiting..." << std::endl;
    return NULL;
  }

  TH1D *hjetpt[fPtbinsDN];
  for (int j=0; j<fptbinsDN; j++) {
      hjetpt[j] = new TH1D(Form("hjetpt_%d",j),"hjetpt",100,0,100);
      hjetpt[j]->Sumw2();
  }

  for (int k=0; k<tree->GetEntries(); k++) {
    tree->GetEntry(k);
    if (brJet->fEta < -fJetEta || brJet->fEta > fJetEta) continue;
    if(!isPrompt){
      if(brD->fPartonType != 5) continue;
    }
    else if(brD->fPartonType != 4) continue;
    if(brD->fAncestorPDG == 2212) continue; // check if not coming from proton

    for (int j=0; j<fptbinsDN; j++) {
        if (brD->fPt < fptbinsDA[j] || brD->fPt >= fptbinsDA[j+1]) continue;
        hjetpt[j]->Fill(brJet->fPt);
    }//end of D-meson pT bin loop

  }

  TH1D*  hPt = new TH1D("hPt","hjetpt",100,0,100);
  for (int j=0; j<fptbinsDN; j++){
    double effC=1, effB=1, eff=1;
    double pt = (fptbinsDA[j]+fptbinsDA[j+1])/2.;
    if(isEffScale)  {
      effC = GetContent(fHEff_prompt,pt);
      effB = GetContent(fHEff_nprompt,pt);
      eff = effB/effC;
    }
    else eff = 1;
    if (!j){
        hPt = (TH1D*)hjetpt[j]->Clone("hPt");
        hPt->Scale(eff);
    }
    else hPt->Add(hjetpt[j],eff);
  }

  hPt->Scale(scaling);

  return hPt;

}

TH1* getPowhegDSpectra(TString data, Bool_t isPrompt=1){

  TFile *fileInput = new TFile(data,"read");
  if(!fileInput){
    std::cout << "File " << fileInput << " cannot be opened! check your file path!" << std::endl; return kFALSE;
  }
  TList* dir=(TList*)fileInput->Get("AliAnalysisTaskDmesonJets_histos");
  if(!dir) {
    std::cout << "Error in getting dir! Exiting..." << std::endl;
    return kFALSE;
  }

  TH1D *hxsection = (TH1D*)dir->FindObject("fHistXsection");
  //TH1D *hxsection = (TH1D*)dir->FindObject("fHistXsectionVsPtHard");
   if(!hxsection) {
    std::cout << "Error in getting x-section hist! Exiting..." << std::endl;
    return kFALSE;
  }
  double xsection = hxsection->GetMean(2);
  double events = (double)hxsection->GetEntries();
  double scaling = xsection/events;

//  cout << "\n\n\n xsection: " << xsection << "\t\t entries: " << events << endl << endl << endl;

  TTree *tree;
  if(!fDmesonSpecie) tree = (TTree*)fileInput->Get("AliAnalysisTaskDmesonJets_D0_MCTruth");
  else tree = (TTree*)fileInput->Get("AliAnalysisTaskDmesonJets_DStar_MCTruth");
  AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary *brD = 0;
  AliAnalysisTaskDmesonJets::AliJetInfoSummary *brJet = 0;
  tree->SetBranchAddress("DmesonJet",&brD);
  //tree->SetBranchAddress(Form("Jet_AKTChargedR0%d0_pt_scheme",Rpar),&brJet);
  tree->SetBranchAddress("Jet_AKTChargedR030_pt_scheme",&brJet);
  if(!tree || !brD || !brJet) {
    std::cout << "Error in setting the tree/branch names! Exiting..." << std::endl;
    return NULL;
  }

  TH1D*  hPt = new TH1D("hPt","hDpt",100,0,100);
  for (int k=0; k<tree->GetEntries(); k++) {
    tree->GetEntry(k);
    if (brJet->fEta < -fJetEta || brJet->fEta > fJetEta) continue;
    if(!isPrompt){
      if(brD->fPartonType != 5) continue;
    }
    else if(brD->fPartonType != 4) continue;
    if(brD->fAncestorPDG == 2212) continue; // check if not coming from proton

    hPt->Fill(brD->fPt);
  }

  hPt->Scale(scaling);

  TH1D *hh_tmp = (TH1D*)hPt->Clone("hPt_reb");
  TH1D *hPt_reb = (TH1D*)hh_tmp->Rebin(fPtbinsDN,"hPt_reb",fPtbinsDA);

  return hPt_reb;

}

TH3D* getInvMass(TString data, Bool_t isBkg = 0) {

  const int ND = 10;

  TString histName;
  if(fDmesonSpecie) histName = "histosDStarMBN";
  else histName = "histosD0UpgradeN";
  TDirectoryFile* dir;
  TList *histList;
  THnSparseF *sparse;
  TH3D *hInvMassptD;

  if(fIsHBins && !isBkg){
    TFile *File = new TFile(data,"read");
  }
  else {
    TFile *File = new TFile(data,"read");;
    if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return ;}
    dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");

    for(int i=0;i<ND; i++){
        histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
        sparse = (THnSparseF*)histList->FindObject("hsDphiz");
        sparse->GetAxis(0)->SetRangeUser(fZmin,fZmax);
        sparse->GetAxis(1)->SetRangeUser(fJetptmin,fJetptmax);
        if(i==0) hInvMassptD=(TH3D*)sparse->Projection(3,1,2);
        else hInvMassptD->Add((TH3D*)sparse->Projection(3,1,2));
    }
  }

  if(fDmesonSpecie) hInvMassptD->GetXaxis()->SetTitle(Form("m(%s)(GeV/c^{2})","K#pi#pi-K#pi"));
  else hInvMassptD->GetXaxis()->SetTitle(Form("m(%s)(GeV/c^{2})","K#pi"));
  hInvMassptD->GetYaxis()->SetTitle("p_{T}^{jet}");
  hInvMassptD->GetZaxis()->SetTitle("p_{T}^{D}");

return hInvMassptD;

}


TH1* getEfficiency (TString file) {

  const int NDMC = 15;
  bool recoPt = 0;

  TString histName;
  if(fDmesonSpecie) histName = "histosDStarMBN";
  else histName = "histosD0UpgradeN";

  TH1D *hMCpt, *hMCpt_reco;
  if(fIsHBins){
    TFile *File = new TFile(effFile,"read");
  }
  else {
    TString effFile = file + fSMBBin;
    TFile *File = new TFile(effFile,"read");
    TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
  	for(int i=0; i<NDMC; i++){
      TList *histList =   (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
      THnSparseF *sparseMC = (THnSparseF*)histList->FindObject("ResponseMatrix");
      sparseMC->GetAxis(6)->SetRangeUser(fJetptmin,fJetptmax); // jet pT gen
      TH1D *hMC = (TH1D*)sparseMC->Projection(7); // Dpt gen
      //hMC->SetName(Form("hMC_%d",i));
  		THnSparseF *sparsereco = (THnSparseF*)histList->FindObject("ResponseMatrix");
      if(recoPt) {
        sparsereco->GetAxis(1)->SetRangeUser(fJetptmin,fJetptmax); // jet pT reco
      }
      else {
        sparsereco->GetAxis(6)->SetRangeUser(fJetptmin,fJetptmax); // jet pT gen
    		sparsereco->GetAxis(1)->SetRangeUser(0,100); // jet pT reco
      }
      TH1D *hreco= (TH1D*)sparsereco->Projection(7); // Dpt gen
      hreco->SetName(Form("hreco_%d",i));
  		if (!i){
  			hMCpt = (TH1D*)hMC->Clone("hMCpt");
  			hMCpt_reco = (TH1D*)hreco->Clone("hMCpt_reco");
        hMCpt->Sumw2();
        hMCpt_reco->Sumw2();
  		}
  		else {
  			hMCpt->Add(hMC);
  			hMCpt_reco->Add(hreco);
  		}
  	}
  }

	TH1D* hEff = (TH1D*)hMCpt_reco->Clone("hEff");
	hEff->Divide(hMCpt_reco,hMCpt,1,1,"b");
	TH1D* hpt_mc_reb = (TH1D*)hMCpt->Rebin(fPtbinsDN,"hpt_mc_reb",fPtbinsDA);
	TH1D* hpt_reco_reb = (TH1D*)hMCpt_reco->Rebin(fPtbinsDN,"hpt_reco_reb",fPtbinsDA);

	TH1D* hEff_reb = (TH1D*)hpt_reco_reb->Clone("hEff_reb");
	hEff_reb->Divide(hpt_reco_reb,hpt_mc_reb,1,1,"b");
	//hEff_reb->GetXaxis()->SetRangeUser(ptmin,ptmax);
	//hEff_reb->SetTitle(Form("|#eta_{jet}|<%.1f",0.9-fRpar));
	hEff_reb->GetXaxis()->SetTitle(Form("p_{T,%s} (GeV/$it{c})",fDmesonS.Data()));
	hEff_reb->GetYaxis()->SetTitle("Acc #times Eff");


/*	TFile *outFile = new TFile(Form("%s/DjetEff_%s_jetpt%d_%d.root", outDir.Data(), isPrompt ? "prompt" : "nonPrompt", (int)fJetptmin, (int)fJetptmax ),"RECREATE");
	hEff->Write();
	hEff_reb->Write();
	hMCpt->Write();
	hMCpt_reco->Write();
	hpt_mc_reb->Write();
	hpt_reco_reb->Write();
	outFile->Close();
*/

return hEff_reb;

}

TH1* getMCJetPt(){



}

double GetContent(TH1 *hh, double Dpt){
    return hh->GetBinContent(hh->GetXaxis()->FindBin(Dpt));
}


void setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Size_t size = 0.9, Width_t width=2){

    h->SetMarkerStyle(Mstyle);
    h->SetMarkerColor(color);
    h->SetMarkerSize(size);
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetTitle(0);
    h->GetXaxis()->SetTitle("p_{T,D^{0}}(GeV/c)");

    return;

}

void SaveCanvas(TCanvas *c, TString name = "tmp"){

    c->SaveAs(Form("%s.png",name.Data()));
    c->SaveAs(Form("%s.pdf",name.Data()));

}
