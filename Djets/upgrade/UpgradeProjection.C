//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------

#include "config.h"

using namespace std;

int fraabin = 0;

void Init() {


  fOutFile = new TFile(Form("%s/D0jets.root",OUTDIRECTORY.Data()),"RECREATE");
  fIsHBins = true;
  fIsHBEffWeighting = true; //this is for efficiencies
  fHBrejection = 1.;
  findex = 2; // reject the 2 HB for the eff calculation (introduces flucutations)
  
  fBRaaSys = 0;

  fJetptmin = 0;
  fJetptmax = 100;

  double sigBkg[] = { 5.824E-2, 1.5562E-1, 3.0889E-1, 5.0915E-1, 1.0686, 2.4596, 5, 8, 9, 10, 10, 10, 10, 10,10,10,10,10,10,10,10,10,10,10 };
  //double sigBkg[] = { 5.824E-2, 1.5562E-1, 3.0889E-1, 5.0915E-1, 1.0686, 2.4596, 6.3135, 6.3135, 9.10, 10, 10, 10, 10, 10 };
  fSB = new double[fPtbinsDN];
  for(int i=0;i<fPtbinsDN;i++) fSB[i] = sigBkg[i];


	double raaC[] = { 0.936, 0.847, 0.575, 0.389, 0.256, 0.250, 0.208, 0.211, 0.228, 0.270, 0.339, 0.397 };
	//double raaB[] = { 1.0557873253440786, 1.0929844575854308, 1.042076159575767, 0.9115762787749808, 0.7845157313889322, 0.6111070408940522, 0.5133346125444087, 0.46539617748472195, 0.44675167605723876, 0.4713357870927599, 0.5151051569650134, 0.5898517859530994 };
	double raaB[] = { 1.05579, 1.09298, 0.80102, 0.50297, 0.40265, 0.33480, 0.30200, 0.30200, 0.30163, 0.32003, 0.37617, 0.397 };
	double raaBSys[] = { 1.05579, 1.09298, 0.92563, 0.56527, 0.45412, 0.36730, 0.32909, 0.32909, 0.32330, 0.34711, 0.40867, 0.4295 };

	int raa = 2;
  fRaaC = new double[fPtbinsDN];
  fRaaB = new double[fPtbinsDN];
  fRaaBSys = new double[fPtbinsDN];
  for(int i=0;i<fPtbinsDN;i++){
    fRaaC[i] = raaC[i+raa];
    fRaaB[i] = raaB[i+raa];
    fRaaBSys[i] = raaBSys[i+raa];
  }

  double weightC[5] = { 6.67666666666666675E+00,4.26666666666666694E-01,3.57000000000000026E-02,3.90000000000000025E-03,6.32666666666666685E-04 };
  double weightB[5] = { 6.68333333333333357E+00,4.26666666666666694E-01,3.56999999999999956E-02,3.89666666666666677E-03,6.27000000000000062E-04 };
  for(int i=0; i<5; i++){
    fHBWeightsC[i] = weightC[i];
    fHBWeightsB[i] = weightB[i];
  }

}

void UpgradeProjection() {

    gStyle->SetOptStat(0000);
	gSystem->Exec(Form("mkdir -p %s/plots",OUTDIRECTORY.Data()));

    Init();
    
   // getRM(fFileEff_prompt);
   // return;
    
	getRaa();
	getSB();
	getEfficiencies();
	getPowhegDmesonSpectra();
	getPowhegDJetSpectra();
	 
	fHJetPtSpectrum = (TH1D*)extractSignal(); // rebined spectrum
    // subtract B FD contribution
	//subtractFD();
	FDsys_POWHEG();
	
	//unfold();
	
	//unfoldBkgFluc();
	
  if(fIsHBins) getMCJetPt(fFileSignal_prompt);

  
	fOutFile->Close();
  //runAllEff();
  //return;


return;
}

TH1* extractSignal() {


  

  TCanvas *cmass = new TCanvas("cmass","cmass",1800,1800);
  cmass->Divide(4,3);
  TCanvas *cmassbkg = new TCanvas("cmassbkg","cmassbkg",1800,1800);
  cmassbkg->Divide(4,3);
  TCanvas *cSB = new TCanvas("cSB","cSB",800,600);
  TCanvas *cJetPt = new TCanvas("cJetPt","cJetPt",1200,800);
  TCanvas *cJetPtReb = new TCanvas("cJetPtReb","cJetPtReb",1200,800);
  TCanvas *cJetPtRebUnc = new TCanvas("cJetPtRebUnc","cJetPtRebUnc",1200,600);
  TCanvas *cJetPtRebUncLog = new TCanvas("cJetPtRebUncLog","cJetPtRebUncLog",1200,600);

  const int ncanvas = 8;
  TCanvas *cpt[ncanvas];
  for(int i=0;i<ncanvas;i++) {
    cpt[i] = new TCanvas(Form("cpt_%d",i),Form("cpt_%d",i),1800,1800);
    cpt[i]->Divide(4,3);
  }
  const int nleg = 4;
  TLegend *leg[nleg];
  for(int i=0;i<nleg;i++){
    leg[i] = new TLegend(0.4,0.6,0.89,0.89);
    leg[i]->SetBorderSize(0);
    leg[i]->SetTextSize(0.055);
    leg[i]->SetLineColor(1);
    leg[i]->SetLineStyle(1);
    leg[i]->SetLineWidth(1);
    leg[i]->SetFillColor(0);
    leg[i]->SetFillStyle(1001);
  }


  TH1D *hJetPt, *hJetPtRebin;
  TH1D *hmean = new TH1D("hmean","hmean",fPtbinsDN,fPtbinsDA);
  TH1D *hsigma = new TH1D("hsigma","hsigma",fPtbinsDN,fPtbinsDA);
  TH1D *hSB = new TH1D("hSB","hSB",fPtbinsDN,fPtbinsDA);
  TH1D *h_tmp = new TH1D("h_tmp","h_tmp",fPtbinsDN,fPtbinsDA);
  h_tmp->SetLineColor(0);
  h_tmp->SetMarkerColor(0);
  h_tmp->SetMarkerSize(0);

  setHistoDetails(hSB,kGreen+2,21,0.9);

  TString signalFile = fFileSignal_prompt;
  if(!fIsHBins) signalFile += fSMBBin;
  TString signalFDFile = fFileSignal_nonprompt;
  if(!fIsHBins) signalFDFile += fSMBBin;
  TString signalBkgFile = fFileSignalBkg_prompt;
  if(!fIsHBins) signalBkgFile += fSMBBin;

  double bkgscale = getBkgScaling();

  for(int i=0; i<fPtbinsDN; i++){

      double pt = (fPtbinsDA[i]+fPtbinsDA[i+1])/2.;
      double effP = GetContent(fHEff_prompt,pt);
      double effNP = GetContent(fHEff_nprompt,pt);
      double sigP = GetContent(fHPowhegDPrompt,pt) * fSimScalingC * fRaaC[i+fraabin] * effP * fDataEv;
      double sigNP = GetContent(fHPowhegDNonPrompt,pt) * fSimScalingB * fRaaB[i+fraabin] * effNP * fDataEv;
	  double signalScale = sigP+sigNP;
	  //signalScale = sigP;
      // inv. mass
      TH1D *hhsignal = (TH1D*)getInvMass(signalFile,fPtbinsDA[i],fPtbinsDA[i+1],1);
      hhsignal->Rebin(fRebinMass);
      hhsignal->GetXaxis()->SetRangeUser(1.71,2.1);
      hhsignal->Scale(1./hhsignal->Integral());
      hhsignal->Scale(sigP);
      
      setHistoDetails(hhsignal,kRed+2,20,0.6);
      hhsignal->SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));
      for(int k=0; k<hhsignal->GetNbinsX();k++){
        hhsignal->SetBinError(k+1,TMath::Sqrt(hhsignal->GetBinContent(k+1)));
      }
      hhsignal->SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));

      // gaussian fit, just a x-check
      TF1 *fitG = new TF1("fitG","gaus",1.78,1.92);
      hhsignal->Fit("fitG","0RM");
      fitG->SetLineColor(2);
      double sigma = fitG->GetParameter(2);
      double mean = fitG->GetParameter(1);
      cmass->cd(i+1);
      TVirtualPad *pad = (TVirtualPad*)cmass->GetPad(i+1);
      hhsignal->DrawCopy();
      //fitG->Draw("same");

      //---------------- extract D-jet pT spectra
      double massmin = 1.78; // mean - 3*sigma;
      double massmax = 1.92; // mean + 3*sigma;
      double sbmin1 = 1.71, sbmax1 = 1.77;
      double sbmin2 = 1.93, sbmax2 = 2.1;

      // get jet pT spectra in a given D pT range, then scale to the POWHEG prediction for the yields (prompt+non-prompt)
      TH1D *hhPt = (TH1D*) getJetPt(signalFile,fPtbinsDA[i],fPtbinsDA[i+1],massmin,massmax,1);
      setHistoDetails(hhPt,4,20,0.7);
      hhPt->GetXaxis()->SetRangeUser(0,100);
      hhPt->GetXaxis()->SetTitle("p_{T, ch. jet}");
      hhPt->SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));
      cpt[0]->cd(i+1);
      TVirtualPad *pad = (TVirtualPad*)cpt[0]->GetPad(i+1);
      pad->SetLogy();
      hhPt->DrawCopy();

	  hhPt->Scale(1./hhPt->Integral());
      hhPt->Scale(sigP);
      
	  TH1D *hhPt_FD = (TH1D*) getJetPt(signalFDFile,fPtbinsDA[i],fPtbinsDA[i+1],massmin,massmax,1);
      setHistoDetails(hhPt_FD,4,20,0.7);
      hhPt_FD->GetXaxis()->SetRangeUser(0,100);
      hhPt_FD->Scale(1./hhPt_FD->Integral());
      hhPt_FD->Scale(sigNP);
      
      hhPt->Add(hhPt_FD);
      // scale to the expected yields
      TH1D *hPt = (TH1D*)hhPt->Clone("hPt");
      hPt->Sumw2();
      
      //hPt->Scale(1./hPt->Integral());
      //hPt->Scale(signalScale);
      
      setHistoDetails(hPt,2,25,0.7);
      hPt->SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));
      hPt->GetYaxis()->SetTitle("dN");
      hPt->SetMinimum(1);
      for(int k=1; k<=hPt->GetNbinsX();k++){
        hPt->SetBinError(k,TMath::Sqrt(hPt->GetBinContent(k)));
      }

      double bkgLevel = getBackground(massmin, massmax,fPtbinsDA[i],fPtbinsDA[i+1],i,cmassbkg);
      double bkg = bkgLevel*bkgscale * fRaaC[i];
      double sb;
      if(pt>6) sb = fSBfun->Eval(pt);
      else sb = GetContent(fHSB,pt);
      //bkg = signalScale/fSB[i];
      bkg = signalScale/sb;
      // get bkg jet pT shape, scale to bkg under the signal
      TH1D *hBkgPt1 = (TH1D*)getJetPt(signalBkgFile,fPtbinsDA[i],fPtbinsDA[i+1],sbmin1,sbmax1,1);
      TH1D *hBkgPt2 = (TH1D*)getJetPt(signalBkgFile,fPtbinsDA[i],fPtbinsDA[i+1],sbmin2,sbmax2,1);
      TH1D *hBkgPt = (TH1D*)hBkgPt1->Clone("hBkgPt");
      hBkgPt->Sumw2();
      hBkgPt->Add(hBkgPt2);
      hBkgPt->Scale(bkg/hBkgPt->Integral());
      setHistoDetails(hBkgPt,8,25,0.7);
      hBkgPt->SetMinimum(1);
      hBkgPt->GetYaxis()->SetTitle("dN");
      hBkgPt->SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));
      hBkgPt->GetXaxis()->SetRangeUser(0,100);
      hBkgPt->GetXaxis()->SetTitle("p_{T, ch. jet}");
      hBkgPt->SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));

      for(int k=1; k<=hhPt->GetNbinsX();k++){
        hBkgPt->SetBinError(k,TMath::Sqrt(hBkgPt->GetBinContent(k+1)));
      }
      hSB->SetBinContent(i+1,signalScale/bkg);

      TH1D *hTotalPt = (TH1D*)hPt->Clone("hTotalPt");
      hTotalPt->Add(hBkgPt);
      setHistoDetails(hTotalPt,4,25,0.7);
      hTotalPt->SetMinimum(1);
      hTotalPt->GetYaxis()->SetTitle("dN");
      hTotalPt->SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));
      for(int k=1; k<=hTotalPt->GetNbinsX();k++){
        hTotalPt->SetBinError(k,TMath::Sqrt(hPt->GetBinContent(k)+hBkgPt->GetBinContent(k)));
      }

      cpt[1]->cd(i+1);
      TVirtualPad *pad = (TVirtualPad*)cpt[1]->GetPad(i+1);
      pad->SetLogy();
      hTotalPt->DrawCopy();
      hBkgPt->DrawCopy("same");
      hPt->DrawCopy("same");

      TH1D *hSignalPt = (TH1D*)hPt->Clone("hSignalPt");
      // increase the stat. uncertanties according to the bkg level under the signal (bkg subtracted jet pT distribution)
      for(int k=1; k<=hSignalPt->GetNbinsX();k++){
        double err1 = TMath::Sqrt(hPt->GetBinContent(k)+hBkgPt->GetBinContent(k));
        double err2 = TMath::Sqrt(hBkgPt->GetBinContent(k));
        hSignalPt->SetBinError(k,TMath::Sqrt(err1*err1+err2*err2));
      }
      cpt[2]->cd(i+1);
      TVirtualPad *pad = (TVirtualPad*)cpt[2]->GetPad(i+1);
      pad->SetLogy();
      hTotalPt->DrawCopy();
      hBkgPt->DrawCopy("same");
      hSignalPt->DrawCopy("same");

      if(!i) {
        leg[0]->AddEntry(h_tmp," ","p");
        leg[0]->AddEntry(hTotalPt,"signal+bkg","p");
        leg[0]->AddEntry(hBkgPt,"bkg","p");
        leg[0]->AddEntry(hSignalPt,"signal","p");
        leg[0]->Draw("same");
      }

      // scale by the prompt efficiency
      TH1D *hSignalPtEff = (TH1D*)hSignalPt->Clone("hSignalPtEff");
      hSignalPtEff->Sumw2();
      setHistoDetails(hSignalPtEff,2,20,0.7);
      hSignalPtEff->SetMinimum(1);
      hSignalPtEff->Scale(1./effP);
      cpt[3]->cd(i+1);
      TVirtualPad *pad = (TVirtualPad*)cpt[3]->GetPad(i+1);
      pad->SetLogy();
      hSignalPtEff->DrawCopy();

      if(!i) hJetPt = (TH1D*)hSignalPtEff->Clone("hJetPt");
      else hJetPt->Add(hSignalPtEff);

      // ===== rebined histograms =====

      TH1D *hPtRebin = (TH1D*)hPt->Rebin(fPtbinsJetMeasN,"hPtRebin",fPtbinsJetMeasA);
      TH1D *hBkgPtRebin = (TH1D*)hBkgPt->Rebin(fPtbinsJetMeasN,"hBkgPtRebin",fPtbinsJetMeasA);

      TH1D *hTotalPtRebin = (TH1D*)hPtRebin->Clone("hTotalPtRebin");
      hTotalPtRebin->Add(hBkgPtRebin);
      setHistoDetails(hTotalPtRebin,4,25,0.7);
      hTotalPtRebin->SetMinimum(1);
      hTotalPtRebin->GetYaxis()->SetTitle("dN");
      hTotalPtRebin->SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));
      for(int k=1; k<=hTotalPtRebin->GetNbinsX();k++){
        hTotalPtRebin->SetBinError(k,TMath::Sqrt(hPtRebin->GetBinContent(k)+hBkgPtRebin->GetBinContent(k)));
      }

      cpt[4]->cd(i+1);
      TVirtualPad *pad = (TVirtualPad*)cpt[4]->GetPad(i+1);
      pad->SetLogy();
      hTotalPtRebin->DrawCopy();
      hBkgPtRebin->DrawCopy("same");
      hPtRebin->DrawCopy("same");

      TH1D *hSignalPtRebin = (TH1D*)hPtRebin->Clone("hSignalPt");
      // increase the stat. uncertanties according to the bkg level under the signal (bkg subtracted jet pT distribution)
      for(int k=1; k<=hSignalPtRebin->GetNbinsX();k++){
        double err1 = TMath::Sqrt(hPtRebin->GetBinContent(k)+hBkgPtRebin->GetBinContent(k));
        double err2 = TMath::Sqrt(hBkgPtRebin->GetBinContent(k));
        hSignalPtRebin->SetBinError(k,TMath::Sqrt(err1*err1+err2*err2));
      }
      cpt[5]->cd(i+1);
      TVirtualPad *pad = (TVirtualPad*)cpt[5]->GetPad(i+1);
      pad->SetLogy();
      hTotalPtRebin->DrawCopy();
      hBkgPtRebin->DrawCopy("same");
      hSignalPtRebin->DrawCopy("same");

      if(!i) {
        leg[1]->AddEntry(h_tmp," ","p");
        leg[1]->AddEntry(hTotalPtRebin,"signal+bkg","p");
        leg[1]->AddEntry(hBkgPtRebin,"bkg","p");
        leg[1]->AddEntry(hSignalPtRebin,"signal","p");
        leg[1]->Draw("same");
      }

      TH1D *hTotalPtRebinEff = (TH1D*)hTotalPtRebin->Clone("hTotalPtRebinEff");
      hTotalPtRebinEff->Scale(1./effP);
      hTotalPtRebinEff->SetMinimum(1);
      TH1D *hBkgPtRebinEff = (TH1D*)hBkgPtRebin->Clone("hBkgPtRebinEff");
      hBkgPtRebinEff->Scale(1./effP);
      TH1D *hSignalPtRebinEff = (TH1D*)hSignalPtRebin->Clone("hSignalPtRebinEff");
      hSignalPtRebinEff->Scale(1./effP);
      cpt[6]->cd(i+1);
      TVirtualPad *pad = (TVirtualPad*)cpt[6]->GetPad(i+1);
      pad->SetLogy();
      hTotalPtRebinEff->DrawCopy();
      hBkgPtRebinEff->DrawCopy("same");
      hSignalPtRebinEff->DrawCopy("same");

      if(!i) {
        leg[2]->AddEntry(h_tmp,"eff scaled","p");
        leg[2]->AddEntry(hTotalPtRebin,"signal+bkg","p");
        leg[2]->AddEntry(hBkgPtRebin,"bkg","p");
        leg[2]->AddEntry(hSignalPtRebin,"signal","p");
        leg[2]->Draw("same");
      }

      // scale by the prompt efficiency
      TH1D *hSignalPtEffRebin = (TH1D*)hSignalPtRebin->Clone("hSignalPtEffRebin");
      hSignalPtEffRebin->Sumw2();
      hSignalPtEffRebin->Scale(1./effP);
      setHistoDetails(hSignalPtEffRebin,2,20,0.7);
      hSignalPtEffRebin->GetYaxis()->SetTitle("dN");
      hSignalPtEffRebin->SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));

      if(!i) hJetPtRebin = (TH1D*)hSignalPtEffRebin->Clone("hJetPtRebin");
      else hJetPtRebin->Add(hSignalPtEffRebin);

      cpt[7]->cd(i+1);
      TVirtualPad *pad = (TVirtualPad*)cpt[7]->GetPad(i+1);
      pad->SetLogy();
      hSignalPtEffRebin->DrawCopy();

      if(!i) {
        leg[3]->AddEntry(hSignalPtEffRebin,"eff scaled signal","p");
        leg[3]->Draw("same");
      }
  }

  cSB->cd();
  cSB->SetLogy();
  hSB->GetXaxis()->SetTitle("p_{T,D^{0}}");
  hSB->GetYaxis()->SetTitle("S/B");
  hSB->Draw();

  cJetPt->cd();
  cJetPt->SetLogy();
  hJetPt->SetTitle();
  hJetPt->GetYaxis()->SetTitle("dN");
  hJetPt->Draw();
  
  TH1D* hPromptSim = (TH1D*)fHPowhegPromptRaa->Clone("hPromptSim");
  hPromptSim->Sumw2();
  hPromptSim->Scale(fDataEv*fSimScalingC);
  setHistoDetails(hPromptSim,kGreen+1,26,1.2);
  //hPromptSim->Draw("same");

  setHistoDetails(hJetPtRebin,2,21,1.2);
  hJetPtRebin->SetTitle();
  TH1D *hJetPtRebinScale = (TH1D*)hJetPtRebin->Clone("hJetPtRebinScale");
  hJetPtRebinScale->Scale(1,"width");
  //hJetPtRebinScale->SetMinimum(1);
  hJetPtRebinScale->GetYaxis()->SetTitle("dN/dp_{T}");

	TH1D* hPromptSimScale = (TH1D*)hPromptSim->Rebin(fPtbinsJetMeasN,"hPromptSimScale",fPtbinsJetMeasA);
	hPromptSimScale->Scale(1,"width");

  cJetPtReb->cd();
  cJetPtReb->SetLogy();
  hJetPtRebinScale->Draw();
  //hPromptSimScale->DrawCopy("same");
  
   TH1D *hjetptspectrumRebUnc = (TH1D*)hJetPtRebin->Clone("hjetptspectrumRebUnc");
   hjetptspectrumRebUnc->GetYaxis()->SetTitle("Rel. unc.");

      for(int mm=1; mm<=hJetPtRebin->GetNbinsX();mm++){
                  double err;
                  if(hJetPtRebin->GetBinContent(mm)) err = hJetPtRebin->GetBinError(mm)/hJetPtRebin->GetBinContent(mm);
                  else err = 0;
                  hjetptspectrumRebUnc->SetBinContent(mm,err);
                  hjetptspectrumRebUnc->SetBinError(mm,0);
      }
      
      
      fHJetPtSpectrumUnc = (TH1D*)hjetptspectrumRebUnc->Clone("fHJetPtSpectrumUnc");
      
    hjetptspectrumRebUnc->SetMarkerColor(kBlue+2);
    hjetptspectrumRebUnc->SetLineColor(kBlue+2);
    cJetPtRebUnc->cd();
	hjetptspectrumRebUnc->DrawCopy();
	
    hjetptspectrumRebUnc->SetMaximum(1);
    hjetptspectrumRebUnc->SetMarkerColor(kBlue+2);
    hjetptspectrumRebUnc->SetLineColor(kBlue+2);
    cJetPtRebUncLog->cd();
	cJetPtRebUncLog->SetLogy();
	hjetptspectrumRebUnc->DrawCopy();
	

  SaveCanvas(cmass, Form("%s/invMass_signal",OUTDIRECTORYPLOTS.Data()));
  SaveCanvas(cmassbkg, Form("%s/invMass_bkg",OUTDIRECTORYPLOTS.Data()));
  SaveCanvas(cSB, Form("%s/SigOverBkg",OUTDIRECTORYPLOTS.Data()));
  SaveCanvas(cJetPt, Form("%s/JetSpectrum_eff",OUTDIRECTORYPLOTS.Data()));
  SaveCanvas(cJetPtReb, Form("%s/JetSpectrum",OUTDIRECTORYPLOTS.Data()));
  SaveCanvas(cJetPtRebUnc, Form("%s/JetSpectrumUnc",OUTDIRECTORYPLOTS.Data()));
  SaveCanvas(cJetPtRebUncLog, Form("%s/JetSpectrumUncLog",OUTDIRECTORYPLOTS.Data()));
  TString canvasname[] = { "POWHEGspectrum", "jetSpectra_DptBins0", "jetSpectra_DptBins", "jetSpectraSignal_DptBins","jetSpectra_DptBins0_reb","jetSpectra_DptBins_reb","jetSpectra_DptBins_reb_eff","jetSpectraSignal_DptBins_reb_eff"}
  for(int i=0;i<ncanvas;i++) {
    SaveCanvas(cpt[i], Form("%s/%s",OUTDIRECTORYPLOTS.Data(),canvasname[i].Data()));
  }


	fOutFile->cd();
	hJetPtRebin->Write();
	hJetPtRebinScale->Write();
	hjetptspectrumRebUnc->Write();
	
  return hJetPtRebin;
}


void subtractFD() {
	

	
	TH1D* hFD = (TH1D*)fHPowhegNonPromptEffScaledRaa->Rebin(fPtbinsJetMeasN,"hFD",fPtbinsJetMeasA);
	hFD->Sumw2();
	hFD->Scale(fDataEv*fSimScalingB);
	
	TH1D* hPromptSim = (TH1D*)fHPowhegNonPromptRaa->Rebin(fPtbinsJetMeasN,"hPromptSim",fPtbinsJetMeasA);
	hPromptSim->Sumw2();
	hPromptSim->Scale(fDataEv*fSimScalingC);
	TH1D *hPromptSimScale = (TH1D*)hPromptSim->Clone("hPromptSimScale");
	hPromptSimScale->Scale(1,"width");
	
	fHJetPtSpectrumPrompt = (TH1D*)fHJetPtSpectrum->Clone("fHJetPtSpectrumPrompt");
	fHJetPtSpectrumPrompt->Sumw2();
	fHJetPtSpectrumPrompt->Add(hFD,-1);
	
	setHistoDetails(fHJetPtSpectrumPrompt,kRed+1,20,1.2);
	
	TH1D *hFDScaled = (TH1D*)hFD->Clone("hFD");
	hFDScaled->Scale(1,"width");
	TH1D *hPt = (TH1D*)fHJetPtSpectrum->Clone("hPt");
	hPt->Scale(1,"width");
	TH1D *hPtPrompt = (TH1D*)fHJetPtSpectrumPrompt->Clone("hPtPrompt");
	hPtPrompt->Scale(1,"width");
	
	setHistoDetails(hFDScaled,kBlue+1,21,1.2);
	setHistoDetails(hPromptSimScale,kGreen,26,1.2);
	setHistoDetails(hPt,kGreen+2,25,1.2);
	
	TCanvas *cspectrum = new TCanvas("cspectrum","cspectrum",1200,800);
	cspectrum->SetLogy();
	hPt->DrawCopy();
	hPtPrompt->DrawCopy("same");
	hFDScaled->DrawCopy("same");
	hPromptSimScale->DrawCopy("same");
	
	SaveCanvas(cspectrum, Form("%s/jetSpectrumPrompt",OUTDIRECTORYPLOTS.Data()));

return;	
}

TH1* FDsys_POWHEG() {
	
  TLegend *leg2 = new TLegend(0.6,0.55,0.85,0.75,"POWHEG,eff-Raa scaled");
  leg2->SetTextSize(0.045);
  leg2->SetBorderSize(0);
  leg2->AddEntry(fHPowhegPromptEffScaledRaa,"prompt","p");
  leg2->AddEntry(fHPowhegNonPromptEffScaledRaa,"non-prompt","p");
  TCanvas *cPow_scaled = new TCanvas("cPow_scaled","cPow_scaled",1200,800);
  cPow_scaled->SetLogy();
  fHPowhegPromptEffScaledRaa->DrawCopy();
  fHPowhegNonPromptEffScaledRaa->DrawCopy("same");
  leg2->Draw("same");
  
  SaveCanvas(cPow_scaled, Form("%s/POWHEG_Djet_effRaaScaled",OUTDIRECTORYPLOTS.Data()));
  
  TH1D *hPowheg = (TH1D*)fHPowhegPromptEffScaledRaa->Clone("hPowheg");
  hPowheg->Sumw2();
  hPowheg->Add(fHPowhegNonPromptEffScaledRaa);
  TH1D *hFDratio = (TH1D*)fHPowhegNonPromptEffScaledRaa->Clone("hFDratio");
  hFDratio->Sumw2();
  hFDratio->Divide(hPowheg);
  
  hFDratio->GetYaxis()->SetRangeUser(0,0.5);
  hFDratio->GetYaxis()->SetTitle("POWHEG, FD fraction");
  setHistoDetails(hFDratio,kBlue-1,20,1.2);
  
  TCanvas *cRatio = new TCanvas("cRatio","cRatio",1200,800);
  hFDratio->DrawCopy();
  
  SaveCanvas(cRatio, Form("%s/POWHEG_Djet_effRaaScaled_ratio",OUTDIRECTORYPLOTS.Data()));
  
  //rebin
    
  TH1D *hPowhegFD_reb = (TH1D*)fHPowhegNonPromptEffScaledRaa->Rebin(fPtbinsJetMeasN,"hPowhegFD_reb",fPtbinsJetMeasA);
  hPowhegFD_reb->Sumw2();
  TH1D *hPowheg_reb = (TH1D*)fHPowhegPromptEffScaledRaa->Rebin(fPtbinsJetMeasN,"hPowheg_reb",fPtbinsJetMeasA);
  hPowheg_reb->Sumw2();
  hPowheg_reb->Add(hPowhegFD_reb);
  TH1D *hFDratio_reb = (TH1D*)hPowhegFD_reb->Clone("hFDratio_reb");
  hFDratio_reb->Sumw2();
  hFDratio_reb->Divide(hPowheg_reb);
  
  hFDratio_reb->GetYaxis()->SetRangeUser(0,0.5);
  hFDratio_reb->GetYaxis()->SetTitle("POWHEG, FD fraction");
  setHistoDetails(hFDratio_reb,kBlue+1,20,1.2);
  
  TCanvas *cRatio_reb = new TCanvas("cRatio_reb","cRatio_reb",1200,800);
  hFDratio_reb->DrawCopy("H");
  
  
   //systematics
  
  //fHPowhegNonPromptEffScaledRaa_powhegSys
  //fHPowhegNonPromptEffScaledRaa_RaaSys
  
  TH1D *hPowhegFD_reb_powhegSys = (TH1D*)fHPowhegNonPromptEffScaledRaa_powhegSys->Rebin(fPtbinsJetMeasN,"hPowhegFD_reb_powhegSys",fPtbinsJetMeasA);
  hPowhegFD_reb_powhegSys->Scale(1.4);
  hPowhegFD_reb_powhegSys->Sumw2();
  TH1D *hPowheg_reb_powhegSys = (TH1D*)fHPowhegPromptEffScaledRaa->Rebin(fPtbinsJetMeasN,"hPowheg_reb_powhegSys",fPtbinsJetMeasA);
  hPowheg_reb_powhegSys->Sumw2();
  hPowheg_reb_powhegSys->Add(hPowhegFD_reb_powhegSys);
  TH1D *hFDratio_reb_powhegSys = (TH1D*)hPowhegFD_reb_powhegSys->Clone("hFDratio_reb_powhegSys");
  hFDratio_reb_powhegSys->Sumw2();
  hFDratio_reb_powhegSys->Divide(hPowheg_reb_powhegSys);
  
  hFDratio_reb_powhegSys->GetYaxis()->SetRangeUser(0,0.5);
  hFDratio_reb_powhegSys->GetYaxis()->SetTitle("POWHEG, FD fraction");
  setHistoDetails(hFDratio_reb_powhegSys,kBlue-1,20,1.2,2,3);
  
  hFDratio_reb_powhegSys->DrawCopy("Hsame");
  
  
  //TH1D *hPowhegFD_reb_RaaSys = (TH1D*)fHPowhegNonPromptEffScaledRaa->Rebin(fPtbinsJetMeasN,"hPowhegFD_reb_RaaSys",fPtbinsJetMeasA);
  //hPowhegFD_reb_RaaSys->Scale(1./3.);
  TH1D *hPowhegFD_reb_RaaSys = (TH1D*)fHPowhegNonPromptEffScaledRaa_RaaSys->Rebin(fPtbinsJetMeasN,"hPowhegFD_reb_RaaSys",fPtbinsJetMeasA);
  hPowhegFD_reb_RaaSys->Sumw2();
  TH1D *hPowheg_reb_RaaSys = (TH1D*)fHPowhegPromptEffScaledRaa->Rebin(fPtbinsJetMeasN,"hPowheg_reb_RaaSys",fPtbinsJetMeasA);
  hPowheg_reb_RaaSys->Sumw2();
  hPowheg_reb_RaaSys->Add(hPowhegFD_reb_RaaSys);
  TH1D *hFDratio_reb_RaaSys = (TH1D*)hPowhegFD_reb_RaaSys->Clone("hFDratio_reb_RaaSys");
  hFDratio_reb_RaaSys->Sumw2();
  hFDratio_reb_RaaSys->Divide(hPowheg_reb_RaaSys);
  
  hFDratio_reb_RaaSys->GetYaxis()->SetRangeUser(0,0.5);
  hFDratio_reb_RaaSys->GetYaxis()->SetTitle("POWHEG, FD fraction");
  setHistoDetails(hFDratio_reb_RaaSys,kBlue-1,20,1.2,3,3);
  
  hFDratio_reb_RaaSys->DrawCopy("Hsame");
  
  TLegend *leg = new TLegend(0.5,0.7,0.7,0.85,"POWHEG,FD fraction");
  leg->SetTextSize(0.035);
  leg->SetBorderSize(0);
  leg->AddEntry(hFDratio_reb,"central","pl");
  leg->AddEntry(hFDratio_reb_powhegSys,"POWHEG syst","pl");
  leg->AddEntry(hFDratio_reb_RaaSys,"R_{AA}^{b}=R_{AA}^{c}","pl");
  leg->Draw("same");
  
  SaveCanvas(cRatio_reb, Form("%s/POWHEG_Djet_effRaaScaled_ratio_reb",OUTDIRECTORYPLOTS.Data()));
  
  
  // final systematics
  
  TH1D *hSys_powheg = (TH1D*)hFDratio_reb_powhegSys->Clone("hSys_powheg");
  TH1D *hSys_Raa = (TH1D*)hFDratio_reb_RaaSys->Clone("hSys_powheg");
  hSys_Raa->SetName("hSys_Raa");
  hSys_powheg->SetName("hSys_powheg");
  
  for(int i=1; i<fPtbinsJetMeasN+1; i++) {
	
	hSys_powheg->SetBinContent(i,TMath::Abs(hFDratio_reb->GetBinContent(i) - hFDratio_reb_powhegSys->GetBinContent(i)) );
	hSys_powheg->SetBinError(i,0);
	
	hSys_Raa->SetBinContent(i,TMath::Abs(hFDratio_reb->GetBinContent(i) - hFDratio_reb_RaaSys->GetBinContent(i)) );
	hSys_Raa->SetBinError(i,0);
  }
  
  hSys_powheg->Scale(100);
  hSys_Raa->Scale(100);
  
  hSys_powheg->GetYaxis()->SetRangeUser(0,16);
  hSys_powheg->GetYaxis()->SetTitle("B feed-down systematic (%)");
  hSys_powheg->GetXaxis()->SetTitle("p_{T,ch. jet} (GeV/c)");
  
  TCanvas *cSys = new TCanvas("cSys","cSys",1200,800);
  hSys_powheg->DrawCopy("H");
  hSys_Raa->DrawCopy("Hsame");
  
  TLegend *legSys = new TLegend(0.35,0.7,0.65,0.85,"POWHEG,FD sys spectrum");
  legSys->SetTextSize(0.035);
  legSys->SetBorderSize(0);
  legSys->AddEntry(hFDratio_reb,"central","pl");
  legSys->AddEntry(hFDratio_reb_powhegSys,"POWHEG syst","pl");
  legSys->AddEntry(hFDratio_reb_RaaSys,"R_{AA}^{b} unc","pl");
  legSys->Draw("same");
  
  
   SaveCanvas(cSys, Form("%s/POWHEG_Djet_sys",OUTDIRECTORYPLOTS.Data()));

	fOutFile->cd();
	hFDratio_reb->Write();
	hFDratio_reb_powhegSys->Write();
	hFDratio_reb_RaaSys->Write();
	hSys_powheg->Write();
	hSys_Raa->Write();
	
	
	return hFDratio_reb;
}

void getPowhegDmesonSpectra() {

  TH1D *hh_tmp = (TH1D*)getPowhegDSpectra(fFilePowhegPrompt,1);  // get D x-section from POWHEG
  fHPowhegDPrompt = (TH1D*)hh_tmp->Rebin(fPtbinsDN,"fHPowhegDPrompt",fPtbinsDA);
  if(!fHPowhegDPrompt) { cout << "\n!!!Prompt D POWHEG spectrum not extracted !!!" << endl; return; }
  setHistoDetails(fHPowhegDPrompt,2,20,1.2);
  hh_tmp = (TH1D*)getPowhegDSpectra(fFilePowhegNonPrompt,0);  // get non-prompt D x-section from POWHEG
  fHPowhegDNonPrompt = (TH1D*)hh_tmp->Rebin(fPtbinsDN,"fHPowhegDNonPrompt",fPtbinsDA);
  if(!fHPowhegDNonPrompt) { cout << "\n!!! Non-Prompt D POWHEG spectrum not extracted !!!" << endl; return; }
  setHistoDetails(fHPowhegDNonPrompt,4,25,1.2);

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
  

}

void getPowhegDJetSpectra(){
	
  fHPowhegPrompt = (TH1D*)getPowhegJetSpectra(fFilePowhegPrompt,1,0,0,0);
  fHPowhegPromptEffScaled = (TH1D*)getPowhegJetSpectra(fFilePowhegPrompt,1,1,0,0);
  fHPowhegPromptEffScaledRaa = (TH1D*)getPowhegJetSpectra(fFilePowhegPrompt,1,1,0,1);
  fHPowhegPromptRaa = (TH1D*)getPowhegJetSpectra(fFilePowhegPrompt,1,0,0,1);
  
  fHPowhegNonPrompt = (TH1D*)getPowhegJetSpectra(fFilePowhegNonPrompt,0,0,0,0);
  fHPowhegNonPromptEffScaled = (TH1D*)getPowhegJetSpectra(fFilePowhegNonPrompt,0,1,0,0);
  fHPowhegNonPromptEffScaledRaa = (TH1D*)getPowhegJetSpectra(fFilePowhegNonPrompt,0,1,0,1);
  fHPowhegNonPromptEffRatioScaledRaa = (TH1D*)getPowhegJetSpectra(fFilePowhegNonPrompt,0,0,1,1);  
  
  fHPowhegNonPromptEffScaledRaa_powhegSys = (TH1D*)getPowhegJetSpectra(fFilePowhegPrompt_sys,0,1,0,1);
  fHPowhegNonPromptEffScaledRaa_RaaSys = (TH1D*)getPowhegJetSpectra(fFilePowhegNonPrompt,0,1,0,1,1);
  
  
  setHistoDetails(fHPowhegPrompt,kRed+1,20,1.2);
  setHistoDetails(fHPowhegPromptEffScaledRaa,kRed+1,20,1.2);
  setHistoDetails(fHPowhegNonPrompt,kBlue+1,21,1.2);
  setHistoDetails(fHPowhegNonPromptEffScaled,kBlue+1,21,1.2);
  setHistoDetails(fHPowhegNonPromptEffScaledRaa,kBlue+1,21,1.2);
  
  fHPowhegPrompt->GetXaxis()->SetTitle("p_{T,ch jet}");
  fHPowhegPromptEffScaled->GetXaxis()->SetTitle("p_{T,ch jet}");
  fHPowhegPromptEffScaledRaa->GetXaxis()->SetTitle("p_{T,ch jet}");
  fHPowhegPromptRaa->GetXaxis()->SetTitle("p_{T,ch jet}");
  fHPowhegNonPrompt->GetXaxis()->SetTitle("p_{T,ch jet}");
  fHPowhegNonPromptEffScaled->GetXaxis()->SetTitle("p_{T,ch jet}");
  fHPowhegNonPromptEffScaledRaa->GetXaxis()->SetTitle("p_{T,ch jet}");
  fHPowhegNonPromptEffRatioScaledRaa->GetXaxis()->SetTitle("p_{T,ch jet}");
  
  
  TLegend *leg = new TLegend(0.6,0.55,0.85,0.75);
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->AddEntry(fHPowhegPrompt,"prompt","p");
  leg->AddEntry(fHPowhegNonPrompt,"non-prompt","p");
  //leg->AddEntry(fHPowhegNonPromptEffScaled,"non-prompt, eff. scaled","p");
  
  TCanvas *cPow = new TCanvas("cPow","cPow",1200,800);
  cPow->SetLogy();
  fHPowhegPrompt->DrawCopy();
  fHPowhegNonPrompt->DrawCopy("same");
  //fHPowhegNonPromptEffScaled->DrawCopy("same");
  //fHPowhegNonPromptEffScaledRaa->DrawCopy("same");
  leg->Draw("same");
  
  SaveCanvas(cPow, Form("%s/POWHEG_Djet",OUTDIRECTORYPLOTS.Data()));
  

  TH1D *hPowheg = (TH1D*)fHPowhegPrompt->Clone("hPowheg");
  hPowheg->Sumw2();
  hPowheg->Add(fHPowhegNonPrompt);
  TH1D *hFDratio = (TH1D*)fHPowhegNonPrompt->Clone("hFDratio");
  hFDratio->Sumw2();
  hFDratio->Divide(hPowheg);
  hFDratio->GetYaxis()->SetRangeUser(0,1);
  
  TCanvas *cRatio = new TCanvas("cRatio","cRatio",1200,800);
  hFDratio->DrawCopy();
  
  SaveCanvas(cRatio, Form("%s/POWHEG_Djet_FDratio",OUTDIRECTORYPLOTS.Data()));
  
 }

TH1* getInvMass(TString data, double ptmin, double ptmax, int isSignal=1){

  TString histName;
  if(fDmesonSpecie) histName = "histosDStarMBN";
  else histName = "histosD0UpgradeN";

  const int ND = 15;
  TH1D *hout;
  if(fIsHBins && isSignal) {
    for(int j=2; j<=fHardBinsN; j++) {
        TFile *File = new TFile(Form("%s_%d.root",data.Data(),j),"read");
        if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
        TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
        TH1D *hh;
        for(int i=0;i<ND; i++){
            histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
            sparse = (THnSparseF*)histList->FindObject("hsDphiz");
            sparse->GetAxis(0)->SetRangeUser(fZmin,fZmax);
            sparse->GetAxis(1)->SetRangeUser(fJetptmin,fHardBinsMax[j-2]*4);
            sparse->GetAxis(2)->SetRangeUser(ptmin,ptmax);    // D meson pT range
            if(i==0) hh=(TH1D*)sparse->Projection(6);
            else hh->Add((TH1D*)sparse->Projection(6));
        }
        TDirectoryFile *dir2 = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
        TList *histList2 = (TList*)dir2->Get("histosD0UpgradeN0MCrec");
        TH1F* hEvents = (TH1F*)histList2->FindObject("hstat");
        double nEvents = hEvents->GetBinContent(1); // number of analyzed events, 2: selected
        hh->Scale(1000000./nEvents);
        //hh->Scale(1./nEvents);
        hh->Scale(fHBWeightsC[j-2]);
        if(j==2) {
          hout = (TH1D*)hh->Clone("hout");
          hout->Sumw2();
        }
        else {
          hout->Add(hh);
        }
        //File->Close();
       delete hh;
    }

  }
  else {
    TFile *File = new TFile(Form("%s",data.Data()),"read");
    if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
    TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
    TH1D *hh;
    for(int i=0;i<ND; i++){
        histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
        sparse = (THnSparseF*)histList->FindObject("hsDphiz");
        sparse->GetAxis(0)->SetRangeUser(fZmin,fZmax);
        sparse->GetAxis(1)->SetRangeUser(fJetptmin,fJetptmax);
        sparse->GetAxis(2)->SetRangeUser(ptmin,ptmax);    // D meson pT range
        if(i==0) hh=(TH1D*)sparse->Projection(3);
        else hh->Add((TH1D*)sparse->Projection(3));
    }
    TDirectoryFile *dir2 = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
    TList *histList2 = (TList*)dir2->Get("histosD0UpgradeN0MCrec");
    TH1F* hEvents = (TH1F*)histList2->FindObject("hstat");
    double nEvents = hEvents->GetBinContent(1); // number of analyzed events, 2: selected
    //hh->Scale(1./nEvents);
    hout = (TH1D*)hh->Clone("hout");
    hout->Sumw2();
   // File->Close();
    delete hh;

  }

return hout;

}

TH1* getJetPt(TString data, double ptmin, double ptmax, double massmin, double massmax, int isPythiaSignal=1){

  TString histName;
  if(fDmesonSpecie) histName = "histosDStarMBN";
  else histName = "histosD0UpgradeN";

  const int ND = 15;
  TH1D *hout;
  if(fIsHBins && isPythiaSignal) {
    for(int j=2; j<=fHardBinsN; j++) {
        TFile *File = new TFile(Form("%s_%d.root",data.Data(),j),"read");
        if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
        TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
        TH1D *hh;
        for(int i=0;i<ND; i++){
            TList *histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
            THnSparseF *sparse = (THnSparseF*)histList->FindObject("hsDphiz");
            sparse->GetAxis(0)->SetRangeUser(fZmin,fZmax);
            sparse->GetAxis(2)->SetRangeUser(ptmin,ptmax);    // D meson pT range
            sparse->GetAxis(3)->SetRangeUser(massmin,massmax);
            if(i==0) hh=(TH1D*)sparse->Projection(1);
            else hh->Add((TH1D*)sparse->Projection(1));

            for(int k=1; k<hh->GetNbinsX()+1;k++) {
              if(hh->GetBinCenter(k)>fHardBinsMax[j-2]*4){
                if(hh->GetBinContent(k)){
                  hh->SetBinContent(k,0);
                  hh->SetBinError(k,0);
                }
              }
            }
        }
        TDirectoryFile *dir2 = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
        TList *histList2 = (TList*)dir2->Get("histosD0UpgradeN0MCrec");
        TH1F* hEvents = (TH1F*)histList2->FindObject("hstat");
        double nEvents = hEvents->GetBinContent(1); // number of analyzed events, 2: selected
      //  hh->Scale(1000000./nEvents); // this 1000.. scaling doesn't really matter
        hh->Scale(1./nEvents); // this 1000.. scaling doesn't really matter
        hh->Scale(fHBWeightsC[j-2]);
        if(j==2) {
          hout = (TH1D*)hh->Clone("hout");
          hout->Sumw2();
        }
        else {
          hout->Add(hh);
        }
        delete hh;
    }
  }
  else {
    TFile *File = new TFile(Form("%s",data.Data()),"read");
    if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
    TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
    TH1D *hh;
    for(int i=0;i<ND; i++){
        TList *histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
        THnSparseF *sparse = (THnSparseF*)histList->FindObject("hsDphiz");
        sparse->GetAxis(0)->SetRangeUser(fZmin,fZmax);
        sparse->GetAxis(1)->SetRangeUser(fJetptmin,fJetptmax);
        sparse->GetAxis(2)->SetRangeUser(ptmin,ptmax);    // D meson pT range
        sparse->GetAxis(3)->SetRangeUser(massmin,massmax);
        if(i==0) hh=(TH1D*)sparse->Projection(3);
        else hh->Add((TH1D*)sparse->Projection(3));
    }
    hout = (TH1D*)hh->Clone("hout");
    hout->Sumw2();
    delete hh;
  }

return hout;

}

void getEfficiencies(){

  TFile *ofile = new TFile(Form("%s/efficiencies_%s%s.root",OUTDIRECTORY.Data(), fIsHBins ? "HBins" : "MB", fOutNamePostfix.Data()),"RECREATE");

  fHEff_prompt = (TH1D*)getEfficiency(fFileEff_prompt,ofile);  // get prompt efficiency
  if(!fHEff_prompt) { cout << "\n!!! Prompt efficincy not extracted !!! " << endl; return; }
  setHistoDetails(fHEff_prompt,2,20,1.2);

  fHEff_nprompt = (TH1D*)getEfficiency(fFileEff_nprompt,ofile,0); // get non-prompt getEfficiency
  if(!fHEff_nprompt) { cout << "\n!!! Non-Prompt efficincy not extracted !!! " << endl; return; }
  setHistoDetails(fHEff_nprompt,4,25,1.2);

  fHEff_prompt->SetMaximum(1.2);
  fHEff_nprompt->SetMaximum(1.2);
  TCanvas* cEffReb = new TCanvas("cEffReb","cEffReb",1200,800);
  cEffReb->SetLogy();
  fHEff_prompt->Draw("ep");
  fHEff_nprompt->Draw("epsame");

	TLegend *leg = new TLegend(0.5,0.22,0.85,0.35);
    leg->SetTextSize(0.045);
    leg->SetBorderSize(0);
    leg->AddEntry(fHEff_prompt,"Prompt","p");
    leg->AddEntry(fHEff_nprompt,"Feed-down","p");
    leg->Draw("same");


  //SaveCanvas(cEffReb, Form("%s/efficiencies_%s%s",OUTDIRECTORYPLOTS.Data(), fIsHBins ? "HBins" : "MB",fOutNamePostfix.Data()));
  SaveCanvas(cEffReb, Form("%s/efficiencies_%s",OUTDIRECTORYPLOTS.Data(), fIsHBins ? "HBins" : "MB"));

  ofile->cd();
  fHEff_prompt->Write();
  fHEff_nprompt->Write();
  ofile->Close();

}

TH1* getEfficiency (TString file, TFile *ofile, bool isprompt=1) {

  const int NDMC = 15;
  bool recoPt = 0;

  TString histName;
  if(fDmesonSpecie) histName = "histosDStarMBN";
  else histName = "histosD0UpgradeN";

  TCanvas *cMC = new TCanvas();
  cMC->SetLogy();
  TCanvas *creco = new TCanvas();
  creco->SetLogy();

  TH1D *hMCpt, *hMCpt_reco;
  TH1D *hpt_mc_reb, *hpt_reco_reb;
	TH1D *hhMC, *hhreco;
  if(fIsHBins){
    for(int j=findex; j<=fHardBinsN; j++) {
    //for(int j=findex; j<=findex; j++) {
        TFile *File = new TFile(Form("%s_%d.root",file.Data(),j),"read");
        if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
        TDirectoryFile* dir = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
        TH1D *hMC, *hreco;
        // loop over D mesons
        for(int i=0;i<NDMC; i++){
            TList *histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
            THnSparseF *sparseMC = (THnSparseF*)histList->FindObject("ResponseMatrix");
            //sparseMC->GetAxis(6)->SetRangeUser(fJetptmin,fJetptmax); // jet pT gen
            sparseMC->GetAxis(6)->SetRangeUser(fJetptmin,fHardBinsMax[j-2]*fHBrejection); // jet pT gen
            hhMC = (TH1D*)sparseMC->Projection(7); // Dpt gen
            THnSparseF *sparsereco = (THnSparseF*)histList->FindObject("ResponseMatrix");
            if(recoPt) {
              sparsereco->GetAxis(1)->SetRangeUser(fJetptmin,fJetptmax); // jet pT reco
              sparsereco->GetAxis(1)->SetRangeUser(fJetptmin,fHardBinsMax[j-2]*fHBrejection); // jet pT reco
            }
            else {
              //sparsereco->GetAxis(6)->SetRangeUser(fJetptmin,fJetptmax); // jet pT gen
              sparsereco->GetAxis(6)->SetRangeUser(fJetptmin,fHardBinsMax[j-2]*fHBrejection); // jet pT gen
              sparsereco->GetAxis(1)->SetRangeUser(0,100); // jet pT reco
            }
            hhreco= (TH1D*)sparsereco->Projection(7); // Dpt gen
			
            if (!i){
              hMC = (TH1D*)hhMC->Clone("hMC");
              hreco = (TH1D*)hhreco->Clone("hreco");
              hMC->Sumw2();
              hreco->Sumw2();
            }
            else {
              hMC->Add(hhMC);
              hreco->Add(hhreco);
            }			
		
		
        }
        // end of loop over D mesons

        TH1D* hMC_reb = (TH1D*)hMC->Rebin(fPtbinsDN,"hMC_reb",fPtbinsDA);
        TH1D* hreco_reb = (TH1D*)hreco->Rebin(fPtbinsDN,"hreco_reb",fPtbinsDA);

        TDirectoryFile *dir2 = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
        TList *histList2 = (TList*)dir2->Get("histosD0UpgradeN0MCrec");
        TH1F* hEvents = (TH1F*)histList2->FindObject("hstat");
        double nEvents = hEvents->GetBinContent(2); // number of analyzed events, 2: selected

        if(fIsHBEffWeighting) {
          hMC->Scale(10000./nEvents);
          if(isprompt) hMC->Scale(fHBWeightsC[j-2]);
          else hMC->Scale(fHBWeightsB[j-2]);
          hreco->Scale(10000./nEvents);
          if(isprompt) hreco->Scale(fHBWeightsC[j-2]);
          else hreco->Scale(fHBWeightsB[j-2]);

          hMC_reb->Scale(10000./nEvents);
          if(isprompt) hMC_reb->Scale(fHBWeightsC[j-2]);
          else hMC_reb->Scale(fHBWeightsB[j-2]);
          hreco_reb->Scale(10000./nEvents);
          if(isprompt) hreco_reb->Scale(fHBWeightsC[j-2]);
          else hreco_reb->Scale(fHBWeightsB[j-2]);
        }
        if(j==findex) {
          hMCpt = (TH1D*)hMC->Clone("hMCpt");
          hMCpt->Sumw2();
          hMCpt_reco = (TH1D*)hreco->Clone("hMCpt_reco");
          hMCpt_reco->Sumw2();

          hpt_mc_reb = (TH1D*)hMC_reb->Clone("hMCpt");
          hpt_mc_reb->Sumw2();
          hpt_reco_reb = (TH1D*)hreco_reb->Clone("hMCpt_reco");
          hpt_reco_reb->Sumw2();
        }
        else {
          hMCpt->Add(hMC);
          hMCpt_reco->Add(hreco);

          hpt_mc_reb->Add(hMC_reb);
          hpt_reco_reb->Add(hreco_reb);
        }

        hMC->SetMarkerColor(fColors[j]);
        hMC->SetLineColor(fColors[j]);
        hreco->SetMarkerColor(fColors[j]);
        hreco->SetLineColor(fColors[j]);
        hMC->GetXaxis()->SetRangeUser(0,50);
        hreco->GetXaxis()->SetRangeUser(0,50);
        hMC->SetMinimum(1);
        hreco->SetMinimum(1);
        cMC->cd();
        if(j==findex) hMC->DrawCopy();
        else hMC->DrawCopy("same");
        creco->cd();
        if(j==findex) hreco->DrawCopy();
        else hreco->DrawCopy("same");
    }
      setHistoDetails(hMCpt,kRed+2,20,0.8);
      setHistoDetails(hMCpt_reco,kBlue+2,20,0.8);
      cMC->cd();
      hMCpt->DrawCopy("same");
      creco->cd();
      hMCpt_reco->DrawCopy("same");
  }
  else {
    TString effFile = file + fSMBBin;
    TFile *File = new TFile(effFile,"read");
    TDirectoryFile* dir = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
  	for(int i=0; i<NDMC; i++){
      TList *histList =   (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
      THnSparseF *sparseMC = (THnSparseF*)histList->FindObject("ResponseMatrix");
      sparseMC->GetAxis(6)->SetRangeUser(fJetptmin,fJetptmax); // jet pT gen
      TH1D *hMC = (TH1D*)sparseMC->Projection(7); // Dpt gen
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

      TH1D* hMC_reb = (TH1D*)hMC->Rebin(fPtbinsDN,"hMC_reb",fPtbinsDA);
      TH1D* hreco_reb = (TH1D*)hreco->Rebin(fPtbinsDN,"hreco_reb",fPtbinsDA);

  		if (!i){
  			hMCpt = (TH1D*)hMC->Clone("hMCpt");
  			hMCpt_reco = (TH1D*)hreco->Clone("hMCpt_reco");
        hMCpt->Sumw2();
        hMCpt_reco->Sumw2();

        hpt_mc_reb = (TH1D*)hMC_reb->Clone("hMCpt");
        hpt_mc_reb->Sumw2();
        hpt_reco_reb = (TH1D*)hreco_reb->Clone("hMCpt_reco");
        hpt_reco_reb->Sumw2();
  		}
  		else {
  			hMCpt->Add(hMC);
  			hMCpt_reco->Add(hreco);

        hpt_mc_reb->Add(hMC_reb);
        hpt_reco_reb->Add(hreco_reb);
  		}

  	}

  }
  setHistoDetails(hMCpt,kRed+2,20,0.8);
  setHistoDetails(hMCpt_reco,kBlue+2,20,0.8);
	TH1D* hEff = (TH1D*)hMCpt_reco->Clone("hEff");
	hEff->Divide(hMCpt_reco,hMCpt,1,1,"b");

  setHistoDetails(hpt_mc_reb,kRed+2,20,0.8);
  setHistoDetails(hpt_reco_reb,kBlue+2,20,0.8);
	TH1D* hEff_reb = (TH1D*)hpt_reco_reb->Clone("hEff_reb");
	hEff_reb->Divide(hpt_reco_reb,hpt_mc_reb,1,1,"b");
	//hEff_reb->GetXaxis()->SetRangeUser(ptmin,ptmax);
	//hEff_reb->SetTitle(Form("|#eta_{jet}|<%.1f",0.9-fRpar));
	hEff_reb->GetXaxis()->SetTitle(Form("p_{T,%s} (GeV/#it{c})",fDmesonS.Data()));
	hEff_reb->GetYaxis()->SetTitle("Acc #times Eff");

  TH1D *hMCpt_int = (TH1D*)hMCpt->Clone("hMCpt_int");
  hMCpt_int->Scale(1./hMCpt_int->Integral());
  TH1D *hMCpt_reco_int = (TH1D*)hMCpt_reco->Clone("hMCpt_reco_int");
  hMCpt_reco_int->Scale(1./hMCpt_reco_int->Integral());
  TH1D *hpt_mc_reb_int = (TH1D*)hpt_mc_reb->Clone("hpt_mc_reb_int");
  hpt_mc_reb_int->Scale(1./hpt_mc_reb_int->Integral());
  TH1D *hpt_reco_reb_int = (TH1D*)hpt_reco_reb->Clone("hpt_reco_reb_int");
  hpt_reco_reb_int->Scale(1./hpt_reco_reb_int->Integral());
  if(isprompt) {
    hMCpt->SetName("mcPt_prompt");
    hMCpt_reco->SetName("recoPt_prompt");
    hpt_mc_reb->SetName("mcPt_reb_prompt");
    hpt_reco_reb->SetName("recoPt_reb_prompt");
    hMCpt_int->SetName("mcPt_prompt_int");
    hMCpt_reco_int->SetName("recoPt_prompt_int");
    hpt_mc_reb_int->SetName("mcPt_reb_prompt_int");
    hpt_reco_reb_int->SetName("recoPt_reb_prompt_int");
    hEff->SetName("eff_prompt");
    hEff_reb->SetName("eff_reb_prompt");
  }
  else {
    hMCpt->SetName("mcPt_nonprompt");
    hMCpt_reco->SetName("recoPt_nonprompt");
    hpt_mc_reb->SetName("mcPt_reb_nonprompt");
    hpt_reco_reb->SetName("recoPt_reb_nonprompt");
    hMCpt_int->SetName("mcPt_nonprompt_int");
    hMCpt_reco_int->SetName("recoPt_nonprompt_int");
    hpt_mc_reb_int->SetName("mcPt_reb_nonprompt_int");
    hpt_reco_reb_int->SetName("recoPt_reb_nonprompt_int");
    hEff->SetName("eff_nonprompt");
    hEff_reb->SetName("eff_reb_nonprompt");
  }



  ofile->cd();
	hMCpt->Write();
	hMCpt_reco->Write();
	hpt_mc_reb->Write();
	hpt_reco_reb->Write();
  hMCpt_int->Write();
	hMCpt_reco_int->Write();
	hpt_mc_reb_int->Write();
	hpt_reco_reb_int->Write();

return hEff_reb;

}


void unfoldBkgFluc() {
	
	
	gSystem->Load("$HOME/Work/alice/RooUnfold-1.1.1/libRooUnfold.so");
	//int kreg = 10;
	
	TString histName = "histosD0UpgradeN";
	TFile *File = new TFile(Form("%s",fFileBkg.Data()),"read");
    if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
    TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
    TH1D *hh;
    for(int i=0;i<5; i++){
        TList *histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
        THnSparseF *sparse = (THnSparseF*)histList->FindObject("hsDphiz");
        //sparse->GetAxis(5)->SetRangeUser(fZmin,fZmax);
        //sparse->GetAxis(2)->SetRangeUser(ptmin,ptmax);    // D meson pT range
        if(i==0) hh=(TH1D*)sparse->Projection(1);
        else hh->Add((TH1D*)sparse->Projection(1));
      }
	
	TH1D*  hinput = (TH1D*) fHPowhegPrompt->Clone("hinput");
	TH1D* priorhisto = (TH1D*) fHPowhegPrompt->Clone("priorhisto");
	//priorhisto->Rebin(2);
	TH1D* priorhisto_reco = (TH1D*) hh->Clone("priorhisto_reco");
	
	TF1 *PriorFunction = new TF1("PriorFunction","[0]* pow(x,-[1]) * exp(-[1]*[2]/x)",3,100);
	priorhisto->Fit(PriorFunction, "MRN","",3,100);
	PriorFunction->SetLineColor(2);

	
	RooUnfoldResponse response; // (120, -20.0, 100.0);
	response.Setup (80, 0, 80, 100, -20, 80);
	RooUnfoldResponse response_trans; // (120, -20.0, 100.0);
	response_trans.Setup (100, -20, 80, 80, 0, 80);
	
	TH2D *hRM = new TH2D("hRM","hRM",120,-20,100,120,-20,100);
	
	Int_t sampling = 10000000;
	for (Int_t i= 0; i<sampling; i++) {
		//Double_t xt = gRandom->BreitWigner (0.3, 2.5);
		Double_t xt = PriorFunction->GetRandom();
		Double_t xsmear = gRandom->Gaus(0.664,7.6912);
		Double_t x = xt + xsmear;
		response_trans.Fill(x,xt);
		response.Fill(xt,x);
		hRM->Fill(xt,x);
		
	}
	hRM->Scale(1./(sampling));
	
	TH1D *hmeas = (TH1D*)hRM->ProjectionX();
	hmeas->SetMarkerColor(4);
	hmeas->SetLineColor(4);
	
	TH1D *htrue = (TH1D*)hRM->ProjectionY();
	htrue->SetMarkerColor(2);
	htrue->SetLineColor(2);
	
	TCanvas *c = new TCanvas;
	c->SetLogy();
	//priorhisto->Draw();
	htrue->Draw();
	PriorFunction->Draw("same");
	hmeas->Draw("same");
	
	TCanvas *c2 = new TCanvas;
	c2->SetLogz();
	hRM->Draw("colz2");
	
	
	 RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
	 Int_t kReg = 5;
	 
	 //RooUnfoldBayes   unfold (&response, fHPowhegPrompt, kReg);
	 RooUnfoldSvd   unfold (&response, fHPowhegPrompt, kReg);
	 TH1D* hReco = (TH1D*) unfold.Hreco();
	 hReco->SetMarkerColor(kGreen+1);
	 hReco->SetLineColor(kGreen+1);
	 
	 TH1D *folded = (TH1D*)hReco->Clone("folded");
	 //RooUnfoldBayes   refold (&response_trans, folded, kReg);
	 RooUnfoldSvd   refold (&response_trans, folded, kReg);
	 TH1D* hRefolded = (TH1D*) refold.Hreco();
	 hRefolded->SetMarkerColor(kBlue+1);
	 hRefolded->SetLineColor(kBlue+1);
	 
	
	 
	TH2D *fPearsonCoeffs = getPearsonCoeffs( unfold.Ereco(RooUnfold::kCovariance) );
	fPearsonCoeffs->SetName(Form("PearsonCoeffs%d",kReg));
	fPearsonCoeffs->GetZaxis()->SetRangeUser(-1,1);
	
	TH2D *fPearsonCoeffs_ref = getPearsonCoeffs( refold.Ereco(RooUnfold::kCovariance) );
	fPearsonCoeffs_ref->SetName(Form("PearsonCoeffs%d",kReg));
	fPearsonCoeffs_ref->GetZaxis()->SetRangeUser(-1,1);
	
	
	
	
	TH1D *hUncTrue = (TH1D*)fHPowhegPrompt->Clone("hUncTrue");
	hUncTrue->GetYaxis()->SetTitle("Rel. unc.");

      for(int mm=1; mm<=hUncTrue->GetNbinsX();mm++){
                  double err;
                  if(fHPowhegPrompt->GetBinContent(mm)) err = fHPowhegPrompt->GetBinError(mm)/fHPowhegPrompt->GetBinContent(mm);
                  else err = 0;
                  hUncTrue->SetBinContent(mm,err);
                  hUncTrue->SetBinError(mm,0);
      }
      
      TH1D *hUncRef = (TH1D*)hRefolded->Clone("hUncTrue");
	  hUncRef->GetYaxis()->SetTitle("Rel. unc.");

      for(int mm=1; mm<=hUncRef->GetNbinsX();mm++){
                  double err;
                  if(hRefolded->GetBinContent(mm)) err = hRefolded->GetBinError(mm)/hRefolded->GetBinContent(mm);
                  else err = 0;
                  hUncRef->SetBinContent(mm,err);
                  hUncRef->SetBinError(mm,0);
      }
      
      TH1D *hUncFold = (TH1D*)hReco->Clone("hUncFold");
	  hUncFold->GetYaxis()->SetTitle("Rel. unc.");

      for(int mm=1; mm<=hUncFold->GetNbinsX();mm++){
                  double err;
                  if(hReco->GetBinContent(mm)) err = hReco->GetBinError(mm)/hReco->GetBinContent(mm);
                  else err = 0;
                  hUncFold->SetBinContent(mm,err);
                  hUncFold->SetBinError(mm,0);
      }
      
     
       hReco->GetXaxis()->SetTitle("p_{T, ch. jet}");
	 hReco->SetMaximum(fHPowhegPrompt->GetMaximum()*2);
	 hinput->SetMarkerSize(0.8);
	 
	TLegend *leg = new TLegend(0.2,0.3,0.4,0.5);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->AddEntry(hinput, "input","l");
	leg->AddEntry(hReco, "folded","l");
	leg->AddEntry(hRefolded, "re-folded","l");
	 
       TCanvas *cun = new TCanvas();
	 cun->SetLogy();
	
	 hReco->Draw();
	 hinput->Draw("same");
	 hRefolded->Draw("same");
	 leg->Draw("same");
	 
	 TCanvas *c = new TCanvas;
	fPearsonCoeffs->Draw("colz2");
	
	TCanvas *cref = new TCanvas;
	fPearsonCoeffs_ref->Draw("colz2");
	
		hUncTrue->GetXaxis()->SetRangeUser(0,80);
		hUncTrue->GetXaxis()->SetTitle("p_{T, ch. jet}");
      TCanvas *cunc = new TCanvas;
      hUncTrue->Draw();
      hUncRef->Draw("same");
      hUncFold->Draw("same");
      leg->Draw("same");
      
      
      SaveCanvas(cun, Form("%s/bkgUnfoldingSVD_%d",OUTDIRECTORYPLOTS.Data(),kReg));
      SaveCanvas(cunc, Form("%s/bkgUnfoldingSVD_unc_%d",OUTDIRECTORYPLOTS.Data(),kReg));

	
	/*TString filename = "BkgFluctuationMtx.root";
	TFile *File = new TFile(filename,"read");

	TH2D *fMatrixDeltaPt = (TH2D*)File->Get("hBkgM");
	TH2D *fMatrixDeltaPt_trans = (TH2D*)File->Get("hBkgM_trans");
	if (!fMatrixDeltaPt) {
		
		return;
	}
	
	TString histName = "histosD0UpgradeN";
	TFile *File = new TFile(Form("%s",fFileBkg.Data()),"read");
    if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
    TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
    TH1D *hh;
    for(int i=0;i<5; i++){
        TList *histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
        THnSparseF *sparse = (THnSparseF*)histList->FindObject("hsDphiz");
        //sparse->GetAxis(5)->SetRangeUser(fZmin,fZmax);
        //sparse->GetAxis(2)->SetRangeUser(ptmin,ptmax);    // D meson pT range
        if(i==0) hh=(TH1D*)sparse->Projection(1);
        else hh->Add((TH1D*)sparse->Projection(1));
      }
	
	TH1D* priorhisto = (TH1D*) fHPowhegPrompt->Clone("priorhisto");
	TH1D* priorhisto_reco = (TH1D*) hh->Clone("priorhisto_reco");
	
	 TH1D* hNormY = (TH1D*)fMatrixDeltaPt->ProjectionY("hNormY");
	 hNormY->Divide(priorhisto);
	 //WeightMatrixY(fMatrixDeltaPt,hNormY,1);
	 
	 TH1D* hNormY_trans = (TH1D*)fMatrixDeltaPt_trans->ProjectionY("hNormY_trans");
	 hNormY_trans->Divide(priorhisto_reco);
	 //WeightMatrixY(fMatrixDeltaPt_trans,hNormY_trans,1);

	
	const int       PtbinsJetTrueN = 15;
    double          PtbinsJetTrueA[PtbinsJetTrueN+1] = { 3,5,6,8,10,12,14,16,20,25,30,40,50,60,80,100 };
    const int       PtbinsJetMeasN = 17;
    double          PtbinsJetMeasA[PtbinsJetMeasN+1] = { -20,-10,-5,0,3,5,6,8,10,12,14,16,20,25,30,40,50,60 };
    
    
	TH2D *fBkgMatrix = Rebin2D("fBkgMatrix", fMatrixDeltaPt, PtbinsJetMeasN, PtbinsJetMeasA, PtbinsJetTrueN, PtbinsJetTrueA,0);
	TH2D *fBkgMatrix_trans = Rebin2D("fBkgMatrix_trans", fMatrixDeltaPt_trans, PtbinsJetTrueN, PtbinsJetTrueA, PtbinsJetMeasN, PtbinsJetMeasA,0);
	
	fBkgMatrix->Scale(1.,"width");
	fBkgMatrix_trans->Scale(1.,"width");
	
	TH1D* hProjYeff = (TH1D*)fBkgMatrix->ProjectionY("hProjYeff");
    TH1D* hProjXeff = (TH1D*)fBkgMatrix->ProjectionX("hProjXeff");
    
    TH1D* hProjYeff_trans = (TH1D*)fBkgMatrix_trans->ProjectionY("hProjYeff_trans");
    TH1D* hProjXeff_trans = (TH1D*)fBkgMatrix_trans->ProjectionX("hProjXeff_trans");
    
    
    TH1D* fRawRebin = (TH1D*)fHJetPtSpectrum->Clone("fRawRebin");
    
   
	RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    //RooUnfoldResponse response_trans(hProjYeff_trans,hProjXeff_trans, fBkgMatrix_trans, "response_trans","response_trans");
    RooUnfoldResponse response_trans(0,0, fBkgMatrix_trans, "response_trans","response_trans");
    response_trans.UseOverflow(1);
    
   
    RooUnfoldBayes unfold_trans (&response_trans, fRawRebin, kreg);
	TH1D *fUnfoldedBayes_trans = (TH1D*)unfold_trans.Hreco();
	TH1D *folded = (TH1D*)fUnfoldedBayes_trans->Clone("folded");
	
    //RooUnfoldResponse response(hProjYeff,hProjXeff, fBkgMatrix, "response","response");
    RooUnfoldResponse response(0,0, fBkgMatrix, "response","response");
    response.UseOverflow(1);
    
   
    RooUnfoldBayes unfold (&response, folded, kreg);
	TH1D *fUnfoldedBayes = (TH1D*)unfold.Hreco();
	//TH1D *folded = (TH1D*)response.ApplyToTruth(fUnfoldedBayes);
	
	*/
	
	
	
	return;
	// ------------ Get Pearson coefficient ------------
	TH2D *fPearsonCoeffs = getPearsonCoeffs( unfold.Ereco(RooUnfold::kCovariance) );
	fPearsonCoeffs->SetName(Form("PearsonCoeffs%d",kreg));
	fPearsonCoeffs->GetZaxis()->SetRangeUser(-1,1);
	
	fUnfoldedBayes->GetXaxis()->SetTitle("p_{T,ch. jet} (GeV/#it{c})");
	fUnfoldedBayes->GetYaxis()->SetTitle("dN/dp_{T}");
	//fUnfoldedBayes->SetTitle("Bayes unfolding");
	fUnfoldedBayes->SetTitle();
	
	setHistoDetails(fUnfoldedBayes,kRed+2,20,0.8);
	setHistoDetails(fRawRebin,kBlue+2,20,0.8);
	setHistoDetails(folded,kGreen+1,20,0.8);

	TH1D *hjetptspectrumRebUnc = (TH1D*)fUnfoldedBayes->Clone("hjetptspectrumRebUnc");
	hjetptspectrumRebUnc->GetYaxis()->SetTitle("Rel. unc.");

      for(int mm=1; mm<=fUnfoldedBayes->GetNbinsX();mm++){
                  double err;
                  if(fUnfoldedBayes->GetBinContent(mm)) err = fUnfoldedBayes->GetBinError(mm)/fUnfoldedBayes->GetBinContent(mm);
                  else err = 0;
                  hjetptspectrumRebUnc->SetBinContent(mm,err);
                  hjetptspectrumRebUnc->SetBinError(mm,0);
      }
      
    
    hjetptspectrumRebUnc->SetMarkerColor(kGreen+2);
    hjetptspectrumRebUnc->SetLineColor(kGreen+2);
    TCanvas *cJetPtRebUnc = new TCanvas("cJetPtRebUnc","cJetPtRebUnc",1200,600);
    cJetPtRebUnc->cd();
    fHJetPtSpectrumUnc->DrawCopy();
	hjetptspectrumRebUnc->DrawCopy("same");
	
	
	fHJetPtSpectrumUncUnf = (TH1D*)hjetptspectrumRebUnc->Clone("fHJetPtSpectrumUncUnf");


	TCanvas *cUnf = new TCanvas("cUnf","cUnf",1200,800);
	cUnf->SetLogy();
	fRawRebin->DrawCopy();
	fUnfoldedBayes->DrawCopy("same");
	folded->DrawCopy("same");
	
	TH1D *hUnfRatio = (TH1D*) fUnfoldedBayes->Clone("hUnfRatio");
	hUnfRatio->Sumw2();
	hUnfRatio->Divide(fRawRebin);
	
	TCanvas *cPear = new TCanvas("cPear","cPear",1200,800);
	fPearsonCoeffs->Draw("colz2");
	
	TCanvas *cRatio = new TCanvas("cRatio","cRatio",1200,600);
	hUnfRatio->Draw();


	TCanvas *cBkgMtx_trans = new TCanvas("cBkgMtx_trans","cBkgMtx_trans",1200,600);
	cBkgMtx_trans->SetLogz();
	fBkgMatrix_trans->Draw("colz2");
	
	TCanvas *cBkgMtx = new TCanvas("cBkgMtx","cBkgMtx",1200,600);
	cBkgMtx->SetLogz();
	fBkgMatrix->Draw("colz2");
	
	
	 SaveCanvas(cPear, Form("%s/PearsonCoef_bkgFluc",OUTDIRECTORYPLOTS.Data()));
	 SaveCanvas(cUnf, Form("%s/JetSpectrumUnfolded_bkgFluc",OUTDIRECTORYPLOTS.Data()));
	 SaveCanvas(cRatio, Form("%s/JetSpectrumUnfoldedRatio_bkgFluc",OUTDIRECTORYPLOTS.Data()));
	 
	 SaveCanvas(cJetPtRebUnc, Form("%s/JetSpectrumUncUnfolded_bkgFluc",OUTDIRECTORYPLOTS.Data()));
	 
	 SaveCanvas(cBkgMtx, Form("%s/RM_bkg",OUTDIRECTORYPLOTS.Data()));
	 SaveCanvas(cBkgMtx_trans, Form("%s/RM_bkgtrans",OUTDIRECTORYPLOTS.Data()));
	

	
}

void unfold() {
	
	gSystem->Load("$HOME/Work/alice/RooUnfold-1.1.1/libRooUnfold.so");
	
	int kreg = 5; 
	
	fHDettMatrix = (TH2D*)getRM(fFileEff_prompt);
	
	TH1D* hProjYeff = (TH1D*)fHDettMatrix->ProjectionY("hProjYeff");
    TH1D* hProjXeff = (TH1D*)fHDettMatrix->ProjectionX("hProjXeff");
    TH1D* fRawRebin = (TH1D*)fHJetPtSpectrum->Clone("fRawRebin");
    
	RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    RooUnfoldResponse response(hProjXeff,hProjYeff, fHDettMatrix, "response","response");
    //RooUnfoldResponse response(0,0, fHDettMatrix, "response","response");
    response.UseOverflow(1);
    
    RooUnfoldBayes unfold (&response, fRawRebin, kreg);
	TH1D *fUnfoldedBayes = (TH1D*)unfold.Hreco();
	TH1D *folded = (TH1D*)response.ApplyToTruth(fUnfoldedBayes);

	// ------------ Get Pearson coefficient ------------
	TH2D *fPearsonCoeffs = getPearsonCoeffs( unfold.Ereco(RooUnfold::kCovariance) );
	fPearsonCoeffs->SetName(Form("PearsonCoeffs%d",kreg));
	fPearsonCoeffs->GetZaxis()->SetRangeUser(-1,1);
	
	fUnfoldedBayes->GetXaxis()->SetTitle("p_{T,ch. jet} (GeV/#it{c})");
	fUnfoldedBayes->GetYaxis()->SetTitle("dN/dp_{T}");
	//fUnfoldedBayes->SetTitle("Bayes unfolding");
	fUnfoldedBayes->SetTitle();
	
	setHistoDetails(fUnfoldedBayes,kRed+2,20,0.8);
	setHistoDetails(fRawRebin,kBlue+2,20,0.8);
	setHistoDetails(folded,kGreen+1,20,0.8);

	TH1D *hjetptspectrumRebUnc = (TH1D*)fUnfoldedBayes->Clone("hjetptspectrumRebUnc");
	hjetptspectrumRebUnc->GetYaxis()->SetTitle("Rel. unc.");

      for(int mm=1; mm<=fUnfoldedBayes->GetNbinsX();mm++){
                  double err;
                  if(fUnfoldedBayes->GetBinContent(mm)) err = fUnfoldedBayes->GetBinError(mm)/fUnfoldedBayes->GetBinContent(mm);
                  else err = 0;
                  hjetptspectrumRebUnc->SetBinContent(mm,err);
                  hjetptspectrumRebUnc->SetBinError(mm,0);
      }
      
    
    hjetptspectrumRebUnc->SetMarkerColor(kGreen+2);
    hjetptspectrumRebUnc->SetLineColor(kGreen+2);
    TCanvas *cJetPtRebUnc = new TCanvas("cJetPtRebUnc","cJetPtRebUnc",1200,600);
    cJetPtRebUnc->cd();
    fHJetPtSpectrumUnc->DrawCopy();
	hjetptspectrumRebUnc->DrawCopy("same");
	
	
	fHJetPtSpectrumUncUnf = (TH1D*)hjetptspectrumRebUnc->Clone("fHJetPtSpectrumUncUnf");


	TCanvas *cUnf = new TCanvas("cUnf","cUnf",1200,800);
	cUnf->SetLogy();
	fRawRebin->DrawCopy();
	fUnfoldedBayes->DrawCopy("same");
	folded->DrawCopy("same");
	
	TH1D *hUnfRatio = (TH1D*) fUnfoldedBayes->Clone("hUnfRatio");
	hUnfRatio->Sumw2();
	hUnfRatio->Divide(fRawRebin);
	
	TCanvas *cPear = new TCanvas("cPear","cPear",1200,800);
	fPearsonCoeffs->Draw("colz2");
	
	TCanvas *cRatio = new TCanvas("cRatio","cRatio",1200,600);
	hUnfRatio->Draw();

	 SaveCanvas(cPear, Form("%s/PearsonCoef",OUTDIRECTORYPLOTS.Data()));
	 SaveCanvas(cUnf, Form("%s/JetSpectrumUnfolded",OUTDIRECTORYPLOTS.Data()));
	 SaveCanvas(cRatio, Form("%s/JetSpectrumUnfoldedRatio",OUTDIRECTORYPLOTS.Data()));
	 
	 SaveCanvas(cJetPtRebUnc, Form("%s/JetSpectrumUncUnfolded",OUTDIRECTORYPLOTS.Data()));
	
}

double getBackground(double massmin, double massmax, double ptmin, double ptmax, int iter, TCanvas *cmassbkg){

      double bkgscale = getBkgScaling();
      TH1D *hhbkg = (TH1D*)getInvMass(fFileBkg,ptmin,ptmax,0);
      hhbkg->Rebin(fRebinMass);
      hhbkg->GetXaxis()->SetRangeUser(minf,maxf);
      hhbkg->SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf",ptmin,fDmesonS.Data(),ptmax));
      setHistoDetails(hhbkg,kBlue+2,20,0.6)
      //hhbkg->Scale(1./hhbkg->Integral());
      //hhbkg->Scale(bkgscale);
      for(int k=0; k<hhbkg->GetNbinsX();k++){
        hhbkg->SetBinError(k+1,TMath::Sqrt(hhbkg->GetBinContent(k+1)));
      }

      TF1 *fitBkg = new TF1("fitBkg","[0]*[1]/(TMath::Exp([1]*2.1)-TMath::Exp([1]*1.71))*TMath::Exp([1]*x)",1.71,2.1);
      fitBkg->SetParameters(10000,-2);
      hhbkg->Fit("fitBkg","0RM");
      fitBkg->SetLineColor(4);

      cmassbkg->cd(iter+1);
      hhbkg->Draw();
      fitBkg->Draw("same");

      double bkgLevel = fitBkg->Integral(massmin,massmax)/hhbkg->GetBinWidth(1);

return bkgLevel;

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



TH2* getRM(TString file, Bool_t isprompt=1) {
	
	const int NDMC = 15;
	
    TString histName = "histosD0UpgradeN";

    float jetmin = 0, jetmax = 120;
    float Dptmin = fPtbinsDA[0], Dptmax = 200; // = fPtbinsDA[fPtbinsDN];

    TH2D *hPtJet2d, *hPtJet;
    TH1D *hPtJetMC, *hPtJetreco, *hreco1d;

	//findex = 5;
	if(fIsHBins){
    for(int j=findex; j<=fHardBinsN-1; j++) {
        TFile *File = new TFile(Form("%s_%d.root",file.Data(),j),"read");
        if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
        TDirectoryFile* dir = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
        TH2D *hMC;
        TH1D *hMC1d, *hreco1d;
        // loop over D mesons
        for(int i=0;i<NDMC; i++){

            TList *histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
            THnSparseF *sparseMC = (THnSparseF*)histList->FindObject("ResponseMatrix");

            sparseMC->GetAxis(6)->SetRangeUser(jetmin,jetmax); // jet pT gen
            sparseMC->GetAxis(7)->SetRangeUser(Dptmin,Dptmax); // jet pT gen

			sparseMC->GetAxis(1)->SetRangeUser(jetmin,jetmax); // jet pT gen
            sparseMC->GetAxis(2)->SetRangeUser(Dptmin,Dptmax); // jet pT gen

			hPtJet = (TH2D*)sparseMC->Projection(6,1,"E");  
			hPtJet->Sumw2();
		
			TH1D* hmc11 = (TH1D*)sparseMC->Projection(6);
			TH1D* hreco11 = (TH1D*)sparseMC->Projection(1);
		
			//fHardBinsMax[j-2]*fHBrejection
			
                     
            if (!i){
              hMC = (TH2D*)hPtJet->Clone("hMC");
              hMC->Sumw2();
              
              hMC1d = (TH1D*)hmc11->Clone("hMC1d");
              hMC1d->Sumw2();
              hreco1d = (TH1D*)hreco11->Clone("hreco1d");
              hreco1d->Sumw2();
            }
            else {
              hMC->Add(hPtJet);
              
              hMC1d->Add(hmc11);
              hreco1d->Add(hreco11);
            }			
		
		
        }
        // end of loop over D mesons
        TDirectoryFile *dir2 = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
        TList *histList2 = (TList*)dir2->Get("histosD0UpgradeN0MCrec");
        TH1F* hEvents = (TH1F*)histList2->FindObject("hstat");
        double nEvents = hEvents->GetBinContent(2); // number of analyzed events, 2: selected

		for(int mm=1; mm<=hMC->GetNbinsX();mm++){
				for(int nn=1; nn<=hMC->GetNbinsY();nn++){
					double cont = hMC->GetBinContent(mm,nn);
					//double cent = hMC->GetBinCenter(mm,nn);
					//if(cent > fHardBinsMax[j-2]*fHBrejection) hMC->SetBinContent(mm,nn,0);
					//if(cont < 5) hMC->SetBinContent(mm,nn,0);
					//else hMC->SetBinContent(mm,nn,hMC->GetBinContent(mm,nn)*fHBWeightsC[j-2]);
				}
		}

		/*for(int j=1; j<=Mtx->GetNbinsY()+1; j++) {

			double value = Mtx->GetYaxis()->GetBinCenter(j);
			double c = h->GetBinContent(h->GetXaxis()->FindBin(value));
			if (divide && c)  c = 1./c;
			//else c = 1.;
			for(int i=1; i<=Mtx->GetNbinsX()+1; i++) {
				Mtx->SetBinContent(i, j, Mtx->GetBinContent(i,j)*c);
				//Mtx->SetBinError(i, j, Mtx->GetBinError(i,j)*c);
			}
		}*/

		 hMC->Scale(100000./nEvents);
         if(isprompt) hMC->Scale(fHBWeightsC[j-2]);
	     else hMC->Scale(fHBWeightsB[j-2]);
	     
	     hMC1d->Scale(100000./nEvents);
         if(isprompt) hMC1d->Scale(fHBWeightsC[j-2]);
	     else hMC1d->Scale(fHBWeightsB[j-2]);
	     
	     hreco1d->Scale(100000./nEvents);
         if(isprompt) hreco1d->Scale(fHBWeightsC[j-2]);
	     else hreco1d->Scale(fHBWeightsB[j-2]);
         
        if(j==findex) {
          hPtJet2d = (TH2D*)hMC->Clone("hPtJet2d");
          hPtJet2d->Sumw2();
          
          hPtJetMC = (TH1D*)hMC1d->Clone("hPtJetMC");
          hPtJetMC->Sumw2();
          
          hPtJetreco = (TH1D*)hreco1d->Clone("hPtJetreco");
          hPtJetreco->Sumw2();
         
        }
        else {
          hPtJet2d->Add(hMC);
          
          hPtJetMC->Add(hMC1d);
          hPtJetreco->Add(hreco1d);
          
        } 
    }
      
  }
  else {
	
   TString datafile = file + fSMBBin;
    TFile *File = new TFile(datafile,"read");
    TDirectoryFile* dir = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
  	for(int i=0; i<NDMC; i++){
    
      TList *histList =   (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
      THnSparseF *sparseMC = (THnSparseF*)histList->FindObject("ResponseMatrix");
      
			sparseMC->GetAxis(6)->SetRangeUser(fJetptmin,jetmax); // jet pT gen
            sparseMC->GetAxis(7)->SetRangeUser(Dptmin,Dptmax); // jet pT gen

			sparseMC->GetAxis(1)->SetRangeUser(fJetptmin,jetmax); // jet pT gen
            sparseMC->GetAxis(2)->SetRangeUser(Dptmin,Dptmax); // jet pT gen

			hPtJet = (TH2D*)sparseMC->Projection(6,1,"E");  
			hPtJet->Sumw2();
		
                     
            if (!i){
              hPtJet2d = (TH2D*)hPtJet->Clone("hMC");
              hPtJet2d->Sumw2();
            }
            else {
              hPtJet2d->Add(hPtJet);
            }

  	}

  }
  
    hPtJet2d->SetTitle();
    hPtJet2d->SetName("hPtJet2d");
    hPtJet2d->GetXaxis()->SetTitle("p_{T,ch jet}^{rec.} (GeV/#it{c})");
    hPtJet2d->GetYaxis()->SetTitle("p_{T,ch jet}^{gen.} (GeV/#it{c})");
	//hPtJet2d->SetMinimum(0.0001);

    TCanvas *cjetPt2d = new TCanvas("cjetPt2d","cjetPt2d",1200,800);
    cjetPt2d->SetLogz();
    hPtJet2d->Draw("colz2");


	SaveCanvas(cjetPt2d, Form("%s/ResMatrix",OUTDIRECTORYPLOTS.Data()));

    //TFile *ofile = new TFile(Form("%s/DetMatrix_Dpt%d_%d.root",outDir, (int)Dptmin, (int)Dptmax),"RECREATE");
  /*  TFile *ofile = new TFile(Form("%s/DetMatrix_%s.root",outDir.Data(), isPrompt ? "prompt" : "nonPrompt" ),"RECREATE");

    hPtJetGen->Write();
    hPtJetRec->Write();
    hPtJet2d->Write();
    ofile->Close();*/
    
    TH2D *detMatrix = (TH2D*)hPtJet2d->Clone("detMatrix");
    detMatrix->SetName("detMatrix");
    //detMatrix->SetMinimum(0.0001);
    
    for(int i=0; i<=detMatrix->GetNbinsX()+1;i++){
        for(int j=0; j<=detMatrix->GetNbinsY()+1;j++){

            double cont = detMatrix->GetBinContent(i,j);
            if(i==0 && j==0)detMatrix->SetBinContent(i,j,0);
            else if(i==detMatrix->GetNbinsX()+1 && j==detMatrix->GetNbinsY()+1)detMatrix->SetBinContent(i,j,0);
            else detMatrix->SetBinContent(i,j,cont);

        }
    }
    
    
	TH2D *fDetMatrix = Rebin2D("fDetMatrix", detMatrix, fPtbinsJetMeasN, fPtbinsJetMeasA, fPtbinsJetTrueN, fPtbinsJetTrueA,0);
	fDetMatrix->SetMinimum(0.001);
	fDetMatrix->SetTitle();
    fDetMatrix->SetName("fDetMatrix");
    fDetMatrix->GetXaxis()->SetTitle("p_{T,ch jet}^{rec.} (GeV/#it{c})");
    fDetMatrix->GetYaxis()->SetTitle("p_{T,ch jet}^{gen.} (GeV/#it{c})");

	TCanvas *cmatrix = new TCanvas("cmatrix","cmatrix",1200,800);
    cmatrix->SetLogz();
    fDetMatrix->Draw("colz2");


	TH1D *hmeas = (TH1D*)detMatrix->ProjectionX("hmeas");
	TH1D *htrue = (TH1D*)detMatrix->ProjectionY("htrue");

	hmeas->SetLineColor(4);
	htrue->SetLineColor(2);
	
	
	hPtJetreco->SetLineColor(4);
	hPtJetMC->SetLineColor(2);

	TCanvas *cPt = new TCanvas("cPt","cPt",800,600);
	cPt->SetLogy();
	hPtJetreco->Draw();
	hPtJetMC->Draw("same");
	//hmeas->Draw();
	//htrue->Draw("same");


	TH1D* hratio = (TH1D*)hPtJetreco->Clone("hratio");
	hratio->Divide(hPtJetMC);
	TCanvas *cRatio = new TCanvas("cRatio","cRatio",1200,600);
	hratio->Draw();

	cout << "\n\n reco jets: " << hPtJetreco->Integral() << "\t MC jets: " << hPtJetMC->Integral() << endl;

	SaveCanvas(cmatrix, Form("%s/ResDetMatrix_rebin",OUTDIRECTORYPLOTS.Data()));
	SaveCanvas(cPt, Form("%s/ResDetMatrixProj",OUTDIRECTORYPLOTS.Data()));
	SaveCanvas(cRatio, Form("%s/ResDetMatrixProjRatio",OUTDIRECTORYPLOTS.Data()));

   return fDetMatrix;

}


TH1* getPowhegJetSpectra(TString data, Bool_t isPrompt=0, Bool_t isEffScale=0, Bool_t isEffRatioScale=0, Bool_t isHIScale=0, Bool_t isRaaSys=0, Bool_t recal=0){
  
  
  TString outname = data;
  outname += "_JetPt_DpT";
  outname += fPtbinsDA[0];
  if(isPrompt) outname += "_prompt";
  else outname += "_nonprompt";
  if(isEffScale) outname += "_effScaled";
  else if(isEffRatioScale) outname += "_effRatioScaled";
  if(isHIScale) outname += "_RaaScaled"; 
  if(isRaaSys) outname += "_RaaSys";
  outname += ".root";
 
  if(!recal) {
    TFile *file = new TFile(outname,"read");
    if(!file->IsZombie()){
      TH1D *hout = (TH1D*)file->Get("hPt");
      return hout;
    }
  }
  
  TFile *fileInput = new TFile(Form("%s.root",data.Data()),"read");
  if(!fileInput){
    std::cout << "File " << fileInput << " cannot be opened! check your file path!" << std::endl; return kFALSE;
  }
  TList* dir=(TList*)fileInput->Get("AliAnalysisTaskDmesonJets_histos");
  if(!dir) {
    std::cout << "Error in getting dir! Exiting..." << std::endl;
    return kFALSE;
  }

  TH1D *hxsection = (TH1D*)dir->FindObject("fHistXsection");
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
  tree->SetBranchAddress("Jet_AKTChargedR030_pt_scheme",&brJet);
  if(!tree || !brD || !brJet) {
    std::cout << "Error in setting the tree/branch names! Exiting..." << std::endl;
    return NULL;
  }

  TH1D *hjetpt[fPtbinsDN];
  for (int j=0; j<fPtbinsDN; j++) {
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

    for (int j=0; j<fPtbinsDN; j++) {
        if (brD->fPt < fPtbinsDA[j] || brD->fPt >= fPtbinsDA[j+1]) continue;
        hjetpt[j]->Fill(brJet->fPt);
    }//end of D-meson pT bin loop

  }

  TH1D*  hPt = new TH1D("hPt","hPt",100,0,100);
  for (int j=0; j<fPtbinsDN; j++){
    double effC=1, effB=1, eff=1;
    double raa = 1;
    double pt = (fPtbinsDA[j]+fPtbinsDA[j+1])/2.;
    effC = GetContent(fHEff_prompt,pt);
    effB = GetContent(fHEff_nprompt,pt);
    
    if(isEffRatioScale) {
		eff = effB/effC;
	}
	else if(isEffScale) {
			if(isPrompt) eff = effC;
			else eff = effB;
	}
	else eff = 1;
	
    if(isHIScale){
		if( isPrompt || isRaaSys) raa = fRaaC[j+fraabin];
		else if (isRaaSys) raa = fRaaBSys[j_fraabin];
		//if(isPrompt) raa = fRaaC[j+fraabin];
		//else if(fBRaaSys || isRaaSys) raa = fRaaB[j+fraabin]/2.;
		else raa = fRaaB[j+fraabin];
	}
	double scale = eff*raa;
	//if(isEffScale && isHIScale) cout << "\n\n =========== B scaling: " << scale << endl;
    if (!j){
        hPt = (TH1D*)hjetpt[j]->Clone("hPt");
        hPt->Scale(scale);
    }
    else{ 
		hjetpt[j]->Scale(scale);
		hPt->Add(hjetpt[j]);
		
	}
  }

  hPt->Scale(scaling);
  hPt->SetName("hPt");

  TFile *outfile = new TFile(outname,"RECREATE");
  hPt->Write();
  //hPt_reb->Write();
  outfile->Close();

  return hPt;

}

TH1* getPowhegDSpectra(TString data, Bool_t isPrompt=1, Bool_t recal=0){

  TString outname = data;
  if(isPrompt) outname += "_prompt";
  else outname += "_nonprompt";
  outname += ".root";
  
  if(!recal) {
    TFile *file = new TFile(outname,"read");
    if(!file->IsZombie()){
      cout << "!!!!We are here" << endl;
      TH1D *hout = (TH1D*)file->Get("hPt_reb");
      return hout;
    }
  }

  TFile *fileInput = new TFile(Form("%s.root",data.Data()),"read");
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

  TH1D *hh_tmp = (TH1D*)hPt->Clone("hh_tmp");
  TH1D *hPt_reb = (TH1D*)hh_tmp->Rebin(fPtbinsDN,"hPt_reb",fPtbinsDA);

  TFile *outfile = new TFile(outname,"RECREATE");
  hPt->Write();
  hPt_reb->Write();
  outfile->Close();

  return hPt;

}

TH3D* get3DInvMass(TString data, Bool_t isBkg = 0) {

  const int ND = 15;

  TString histName;
  if(fDmesonSpecie) histName = "histosDStarMBN";
  else histName = "histosD0UpgradeN";
  TDirectoryFile* dir;
  TList *histList;
  THnSparseF *sparse;
  TH3D *hInvMassptD;

  if(fIsHBins && !isBkg){

    TCanvas *cJetPt = new TCanvas();
    TH1D *hjetpt;
    for(int j=2; j<=fHardBinsN; j++) {
        TFile *File = new TFile(Form("%s_%d.root",data.Data(),j),"read");
        if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
        dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");

        TH3D *hmass;
        for(int i=0;i<ND; i++){
            histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
            sparse = (THnSparseF*)histList->FindObject("hsDphiz");
            sparse->GetAxis(0)->SetRangeUser(fZmin,fZmax);
            //sparse->GetAxis(1)->SetRangeUser(fJetptmin,fJetptmax);
            sparse->GetAxis(1)->SetRangeUser(fJetptmin,fHardBinsMax[j-2]*4);
            if(i==0) hmass=(TH3D*)sparse->Projection(3,1,2);
            else hmass->Add((TH3D*)sparse->Projection(3,1,2));
        }
        TDirectoryFile *dir2 = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
        TList *histList2 = (TList*)dir2->Get("histosD0UpgradeN0MCrec");
        TH1F* hEvents = (TH1F*)histList2->FindObject("hstat");
        //double nEvSel = hEvents->GetBinContent(2);
        double nEvents = hEvents->GetBinContent(1); // number of analyzed events
        hmass->Scale(1000000./nEvents);
        //hmass->Scale(1./nEvents);
        hmass->Scale(fHBWeightsC[j-2]);

        if(j==2) {
          hInvMassptD = (TH3D*)hmass->Clone("hInvMassptD");
          //hInvMassptD->Sumw2();
        }
        else {
          hInvMassptD->Add(hmass);
        }
        TH1D *hpt = (TH1D*)hInvMassptD->ProjectionY();
        if(j==2) hjetpt = (TH1D*)hpt->Clone("hjetpt");
        else hjetpt->Add(hpt);
        hpt->SetLineColor(j+1);
        hpt->SetMarkerColor(j+1);
        if(j==2) hpt->DrawCopy();
        else hpt->DrawCopy("same");
    }

    hjetpt->SetLineColor(1);
    hjetpt->Draw("same");

  }
  else {
    TFile *File = new TFile(data+fSMBBin,"read");;
    if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
    dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");

    for(int i=0;i<ND; i++){
        histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
        sparse = (THnSparseF*)histList->FindObject("hsDphiz");
        sparse->GetAxis(0)->SetRangeUser(fZmin,fZmax);
        sparse->GetAxis(1)->SetRangeUser(fJetptmin,fJetptmax);
        if(i==0) { hInvMassptD=(TH3D*)sparse->Projection(3,1,2); hInvMassptD->Sumw2(); }
        else hInvMassptD->Add((TH3D*)sparse->Projection(3,1,2));
    }

    if(fDmesonSpecie) hInvMassptD->GetXaxis()->SetTitle(Form("m(%s)(GeV/c^{2})","K#pi#pi-K#pi"));
    else hInvMassptD->GetXaxis()->SetTitle(Form("m(%s)(GeV/c^{2})","K#pi"));
    hInvMassptD->GetYaxis()->SetTitle("p_{T}^{jet}");
    hInvMassptD->GetZaxis()->SetTitle("p_{T}^{D}");
  }



return hInvMassptD;

}


void getSB() {

/*
  1.5289473684210524, 0.059095553466660504
2.515789473684211, 0.15771704865301248
3.489473684210526, 0.3044765420856639
4.489473684210527, 0.49634797527710756
5.489473684210526, 1.072365490319609
7.015789473684211, 2.411472486746194
10.00263157894737, 6.2924564093126225
14.015789473684208, 8.553215338100452 */

const int       ptbinsDN = 8;
double          ptbinsDA[ptbinsDN+1] = { 1,2,3,4,5,6,8,12,16 };
double SBtdr[ptbinsDN] = {
0.059095553466660504,
0.15771704865301248,
0.3044765420856639,
0.49634797527710756,
1.072365490319609,
2.411472486746194,
6.2924564093126225,
8.553215338100452
};

  //fHSB = new TH1D("fHSB","S/B",fPtbinsDN,fPtbinsDA);
  fHSB = new TH1D("fHSB","S/B",ptbinsDN,ptbinsDA);
  for(int i=0; i<ptbinsDN; i++){
    //fHSB->SetBinContent(i+1,fSB[i]);
    fHSB->SetBinContent(i+1,SBtdr[i]);

  }

  //fSBfun = new TF1("fSBfun","[1]*(1-exp(-x/[0]))",1,24);
  TF1 *SBfun2 = new TF1("SBfun2","[1]*(tanh(x/[0])+[2])",5,16);
  SBfun2->SetParameter(0,1);
  SBfun2->SetParameter(1,1);
  SBfun2->SetParameter(2,1);
  //fSBfun->FixParameter(0,4.36070e+03); // fitted
  //fSBfun->FixParameter(0,10e+03);
  //fSBfun->FixParameter(1,1.14255e+03);
  fHSB->Fit(SBfun2,"MR0");
  setHistoDetails(fHSB,kMagenta+2,20,1.2);
  SBfun2->SetLineColor(kRed+1);
  SBfun2->SetLineWidth(3);
  SBfun2->SetLineStyle(2);

  fSBfun = new TF1("fSBfun","[1]*(tanh(x/[0])+[2])",5,50);
  fSBfun->SetParameter(0,SBfun2->GetParameter(0));
  fSBfun->SetParameter(1,SBfun2->GetParameter(1));
  fSBfun->SetParameter(2,SBfun2->GetParameter(2));
  fSBfun->SetLineColor(kRed+1);
  fSBfun->SetLineWidth(3);
  fSBfun->SetLineStyle(2);
  fHSB->GetXaxis()->SetRangeUser(1,50);

TH1D *htmp = new TH1D("htmp","htmp",fPtbinsDN,fPtbinsDA);
htmp->SetTitle();
htmp->GetXaxis()->SetTitle("p_{T,D^{0}}");
htmp->GetYaxis()->SetTitle("S/B");
htmp->SetMinimum(0.01);
  TCanvas *cSB = new TCanvas("cSB","cSB",1200,800);
  cSB->SetLogy();
  fHSB->GetXaxis()->SetTitle("p_{T, D^{0}}");
  fHSB->GetYaxis()->SetTitle("S/B");
  htmp->Draw();
  fSBfun->Draw("same");
  fHSB->Draw("psame");
  //fSBfun->Draw("same");


    SaveCanvas(cSB, Form("%s/SigOverBkg_fit",OUTDIRECTORYPLOTS.Data()));

return;
}



void getRaa() {


  TH1D *hRaaC = new TH1D("hRaaC","RAA",fPtbinsDN,fPtbinsDA);
  TH1D *hRaaB = new TH1D("hRaaB","RAA",fPtbinsDN,fPtbinsDA);

  for(int i=0; i<fPtbinsDN; i++){
    hRaaC->SetBinContent(i+1,fRaaC[i]);
    hRaaB->SetBinContent(i+1,fRaaB[i]);
  }

	hRaaC->GetYaxis()->SetRangeUser(0,1.2);
	hRaaC->GetXaxis()->SetTitle("p_{T,D^{0}}");
	hRaaC->GetYaxis()->SetTitle("R_{AA}");
	
	setHistoDetails(hRaaC,kRed+1,20,1.4);
	setHistoDetails(hRaaB,kBlue+1,21,1.4);
  
  	TLegend *leg = new TLegend(0.5,0.7,0.85,0.8);
    leg->SetTextSize(0.045);
    leg->SetBorderSize(0);
    leg->AddEntry(hRaaC,"Prompt","p");
    leg->AddEntry(hRaaB,"Feed-down","p");
  TCanvas *cRaa = new TCanvas("cRaa","cRaa",1200,800);
//  cRaa->SetLogy();
cRaa->cd();
  hRaaC->Draw("p");
  hRaaB->Draw("psame");

    leg->Draw("same");

    SaveCanvas(cRaa, Form("%s/Raa",OUTDIRECTORYPLOTS.Data()));

return;
}


void getMCJetPt(TString data){

  TString histName;
  if(fDmesonSpecie) histName = "histosDStarMBN";
  else histName = "histosD0UpgradeN";

  const int ND = 15;
  TCanvas *cJetPt = new TCanvas("cJetPt","cJetPt",1200,800);
  cJetPt->SetLogy();
  TH1D *hjetpt;

  for(int j=2; j<=fHardBinsN; j++) {
      TFile *File = new TFile(Form("%s_%d.root",data.Data(),j),"read");
      if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
      TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");

      TH1D *hjet;
      for(int i=0;i<ND; i++){
          TList *histList =   (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
          THnSparseF *sparseMC = (THnSparseF*)histList->FindObject("ResponseMatrix");
          TH1D *h = (TH1D*)sparseMC->Projection(6);
          if(!i) hjet = (TH1D*)h->Clone("hjet");
          else hjet->Add(h);

      }
      TDirectoryFile *dir2 = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
      TList *histList2 = (TList*)dir2->Get("histosD0UpgradeN0MCrec");
      TH1F* hEvents = (TH1F*)histList2->FindObject("hstat");
      //double nEvSel = hEvents->GetBinContent(2);
      double nEvents = hEvents->GetBinContent(1); // number of analyzed events
      hjet->Scale(1./nEvents);
      hjet->Scale(fHBWeightsC[j-2]);
      hjet->SetMinimum(0.00000001);
      hjet->GetXaxis()->SetRangeUser(0,150);
      hjet->GetXaxis()->SetTitle("jet p_{T}");
      hjet->SetTitle();
      if(j==2) {
        hjetpt = (TH1D*)hjet->Clone("hjetpt");
        hjetpt->Sumw2();
      }
      else {
        hjetpt->Add(hjet);
      }

      hjet->SetLineColor(colors2[j-2]);
      hjet->SetMarkerColor(colors2[j-2]);
      if(j==2) hjet->DrawCopy();
      else hjet->DrawCopy("same");
  }

  hjetpt->SetLineColor(1);
  hjetpt->SetMarkerStyle(20);
  hjetpt->SetMarkerSize(0.8);
  hjetpt->Draw("same");

  TFile *File = new TFile(Form("%s_%d.root",data.Data(),1),"read");
  if(!File) { cout << "==== WRONG FILE WITH DATA =====" << endl; return;}
  TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");

  TH1D *hjet2;
  for(int i=0;i<ND; i++){
      TList *histList =   (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
      THnSparseF *sparseMC = (THnSparseF*)histList->FindObject("ResponseMatrix");
      TH1D *h = (TH1D*)sparseMC->Projection(6);
      if(!i) hjet2 = (TH1D*)h->Clone("hjet");
      else hjet2->Add(h);

  }
  TDirectoryFile *dir2 = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
  TList *histList2 = (TList*)dir2->Get("histosD0UpgradeN0MCrec");
  TH1F* hEvents = (TH1F*)histList2->FindObject("hstat");
  //double nEvSel = hEvents->GetBinContent(2);
  double nEvents = hEvents->GetBinContent(1); // number of analyzed events
  //hmass->Scale(1000000./nEvents);
  hjet2->Scale(1./nEvents);
  hjet2->Scale(1.54866666666666664E+01);
  hjet2->SetMinimum(0.00000001);
  hjet2->SetLineColor(kRed+2);
  hjet2->SetMarkerColor(kRed+2);
  hjet2->DrawCopy("same");
  
  SaveCanvas(cJetPt, Form("%s/DjetpT_POWHEG",OUTDIRECTORYPLOTS.Data()));

}



/// rebin in 2d variable size - no such routine in Root
TH2D * Rebin2D(const char* name, TH2D *h, int nx, const double *binx, int ny, const double *biny, bool crop) {
	if (!h) {
		cerr << "Warning in <AliHeavyUnfoldTools::Rebin2D> : h==0." << endl;
		return 0;
	}

	TAxis *xaxis = h->GetXaxis();
    TAxis *yaxis = h->GetYaxis();

	TH2D * hre = new TH2D(name,name,nx,binx,ny,biny);
    hre->Sumw2();
    for (int i=1; i<=xaxis->GetNbins();i++) {
        for (int j=1; j<=yaxis->GetNbins();j++) {
            hre->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),h->GetBinContent(i,j));
        }
    }

    //for(int i=0;i<=hre->GetNbinsX();i++){
        for(int j=0;j<=hre->GetNbinsY()+1;j++){
            hre->SetBinContent(0,j,0);
            hre->SetBinError(0,j,0);
			hre->SetBinContent(hre->GetNbinsX()+1,j,0);
            hre->SetBinError(hre->GetNbinsX()+1,j,0);

    }
    //}


	return hre;
}

/// get pearson coeffs from covariance matrix
TH2D * getPearsonCoeffs(const TMatrixD &covMatrix) {

	Int_t nrows = covMatrix.GetNrows();
	Int_t ncols = covMatrix.GetNcols();

	TH2D* PearsonCoeffs = new TH2D("PearsonCoeffs","Pearson Coefficients", nrows, 0, nrows, ncols, 0, ncols);
	for(Int_t row = 0; row<nrows; row++) {
		for(Int_t col = 0; col<ncols; col++) {
			Double_t pearson = 0.;
			if(covMatrix(row,row)!=0. && covMatrix(col,col)!=0.)
				pearson = covMatrix(row,col)/TMath::Sqrt(covMatrix(row,row)*covMatrix(col,col));
			PearsonCoeffs->SetBinContent(row+1,col+1, pearson);
		}
	}

	return PearsonCoeffs;
}

/// Weight matrix along y axis by histo values
void WeightMatrixY(TH2D * Mtx, TH1D * h, bool divide) {
	if (!Mtx) {
		cerr << "Warning in <AliHeavyUnfoldTools::WeightMatrixY> : Mtx==0." << endl;
		return;
	}
	if (!h) {
		cerr << "Warning in <AliHeavyUnfoldTools::WeightMatrixY> : h==0." << endl;
		return;
	}

	for(int j=1; j<=Mtx->GetNbinsY()+1; j++) {

		double value = Mtx->GetYaxis()->GetBinCenter(j);
		double c = h->GetBinContent(h->GetXaxis()->FindBin(value));
		if (divide && c)  c = 1./c;
		//else c = 1.;
		for(int i=1; i<=Mtx->GetNbinsX()+1; i++) {
			Mtx->SetBinContent(i, j, Mtx->GetBinContent(i,j)*c);
			//Mtx->SetBinError(i, j, Mtx->GetBinError(i,j)*c);
		}
	}
}

/// Weight matrix along y axis by histo values
void WeightMatrixX(TH2D * Mtx, TH1D * h, bool divide) {
	if (!Mtx) {
		cerr << "Warning in <AliHeavyUnfoldTools::WeightMatrixY> : Mtx==0." << endl;
		return;
	}
	if (!h) {
		cerr << "Warning in <AliHeavyUnfoldTools::WeightMatrixY> : h==0." << endl;
		return;
	}

	for(int j=1; j<=Mtx->GetNbinsX()+1; j++) {

		double value = Mtx->GetXaxis()->GetBinCenter(j);
		double c = h->GetBinContent(h->GetYaxis()->FindBin(value));
		if (divide && c)  c = 1./c;
		//else c = 1.;
		for(int i=1; i<=Mtx->GetNbinsY()+1; i++) {
			Mtx->SetBinContent(i, j, Mtx->GetBinContent(i,j)*c);
			//Mtx->SetBinError(i, j, Mtx->GetBinError(i,j)*c);
		}
	}
}

double GetContent(TH1 *hh, double bin){
    return hh->GetBinContent(hh->GetXaxis()->FindBin(bin));
}


void setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Size_t size = 0.9, Style_t Lstyle=1, Width_t width=2){

    h->SetMarkerStyle(Mstyle);
    h->SetMarkerColor(color);
    h->SetMarkerSize(size);
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetLineStyle(Lstyle);
    h->SetTitle(0);

    return;

}

void SaveCanvas(TCanvas *c, TString name = "tmp"){

    c->SaveAs(Form("%s.png",name.Data()));
    c->SaveAs(Form("%s.pdf",name.Data()));

}



void runAllEff(){


  TString posbin = "_moreBins";

  TString posf = "";
  TString fFileEff_prompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_charm";
  TString fFileEff_nprompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_beauty";
  fFileEff_prompt += posf;
  fFileEff_nprompt += posf;

  fIsHBins = 1;
  fIsHBEffWeighting = 1;
  fOutNamePostfix =  "HB";
  fOutNamePostfix += posf;
  fOutNamePostfix += "_weighted";
  fOutNamePostfix += posbin;
  getEfficiencies();

  fIsHBins = 1;
  fIsHBEffWeighting = 0;
  fOutNamePostfix =  "HB";
  fOutNamePostfix += posf;
  fOutNamePostfix += "_notweighted";
  fOutNamePostfix += posbin;
  getEfficiencies();

  fIsHBins = 0;
  fIsHBEffWeighting = 0;
  fOutNamePostfix = posf;
  fOutNamePostfix += posbin;
  getEfficiencies();


  TString fFileEff_prompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_charm";
  TString fFileEff_nprompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_beauty";
  TString posf = "Bis";
  fFileEff_prompt += posf;
  fFileEff_nprompt += posf;

  fIsHBins = 1;
  fIsHBEffWeighting = 1;
  fOutNamePostfix =  "HB";
  fOutNamePostfix += posf;
  fOutNamePostfix += "_weighted";
  fOutNamePostfix += posbin;
  getEfficiencies();

  fIsHBins = 1;
  fIsHBEffWeighting = 0;
  fOutNamePostfix =  "HB";
  fOutNamePostfix += posf;
  fOutNamePostfix += "_notweighted";
  fOutNamePostfix += posbin;
  getEfficiencies();

  fIsHBins = 0;
  fIsHBEffWeighting = 0;
  fOutNamePostfix = posf;
  fOutNamePostfix += posbin;
  getEfficiencies();

  TString fFileEff_prompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_charm";
  TString fFileEff_nprompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_beauty";
  TString posf = "Tres";
  fFileEff_prompt += posf;
  fFileEff_nprompt += posf;

  fIsHBins = 1;
  fIsHBEffWeighting = 1;
  fOutNamePostfix =  "HB";
  fOutNamePostfix += posf;
  fOutNamePostfix += "_weighted";
  fOutNamePostfix += posbin;
  getEfficiencies();

  fIsHBins = 1;
  fIsHBEffWeighting = 0;
  fOutNamePostfix =  "HB";
  fOutNamePostfix += posf;
  fOutNamePostfix += "_notweighted";
  fOutNamePostfix += posbin;
  getEfficiencies();

  fIsHBins = 0;
  fIsHBEffWeighting = 0;
  fOutNamePostfix = posf;
  fOutNamePostfix += posbin;
  getEfficiencies();


/*
  TString posbin = "_stdBins";
  //TString posbin = "_moreBins";

  TString posf = "";
  TString fFileEff_prompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_charm";
  TString fFileEff_nprompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_beauty";
  fFileEff_prompt += posf;
  fFileEff_nprompt += posf;

  for (int i=2; i<7;i++) {
    fIsHBins = 1;
    fIsHBEffWeighting = 1;
    findex = i;
    fOutNamePostfix =  "HB";
    fOutNamePostfix += posf;
    fOutNamePostfix += findex;
    fOutNamePostfix += "_weighted";
    fOutNamePostfix += posbin;
    getEfficiencies();
  }

  for (int i=2; i<7;i++) {
    fIsHBins = 1;
    fIsHBEffWeighting = 0;
    findex = i;
    fOutNamePostfix =  "HB";
    fOutNamePostfix += posf;
    fOutNamePostfix += findex;
    fOutNamePostfix += "_notweighted";
    fOutNamePostfix += posbin;
    getEfficiencies();
  }

  TString posf = "Bis";
  TString fFileEff_prompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_charm";
  TString fFileEff_nprompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_beauty";
  fFileEff_prompt += posf;
  fFileEff_nprompt += posf;

  for (int i=2; i<7;i++) {
    fIsHBins = 1;
    fIsHBEffWeighting = 1;
    findex = i;
    fOutNamePostfix =  "HB";
    fOutNamePostfix += posf;
    fOutNamePostfix += findex;
    fOutNamePostfix += "_weighted";
    fOutNamePostfix += posbin;
    getEfficiencies();
  }

  for (int i=2; i<7;i++) {
    fIsHBins = 1;
    fIsHBEffWeighting = 0;
    findex = i;
    fOutNamePostfix =  "HB";
    fOutNamePostfix += posf;
    fOutNamePostfix += findex;
    fOutNamePostfix += "_notweighted";
    fOutNamePostfix += posbin;
    getEfficiencies();
  }

  TString posf = "Tres";
  TString fFileEff_prompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_charm";
  TString fFileEff_nprompt = "/home/basia/Work/alice/analysis/upgradeProjections/out/AnalysisResults_beauty";
  fFileEff_prompt += posf;
  fFileEff_nprompt += posf;

  for (int i=2; i<7;i++) {
    fIsHBins = 1;
    fIsHBEffWeighting = 1;
    findex = i;
    fOutNamePostfix =  "HB";
    fOutNamePostfix += posf;
    fOutNamePostfix += findex;
    fOutNamePostfix += "_weighted";
    fOutNamePostfix += posbin;
    getEfficiencies();
  }

  for (int i=2; i<7;i++) {
    fIsHBins = 1;
    fIsHBEffWeighting = 0;
    findex = i;
    fOutNamePostfix =  "HB";
    fOutNamePostfix += posf;
    fOutNamePostfix += findex;
    fOutNamePostfix += "_notweighted";
    fOutNamePostfix += posbin;
    getEfficiencies();
  }
*/

  return;
}


