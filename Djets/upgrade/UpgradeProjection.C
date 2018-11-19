//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------

#include "config.h"

using namespace std;

int fraabin = 0;

void Init() {

  fIsHBins = true;
  fIsHBEffWeighting = true; //this is for efficiencies
  fHBrejection = 1.;
  findex = 2; // reject the 2 HB for the eff calculation (introduces flucutations)

  fJetptmin = 0;
  fJetptmax = 100;

  double sigBkg[] = { 5.824E-2, 1.5562E-1, 3.0889E-1, 5.0915E-1, 1.0686, 2.4596, 5, 8, 9, 10, 10, 10, 10, 10,10,10,10,10,10,10,10,10,10,10 };
  //double sigBkg[] = { 5.824E-2, 1.5562E-1, 3.0889E-1, 5.0915E-1, 1.0686, 2.4596, 6.3135, 6.3135, 9.10, 10, 10, 10, 10, 10 };
  fSB = new double[fPtbinsDN];
  for(int i=0;i<fPtbinsDN;i++) fSB[i] = sigBkg[i];

/* Raa_C
1.600710299340438, 0.9361111111111113
2.5680703534584826, 0.8472222222222225
3.5214780991036694, 0.5750000000000004
4.542322002367664, 0.38888888888888906
5.445416878065279, 0.255555555555556
6.478521900896329, 0.23055555555555607
7.3913411128023, 0.2250000000000001
8.973025536952477, 0.20833333333333393
11.043252156265854, 0.21111111111111125
14.149543378995428, 0.22777777777777786
20.058346017250123, 0.2694444444444448
29.98752748181971, 0.33888888888888924
43.14265178420428, 0.39722222222222237
*/
/*
1.490693739424704, 1.0557873253440786
2.4720812182741128, 1.0929844575854308
3.506627185561195, 1.042076159575767
4.494829855235947, 0.9115762787749808
5.5171084790374145, 0.7845157313889322
7.050526414739613, 0.6111070408940522
8.992855799962399, 0.5133346125444087
11.003337093438615, 0.46539617748472195
14.002021056589586, 0.44675167605723876
20.20384470765181, 0.4713357870927599
29.745111863132166, 0.5151051569650134
39.86567023876668, 0.5898517859530994
*/


double raaC[] = { 0.936, 0.847, 0.575, 0.389, 0.256, 0.250, 0.208, 0.211, 0.228, 0.270, 0.339, 0.397 };
double raaB[] = { 1.0557873253440786,
1.0929844575854308,
1.042076159575767,
0.9115762787749808,
0.7845157313889322,
0.6111070408940522,
0.5133346125444087,
0.46539617748472195,
0.44675167605723876,
0.4713357870927599,
0.5151051569650134,
0.5898517859530994 };


 
  fRaaC = new double[fPtbinsDN];
  fRaaB = new double[fPtbinsDN];
  for(int i=0;i<fPtbinsDN;i++){
    fRaaC[i] = raaC[i];
    fRaaB[i] = raaB[i];
  }

  double weightC[5] = { 6.67666666666666675E+00,4.26666666666666694E-01,3.57000000000000026E-02,3.90000000000000025E-03,6.32666666666666685E-04 };
  double weightB[5] = { 6.68333333333333357E+00,4.26666666666666694E-01,3.56999999999999956E-02,3.89666666666666677E-03,6.27000000000000062E-04};
  for(int i=0; i<5; i++){
    fHBWeightsC[i] = weightC[i];
    fHBWeightsB[i] = weightB[i];
  }

}

void UpgradeProjection() {

  gStyle->SetOptStat(0000); //Mean and RMS shown
	gSystem->Exec(Form("mkdir -p %s/plots",OUTDIRECTORY.Data()));

    Init();

  //runAllEff();
  //return;

	getRaa();
	getSB();
  
	getEfficiencies();
  
	getPowhegDmesonSpectra();
	fHJetPtSpecrum = (TH1D*)extractSignal();

  if(fIsHBins) getMCJetPt(fFileSignal_prompt);

return;

  // get POWHEG D-jet x-sections
  fHPowhegPrompt = getPowhegJetSpectra(fFilePowhegPrompt);
  fHPowhegNonPrompt = getPowhegJetSpectra(fFilePowhegNonPrompt,0);
  TH1D *hPowhegNonPromptScaled = getPowhegJetSpectra(fFilePowhegNonPrompt,0,1);

return;
}

TH1* extractSignal() {


	 TFile *ofile = new TFile(Form("%s/D0jets.root",OUTDIRECTORY.Data()),"RECREATE");

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

  setHistoDetails(hSB,kGreen+2,21,0.9);

  TString signalFile = fFileSignal_prompt;
  if(!fIsHBins) signalFile += fSMBBin;
  TString signalBkgFile = fFileSignalBkg_prompt;
  if(!fIsHBins) signalBkgFile += fSMBBin;

  double bkgscale = getBkgScaling();

  for(int i=0; i<fPtbinsDN; i++){

      double pt = (fPtbinsDA[i]+fPtbinsDA[i+1])/2.;
      double effP = GetContent(fHEff_prompt,pt);
      double effNP = GetContent(fHEff_nprompt,pt);
      double sigP = GetContent(fHPowhegDPrompt,pt) * fSimScalingC * fRaaC[i+fraabin] * effP * fDataEv;
      double sigNP = GetContent(fHPowhegDNonPrompt,pt) * fSimScalingB * fRaaB[i+fraabin] * effNP * fDataEv;

      // inv. mass
      TH1D *hhsignal = (TH1D*)getInvMass(signalFile,fPtbinsDA[i],fPtbinsDA[i+1],1);
      hhsignal->Rebin(fRebinMass);
      hhsignal->GetXaxis()->SetRangeUser(1.71,2.1);
      setHistoDetails(hhsignal,kRed+2,20,0.6);
      hhsignal->SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));
      // scale to the POWHEG predictions (prompt + non-prompt expected yields)
      hhsignal->Scale(1./hhsignal->Integral());
      hhsignal->Scale(sigP+sigNP);
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

      // scale to the expected yields
      TH1D *hPt = (TH1D*)hhPt->Clone("hPt");
      hPt->Sumw2();
      setHistoDetails(hPt,2,25,0.7);
      hPt->SetTitle(Form("%.1lf < p_{T}^{%s} < %.1lf",fPtbinsDA[i],fDmesonS.Data(),fPtbinsDA[i+1]));
      hPt->GetYaxis()->SetTitle("dN");
      hPt->Scale(1./hPt->Integral());
      hPt->Scale(sigP+sigNP);
      hPt->SetMinimum(1);
      for(int k=1; k<=hPt->GetNbinsX();k++){
        hPt->SetBinError(k,TMath::Sqrt(hPt->GetBinContent(k)));
      }

      double bkgLevel = getBackground(massmin, massmax,fPtbinsDA[i],fPtbinsDA[i+1],i,cmassbkg);
      double bkg = bkgLevel*bkgscale * fRaaC[i];
      double sb;
      if(pt>6) sb = fSBfun->Eval(pt);
      else sb = GetContent(fHSB,pt);
      //bkg = (sigP+sigNP)/fSB[i];
      bkg = (sigP+sigNP)/sb;
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
      hSB->SetBinContent(i+1,(sigP+sigNP)/bkg);

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
  hSB->GetXaxis()->SetTitle("p_{T, D^{0}}");
  hSB->GetYaxis()->SetTitle("S/B");
  hSB->Draw();

  cJetPt->cd();
  cJetPt->SetLogy();
  hJetPt->SetTitle();
  hJetPt->GetYaxis()->SetTitle("dN");
  hJetPt->Draw();

  setHistoDetails(hJetPtRebin,2,21,1.2);
  hJetPtRebin->SetTitle();
  TH1D *hJetPtRebinScale = (TH1D*)hJetPtRebin->Clone("hJetPtRebinScale");
  hJetPtRebinScale->Scale(1,"width");
  hJetPtRebinScale->SetMinimum(1);
  hJetPtRebinScale->GetYaxis()->SetTitle("dN/dp_{T}");

  cJetPtReb->cd();
  cJetPtReb->SetLogy();
  hJetPtRebinScale->Draw();
  
   TH1D *hjetptspectrumRebUnc = (TH1D*)hJetPtRebin->Clone("hjetptspectrumRebUnc");
   hjetptspectrumRebUnc->GetYaxis()->SetTitle("Rel. unc.");

      for(int mm=1; mm<=hJetPtRebin->GetNbinsX();mm++){
                  double err;
                  if(hJetPtRebinScale->GetBinContent(mm)) err = hJetPtRebin->GetBinError(mm)/hJetPtRebin->GetBinContent(mm);
                  else err = 0;
                  hjetptspectrumRebUnc->SetBinContent(mm,err);
                  hjetptspectrumRebUnc->SetBinError(mm,0);
      }
      
      
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

	ofile->cd();
	hJetPtRebin->Write();
	hJetPtRebinScale->Write();
	hjetptspectrumRebUnc->Write();
	ofile->Close();

  SaveCanvas(cmass, Form("%s/invMass_signal",OUTDIRECTORYPLOTS.Data()));
  SaveCanvas(cmassbkg, Form("%s/invMass_bkg",OUTDIRECTORYPLOTS.Data()));
  SaveCanvas(cSB, Form("%s/SigOverBkg",OUTDIRECTORYPLOTS.Data()));
  SaveCanvas(cJetPt, Form("%s/JetSpectra_eff",OUTDIRECTORYPLOTS.Data()));
  SaveCanvas(cJetPtReb, Form("%s/JetSpectra",OUTDIRECTORYPLOTS.Data()));
  SaveCanvas(cJetPtRebUnc, Form("%s/JetSpectraUnc",OUTDIRECTORYPLOTS.Data()));
  SaveCanvas(cJetPtRebUncLog, Form("%s/JetSpectraUncLog",OUTDIRECTORYPLOTS.Data()));
  TString canvasname[] = { "POWHEGspectrum", "jetSpectra_DptBins0", "jetSpectra_DptBins", "jetSpectraSignal_DptBins","jetSpectra_DptBins0_reb","jetSpectra_DptBins_reb","jetSpectra_DptBins_reb_eff","jetSpectraSignal_DptBins_reb_eff"}
  for(int i=0;i<ncanvas;i++) {
    SaveCanvas(cpt[i], Form("%s/%s",OUTDIRECTORYPLOTS.Data(),canvasname[i].Data()));
  }

  return hJetPtRebin;
}

void getPowhegDmesonSpectra() {

  fHPowhegDPrompt = (TH1D*)getPowhegDSpectra(fFilePowhegPrompt);  // get D x-section from POWHEG
  if(!fHPowhegDPrompt) { cout << "\n!!!Prompt D POWHEG spectrum not extracted !!!" << endl; return; }
  setHistoDetails(fHPowhegDPrompt,2,20,1.2);
  fHPowhegDNonPrompt = (TH1D*)getPowhegDSpectra(fFilePowhegNonPrompt,0);  // get non-prompt D x-section from POWHEG
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
          
/*	
            for(int k=1; k<hMC->GetNbinsX()+1;k++) {
              //if(hMC->GetBinCenter(k)>fHardBinsMax[j-2]*4){
              if(hMC->GetBinContent(k)<3){
                if(hMC->GetBinContent(k)){
                  hMC->SetBinContent(k,0);
                  hMC->SetBinError(k,0);
                }
              }
            }
            for(int k=1; k<hreco->GetNbinsX()+1;k++) {
              //if(hreco->GetBinCenter(k)>fHardBinsMax[j-2]*4){
              if(hreco->GetBinContent(k)<3){
                if(hreco->GetBinContent(k)){
                  hreco->SetBinContent(k,0);
                  hreco->SetBinError(k,0);
                }
              }
            }
*/
				
		
		
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

TH1* getPowhegDSpectra(TString data, Bool_t isPrompt=1, Bool_t recal=0){

  TString outname;
  if(isPrompt) outname = "Powheg_DpT_prompt.root";
  else outname = "Powheg_DpT_nonprompt.root";
  if(!recal) {
    TFile *file = new TFile(outname,"read");
    if(!file->IsZombie()){
      cout << "!!!!We are here" << endl;
      TH1D *hout = (TH1D*)file->Get("hPt_reb");
      return hout;
    }
  }

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

  TH1D *hh_tmp = (TH1D*)hPt->Clone("hh_tmp");
  TH1D *hPt_reb = (TH1D*)hh_tmp->Rebin(fPtbinsDN,"hPt_reb",fPtbinsDA);

  TFile *outfile = new TFile(outname,"RECREATE");
  hPt->Write();
  hPt_reb->Write();
  outfile->Close();

  return hPt_reb;

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

double GetContent(TH1 *hh, double bin){
    return hh->GetBinContent(hh->GetXaxis()->FindBin(bin));
}


void setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Size_t size = 0.9, Width_t width=2){

    h->SetMarkerStyle(Mstyle);
    h->SetMarkerColor(color);
    h->SetMarkerSize(size);
    h->SetLineColor(color);
    h->SetLineWidth(width);
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


