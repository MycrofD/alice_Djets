/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
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

 // Extraction of raw (or efficiency-corrected) D-jet pT spectrum
 // inv. mass method: side-band
//
//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------


#include "signalExtraction.h"

void signalExtraction_SB(
  TString dataFile = "/home/basia/Work/alice/analysis/pp13tev/outData/losser_15Nov/AnalysisResults",
  TString lhcprod = "LHC16kl", // if one file: e.g. LHC16k, LHC16kl ... ; for more than one file: LHC16
  bool isMoreFiles = 0, TString prod = "kl",    // for more than 1 file, for one file leave it empty
  bool isEff = 0, TString efffile = "../efficiency/DjetEff_prompt.root",
  bool isRef = 0, TString refFile = "",
  bool save = 1,
  bool postfix = 0, TString listName = "Cut",
  TString out = "signalExtraction")

{

    fUseRefl = isRef;
    if(fUseRefl) fReflFilename = refFile;
    savePlots = save;
    bEff = isEff;
    if(bEff)plotsDir="/plots";
    else plotsDir = "/plotsNoEff";
    TString outdir = OUTDIRECTORY + out;
    gSystem->Exec(Form("mkdir %s",outdir.Data()));
    gSystem->Exec(Form("mkdir %s%s",outdir.Data(),plotsDir.Data()));
    setStyle();

    if(!isMoreFiles) prod="";
    int nFiles = (int)prod.Length();

    TString histName;
    if(fDmesonSpecie) histName = "histosDStarMBN";
    else histName = "histosD0MBN";

    // get analysis output file
    TString datafile;
    TFile *File;
    TDirectoryFile* dir;
    TList *histList;
    THnSparseF *sparse;

    if(!nFiles) {
      datafile = dataFile;
      datafile += "_";
      datafile += lhcprod;
      datafile += ".root";
      File = new TFile(datafile,"read");
      dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");

      for(int i=0;i<ND; i++){
          if(postfix) histList =  (TList*)dir->Get(Form("%s%d%s",histName.Data(),i,listName.Data()));
          else histList =  (TList*)dir->Get(Form("%s%d",histName.Data(),i));
          sparse = (THnSparseF*)histList->FindObject("hsDphiz");
          sparse->GetAxis(0)->SetRangeUser(zmin,zmax);
          sparse->GetAxis(1)->SetRangeUser(jetmin,jetmax);
          if(isEta) sparse->GetAxis(5)->SetRangeUser(-jetEta,jetEta);
          if(i==0) hInvMassptD=(TH3D*)sparse->Projection(3,1,2);
          else hInvMassptD->Add((TH3D*)sparse->Projection(3,1,2));
      }
    }
    else {
      for (int j=0;j<nFiles;j++){
          datafile = dataFile;
          datafile += "_";
          datafile += lhcprod;
          datafile += prod.Data()[j];
          datafile += ".root";
          File = new TFile(datafile,"read");
          dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");

          for(int i=0;i<ND; i++){
              if(postfix) histList =  (TList*)dir->Get(Form("%s%d%s",histName.Data(),i,listName.Data()));
              else histList =  (TList*)dir->Get(Form("%s%d",histName.Data(),i));
              sparse = (THnSparseF*)histList->FindObject("hsDphiz");
              sparse->GetAxis(0)->SetRangeUser(zmin,zmax);
              sparse->GetAxis(1)->SetRangeUser(jetmin,jetmax);
              if(isEta) sparse->GetAxis(5)->SetRangeUser(-jetEta,jetEta);
              if(j==0 && i==0) hInvMassptD=(TH3D*)sparse->Projection(3,1,2);
              else hInvMassptD->Add((TH3D*)sparse->Projection(3,1,2));
          }
      }
    }

    efficiency = new double[fptbinsDN];
    if(bEff){
        TFile *FileEff = new TFile(efffile.Data(),"read");
        TH1F *hEff = (TH1F*)FileEff->Get("hEff_reb");
        for(int i=0;i<fptbinsDN;i++){
            double pt = (fptbinsDA[i]+fptbinsDA[i+1]) / 2.;
            efficiency[i] = hEff->GetBinContent(hEff->GetXaxis()->FindBin(pt));
        }
    }
    else {
        for(int i=0;i<fptbinsDN;i++){
            efficiency[i] = 1;
        }
    }

    // --------------------------------------------------------
    // fit inv. mass in D pT bins and get raw jet pT spectra in the signal and side-bands regions
     Bool_t okSignalExt = rawJetSpectra(outdir,lhcprod+prod);
     if(!okSignalExt) { std::cout << "!!!!!! Something wrong in the raw signal extraction !!!!" << endl; return; }
    // --------------------------------------------------------

    // --------------------------------------------------------
   //------------------- draw output histos -------------------

    if(savePlots){
      saveFitParams(outdir,lhcprod+prod);
      saveSpectraPlots(outdir,lhcprod+prod);
    }

    // --------------------------------------------------------

    // --------------------------------------------------------
    // ----------- write to output file
    TFile *ofile = new TFile(Form("%s/JetPtSpectra_SB_%s.root",outdir.Data(), bEff ? "eff" : "noEff"),"RECREATE");
    hmean->Write();
    hsigma->Write();
    hsign->Write();
    hsb->Write();
    hSignal->Write();
    hrelErr->Write();

    hjetptspectrum->Write();
    hjetptspectrumReb->Write();
    hjetptspectrumRebUnc->Write();

    for(int i=0; i<fptbinsDN; i++){
        if(hjetpt[i]) hjetpt[i]->Write();
        if(hjetpt_s[i]) hjetpt_s[i]->Write();
        if(hjetptsub[i]) hjetptsub[i]->Write();
        if(hmass[i]) hmass[i]->Write();
        if(hmass_l[i]) hmass_l[i]->Write();
        if(hmass_u[i]) hmass_u[i]->Write();
        if(hmass_c[i]) hmass_c[i]->Write();
        if(fullfit[i]) fullfit[i]->Write();
        if(massfit[i]) massfit[i]->Write();
        if(bkgfit[i]) bkgfit[i]->Write();
        if(bkgRfit[i] && fUseRefl && fDmesonSpecie == 0) bkgRfit[i]->Write();
        if(hjetptcorr[i]) hjetptcorr[i]->Write();
    }

    ofile->Close();
    // --------------------------------------------------------

}


Bool_t rawJetSpectra(TString outdir, TString prod){

    hmean = new TH1F("hmean","hmean",fptbinsDN,fptbinsDA);
    hsigma = new TH1F("hsigma","hsigma",fptbinsDN,fptbinsDA);
    hrelErr = new TH1F("hrelErr","hrelErr",fptbinsDN,fptbinsDA);
    hsign = new TH1F("hsign","hsign",fptbinsDN,fptbinsDA);
    hsb = new TH1F("hsb","hsb",fptbinsDN,fptbinsDA);
    hSignal = new TH1F("hSignal","hSignal",fptbinsDN,fptbinsDA);
    hSignal->Sumw2();
    hReflRS = new TH1F("hReflRS","hReflRS",fptbinsDN,fptbinsDA);

    if(fDmesonSpecie) hInvMassptD->GetXaxis()->SetTitle(Form("m(%s)(GeV/c^{2})","K#pi#pi-K#pi"));
    else hInvMassptD->GetXaxis()->SetTitle(Form("m(%s)(GeV/c^{2})","K#pi"));
    hInvMassptD->GetXaxis()->SetTitleSize(0.06);
    hInvMassptD->GetXaxis()->SetTitleOffset(0.8);
    hInvMassptD->GetYaxis()->SetTitle("p_{T}^{jet}");
    hInvMassptD->SetTitle();

    TPaveText *pvProd = new TPaveText(0.8,0.25,0.98,0.3,"brNDC");
    pvProd->SetFillStyle(0);
    pvProd->SetBorderSize(0);
    pvProd->AddText(Form("%s",prod.Data()));

    TPaveText *pvCuts = new TPaveText(0.8,0.3,0.98,0.35,"brNDC");
    pvCuts->SetFillStyle(0);
    pvCuts->SetBorderSize(0);
    pvCuts->AddText(Form("%s",outdir.Data()));

    int xnx = 3, xny=4;
    if(fptbinsDN>4 && fptbinsDN<7) { xnx = 2; xny=3; }
    else if(fptbinsDN>6 && fptbinsDN<10) { xnx = 3; xny=3; }
    else if(fptbinsDN>9 && fptbinsDN<13) { xnx = 3; xny=4; }
    else { xnx = 4; xny=4; }

    TCanvas *c2 = new TCanvas("c2","c2",1200,1200);
    c2->Divide(xnx,xny);
    TCanvas *c2jet = new TCanvas("c2jet","c2jet",1200,1200);
    c2jet->Divide(xnx,xny);
    TCanvas *c2jetcorr = new TCanvas("c2jetcorr","c2jetcorr",1200,1200);
    c2jetcorr->Divide(xnx,xny);

    int firstPtBin = 0;
    if(fptbinsDA[0] == 2) firstPtBin = 3;
    else if(fptbinsDA[0] == 3) firstPtBin = 4;
    else if(fptbinsDA[0] == 4) firstPtBin = 5;
    if(!firstPtBin) { std::cout << "==== Wrong first value of the D pT (should be 2,3 or 4) === \n"; return kFALSE; }
    Float_t RS = 0;

    for(int i=0; i<fptbinsDN; i++){

        TH1F *hh=(TH1F*)hInvMassptD->ProjectionX(Form("hh_%d",i),hInvMassptD->GetYaxis()->FindBin(jetmin), hInvMassptD->GetYaxis()->FindBin(jetmax)-1,hInvMassptD->GetZaxis()->FindBin(fptbinsDA[i]), hInvMassptD->GetZaxis()->FindBin(fptbinsDA[i+1])-1);
        hh->Rebin(fRebinMass);

        hh->GetXaxis()->SetRangeUser(minf,maxf);
        hh->SetTitle(Form("%.1lf < pt^{%s} < %.1lf",fptbinsDA[i],fDmesonS.Data(),fptbinsDA[i+1]));

        TH1F *hmassfit = (TH1F*)hh->Clone("hmassfit");
        if(fDmesonSpecie) hmassfit->SetMaximum(hmassfit->GetMaximum()*1.3);

        float hmin = TMath::Max(minf,hmassfit->GetBinLowEdge(2));
        float hmax = TMath::Min(maxf,hmassfit->GetBinLowEdge(hmassfit->GetNbinsX()));
       // AliHFMassFitter* fitterp=new AliHFMassFitter((TH1F*)hmassfit,hmin,hmax,1,fbkgtype,0);
        AliHFInvMassFitter* fitterp = new AliHFInvMassFitter((TH1F*)hmassfit,hmin,hmax,fbkgtype,0);
        fitterp->SetInitialGaussianMean(fDmass);
        fitterp->SetInitialGaussianSigma(fDsigma);

        if(fUseRefl && fDmesonSpecie == 0) {
          SetReflection(fitterp,i+firstPtBin,hmin,hmax,RS);
        }

        fitterp->MassFitter(kFALSE);
        //fitterp->PrintFunctions();

        TH1F* h=fitterp->GetHistoClone();
        massfit[i] = fitterp->GetMassFunc();
        massfit[i]->SetRange(hmin,hmax);
        massfit[i]->SetLineColor(4);
        fullfit[i] = h->GetFunction("funcmass");
        if(fullfit[i]) fullfit[i]->SetName(Form("fullfit_%d",i));
        hmass[i] = (TH1F*)h->Clone(Form("hmass_%d",i));
        hmass[i]->SetName(Form("hmass_%d",i));
        hmass[i]->GetYaxis()->SetTitle("Entries");
        hmass[i]->GetYaxis()->SetTitleOffset(1.4);
        hmass[i]->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
        bkgfit[i] = fitterp->GetBackgroundRecalcFunc();
        bkgfit[i]->SetRange(hmin,hmax);
        bkgfit[i]->SetLineColor(2);
        bkgfit[i]->SetName("bkgFit");

        if(fUseRefl && fDmesonSpecie == 0) {
          bkgRfit[i] = fitterp->GetBkgPlusReflFunc();
          bkgRfit[i]->SetName("bkgFitWRef");
          bkgRfit[i]->SetRange(hmin,hmax);
          bkgRfit[i]->SetLineColor(15);
          hReflRS->SetBinContent(i+1,RS);
          hReflRS->SetBinError(i+1,0);
        }

        TVirtualPad *pad = (TVirtualPad*)c2->GetPad(i+1);
        fitterp->DrawHere(pad,3,0);

        float Dsigma = 0, Dmean = 0, DmeanUnc = 0, DsigmaUnc = 0;

        if(!fullfit[i]) { std::cout << "======= Fit failed for bin: " << i << endl; continue; }

        Dsigma = fitterp->GetSigma();
        DsigmaUnc  = fitterp->GetSigmaUncertainty();
        Dmean = fitterp->GetMean();
        DmeanUnc = fitterp->GetMeanUncertainty();

        float signal_c_min = Dmean-fsigmaSignal*Dsigma;
        float signal_c_max = Dmean+fsigmaSignal*Dsigma;
        float signal_l_min = Dmean+fsigmaBkg[0]*Dsigma;
        float signal_l_max = Dmean+fsigmaBkg[1]*Dsigma;
        float signal_u_min = Dmean+fsigmaBkg[2]*Dsigma;
        float signal_u_max = Dmean+fsigmaBkg[3]*Dsigma;

        int binmin = hmass[i]->GetXaxis()->FindBin(signal_c_min);
        int binmax = hmass[i]->GetXaxis()->FindBin(signal_c_max);
        double binwidth = hmass[i]->GetXaxis()->GetBinWidth(1)*0.5;

        Double_t s=0,serr=0,srelerr=0,b=0,berr=0,signf=0,signferr=0,sob=0,soberr=0;
        Double_t br=0; //background with reflections

        Double_t min = hmass[i]->GetXaxis()->GetBinCenter(binmin)-binwidth;
        Double_t max = hmass[i]->GetXaxis()->GetBinCenter(binmax-1)+binwidth;
        fitterp->Signal(fsigmaSignal,s,serr);
        fitterp->Background(min, max ,b,berr);
        //fitterp->Background(fsigmaSignal,b,berr);
        fitterp->Significance(fsigmaSignal,signf,signferr);
        if(s) srelerr = serr/s;
        if(b) sob = s/b;
        else sob = s;
        if(b && berr) soberr = TMath::Sqrt((serr/b)*(serr/b) + (s/b/b*berr)*(s/b/b*berr));
        else soberr = serr;
        if(fUseRefl && fDmesonSpecie == 0) {
          br = bkgRfit[i]->Integral(min,max)/(Double_t)hmass[i]->GetBinWidth(1);
        }
        else br = b;

        TPaveText *pvSig;
        if(fDmesonSpecie) pvSig = new TPaveText(0.55,0.47,0.95,0.75,"brNDC");
        else pvSig = new TPaveText(0.15,0.55,0.47,0.9,"brNDC");
        pvSig->SetFillStyle(0);
        pvSig->SetBorderSize(0);
        Bool_t twodigits=kTRUE;
        if(soberr*100. > 35.) twodigits=kFALSE;
        //if(twodigits) pv->AddText(Form("S/B (3#sigma) = (%.2f #pm %.2f)", sob,soberr));
        //else pv->AddText(Form("S/B (3#sigma) = (%.1f #pm %.1f)", sob,soberr));
        if(twodigits) pvSig->AddText(Form("S (3#sigma) = %.2f #pm %.2f", s,serr));
        else pvSig->AddText(Form("S (3#sigma) = %.1f #pm %.1f", s,serr));
        if(twodigits) pvSig->AddText(Form("B (3#sigma) = %.2f #pm %.2f", b,berr));
        else pvSig->AddText(Form("B (3#sigma) = %.1f #pm %.1f", b,berr));
        pvSig->AddText(Form("Signif.(3#sigma) = %.1f #pm %.1f", signf,signferr));
        pvSig->AddText(Form("S/B(3#sigma) = %.2f #pm %.2f", sob,soberr));
        if(fUseRefl && fDmesonSpecie == 0) pvSig->AddText(Form("R/S = %.2f", RS));
        pvSig->Draw("same");
        //if(isdetails) pvProd->Draw("same");
        //if(isdetails) pvCuts->Draw("same");

        TPaveText *pv;
        if(fDmesonSpecie) pv = new TPaveText(0.55,0.77,0.95,0.9,"brNDC");
        else pv = new TPaveText(0.57,0.77,0.95,0.9,"brNDC");
        pv->SetFillStyle(0);
        pv->SetBorderSize(0);
        if(fDmesonSpecie) pv->AddText(Form("#mu = (%.2f #pm %.2f) MeV/#it{c}^{2}", Dmean*1000,DmeanUnc*1000));
        else pv->AddText(Form("#mu = (%.2f #pm %.2f) GeV/#it{c}^{2}", Dmean,DmeanUnc));
        pv->AddText(Form("#sigma = (%.2f #pm %.2f) MeV/#it{c}^{2}", Dsigma*1000,DsigmaUnc*1000));
        pv->Draw("same");

        // ---------------- fitting results
        if(fDmesonSpecie) {
          hmean->SetBinContent(i+1,Dmean*1000);
          hmean->SetBinError(i+1,DmeanUnc*1000);
        }
        else {
          hmean->SetBinContent(i+1,Dmean);
          hmean->SetBinError(i+1,DmeanUnc);
        }
        hsigma->SetBinContent(i+1,Dsigma*1000);
        hsigma->SetBinError(i+1,DsigmaUnc*1000);

        hrelErr->SetBinContent(i+1,srelerr);
        hsign->SetBinContent(i+1,signf);
        hsign->SetBinError(i+1,signferr);
        hsb->SetBinContent(i+1,sob);
        hsb->SetBinError(i+1,soberr);
        hSignal->SetBinContent(i+1,s);
        hSignal->SetBinError(i+1,serr);

        // ---------------- side-band drawing
        hmass_l[i] = (TH1F*)hmass[i]->Clone("hmass_l");
        hmass_l[i]->GetXaxis()->SetRangeUser(signal_l_min,signal_l_max);
        hmass_l[i]->SetName(Form("hmass_l_%d",i));
        hmass_u[i] = (TH1F*)hmass[i]->Clone("hmass_u");
        hmass_u[i]->GetXaxis()->SetRangeUser(signal_u_min,signal_u_max);
        hmass_u[i]->SetName(Form("hmass_u_%d",i));
        hmass_c[i] = (TH1F*)hmass[i]->Clone("hmass_c");
        hmass_c[i]->GetXaxis()->SetRangeUser(signal_c_min,signal_c_max);
        hmass_c[i]->SetName(Form("hmass_c_%d",i));

        hmass_l[i]->SetFillColor(kBlue+2);
        hmass_u[i]->SetFillColor(kBlue+2);
        hmass_c[i]->SetFillColor(kRed+2);
        hmass_l[i]->SetFillStyle(3004);
        hmass_u[i]->SetFillStyle(3004);
        hmass_c[i]->SetFillStyle(3005);

        hmass_l[i]->Draw("hsame");
        hmass_u[i]->Draw("hsame");
        hmass_c[i]->Draw("hsame");

        //----------------------------------
        //-------- jet pt spectrum - signal
        hjetpt[i]=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_%d",i),hInvMassptD->GetXaxis()->FindBin(signal_c_min), hInvMassptD->GetXaxis()->FindBin(signal_c_max)-1,hInvMassptD->GetZaxis()->FindBin(fptbinsDA[i]), hInvMassptD->GetZaxis()->FindBin(fptbinsDA[i+1])-1);
        hjetpt[i]->GetXaxis()->SetRangeUser(jetplotmin,jetplotmax);
        hjetpt[i]->SetTitle(Form("%.1lf < pt^{%s} < %.1lf",fptbinsDA[i],fDmesonS.Data(),fptbinsDA[i+1]));
        hjetpt[i]->SetMarkerColor(kRed+2);
        hjetpt[i]->SetLineColor(kRed+2);

        //------ jet pt spectrum - side bands
        TH1F* hjetpt_s1=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_s1%d",i),hInvMassptD->GetXaxis()->FindBin(signal_l_min), hInvMassptD->GetXaxis()->FindBin(signal_l_max)-1,hInvMassptD->GetZaxis()->FindBin(fptbinsDA[i]), hInvMassptD->GetZaxis()->FindBin(fptbinsDA[i+1])-1);
        hjetpt_s1->SetTitle(Form("%.1lf < pt^{%s} < %.1lf",fptbinsDA[i],fDmesonS.Data(),fptbinsDA[i+1]));
        hjetpt_s1->SetMarkerColor(kBlue+2);
        hjetpt_s1->SetLineColor(kBlue+2);
        TH1F* hjetpt_s2=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_s2%d",i),hInvMassptD->GetXaxis()->FindBin(signal_u_min), hInvMassptD->GetXaxis()->FindBin(signal_u_max)-1,hInvMassptD->GetZaxis()->FindBin(fptbinsDA[i]), hInvMassptD->GetZaxis()->FindBin(fptbinsDA[i+1])-1);
        hjetpt_s2->SetTitle(Form("%.1lf < pt^{%s} < %.1lf",fptbinsDA[i],fDmesonS.Data(),fptbinsDA[i+1]));
        hjetpt_s2->SetMarkerColor(kBlue+2);
        hjetpt_s2->SetLineColor(kBlue+2);
        hjetpt_s[i] = (TH1F*)hjetpt_s1->Clone(Form("hjetpt_s_%d",i));
        hjetpt_s[i]->Add(hjetpt_s2);

        // scale background from side bands to the background under the peak
        double scaling = br / hjetpt_s[i]->Integral();
        hjetpt_s[i]->Scale(scaling);

        c2jet->cd(i+1);
        gPad->SetLogy();

        hjetpt[i]->GetXaxis()->SetRangeUser(jetplotmin,jetplotmax);
        hjetpt_s[i]->GetXaxis()->SetRangeUser(jetplotmin,jetplotmax);
        hjetpt[i]->GetXaxis()->SetTitle("p_{T,ch jet} (GeV/c)");
        hjetpt[i]->GetXaxis()->SetTitleSize(0.05);
        hjetpt[i]->GetXaxis()->SetTitleOffset(0.8);

        hjetpt[i]->Draw("ep");
        hjetpt_s[i]->Draw("epsame");
        TLegend *l1 = new TLegend(0.75,0.6,0.9,0.8);
        l1->AddEntry(hjetpt[i],"signal","l");
        l1->AddEntry(hjetpt_s[i],"SB","l");
        l1->Draw("same");

        //------- subtract background from signal jet
        hjetptsub[i] = (TH1F*)hjetpt[i]->Clone(Form("hjetptsub_%d",i));
        hjetptsub[i]->Add(hjetpt_s[i],-1);
        hjetptsub[i]->SetMarkerColor(kGreen+3);
        hjetptsub[i]->SetLineColor(kGreen+3);
        hjetptsub[i]->Draw("epsame");
        if(isdetails) pvProd->Draw("same");
        if(isdetails) pvCuts->Draw("same");

        TLegend *l2 = new TLegend(0.6,0.75,0.9,0.9);
        l2->AddEntry(hjetptsub[i],"bkg subtracted","l");
        l1->AddEntry(hjetptsub[i],"sig-SB","l");
        //l2->Draw("same");
        if(!hrawjetptspectrum) hrawjetptspectrum = (TH1F*)hjetptsub[i]->Clone("hrawjetptspectrum");
        else hrawjetptspectrum->Add(hjetptsub[i]);

        //------- correct for D* efficiency
        hjetptcorr[i] = (TH1F*)hjetptsub[i]->Clone(Form("hjetptcorr_%d",i));
        if(bEff) hjetptcorr[i]->Scale(1./efficiency[i]); // D efficiency
        hjetptcorr[i]->SetMarkerColor(kBlue+3);
        hjetptcorr[i]->SetLineColor(kBlue+3);
        hjetptcorr[i]->GetXaxis()->SetTitleSize(0.04);
        hjetptcorr[i]->GetXaxis()->SetTitleOffset(0.9);
        //hjetptcorr->GetXaxis()->SetRangeUser(jetplotmin,jetplotmax);

        c2jetcorr->cd(i+1);
        hjetptcorr[i]->Draw("ep");
        if(isdetails) pvProd->Draw("same");
        if(isdetails) pvCuts->Draw("same");
        TLegend *l3 = new TLegend(0.6,0.7,0.9,0.9);
        l3->AddEntry(hjetptcorr[i],"#splitline{corrected jet p_{T}}{spectrum}","l");
        l3->Draw("same");

        //----- add corrected jet pt spectrum distributions in each D pt bin into a one distribution
        if(!hjetptspectrum) hjetptspectrum = (TH1F*)hjetptcorr[i]->Clone("hjetptspectrum");
        else hjetptspectrum->Add(hjetptcorr[i]);

    }

    if(savePlots) SaveCanvas(c2,outdir+plotsDir+"/invMass_"+prod);
    if(savePlots) SaveCanvas(c2jet,outdir+plotsDir+"/jetRawSpectrum"+prod);

    return kTRUE;

}

Bool_t SetReflection(AliHFInvMassFitter* &fitter, Int_t iBin, Float_t fLeftFitRange, Float_t fRightFitRange, Float_t &RS) {

  TFile *fileRefl = TFile::Open(fReflFilename.Data());
  if(!fileRefl){
    std::cout << "File " << fReflFilename << " (reflection templates) cannot be opened! check your file path!"; return kFALSE;
  }

  TString fHistnameRefl = "histRflFittedDoubleGaus_ptBin";
  TString fHistnameSign = "histSgn_";
  TH1F *histRefl = (TH1F*)fileRefl->Get(Form("%s%d",fHistnameRefl.Data(),iBin));
  TH1F *histSign = (TH1F*)fileRefl->Get(Form("%s%d",fHistnameSign.Data(),iBin));
  if(!histRefl || !histSign){
    std::cout << "Error in loading the template/signal histrograms! Exiting..." << endl; return kFALSE;
  }

  fitter->SetTemplateReflections(histRefl,"template",fLeftFitRange,fRightFitRange);
  Double_t RoverS = histRefl->Integral(histRefl->FindBin(fLeftFitRange),histRefl->FindBin(fRightFitRange))/histSign->Integral(histSign->FindBin(fLeftFitRange),histSign->FindBin(fRightFitRange));
  printf("R/S ratio in fit range for bin %d = %1.3f\n",iBin,RoverS);
  fitter->SetFixReflOverS(RoverS);

  RS = (Float_t)RoverS;
  return kTRUE;

}

void  saveSpectraPlots(TString outdir,TString prod){

      TPaveText *pvProd = new TPaveText(0.75,0.65,0.9,0.7,"brNDC");
      pvProd->SetFillStyle(0);
      pvProd->SetBorderSize(0);
      pvProd->AddText(Form("%s",prod.Data()));

      TPaveText *pvCuts = new TPaveText(0.75,0.6,0.9,0.65,"brNDC");
      pvCuts->SetFillStyle(0);
      pvCuts->SetBorderSize(0);
      pvCuts->AddText(Form("%s",outdir.Data()));

       // before eff. correction
      TCanvas *cRawSpectrum = new TCanvas("cRawSpectrum","cRawSpectrum",800,600);
      setHistoDetails(hrawjetptspectrum,0,kBlue+1,20,1.2);
      hrawjetptspectrum->GetXaxis()->SetTitle("p_{T}^{ch,jet} (GeV/c)");
      hrawjetptspectrum->GetYaxis()->SetTitle("Raw dN/dp_{T}");
      hrawjetptspectrum->GetXaxis()->SetRangeUser(jetplotmin,jetplotmax);
      hrawjetptspectrum->Draw();
      if(isdetails) pvProd->Draw("same");
      if(isdetails) pvCuts->Draw("same");

      TPaveText *pv3 = new TPaveText(0.68,0.65,0.95,0.85,"brNDC");
      pv3->SetFillStyle(0);
      pv3->SetBorderSize(0);
      if(!isptcut) pv3->AddText(Form("charged jets,  R=0.%d",Rpar));
      else pv3->AddText(Form("#splitline{charged jets, R=0.%d}{p_{T,%s}>%dGeV/c}",Rpar,fDmesonS.Data(),(Int_t)fptbinsDA[0]));

      hjetptspectrum->GetYaxis()->SetTitle("dN/dp_{T}");
      hjetptspectrum->GetXaxis()->SetRangeUser(jetplotmin,jetplotmax);

      TCanvas *cSpectrum = new TCanvas("cSpectrum","cSpectrum",800,600);
      cSpectrum->SetLogy();
      setHistoDetails(hjetptspectrum,0,kRed,20,1.2);
      hjetptspectrum->GetXaxis()->SetTitle("p_{T,ch jet} (GeV/c)");
      hjetptspectrum->Draw();
      pv3->Draw("same");
      if(isdetails) pvProd->Draw("same");
      if(isdetails) pvCuts->Draw("same");
      if(savePlots) SaveCanvas(cSpectrum,outdir+plotsDir+"/jetPtSpectrum_SB_"+prod);

      TH1F *hjetptspectrumReb_tmp = (TH1F*)hjetptspectrum->Clone("hjetptspectrumReb_tmp");
      hjetptspectrumReb = (TH1F*)hjetptspectrumReb_tmp->Rebin(fptbinsJetMeasN,"hjetptspectrumReb",fptbinsJetMeasA);
      TH1F* hjetptspectrumReb2 = (TH1F*)hjetptspectrumReb_tmp->Rebin(fptbinsJetMeasN,"hjetptspectrumReb",fptbinsJetMeasA);
      setHistoDetails(hjetptspectrumReb,0,kBlue,20,1.2); // with bin width scaling
      setHistoDetails(hjetptspectrumReb2,1,kBlue,20,1.2); // with bin width scaling
      hjetptspectrumReb->GetXaxis()->SetTitle("p_{T,ch jet} (GeV/c)");
      hjetptspectrumReb2->GetXaxis()->SetTitle("p_{T,ch jet} (GeV/c)");
      TCanvas *cSpectrumRebin = new TCanvas("cSpectrumRebin","cSpectrumRebin",800,600);
      cSpectrumRebin->SetLogy();
      hjetptspectrumReb2->Draw();
      if(isdetails) pvProd->Draw("same");
      if(isdetails) pvCuts->Draw("same");
      pv3->Draw("same");

      SaveCanvas(cSpectrumRebin,outdir+plotsDir+"/jetPtSpectrum_SB_Rebin_"+prod);

      hjetptspectrumRebUnc = (TH1F*)hjetptspectrumReb->Clone("hjetptspectrumRebUnc");
      hjetptspectrumRebUnc->GetYaxis()->SetTitle("Rel. unc.");

      for(int j=1; j<=hjetptspectrumReb->GetNbinsX();j++){
                  double err;
                  if(hjetptspectrumReb->GetBinContent(j)) err = hjetptspectrumReb->GetBinError(j)/hjetptspectrumReb->GetBinContent(j);
                  else err = 0;
                  hjetptspectrumRebUnc->SetBinContent(j,err);
                  hjetptspectrumRebUnc->SetBinError(j,0);
      }

      TCanvas *cSpectrumRebinUnc = new TCanvas("cSpectrumRebinUnc","cSpectrumRebinUnc",800,600);

      hjetptspectrumRebUnc->GetXaxis()->SetRangeUser(jetplotmin,jetplotmax);
      hjetptspectrumRebUnc->Draw();
      if(isdetails) pvProd->Draw("same");
      if(isdetails) pvCuts->Draw("same");
      pv3->Draw("same");

      SaveCanvas(cSpectrumRebinUnc,outdir+plotsDir+"/jetPtSpectrumUnc_SB_Rebin_"+prod);

}


void  saveFitParams(TString outdir,TString prod){

    TPaveText *pvProd = new TPaveText(0.78,0.75,0.9,0.8,"brNDC");
    pvProd->SetFillStyle(0);
    pvProd->SetBorderSize(0);
    pvProd->AddText(Form("%s",prod.Data()));

    TPaveText *pvCuts = new TPaveText(0.78,0.8,0.9,0.85,"brNDC");
    pvCuts->SetFillStyle(0);
    pvCuts->SetBorderSize(0);
    pvCuts->AddText(Form("%s",outdir.Data()));

    setHistoDetails(hmean,0,2,20);
    if(fDmesonSpecie) hmean->GetYaxis()->SetTitle("signal #mu (MeV/c^{2})");
    else hmean->GetYaxis()->SetTitle("signal #mu (GeV/c^{2})");
    setHistoDetails(hsigma,0,4,20);
    hsigma->GetYaxis()->SetTitle("signal #sigma (MeV/c^{2})");
    setHistoDetails(hsign,0,4,20);
    hsign->GetYaxis()->SetTitle("significance");
    setHistoDetails(hsb,0,8,20);
    hsb->GetYaxis()->SetTitle("S/B");
    setHistoDetails(hrelErr,0,6,20);
    hrelErr->GetYaxis()->SetTitle("rel. unc.");

    hmean->GetYaxis()->SetTitleOffset(1.8);
    hmean->GetYaxis()->SetTitleSize(0.04);
    hmean->GetYaxis()->SetLabelSize(0.03);

    hsigma->GetYaxis()->SetTitleOffset(1.4);
    hsigma->GetYaxis()->SetTitleSize(0.04);
    hsigma->GetYaxis()->SetLabelSize(0.03);

    if(fDmesonSpecie) {
      hmean->GetYaxis()->SetRangeUser(144.8,146);
      hsigma->GetYaxis()->SetRangeUser(0.35,0.8);
    }
    else {
      hmean->GetYaxis()->SetRangeUser(1.855,1.885);
      hsigma->GetYaxis()->SetRangeUser(9,20);
    }

    hSignal->SetName("hSignal");
    setHistoDetails(hSignal,0,2,20);
    hSignal->GetYaxis()->SetTitle("yield");
    hSignal->GetYaxis()->SetTitleOffset(1.2);
    hSignal->GetYaxis()->SetTitleSize(0.04);
    hSignal->GetYaxis()->SetLabelSize(0.03);
    hsign->GetYaxis()->SetTitleOffset(1.2);
    hsign->GetYaxis()->SetTitleSize(0.04);
    hsign->GetYaxis()->SetLabelSize(0.03);
    hrelErr->GetYaxis()->SetTitleOffset(1.2);
    hrelErr->GetYaxis()->SetTitleSize(0.04);
    hrelErr->GetYaxis()->SetLabelSize(0.03);


    double mean = 145.421;
    TLine *lm = new TLine(fptbinsDA[0], mean, fptbinsDA[fptbinsDN], mean);
    lm->SetLineStyle(2);

    TCanvas *cMassFit = new TCanvas("cMassFit","cMassFit",1200,600);
    cMassFit->Divide(2,1);
    cMassFit->cd(1);
    gPad->SetLeftMargin(0.15);
    hmean->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    lm->Draw("same");
    cMassFit->cd(2);
    hsigma->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");

     if(savePlots)  SaveCanvas(cMassFit,outdir+plotsDir+"/gaussianParams_"+prod);

    TCanvas *cSignal = new TCanvas("cSignal","cSignal",1200,1200);
    cSignal->Divide(2,2);

    cSignal->cd(1);
    hSignal->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    cSignal->cd(2);
    hsign->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    cSignal->cd(3);
    hrelErr->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    cSignal->cd(4);
    hsb->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");

    if(savePlots) SaveCanvas(cSignal,outdir+plotsDir+"/signalParams_"+prod);

    if(hReflRS && fUseRefl && fDmesonSpecie == 0) {
      setHistoDetails(hReflRS,0,kGreen+2,20);
      hReflRS->GetYaxis()->SetTitle("R/S");
      hReflRS->GetYaxis()->SetTitleOffset(1.2);
      hReflRS->GetYaxis()->SetTitleSize(0.04);
      hReflRS->GetYaxis()->SetLabelSize(0.03);
      //hReflRS->GetYaxis()->SetRangeUser(0,0.4);
      TCanvas *cRS = new TCanvas("cRS","cRS",800,600);
      cRS->cd();
      hReflRS->Draw();
      if(savePlots) SaveCanvas(cRS,outdir+plotsDir+"/RefOverS_"+prod);
    }

}
void setHistoDetails(TH1 *h, int scale, Color_t color, Style_t Mstyle, Size_t size = 0.9, Width_t width=2){

    if(scale)h->Scale(1,"width");
    h->SetMarkerStyle(Mstyle);
    h->SetMarkerColor(color);
    h->SetMarkerSize(size);
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetTitle(0);
    h->GetXaxis()->SetTitle(Form("p_{T,%s}(GeV/c)",fDmesonS.Data()));

    return;
}

void SaveCanvas(TCanvas *c, TString name = "tmp"){

    c->SaveAs(Form("%s_pTD%d.png",name.Data(),(int)fptbinsDA[0]));
    c->SaveAs(Form("%s_pTD%d.pdf",name.Data(),(int)fptbinsDA[0]));
}

void setStyle(){

    gStyle->SetOptStat(000);
    gStyle->SetLegendFont(42);
    gStyle->SetTextFont(22) ;
    //gStyle->SetLegendTextSize(0.05);
    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetPadRightMargin(0.02);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetTitleOffset(1.,"x");
    gStyle->SetTitleOffset(0.9.,"y");
    gStyle->SetTitleSize(0.04,"xyz");
    gStyle->SetLabelSize(0.03,"xyz");
    //style();

}
