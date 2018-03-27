
#include <string>
#include <sstream>
#include <iostream>

#include "config.h"
#include "sys.h"


    Int_t colors[] = {2,4,6,1,kOrange-1,kGray+1,kCyan+1,kMagenta+2,kGreen+3,kViolet+5,kYellow+2,8};
    Int_t markers[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};
    Int_t linestyle[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};

    const int nFiles = 8;
    TString desc[nFiles] = { "central", "looser 1", "looser 2", "looser 3", "looser 4", "tighter 1", "tighter 2", "tighter 3" }

    //const int nFiles = 5;
    //TString desc[nFiles] = { "central", "looser 1", "looser 2", "tighter 1", "tighter 2" }

    TString dirBase;
    TString dirSignal[nFiles] = { "CutVarBase","CutL0","CutL1","CutL2","CutL3","CutT0","CutT1","CutT2" };
  //  TString dirSignal[nFiles] = { "CutVarBase","CutL0","CutL1","CutT0","CutT2" };
    TString inDir[nFiles];
    TString outDir;

TString outPlotName = "CutVariationSyst_";


void cutsSystematics(int reg = 3, TString indirbase = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_", TString input = "Default_jetMeas3_50_jetTrue3_50", bool isRaw = 0, bool unf = 1, bool isChain = 0, int measmin=3, int measmax=50, int truemin=5, int truemax=50 )
{


  dirBase = indirbase;
  if(!unf) reg=0;

  if(!isChain) {
    if (truemin == 3){
      nJetBins = 9;
      ptJetbins = new double[nJetBins+1];
      for(int i=0;i<nJetBins+1;i++) ptJetbins[i] = ptbins[i];
    }
    else if (truemin == 4){
      nJetBins = 8;
      ptJetbins = new double[nJetBins+1];
      for(int i=0;i<nJetBins+1;i++) ptJetbins[i] = ptbins[i+1];
    }
    else if (truemin == 5){
      nJetBins = 7;
      ptJetbins = new double[nJetBins+1];
      for(int i=0;i<nJetBins+1;i++) ptJetbins[i] = ptbins[i+2];
    }
    else {
      cout << "WRONG true minimum pT !!!!" << endl;
      return;
    }
  }
  else {
    nJetBins = fptbinsJetTrueN;
    ptJetbins = new double[nJetBins+1];
    for(int i=0;i<nJetBins+1;i++) ptJetbins[i] = fptbinsJetTrueA[i];
  }


  for(int i=0; i<nFiles; i++) {
    inDir[i] = dirBase;
    inDir[i] += dirSignal[i];
    inDir[i] += "/";
    inDir[i] += input;

  }

  outDir = dirBase;
  outDir += "BaseCuts/";
  outDir += input;
  outDir += "/systematics";

  outPlotName += "reg";
  outPlotName += reg;

  pvEn = new TPaveText(0.25,0.80,0.8,0.85,"brNDC");
  pvEn->SetFillStyle(0);
  pvEn->SetBorderSize(0);
  pvEn->SetTextFont(42);
  pvEn->SetTextSize(0.045);
  pvEn->SetTextAlign(11);
  pvEn->AddText(Form("%s",fSystemS.Data()));

  double shift = 0.35;
  double dshift = 0.05;
  pvJet = new TPaveText(0.12,0.65-shift,0.9,0.7-shift,"brNDC");
  pvJet->SetFillStyle(0);
  pvJet->SetBorderSize(0);
  pvJet->SetTextFont(42);
  pvJet->SetTextSize(0.03);
  pvJet->SetTextAlign(11);
  pvJet->AddText(Form("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.%d",Rpar));

  shift+=dshift;
  pvD = new TPaveText(0.12,0.65-shift,0.9,0.7-shift,"brNDC");
  pvD->SetFillStyle(0);
  pvD->SetBorderSize(0);
  pvD->SetTextFont(42);
  pvD->SetTextSize(0.03);
  pvD->SetTextAlign(11);
  if(fDmesonSpecie) pvD->AddText("with D^{*+} #rightarrow D^{0}#pi^{+} and charge conj.");
  else pvD->AddText("with D^{0} #rightarrow K^{-}#pi^{+} and charge conj.");

  shift+=dshift;
  pvEta = new TPaveText(0.12,0.65-shift,0.9,0.7-shift,"brNDC");
  pvEta->SetFillStyle(0);
  pvEta->SetBorderSize(0);
  pvEta->SetTextFont(42);
  pvEta->SetTextSize(0.03);
  pvEta->SetTextAlign(11);
  pvEta->AddText(Form("|#it{#eta}_{jet}| < 0.%d",9-Rpar));

  shift+=dshift;
  pv3 = new TPaveText(0.12,0.65-shift,0.9,0.7-shift,"brNDC");
  pv3->SetFillStyle(0);
  pv3->SetBorderSize(0);
  pv3->SetTextFont(42);
  pv3->SetTextSize(0.03);
  pv3->SetTextAlign(11);
  pv3->AddText(Form("%d < p_{T,%s} < %d GeV/#it{c}",(Int_t)fptbinsDA[0],fDmesonS.Data(),(Int_t)fptbinsDA[fptbinsDN]));

  TPaveText *pvEn= new TPaveText(0.15,0.80,0.8,0.85,"brNDC");
  pvEn->SetFillStyle(0);
  pvEn->SetBorderSize(0);
  pvEn->SetTextFont(42);
  pvEn->SetTextSize(0.045);
  pvEn->SetTextAlign(11);
  pvEn->AddText(Form("%s",fSystemS.Data()));
 // pvEn->AddText("PYTHIA6+HIJING, p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");

  shift = 0.3;
  pvD2 = new TPaveText(0.2,0.66-shift,0.9,0.7-shift,"brNDC");
  pvD2->SetFillStyle(0);
  pvD2->SetBorderSize(0);
  pvD2->SetTextFont(42);
  pvD2->SetTextSize(0.03);
  pvD2->SetTextAlign(11);
  if(fDmesonSpecie) pvD2->AddText("D^{*+} #rightarrow D^{0}#pi^{+} and charge conj.");
  else pvD2->AddText("D^{0} #rightarrow K^{-}#pi^{+} and charge conj.");

  pvJet2 = new TPaveText(0.2,0.61-shift,0.9,0.65-shift,"brNDC");
  pvJet2->SetFillStyle(0);
  pvJet2->SetBorderSize(0);
  pvJet2->SetTextFont(42);
  pvJet2->SetTextSize(0.03);
  pvJet2->SetTextAlign(11);
  pvJet2->AddText(Form("in Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.%d",Rpar));

  pvEta2 = new TPaveText(0.2,0.56-shift,0.8,0.6-shift,"brNDC");
  pvEta2->SetFillStyle(0);
  pvEta2->SetBorderSize(0);
  pvEta2->SetTextFont(42);
  pvEta2->SetTextSize(0.03);
  pvEta2->SetTextAlign(11);
  pvEta2->AddText(Form("|#it{#eta}_{jet}| < 0.%d",9-Rpar));

  pvJetPt2 = new TPaveText(0.2,0.48-shift,0.8,0.52-shift,"brNDC");
  pvJetPt2->SetFillStyle(0);
  pvJetPt2->SetBorderSize(0);
  pvJetPt2->SetTextFont(42);
  pvJetPt2->SetTextSize(0.03);
  pvJetPt2->SetTextAlign(11);
  pvJetPt2->AddText(Form("%.0f < p_{T.ch jet} < %.0f GeV/#it{c}",5.,50.));

  compareEfficiencies(1);
  compareEfficiencies(0);
  if(isRaw) compareRawSpectra();
  compareCorrSpectra(reg,unf);
  compareFD();

return;

}

void compareEfficiencies(int isprompt=1)
{
    gStyle->SetOptStat(000);

    TString fileName;
    if(isprompt) fileName = "DjetEff_prompt_jetpt5_50.root";
    else fileName = "DjetEff_nonPrompt_jetpt5_50.root";
    TH1D *hSpectrum[nFiles];

    TLegend  *leg;
    if(isprompt) leg = new TLegend(0.6,0.15,0.85,0.55,"Prompt efficiency");
    else leg = new TLegend(0.6,0.15,0.85,0.55,"Non-Prompt efficiency");
    TCanvas *cJetPt = new TCanvas("cJetPt","cJetPt",800,600);
    cJetPt->SetLogy();
    for(int i=0;i<nFiles;i++){
        TFile *fileIn = new TFile(Form("%s/efficiency/%s", inDir[i].Data(),fileName.Data()) );
        hSpectrum[i] = (TH1D*)fileIn->Get("hEff_reb");
        hSpectrum[i]->SetLineColor(colors[i]);
        hSpectrum[i]->SetMarkerColor(colors[i]);
        hSpectrum[i]->SetMarkerSize(0.9);
        hSpectrum[i]->SetMarkerStyle(markers[i]);

        hSpectrum[i]->GetXaxis()->SetRangeUser(3,36);
        //hSpectrum[i]->GetYaxis()->SetTitle("efficiency");
        if(!i) hSpectrum[i]->Draw();
        else hSpectrum[i]->Draw("same");
        leg->AddEntry(hSpectrum[i],Form("%s",desc[i].Data()),"p");
    }
     leg->Draw("same");
     pvEn->Draw("same");
     pvJet2->Draw("same");
     pvD2->Draw("same");
     pvEta2->Draw("same");
     pvJetPt2->Draw("same");

     if (isprompt){
          cJetPt->SaveAs(Form("%s/%s_PromptEfficiencies.pdf",outDir.Data(),outPlotName.Data()));
          cJetPt->SaveAs(Form("%s/%s_PromptEfficiencies.png",outDir.Data(),outPlotName.Data()));
    }
    else {
          cJetPt->SaveAs(Form("%s/%s_NonPromptEfficiencies.pdf",outDir.Data(),outPlotName.Data()));
          cJetPt->SaveAs(Form("%s/%s_NonPromptEfficiencies.png",outDir.Data(),outPlotName.Data()));
    }

     TCanvas *cRatio = new TCanvas("cRatio","cRatio",800,400);

     TH1D *hcentral = (TH1D*) hSpectrum[0]->Clone("hcentral");
     TH1D *hratios[nFiles-1];
     //TF1 *fline;
     TF1 *fline[nFiles-1];

     for(int i=0; i<nFiles-1; i++){
            hratios[i] = (TH1D*)hSpectrum[i+1]->Clone(Form("hratios_%d",i));
            hratios[i]->Divide(hcentral);
            hratios[i]->GetYaxis()->SetTitle("ratio");
            hratios[i]->GetYaxis()->SetRangeUser(0.5,1.8);
            if(!i)hratios[i]->Draw();
            else hratios[i]->Draw("same");

    }
    TLine *l=new TLine(3,1,36,1);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw();

    if (isprompt){
        cRatio->SaveAs(Form("%s/%s_PromptEfficiencies_ratio.pdf",outDir.Data(),outPlotName.Data()));
        cRatio->SaveAs(Form("%s/%s_PromptEfficiencies_ratio.png",outDir.Data(),outPlotName.Data()));
   }
   else {
        cRatio->SaveAs(Form("%s/%s_NonPromptEfficiencies_ratio.pdf",outDir.Data(),outPlotName.Data()));
        cRatio->SaveAs(Form("%s/%s_NonPromptEfficiencies_ratio.png",outDir.Data(),outPlotName.Data()));
   }

}


void compareRawSpectra()
{
    gStyle->SetOptStat(000);

    TString fileName;
    fileName = "JetPtSpectra_SB_noEff.root";
    TH1D *hSpectrum[nFiles];
    double events = 3.87988E8;

    TLegend  *leg = new TLegend(0.65,0.45,0.85,0.85,"raw yields");
    TCanvas *cJetPt = new TCanvas("cJetPt","cJetPt",800,600);
    cJetPt->SetLogy();
    for(int i=0;i<nFiles;i++){
        TFile *fileIn = new TFile(Form("%s/signalExtraction/%s", inDir[i].Data(),fileName.Data()) );
        hSpectrum[i] = (TH1D*)fileIn->Get("hjetptspectrumReb");
        hSpectrum[i]->Scale(1./events);
        hSpectrum[i]->Scale(1,"width");

        hSpectrum[i]->SetLineColor(colors[i]);
        hSpectrum[i]->SetMarkerColor(colors[i]);
        hSpectrum[i]->SetMarkerSize(0.9);
        hSpectrum[i]->SetMarkerStyle(markers[i]);
        //hSpectrum[i]->GetXaxis()->SetTitle("#it{p}_{T,D*} (GeV/#it{c})");
        hSpectrum[i]->GetYaxis()->SetTitle("dN/dp_{T} / events");
        if(!i) hSpectrum[i]->Draw();
        else hSpectrum[i]->Draw("same");
        leg->AddEntry(hSpectrum[i],Form("%s",desc[i].Data()),"p");
    }
     leg->Draw("same");
     pv3->Draw("same");
     pvEn->Draw("same");
     pvD->Draw("same");
     pvJet->Draw("same");
     pvEta->Draw("same");

          cJetPt->SaveAs(Form("%s/%s_RawSpectra.pdf",outDir.Data(),outPlotName.Data()));
          cJetPt->SaveAs(Form("%s/%s_RawSpectra.png",outDir.Data(),outPlotName.Data()));


     TCanvas *cRatio = new TCanvas("cRatio","cRatio",800,400);

     TH1D *hcentral = (TH1D*) hSpectrum[0]->Clone("hcentral");
     TH1D *hratios[nFiles-1];
     //TF1 *fline;
     TF1 *fline[nFiles-1];

     for(int i=0; i<nFiles-1; i++){
            hratios[i] = (TH1D*)hSpectrum[i+1]->Clone(Form("hratios_%d",i));
            hratios[i]->Divide(hcentral);
            hratios[i]->GetYaxis()->SetTitle("ratio");
            hratios[i]->GetYaxis()->SetRangeUser(0.5,1.5);
            if(!i)hratios[i]->Draw();
            else hratios[i]->Draw("same");

    }
    TLine *l=new TLine(3,1,36,1);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw();


        cRatio->SaveAs(Form("%s/%s_RawSpectra_ratio.pdf",outDir.Data(),outPlotName.Data()));
        cRatio->SaveAs(Form("%s/%s_RawSpectra_ratio.png",outDir.Data(),outPlotName.Data()));
}

void compareCorrSpectra(int reg, int unfold)
{
    gStyle->SetOptStat(000);

    TString fileName;
    if(unfold) fileName = "unfoldedSpectrum_unfoldedJetSpectrum.root";
    else fileName = "JetPtSpectrum_FDsub.root";


    TH1D *hSpectrum[nFiles];
    double events = 3.87988E8;

    TLegend  *leg = new TLegend(0.65,0.45,0.85,0.85, "corrected yields");
    TCanvas *cJetPt = new TCanvas("cJetPt","cJetPt",800,600);
    cJetPt->SetLogy();
    for(int i=0;i<nFiles;i++){
        TFile *fileIn;
        if(unfold) {
          fileIn = new TFile(Form("%s/unfolding_Bayes_%d/%s", inDir[i].Data(),reg,fileName.Data()) );
          hSpectrum[i] = (TH1D*)fileIn->Get("unfoldedSpectrum");
        }
        else {
          fileIn = new TFile(Form("%s/FDsubtraction/%s", inDir[i].Data(),fileName.Data()) );
          hSpectrum[i] = (TH1D*)fileIn->Get("hData_binned_sub");
        }
        hSpectrum[i]->Scale(1./events);
        hSpectrum[i]->Scale(1,"width");
        hSpectrum[i]->SetTitle();

        hSpectrum[i]->SetLineColor(colors[i]);
        hSpectrum[i]->SetMarkerColor(colors[i]);
        hSpectrum[i]->SetMarkerSize(0.9);
        hSpectrum[i]->SetMarkerStyle(markers[i]);
        //hSpectrum[i]->GetXaxis()->SetTitle("#it{p}_{T,D*} (GeV/#it{c})");
        hSpectrum[i]->GetYaxis()->SetTitle("dN/dp_{T} / events");
        if(!i) hSpectrum[i]->Draw();
        else hSpectrum[i]->Draw("same");
        leg->AddEntry(hSpectrum[i],Form("%s",desc[i].Data()),"p");
    }
     leg->Draw("same");
     pv3->Draw("same");
     pvEn->Draw("same");
     pvD->Draw("same");
     pvJet->Draw("same");
     pvEta->Draw("same");


          cJetPt->SaveAs(Form("%s/%s_CorrectedSpectra.pdf",outDir.Data(),outPlotName.Data()));
          cJetPt->SaveAs(Form("%s/%s_CorrectedSpectra.png",outDir.Data(),outPlotName.Data()));


     TCanvas *cRatio = new TCanvas("cRatio","cRatio",800,400);

     TH1F *hcentral = (TH1F*) hSpectrum[0]->Clone("hcentral");
     TH1F *hratios[nFiles-1];
     //TF1 *fline;
     TF1 *fline[nFiles-1];

     for(int i=0; i<nFiles-1; i++){
            hratios[i] = (TH1F*)hSpectrum[i+1]->Clone(Form("hratios_%d",i));
            hratios[i]->Divide(hcentral);
            hratios[i]->GetYaxis()->SetTitle("ratio");
            hratios[i]->GetYaxis()->SetRangeUser(0.7,1.5);
            if(!i)hratios[i]->Draw();
            else hratios[i]->Draw("same");

    }
    TLine *l=new TLine(3,1,50,1);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw();

        cRatio->SaveAs(Form("%s/%s_CorrectedSpectra_ratio.pdf",outDir.Data(),outPlotName.Data()));
        cRatio->SaveAs(Form("%s/%s_CorrectedSpectra_ratio.png",outDir.Data(),outPlotName.Data()));

        TH1F *hsys = new TH1F("hsys","syst. rms; p_{T,ch jet};  sys [%] (rms)",nJetBins,ptJetbins);
        TH1F *hmean = (TH1F*)hsys->Clone("hmean");
        getRMS(nFiles,hratios,hmean,hsys);

        hsys->GetYaxis()->SetRangeUser(0,20);
        hsys->SetLineColor(kViolet+2);
        TCanvas *cspecRMS = new TCanvas("cspecRMS","cspecRMS",800,400);
        hsys->Draw("hist");

        cspecRMS->SaveAs(Form("%s/%s_CorrectedSpectra_rms.pdf",outDir.Data(),outPlotName.Data()));
        cspecRMS->SaveAs(Form("%s/%s_CorrectedSpectra_rms.png",outDir.Data(),outPlotName.Data()));

        hmean->GetYaxis()->SetRangeUser(0.93,1.2);
        hmean->SetLineColor(kMagenta+1);
        TCanvas *cspecMean = new TCanvas("cspecMean","cspecMean",800,400);
        hmean->Draw("hist");
        l->Draw("same");

        cspecMean->SaveAs(Form("%s/%s_CorrectedSpectra_mean.pdf",outDir.Data(),outPlotName.Data()));
        cspecMean->SaveAs(Form("%s/%s_CorrectedSpectra_mean.png",outDir.Data(),outPlotName.Data()));

}

void compareFD()
{
    gStyle->SetOptStat(000);

    TString fileName;
    fileName = "JetPtSpectrum_FDsub.root";
    TH1D *hSpectrum[nFiles];
    double events = 3.87988E8;

    TLegend  *leg = new TLegend(0.6,0.12,0.85,0.45);
    TCanvas *cJetPt = new TCanvas("cJetPt","cJetPt",800,600);
    //cJetPt->SetLogy();
    for(int i=0;i<nFiles;i++){
        TFile *fileIn = new TFile(Form("%s/FDsubtraction/%s", inDir[i].Data(),fileName.Data()) );
        hSpectrum[i] = (TH1D*)fileIn->Get("hFD_ratio");
        //hSpectrum[i]->Scale(1./events);
        //hSpectrum[i]->Scale(1,"width");

        hSpectrum[i]->SetLineColor(colors[i]);
        hSpectrum[i]->SetMarkerColor(colors[i]);
        hSpectrum[i]->SetMarkerSize(0.9);
        hSpectrum[i]->SetMarkerStyle(markers[i]);
        //hSpectrum[i]->GetXaxis()->SetTitle("#it{p}_{T,D*} (GeV/#it{c})");
        hSpectrum[i]->GetYaxis()->SetTitle("FD ratio");
        hSpectrum[i]->GetYaxis()->SetRangeUser(0,0.4);
        if(!i) hSpectrum[i]->Draw();
        else hSpectrum[i]->Draw("same");
        leg->AddEntry(hSpectrum[i],Form("%s",desc[i].Data()),"p");
    }
     leg->Draw("same");


          cJetPt->SaveAs(Form("%s/%s_FDFraction.pdf",outDir.Data(),outPlotName.Data()));
          cJetPt->SaveAs(Form("%s/%s_FDFraction.png",outDir.Data(),outPlotName.Data()));


     TCanvas *cRatio = new TCanvas("cRatio","cRatio",800,400);

     TH1D *hcentral = (TH1D*) hSpectrum[0]->Clone("hcentral");
     TH1D *hratios[nFiles-1];
     //TF1 *fline;
     TF1 *fline[nFiles-1];

     for(int i=0; i<nFiles-1; i++){
            hratios[i] = (TH1D*)hSpectrum[i+1]->Clone(Form("hratios_%d",i));
            hratios[i]->Divide(hcentral);
            hratios[i]->GetYaxis()->SetTitle("ratio");
            hratios[i]->GetYaxis()->SetRangeUser(0.7,1.5);
            if(!i)hratios[i]->Draw();
            else hratios[i]->Draw("same");

    }
    TLine *l=new TLine(3,1,36,1);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw();

        cRatio->SaveAs(Form("%s/%s_FDFraction_ratio.pdf",outDir.Data(),outPlotName.Data()));
        cRatio->SaveAs(Form("%s/%s_FDFraction_ratio.png",outDir.Data(),outPlotName.Data()));

}


void getRMS(const int nFiles, TH1F **hratio, TH1F *hmean, TH1F *hsys)
{

  //TH1D *hsys = new TH1D("hsys","syst. rms; p_{T,ch jet};  sys [%] (rms)",nJetBins,ptJetbins);
  hsys->SetTitle();
  hsys->SetLineColor(1);
  hsys->SetLineWidth(2);
  hsys->SetLineStyle(2);
  hsys->SetMaximum(0.12);
  hsys->GetXaxis()->SetLabelSize(0.05);
  hsys->GetXaxis()->SetTitleSize(0.05);
  hsys->GetYaxis()->SetTitleSize(0.05);
  hsys->GetYaxis()->SetLabelSize(0.05);
  hsys->GetYaxis()->SetTitleOffset(0.8);

  //hmean = (TH1F*)hsys->Clone("hmean");
  hmean->GetYaxis()->SetTitle("mean");
//  hmean->GetYaxis()->SetRangeUser(0.95,1.1);
  hmean->SetMarkerStyle(20);
  hmean->SetLineStyle(1);
  hmean->SetTitle();

  double *rms = new double[nJetBins];
  double *mean = new double[nJetBins];
  for(int i=0; i<nJetBins; i++){
      rms[i] = 0;
       mean[i] = 0;
       for (int j=0; j<nFiles-1; j++){
         mean[i] = mean[i]+ ( hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) );
      //double m = ( hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) );
      //rms[i] = rms[i]+ ( 1-hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) ) * ( 1-hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) ) ;
      }
   mean[i] = mean[i]/(double)(nFiles-1);

   for (int j=0; j<nFiles-1; j++){
      //mean[i] = mean[i]+ ( hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) );
      //double m = ( hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) );
    //  rms[i] = rms[i]+ ( mean[i]-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) ) * ( mean[i]-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) ) ;
      rms[i] = rms[i]+ ( 1-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) ) * ( 1-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins[i]+ptJetbins[i+1])/2. )) ) ;


  }
      rms[i] = sqrt(rms[i]/(double)(nFiles-1));
      hsys->SetBinContent(i+1,rms[i]*100);

      hmean->SetBinContent(i+1,mean[i]);
      cout << "RMS pT " << (ptJetbins[i]+ptJetbins[i+1])/2. << " GeV/c:\t" << rms[i]*100 << endl;
      cout << "Mean pT " << (ptJetbins[i]+ptJetbins[i+1])/2. << " GeV/c:\t" << mean[i] << endl;
  }

}
