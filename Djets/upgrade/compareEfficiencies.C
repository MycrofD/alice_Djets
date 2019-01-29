
#include <string>
#include <sstream>
#include <iostream>

#include "config.h"

Int_t colors2[] = {kBlue+1,kGreen+2,kOrange+1,kViolet+1,kGray+2,kMagenta+2,kYellow+2,kCyan+1,4,6,8};
Int_t markers2[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};
Int_t linestyle2[] = {1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

const int nHist = 10;
/*
TString histName[] = {
  "mcPt_prompt_int",
  "recoPt_prompt_int",
  "mcPt_reb_prompt_int",
  "recoPt_reb_prompt_int",
  "mcPt_nonprompt_int",
  "recoPt_nonprompt_int",
  "mcPt_reb_nonprompt_int",
  "recoPt_reb_nonprompt_int",
  "eff_reb_prompt",
  "eff_reb_nonprompt"
}
*/
TString histName[] = {
  "mcPt_prompt",
  "recoPt_prompt",
  "mcPt_reb_prompt",
  "recoPt_reb_prompt",
  "mcPt_nonprompt",
  "recoPt_nonprompt",
  "mcPt_reb_nonprompt",
  "recoPt_reb_nonprompt",
  "eff_reb_prompt",
  "eff_reb_nonprompt"
}

TString titles[] = {
  "MC p_{T}, prompt",
  "MC p_{T}, reco, prompt",
  "MC p_{T}, prompt",
  "MC p_{T}, reco, prompt",
  "MC p_{T}, non-prompt",
  "MC p_{T}, reco, non-prompt",
  "MC p_{T}, non-prompt",
  "MC p_{T}, reco, non-prompt",
  "efficiency, prompt",
  "efficiency, non-prompt"

}

TString inDir = "/home/basia/Work/alice/analysis/upgradeProjections/results_HB4/";

void compareEfficiencies(TString outName = "/home/basia/Work/alice/analysis/upgradeProjections/results_HB4/plots/efficiencies/")
{

  gStyle->SetOptStat(0000); //Mean and RMS shown
  gSystem->Exec(Form("mkdir -p %s",outName.Data()));

  //compare MB and HB

  compareBins(outName,inDir,"Bis",0,0,"Prodcomparison_moreBins");
  compareBins(outName,inDir,"Bis",0,1,"Prodcomparison_stdBins");
  compareBins(outName,inDir,"Tres",0,0,"Prodcomparison_moreBins");
  compareBins(outName,inDir,"Tres",0,1,"Prodcomparison_stdBins");
  compareBins(outName,inDir,"",0,0,"Prodcomparison_moreBins");
  compareBins(outName,inDir,"",0,1,"Prodcomparison_stdBins");

  // compare bis, tres, all
  // MB
  compareProd(outName,inDir,"efficiencies_MB","_stdBins","Prodcomparison_MB_stdBins");
  compareProd(outName,inDir,"efficiencies_MB","_moreBins","Prodcomparison_MB_moreBins");

  // HB
  compareProd(outName,inDir,"efficiencies_HBinsHB","_weighted_stdBins","Prodcomparison_HB_stdBins");
  compareProd(outName,inDir,"efficiencies_HBinsHB","_weighted_moreBins","Prodcomparison_HB_moreBins");

  compareProd(outName,inDir,"efficiencies_HBinsHB","_notweighted_stdBins","Prodcomparison_HB_stdBins");
  compareProd(outName,inDir,"efficiencies_HBinsHB","_notweighted_moreBins","Prodcomparison_HB_moreBins");
return;

  // compare hard bins between each other
   compareHB(outName,inDir,"efficiencies_HBinsHBBis","_weighted_stdBins","HBcomparison_Bis");
   compareHB(outName,inDir,"efficiencies_HBinsHBTres","_weighted_stdBins","HBcomparison_Tres");
   compareHB(outName,inDir,"efficiencies_HBinsHB","_weighted_stdBins","HBcomparison");

   compareHB(outName,inDir,"efficiencies_HBinsHBBis","_notweighted_stdBins","HBcomparison_Bis");
   compareHB(outName,inDir,"efficiencies_HBinsHBTres","_notweighted_stdBins","HBcomparison_Tres");
   compareHB(outName,inDir,"efficiencies_HBinsHB","_notweighted_stdBins","HBcomparison");

}

void compareBins(TString outName, TString inDir, TString postfix, Bool_t isWght = 0, Bool_t isStdBins = 1, TString outHistName)
{

          /*  efficiencies_MB.root
            efficiencies_MBBis.root
            efficiencies_HBinsHBBis_weighted.root
            efficiencies_HBinsHBBis_weighted_stdBins.root */

            const int nHBFiles = 2;
            TString name[] = {"HB w/o weight","HB w/ weight"};

            TString out = outName;

            TCanvas *cEff[nHist];
            TCanvas *cEffRatio[nHist];
            TLegend *leg = new TLegend(0.7,0.7,0.85,0.85, Form("%s",postfix.Data()));
            leg->SetBorderSize(0);
            for (int j=0; j<nHist; j++) {
              cEff[j] = new TCanvas(Form("cEff_%d",j),Form("cEff_%d",j),1200,800);
              cEff[j]->SetLogy();
              cEffRatio[j] = new TCanvas(Form("cEffRatio_%d",j),Form("cEffRatio_%d",j),1200,800);
              //leg[j] = new TLegend(0.65,0.21,0.89,0.65);
            //  leg[j]->SetBorderSize(0);
            }

            TString fileName  = inDir;
            fileName += "efficiencies_MB";
            fileName += postfix;
            if(isStdBins) fileName += "_stdBins";
            else fileName += "_moreBins";
            fileName += ".root";
            TH1D *hist0[nHist];
            TFile *file = new TFile(Form("%s",fileName.Data()),"READ");

            for (int j=0; j<nHist-2; j++) {
              hist0[j] = (TH1D*)file->Get(Form("%s",histName[j].Data()));
              if(!hist0[j]) { cout << "NO HISTOGRAM !!!!!" << histName[j] << endl; return; }
              hist0[j]->Sumw2();
              hist0[j]->SetTitle(Form("%s",titles[j].Data()));
              hist0[j]->SetLineColor(2);
              hist0[j]->SetMarkerColor(2);
              hist0[j]->SetMarkerStyle(29);
              hist0[j]->SetMarkerSize(1.8);
              hist0[j]->Scale(1./hist0[j]->Integral());
              hist0[j]->Scale(1,"width");
              hist0[j]->SetMaximum(hist0[j]->GetMaximum()*1.5);
              hist0[j]->SetMinimum(1);
              hist0[j]->GetXaxis()->SetTitle("p_{T,D}");
              hist0[j]->GetYaxis()->SetTitle("dN/dp_{T}");

              cEff[j]->cd();
              hist0[j]->DrawCopy();

              TH1D *histRatio = (TH1D*)hist0[j]->Clone("histRatio");
              histRatio->Divide(hist0[j],hist0[j],1,1,"b");
              histRatio->GetYaxis()->SetRangeUser(0.4,2);
              cEffRatio[j]->cd();
              histRatio->DrawCopy();

            }

            leg->AddEntry(hist0[0],"MB","p");

            for (int i=0; i<nHBFiles; i++){
                TString fileName  = inDir;
                fileName += "efficiencies_HBinsHB";
                fileName += postfix;
                if(!i) fileName += "_notweighted";
                else fileName += "_weighted";
                if(isStdBins) fileName += "_stdBins";
                else fileName += "_moreBins";
                fileName += ".root";
                TFile *file = new TFile(Form("%s",fileName.Data()),"READ");
                if(file->IsZombie()) { cout << "NO FILE !!!!!" << fileName << endl; return; }
                for (int j=0; j<nHist-2; j++) {
                  TH1D *hist = (TH1D*)file->Get(Form("%s",histName[j].Data()));
                  if(!hist) { cout << "NO HISTOGRAM !!!!!" << histName[j] << endl; return; }
                  hist->Sumw2();
                  hist->SetTitle();
                  hist->SetLineColor(colors2[i]);
                  hist->SetMarkerColor(colors2[i]);
                  hist->SetMarkerStyle(markers2[i]);
                  hist->SetMarkerSize(1.4);
                  hist->Scale(1./hist->Integral());
                  hist->Scale(1,"width");
                  hist->SetMaximum(hist->GetMaximum()*2);
                  hist->GetXaxis()->SetTitle("p_{T,D}");
                  hist->GetYaxis()->SetTitle("dN/dp_{T}");
                  if(!j) leg->AddEntry(hist,Form("%s",name[i].Data()),"p");
                  cEff[j]->cd();
                  hist->DrawCopy("same");

                  TH1D *histRatio = (TH1D*)hist->Clone("histRatio");
                  histRatio->Divide(hist,hist0[j],1,1,"b");
                  histRatio->GetYaxis()->SetRangeUser(0.4,2);
                  cEffRatio[j]->cd();
                  histRatio->DrawCopy("same");
                }

            }

            for (int j=0; j<nHist-2; j++) {
              cEff[j]->cd();
              hist0[j]->DrawCopy("same");
              leg->Draw("same");
              SaveCanvas(cEff[j],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),histName[j].Data()));

              cEffRatio[j]->cd();
              leg->Draw("same");
              SaveCanvas(cEffRatio[j],Form("%s/%s%s_%s_ratio",outName.Data(),outHistName.Data(),postfix.Data(),histName[j].Data()));
            }


            // efficiencies
            TLegend *leg2 = new TLegend(0.55,0.15,0.7,0.3,Form("prompt eff, %s",postfix.Data()));
            leg2->SetBorderSize(0);
            TLegend *leg3 = new TLegend(0.55,0.15,0.7,0.3,Form("non-prompt eff, %s",postfix.Data()));
            leg3->SetBorderSize(0);

            fileName  = inDir;
            fileName += "efficiencies_MB";
            fileName += postfix;
            if(isStdBins) fileName += "_stdBins";
            else fileName += "_moreBins";
            fileName += ".root";
            file = new TFile(Form("%s",fileName.Data()),"READ");
            TH1D *heffPrompt0 = (TH1D*)file->Get(Form("%s",histName[nHist-2].Data()));
            if(!heffPrompt0) { cout << "NO HISTOGRAM !!!!!" << histName[nHist-2] << endl; return; }
            heffPrompt0->Sumw2();
            heffPrompt0->SetTitle();
            heffPrompt0->SetLineColor(2);
            heffPrompt0->SetMarkerColor(2);
            heffPrompt0->SetMarkerStyle(29);
            heffPrompt0->SetMarkerSize(1.8);
            heffPrompt0->SetMaximum(0.8);
            heffPrompt0->SetMinimum(0.01);
            cEff[nHist-2]->cd();
            heffPrompt0->DrawCopy();
            leg2->AddEntry(heffPrompt0,"MB","p");

            TH1D *heffPromptRatio = (TH1D*)heffPrompt0->Clone("heffPromptRatio");
            heffPromptRatio->Divide(heffPrompt0,heffPrompt0,1,1,"b");
            heffPromptRatio->GetYaxis()->SetRangeUser(0.4,1.6);
            cEffRatio[nHist-2]->cd();
            heffPromptRatio->DrawCopy();

            TH1D *heffNonPrompt0 = (TH1D*)file->Get(Form("%s",histName[nHist-1].Data()));
            if(!heffNonPrompt0) { cout << "NO HISTOGRAM !!!!!" << histName[nHist-1] << endl; return; }
            heffNonPrompt0->Sumw2();
            heffNonPrompt0->SetTitle();
            heffNonPrompt0->SetLineColor(2);
            heffNonPrompt0->SetMarkerColor(2);
            heffNonPrompt0->SetMarkerStyle(29);
            heffNonPrompt0->SetMarkerSize(1.8);
            heffNonPrompt0->SetMaximum(0.8);
            heffNonPrompt0->SetMinimum(0.01);
            cEff[nHist-1]->cd();
            heffNonPrompt0->DrawCopy();
            leg3->AddEntry(heffNonPrompt0,"MB","p");

            TH1D *heffNonPromptRatio = (TH1D*)heffNonPrompt0->Clone("heffNonPromptRatio");
            heffNonPromptRatio->Divide(heffNonPrompt0,heffNonPrompt0,1,1,"b");
            heffNonPromptRatio->GetYaxis()->SetRangeUser(0.4,1.6);
            cEffRatio[nHist-1]->cd();
            heffNonPromptRatio->DrawCopy();

            for (int i=0; i<nHBFiles; i++){
                TString fileName  = inDir;
                fileName += "efficiencies_HBinsHB";
                fileName += postfix;
                if(!i) fileName += "_notweighted";
                else fileName += "_weighted";
                if(isStdBins) fileName += "_stdBins";
                else fileName += "_moreBins";
                fileName += ".root";
                TFile *file = new TFile(Form("%s",fileName.Data()),"READ");
                if(file->IsZombie()) { cout << "NO FILE !!!!!" << fileName << endl; return; }

                TH1D *heffPrompt = (TH1D*)file->Get(Form("%s",histName[nHist-2].Data()));
                if(!heffPrompt) { cout << "NO HISTOGRAM !!!!!" << histName[nHist-2] << endl; return; }
                heffPrompt->Sumw2();
                heffPrompt->SetTitle();
                heffPrompt->SetLineColor(colors2[i]);
                heffPrompt->SetMarkerColor(colors2[i]);
                heffPrompt->SetMarkerStyle(markers2[i]);
                heffPrompt->SetMarkerSize(1.8);
                heffPrompt->SetMaximum(0.8);
                cEff[nHist-2]->cd();
                heffPrompt->DrawCopy("same");

                leg2->AddEntry(heffPrompt,Form("%s",name[i].Data()),"p");

                TH1D *heffPromptRatio = (TH1D*)heffPrompt->Clone("heffPromptRatio");
                heffPromptRatio->Divide(heffPrompt,heffPrompt0,1,1,"b");
                heffPromptRatio->GetYaxis()->SetRangeUser(0.4,2);
                cEffRatio[nHist-2]->cd();
                heffPromptRatio->DrawCopy("same");

                TH1D *heffNonPrompt = (TH1D*)file->Get(Form("%s",histName[nHist-1].Data()));
                if(!heffNonPrompt) { cout << "NO HISTOGRAM !!!!!" << histName[nHist-1] << endl; return; }
                heffNonPrompt->Sumw2();
                heffNonPrompt->SetTitle();
                heffNonPrompt->SetLineColor(colors2[i]);
                heffNonPrompt->SetMarkerColor(colors2[i]);
                heffNonPrompt->SetMarkerStyle(markers2[i]);
                heffNonPrompt->SetMarkerSize(1.8);
                heffNonPrompt->SetMaximum(0.5);
                cEff[nHist-1]->cd();
                heffNonPrompt->DrawCopy("same");

                leg3->AddEntry(heffNonPrompt,Form("%s",name[i].Data()),"p");

                TH1D *heffNonPromptRatio = (TH1D*)heffNonPrompt->Clone("heffNonPromptRatio");
                heffNonPromptRatio->Divide(heffNonPrompt,heffNonPrompt0,1,1,"b");
                heffNonPromptRatio->GetYaxis()->SetRangeUser(0.4,2);
                cEffRatio[nHist-1]->cd();
                heffNonPromptRatio->DrawCopy("same");

            }

            cEff[nHist-2]->cd();
            heffPrompt0->DrawCopy("same");
            leg2->Draw("same");
            SaveCanvas(cEff[nHist-2],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),"effPrompt"));
            cEff[nHist-1]->cd();
            heffNonPrompt0->DrawCopy("same");
            leg3->Draw("same");
            SaveCanvas(cEff[nHist-1],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),"effNonPrompt"));

            cEffRatio[nHist-2]->cd();
            leg2->Draw("same");
            SaveCanvas(cEffRatio[nHist-2],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),"effPromptRatio"));
            cEffRatio[nHist-1]->cd();
            leg3->Draw("same");
            SaveCanvas(cEffRatio[nHist-1],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),"effNonPromptRatio"));

            return;

}

void compareProd(TString outName, TString inDir, TString fileN, TString postfix, TString outHistName)
{

            TString name[] = {"Bis","Tres"};
            TString out = outName;

            TCanvas *cEff[nHist];
            TCanvas *cEffRatio[nHist];
            TLegend *leg = new TLegend(0.7,0.7,0.85,0.85);
            leg->SetBorderSize(0);
            for (int j=0; j<nHist; j++) {
              cEff[j] = new TCanvas(Form("cEff_%d",j),Form("cEff_%d",j),1200,800);
              cEff[j]->SetLogy();
              cEffRatio[j] = new TCanvas(Form("cEffRatio_%d",j),Form("cEffRatio_%d",j),1200,800);
              //leg[j] = new TLegend(0.65,0.21,0.89,0.65);
            //  leg[j]->SetBorderSize(0);
            }

            TString fileName  = inDir;
            fileName += fileN;
            fileName += postfix;
            fileName += ".root";
            TH1D *hist0[nHist];
            TFile *file = new TFile(Form("%s",fileName.Data()),"READ");

            for (int j=0; j<nHist-2; j++) {
              hist0[j] = (TH1D*)file->Get(Form("%s",histName[j].Data()));
              if(!hist0[j]) { cout << "NO HISTOGRAM !!!!!" << histName[j] << endl; return; }
              hist0[j]->Sumw2();
              hist0[j]->SetTitle(Form("%s",titles[j].Data()));
              hist0[j]->SetLineColor(2);
              hist0[j]->SetMarkerColor(2);
              hist0[j]->SetMarkerStyle(29);
              hist0[j]->SetMarkerSize(1.8);
              hist0[j]->Scale(1./hist0[j]->Integral());
              hist0[j]->Scale(1,"width");
              hist0[j]->SetMaximum(hist0[j]->GetMaximum()*1.5);
              hist0[j]->SetMinimum(1);
              hist0[j]->GetXaxis()->SetTitle("p_{T,D}");
              hist0[j]->GetYaxis()->SetTitle("dN/dp_{T}");

              cEff[j]->cd();
              hist0[j]->DrawCopy();

              TH1D *histRatio = (TH1D*)hist0[j]->Clone("histRatio");
              histRatio->Divide(hist0[j],hist0[j],1,1,"b");
              histRatio->GetYaxis()->SetRangeUser(0.4,2);
              cEffRatio[j]->cd();
              histRatio->DrawCopy();

            }
            leg->AddEntry(hist0[0],"All","p");

            for (int i=0; i<2; i++){
                TString fileName  = inDir;
                fileName += fileN;
                fileName += name[i];
                fileName += postfix;
                fileName += ".root";
                TFile *file = new TFile(Form("%s",fileName.Data()),"READ");
                if(file->IsZombie()) { cout << "NO FILE !!!!!" << fileName << endl; return; }
                for (int j=0; j<nHist-2; j++) {
                  TH1D *hist = (TH1D*)file->Get(Form("%s",histName[j].Data()));
                  if(!hist) { cout << "NO HISTOGRAM !!!!!" << histName[j] << endl; return; }
                  hist->Sumw2();
                  hist->SetTitle();
                  hist->SetLineColor(colors2[i]);
                  hist->SetMarkerColor(colors2[i]);
                  hist->SetMarkerStyle(markers2[i]);
                  hist->SetMarkerSize(1.4);
                  hist->Scale(1./hist->Integral());
                  hist->Scale(1,"width");
                  hist->SetMaximum(hist->GetMaximum()*2);
                  hist->GetXaxis()->SetTitle("p_{T,D}");
                  hist->GetYaxis()->SetTitle("dN/dp_{T}");
                  if(!j) leg->AddEntry(hist,Form("%s",name[i].Data()),"p");
                  cEff[j]->cd();
                  hist->DrawCopy("same");

                  TH1D *histRatio = (TH1D*)hist->Clone("histRatio");
                  histRatio->Divide(hist,hist0[j],1,1,"b");
                  histRatio->GetYaxis()->SetRangeUser(0.4,2);
                  cEffRatio[j]->cd();
                  histRatio->DrawCopy("same");
                }

            }

            for (int j=0; j<nHist-2; j++) {
              cEff[j]->cd();
              hist0[j]->DrawCopy("same");
              leg->Draw("same");
              SaveCanvas(cEff[j],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),histName[j].Data()));

              cEffRatio[j]->cd();
              leg->Draw("same");
              SaveCanvas(cEffRatio[j],Form("%s/%s%s_%s_ratio",outName.Data(),outHistName.Data(),postfix.Data(),histName[j].Data()));
            }

            // efficiencies
            TLegend *leg2 = new TLegend(0.35,0.15,0.55,0.3,"prompt eff");
            leg2->SetBorderSize(0);
            TLegend *leg3 = new TLegend(0.35,0.15,0.55,0.3,"non-prompt eff");
            leg3->SetBorderSize(0);

            fileName  = inDir;
            fileName += fileN;
            fileName += postfix;
            fileName += ".root";
            file = new TFile(Form("%s",fileName.Data()),"READ");
            TH1D *heffPrompt0 = (TH1D*)file->Get(Form("%s",histName[nHist-2].Data()));
            if(!heffPrompt0) { cout << "NO HISTOGRAM !!!!!" << histName[nHist-2] << endl; return; }
            heffPrompt0->Sumw2();
            heffPrompt0->SetTitle();
            heffPrompt0->SetLineColor(2);
            heffPrompt0->SetMarkerColor(2);
            heffPrompt0->SetMarkerStyle(29);
            heffPrompt0->SetMarkerSize(1.8);
            heffPrompt0->SetMaximum(0.8);
            heffPrompt0->SetMinimum(0.01);
            cEff[nHist-2]->cd();
            heffPrompt0->DrawCopy();
            leg2->AddEntry(heffPrompt0,"All","p");

            TH1D *heffPromptRatio = (TH1D*)heffPrompt0->Clone("heffPromptRatio");
            heffPromptRatio->Divide(heffPrompt0,heffPrompt0,1,1,"b");
            heffPromptRatio->GetYaxis()->SetRangeUser(0.4,1.6);
            cEffRatio[nHist-2]->cd();
            heffPromptRatio->DrawCopy();

            TH1D *heffNonPrompt0 = (TH1D*)file->Get(Form("%s",histName[nHist-1].Data()));
            if(!heffNonPrompt0) { cout << "NO HISTOGRAM !!!!!" << histName[nHist-1] << endl; return; }
            heffNonPrompt0->Sumw2();
            heffNonPrompt0->SetTitle();
            heffNonPrompt0->SetLineColor(2);
            heffNonPrompt0->SetMarkerColor(2);
            heffNonPrompt0->SetMarkerStyle(29);
            heffNonPrompt0->SetMarkerSize(1.8);
            heffNonPrompt0->SetMaximum(0.8);
            heffNonPrompt0->SetMinimum(0.01);
            cEff[nHist-1]->cd();
            heffNonPrompt0->DrawCopy();
            leg3->AddEntry(heffNonPrompt0,"All","p");

            TH1D *heffNonPromptRatio = (TH1D*)heffNonPrompt0->Clone("heffNonPromptRatio");
            heffNonPromptRatio->Divide(heffNonPrompt0,heffNonPrompt0,1,1,"b");
            heffNonPromptRatio->GetYaxis()->SetRangeUser(0.4,1.6);
            cEffRatio[nHist-1]->cd();
            heffNonPromptRatio->DrawCopy();

            for (int i=0; i<2; i++){
                TString fileName  = inDir;
                fileName += fileN;
                fileName += name[i];
                fileName += postfix;
                fileName += ".root";
                TFile *file = new TFile(Form("%s",fileName.Data()),"READ");
                if(file->IsZombie()) { cout << "NO FILE !!!!!" << fileName << endl; return; }

                TH1D *heffPrompt = (TH1D*)file->Get(Form("%s",histName[nHist-2].Data()));
                if(!heffPrompt) { cout << "NO HISTOGRAM !!!!!" << histName[nHist-2] << endl; return; }
                heffPrompt->Sumw2();
                heffPrompt->SetTitle();
                heffPrompt->SetLineColor(colors2[i]);
                heffPrompt->SetMarkerColor(colors2[i]);
                heffPrompt->SetMarkerStyle(markers2[i]);
                heffPrompt->SetMarkerSize(1.8);
                heffPrompt->SetMaximum(0.8);
                cEff[nHist-2]->cd();
                heffPrompt->DrawCopy("same");

                leg2->AddEntry(heffPrompt,Form("%s",name[i].Data()),"p");

                TH1D *heffPromptRatio = (TH1D*)heffPrompt->Clone("heffPromptRatio");
                heffPromptRatio->Divide(heffPrompt,heffPrompt0,1,1,"b");
                heffPromptRatio->GetYaxis()->SetRangeUser(0.4,2);
                cEffRatio[nHist-2]->cd();
                heffPromptRatio->DrawCopy("same");

                TH1D *heffNonPrompt = (TH1D*)file->Get(Form("%s",histName[nHist-1].Data()));
                if(!heffNonPrompt) { cout << "NO HISTOGRAM !!!!!" << histName[nHist-1] << endl; return; }
                heffNonPrompt->Sumw2();
                heffNonPrompt->SetTitle();
                heffNonPrompt->SetLineColor(colors2[i]);
                heffNonPrompt->SetMarkerColor(colors2[i]);
                heffNonPrompt->SetMarkerStyle(markers2[i]);
                heffNonPrompt->SetMarkerSize(1.8);
                heffNonPrompt->SetMaximum(0.5);
                cEff[nHist-1]->cd();
                heffNonPrompt->DrawCopy("same");

                leg3->AddEntry(heffNonPrompt,Form("%s",name[i].Data()),"p");

                TH1D *heffNonPromptRatio = (TH1D*)heffNonPrompt->Clone("heffNonPromptRatio");
                heffNonPromptRatio->Divide(heffNonPrompt,heffNonPrompt0,1,1,"b");
                heffNonPromptRatio->GetYaxis()->SetRangeUser(0.4,2);
                cEffRatio[nHist-1]->cd();
                heffNonPromptRatio->DrawCopy("same");

            }

            cEff[nHist-2]->cd();
            heffPrompt0->DrawCopy("same");
            leg2->Draw("same");
            SaveCanvas(cEff[nHist-2],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),"effPrompt"));
            cEff[nHist-1]->cd();
            heffNonPrompt0->DrawCopy("same");
            leg3->Draw("same");
            SaveCanvas(cEff[nHist-1],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),"effNonPrompt"));

            cEffRatio[nHist-2]->cd();
            leg2->Draw("same");
            SaveCanvas(cEffRatio[nHist-2],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),"effPromptRatio"));
            cEffRatio[nHist-1]->cd();
            leg3->Draw("same");
            SaveCanvas(cEffRatio[nHist-1],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),"effNonPromptRatio"));

            return;
}

void compareHB(TString outName, TString inDir, TString fileN, TString postfix, TString outHistName)
{
            TString out = outName;
            TCanvas *cEff[nHist+2];
            TLegend *leg = new TLegend(0.75,0.6,0.89,0.85);
            leg->SetBorderSize(0);
            for (int j=0; j<nHist+2; j++) {
              cEff[j] = new TCanvas(Form("cEff_%d",j),Form("cEff_%d",j),1200,800);
              if(j<nHist)cEff[j]->SetLogy();
              //leg[j] = new TLegend(0.65,0.21,0.89,0.65);
            //  leg[j]->SetBorderSize(0);
            }

            TString fileName  = inDir;
            fileName += fileN;
            fileName += postfix;
            fileName += ".root";
            TFile *file = new TFile(Form("%s",fileName.Data()),"READ");
            TH1D *hist0[nHist];
            for (int j=0; j<nHist-2; j++) {
              hist0[j] = (TH1D*)file->Get(Form("%s",histName[j].Data()));
              if(!hist0[j]) { cout << "NO HISTOGRAM !!!!!" << histName[j] << endl; return; }
              hist0[j]->Sumw2();
              hist0[j]->SetTitle(Form("%s",titles[j].Data()));
              hist0[j]->SetLineColor(2);
              hist0[j]->SetMarkerColor(2);
              hist0[j]->SetMarkerStyle(29);
              hist0[j]->SetMarkerSize(1.8);
              hist0[j]->Scale(1,"width");
              hist0[j]->SetMaximum(hist0[j]->GetMaximum()*1.5);
              //hist0[j]->SetMinimum(hist0[j]->GetMinimum()*0.8);
              hist0[j]->GetXaxis()->SetTitle("p_{T,D}");
              hist0[j]->GetYaxis()->SetTitle("dN/dp_{T}");

              cEff[j]->cd();
              hist0[j]->DrawCopy();
            }
            leg->AddEntry(hist0[0],"All HB","p");

            for (int i=2; i<7; i++){
                TString fileName  = inDir;
                fileName += fileN;
                fileName += i;
                fileName += postfix;
                fileName += ".root";
                TFile *file = new TFile(Form("%s",fileName.Data()),"READ");
                if(file->IsZombie()) { cout << "NO FILE !!!!!" << fileName << endl; return; }
                for (int j=0; j<nHist-2; j++) {
                  TH1D *hist = (TH1D*)file->Get(Form("%s",histName[j].Data()));
                  if(!hist) { cout << "NO HISTOGRAM !!!!!" << histName[j] << endl; return; }
                  hist->Sumw2();
                  hist->SetTitle();
                  hist->SetLineColor(colors2[i-2]);
                  hist->SetMarkerColor(colors2[i-2]);
                  hist->SetMarkerStyle(markers2[i-2]);
                  hist->SetMarkerSize(1.4);
                  hist->Scale(1,"width");
                  hist->SetMaximum(hist->GetMaximum()*2);
                  hist->GetXaxis()->SetTitle("p_{T,D}");
                  hist->GetYaxis()->SetTitle("dN/dp_{T}");
                  if(!j) {
                    leg->AddEntry(hist,Form("HB_%d",i),"p");
                  }
                  cEff[j]->cd();
                  hist->DrawCopy("same");
                }

            }

            for (int j=0; j<nHist-2; j++) {
              cEff[j]->cd();
              hist0[j]->DrawCopy("same");
              leg->Draw("same");
              SaveCanvas(cEff[j],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),histName[j].Data()));
            }

            // efficiencies
            TLegend *leg2 = new TLegend(0.6,0.15,0.75,0.4,"prompt eff");
            leg2->SetBorderSize(0);
            TLegend *leg3 = new TLegend(0.6,0.15,0.75,0.4,"non-prompt eff");
            leg3->SetBorderSize(0);

            fileName  = inDir;
            fileName += fileN;
            fileName += postfix;
            fileName += ".root";
            file = new TFile(Form("%s",fileName.Data()),"READ");
            TH1D *heffPrompt0 = (TH1D*)file->Get(Form("%s",histName[nHist-2].Data()));
            if(!heffPrompt0) { cout << "NO HISTOGRAM !!!!!" << histName[nHist-2] << endl; return; }
            heffPrompt0->Sumw2();
            heffPrompt0->SetTitle(Form("%s",titles[nHist-2].Data()));
            heffPrompt0->SetLineColor(2);
            heffPrompt0->SetMarkerColor(2);
            heffPrompt0->SetMarkerStyle(29);
            heffPrompt0->SetMarkerSize(1.8);
            heffPrompt0->SetMaximum(0.8);
            heffPrompt0->SetMinimum(0.01);
            cEff[nHist-2]->cd();
            heffPrompt0->DrawCopy();
            leg2->AddEntry(heffPrompt0,"All HB","p");

            TH1D *heffPromptRatio = (TH1D*)heffPrompt0->Clone("heffPromptRatio");
            heffPromptRatio->Divide(heffPrompt0,heffPrompt0,1,1,"b");
            heffPromptRatio->GetYaxis()->SetRangeUser(0.4,1.6);
            cEff[nHist]->cd();
            heffPromptRatio->DrawCopy();

            TH1D *heffNonPrompt0 = (TH1D*)file->Get(Form("%s",histName[nHist-1].Data()));
            if(!heffNonPrompt0) { cout << "NO HISTOGRAM !!!!!" << histName[nHist-1] << endl; return; }
            heffNonPrompt0->Sumw2();
            heffNonPrompt0->SetTitle(Form("%s",titles[nHist-1].Data()));
            heffNonPrompt0->SetLineColor(2);
            heffNonPrompt0->SetMarkerColor(2);
            heffNonPrompt0->SetMarkerStyle(29);
            heffNonPrompt0->SetMarkerSize(1.8);
            heffNonPrompt0->SetMaximum(0.6);
            heffNonPrompt0->SetMinimum(0.01);
            cEff[nHist-1]->cd();
            heffNonPrompt0->DrawCopy();
            leg3->AddEntry(heffNonPrompt0,"All HB","p");

            TH1D *heffNonPromptRatio = (TH1D*)heffNonPrompt0->Clone("heffNonPromptRatio");
            heffNonPromptRatio->Divide(heffNonPrompt0,heffNonPrompt0,1,1,"b");
            heffNonPromptRatio->GetYaxis()->SetRangeUser(0.4,1.6);
            cEff[nHist+1]->cd();
            heffNonPromptRatio->DrawCopy();

            for (int i=2; i<7; i++){
                TString fileName  = inDir;
                fileName += fileN;
                fileName += i;
                fileName += postfix;
                fileName += ".root";
                TFile *file = new TFile(Form("%s",fileName.Data()),"READ");
                if(file->IsZombie()) { cout << "NO FILE !!!!!" << fileName << endl; return; }

                TH1D *heffPrompt = (TH1D*)file->Get(Form("%s",histName[nHist-2].Data()));
                if(!heffPrompt) { cout << "NO HISTOGRAM !!!!!" << histName[nHist-2] << endl; return; }
                heffPrompt->Sumw2();
                heffPrompt->SetTitle();
                heffPrompt->SetLineColor(colors2[i-2]);
                heffPrompt->SetMarkerColor(colors2[i-2]);
                heffPrompt->SetMarkerStyle(markers2[i-2]);
                heffPrompt->SetMarkerSize(1.8);
                heffPrompt->SetMaximum(0.8);
                cEff[nHist-2]->cd();
                heffPrompt->DrawCopy("same");

                leg2->AddEntry(heffPrompt,Form("HB_%d",i),"p");

                TH1D *heffPromptRatio = (TH1D*)heffPrompt->Clone("heffPromptRatio");
                heffPromptRatio->Divide(heffPrompt,heffPrompt0,1,1,"b");
                heffPromptRatio->GetYaxis()->SetRangeUser(0.4,2);
                cEff[nHist]->cd();
                heffPromptRatio->DrawCopy("same");

                TH1D *heffNonPrompt = (TH1D*)file->Get(Form("%s",histName[nHist-1].Data()));
                if(!heffNonPrompt) { cout << "NO HISTOGRAM !!!!!" << histName[nHist-1] << endl; return; }
                heffNonPrompt->Sumw2();
                heffNonPrompt->SetTitle();
                heffNonPrompt->SetLineColor(colors2[i-2]);
                heffNonPrompt->SetMarkerColor(colors2[i-2]);
                heffNonPrompt->SetMarkerStyle(markers2[i-2]);
                heffNonPrompt->SetMarkerSize(1.8);
                heffNonPrompt->SetMaximum(0.5);
                cEff[nHist-1]->cd();
                heffNonPrompt->DrawCopy("same");

                leg3->AddEntry(heffPrompt,Form("HB_%d",i),"p");

                TH1D *heffNonPromptRatio = (TH1D*)heffNonPrompt->Clone("heffNonPromptRatio");
                heffNonPromptRatio->Divide(heffNonPrompt,heffNonPrompt0,1,1,"b");
                heffNonPromptRatio->GetYaxis()->SetRangeUser(0.4,2);
                cEff[nHist+1]->cd();
                heffNonPromptRatio->DrawCopy("same");

            }

            cEff[nHist-2]->cd();
            heffPrompt0->DrawCopy("same");
            leg2->Draw("same");
            SaveCanvas(cEff[nHist-2],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),"effPrompt"));
            cEff[nHist-1]->cd();
            heffNonPrompt0->DrawCopy("same");
            leg3->Draw("same");
            SaveCanvas(cEff[nHist-1],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),"effNonPrompt"));

            cEff[nHist]->cd();
            leg2->Draw("same");
            SaveCanvas(cEff[nHist],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),"effPromptRatio"));
            cEff[nHist+1]->cd();
            leg3->Draw("same");
            SaveCanvas(cEff[nHist+1],Form("%s/%s%s_%s",outName.Data(),outHistName.Data(),postfix.Data(),"effNonPromptRatio"));

            return;

}

void SaveCanvas(TCanvas *c, TString name = "tmp"){

    c->SaveAs(Form("%s.png",name.Data()));
    c->SaveAs(Form("%s.pdf",name.Data()));

}
