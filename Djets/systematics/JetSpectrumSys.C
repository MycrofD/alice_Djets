
#include <string>
#include <sstream>
#include <iostream>

#include "sys.h"
#include "config.h"


void JetSpectrumSys(int reg=3,  TString inDirBase = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts", TString input = "Default_jetMeas3_50_jetTrue3_50", bool isChain = 0, TString int measmin=3, int measmax=50, int truemin=5, int truemax=50)
{

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

  TString inDir = inDirBase;
  inDir += "/";
  inDir += input;
  gSystem->Exec(Form("mkdir %s/systematics",inDir.Data()));


  FDsys(reg,inDir);
  comparePriors(reg,inDir);
  compareBkgM(reg,inDir);
  compareJES(reg,inDir);
//  compareRanges(reg,inDirBase,measmin,measmax,truemin,truemax);


return;
}

void compareJES(int reg = 3 , TString inDir = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50")
{

        gStyle->SetOptStat(0000); //Mean and RMS shown
        gSystem->Exec(Form("mkdir %s/systematics",inDir.Data()));

        const int nFiles = 4;
        TString tab[nFiles-1] = { "96", "95","90"};
        TString dirName[nFiles];
        dirName[0] = inDir;
        dirName[0] += "_JES";
        dirName[0] += "/unfolding_Bayes_";
        dirName[0] += reg;

        dirName[1] = inDir;
        dirName[1] += "_JES96";
        dirName[1] += "/unfolding_Bayes_";
        dirName[1] += reg;

        dirName[2] = inDir;
        dirName[2] += "_JES95";
        dirName[2] += "/unfolding_Bayes_";
        dirName[2] += reg;

        dirName[3] = inDir;
        dirName[3] += "_JES90";
        dirName[3] += "/unfolding_Bayes_";
        dirName[3] += reg;

        TString desc[nFiles] = {"central","inefficiency 4%","inefficiency 5%","inefficiency 10%"};

        double plotmin = ptJetbins[0], plotmax = ptJetbins[nJetBins];
          double plotmin = 5;

        TFile *fproj[nFiles];
        for(int i=0; i<nFiles; i++) fproj[i] = new TFile(Form("%s/unfoldedSpectrum_unfoldedJetSpectrum.root",dirName[i].Data()),"READ");

        TCanvas *cspec = new TCanvas("cspec","cspec",800,600);
        cspec->SetLogy();

        TLegend *leg = new TLegend(0.5,0.65,0.85,0.8);
        leg->SetBorderSize(0);

        TH1F *spec[nFiles];
        for(int i=0; i<nFiles; i++) {
            spec[i] = (TH1F*)fproj[i]->Get("unfoldedSpectrum");
            spec[i]->Sumw2();
            spec[i] -> Scale(1,"width");
            spec[i]->SetTitle();
            spec[i]->SetLineColor(colors[i]);
            spec[i]->SetMarkerColor(colors[i]);
            spec[i]->SetMarkerStyle(markers[i]);
            spec[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
            if(!i) spec[i]->Draw();
            else spec[i]->Draw("same");
            leg->AddEntry(spec[i],desc[i].Data());
        }
        leg->Draw("same");

      cspec->SaveAs(Form("%s/systematics/JES_reg%d.pdf",inDir.Data(),reg));
      cspec->SaveAs(Form("%s/systematics/JES_reg%d.png",inDir.Data(),reg));


        TLegend *leg2 = new TLegend(0.6,0.65,0.85,0.85);
        leg2->SetBorderSize(0);
        TCanvas *cspec2 = new TCanvas("cspec2","cspec2",800,400);
        TH1F *hratio[nFiles-1];
        TH1F *hratiof[nFiles-1];

        TF1 *fr[nFiles-1];
        for(int i=0; i<nFiles-1; i++){
            hratio[i] = (TH1F*)spec[i+1]->Clone( Form("hratio_%d",i));
            //TH1D * hEff = (TH1D*)hMCpt_reco->Clone("hEff");
          	hratio[i] -> Divide(spec[i+1],spec[0],1,1,"b");
            hratio[i]->SetLineStyle(linestyle[i]);
            hratio[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
            hratio[i]->GetYaxis()->SetRangeUser(0.82,1.1);
            hratio[i]->GetYaxis()->SetTitle(Form("ratio to central (%s)",desc[0].Data()));

            fr[i] = new TF1(Form("fr_%d",i),"[0]*x+[1]",4,plotmax);
            fr[i]->SetLineStyle(1);
            fr[i]->SetLineWidth(2);
            fr[i]->SetLineColor(colors[i+1]);
            hratio[i]->Fit(fr[i],"EM0","",6,plotmax);

            //if(!i) hratio[i]->Draw("hist");
            //else hratio[i]->Draw("histsame");
            if(!i) hratio[i]->Draw();
            else hratio[i]->Draw("same");
             fr[i]->Draw("same");
            leg2->AddEntry(hratio[i],desc[i+1].Data());
        }
        leg2->Draw("same");

        TLine *line = new TLine(plotmin,1,plotmax,1);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");


        cspec2->SaveAs(Form("%s/systematics/JES_reg%d_ratio.pdf",inDir.Data(),reg));
        cspec2->SaveAs(Form("%s/systematics/JES_reg%d_ratio.png",inDir.Data(),reg));


    TCanvas *cspecf2 = new TCanvas("cspecf2","cspecf2",800,400);
    for(int i=0; i<nFiles-1; i++){
        double value = 0;
        hratiof[i] = (TH1F*)hratio[i]->Clone(Form("hratiof_%d",i));
        for(int j=0; j<hratio[i]->GetNbinsX();j++){
            value = (1 - fr[i]->Eval(hratio[i]->GetBinCenter(j+1))) *100;
            hratiof[i]->SetBinContent(j+1,value);
            //if(!i)cout << "bin: " << j+1 << "\t\t 4% value in: " << hratio[i]->GetBinCenter(j+1) << ":\t\t" << value << endl;
        }
        hratiof[i]->GetYaxis()->SetRangeUser(-1,20);
        hratiof[i]->GetYaxis()->SetTitle("unc.from fit [%]");
        if(!i) hratiof[i]->Draw();
        else hratiof[i]->Draw("same");

    }
    leg2->Draw("same");

    for(int j=0; j<hratiof[0]->GetNbinsX();j++){
        cout << "bin: " << j+1 << "\t\t 4% value in: " << hratiof[0]->GetBinCenter(j+1) << ":\t\t" << hratiof[0]->GetBinContent(j+1) << endl;
    }

    cspecf2->SaveAs(Form("%s/systematics/JES_reg%d_unc.pdf",inDir.Data(),reg));
    cspecf2->SaveAs(Form("%s/systematics/JES_reg%d_unc.png",inDir.Data(),reg));

        return;

}

void compareRanges(int reg, TString inDirBase, int measmin, int measmax, int truemin, int truemax)
{

            TString out = inDirBase;
            out+= "/systematics";
            gStyle->SetOptStat(0000); //Mean and RMS shown
            gSystem->Exec(Form("mkdir %s",out.Data()));

            const int nFiles = 4;
            TString inDir[nFiles] = {
              "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue5_50",
              "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50",
              "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas4_50_jetTrue5_50",
              "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas5_50_jetTrue5_50"
            }

            TString dirName[nFiles];
            //int regList[nFiles] = {3,4};

            for (int i=0; i<nFiles; i++){
                dirName[i] = inDir[i];
                dirName[i] += "/unfolding_Bayes_";
                dirName[i] += reg;
            }

            TString desc[nFiles] = {
              "meas: 3-50, true: 5-50",
              "meas: 3-50, true: 3-50",
              "meas: 4-50, true: 5-50",
              "meas: 5-50, true: 5-50"
            };
          //  for(int i=0; i<nFiles; i++){
          //    desc[i] = "reg=";
          //    desc[i] += regList[i];

          //  }

            double plotmin = 3, plotmax = 50;
            const int ptbinsJetMeasN = 7;
            double ptbinsJetMeasA[ptbinsJetMeasN+1] = { 5,6,8,10,14,20,30,50 };

            TFile *fproj[nFiles];
            for(int i=0; i<nFiles; i++) fproj[i] = new TFile(Form("%s/unfoldedSpectrum_unfoldedJetSpectrum.root",dirName[i].Data()),"READ");

            TCanvas *cspec = new TCanvas("cspec","cspec",800,600);
            cspec->SetLogy();

            TLegend *leg = new TLegend(0.5,0.6,0.85,0.8);
            leg->SetBorderSize(0);

            TH1F *spec[nFiles];
            TH1F *specReb[nFiles];
            for(int i=0; i<nFiles; i++) {
                spec[i] = (TH1F*)fproj[i]->Get("unfoldedSpectrum");
                spec[i]->Sumw2();
                spec[i] -> Scale(1,"width");
                spec[i]->SetTitle();
                spec[i]->SetLineColor(colors[i]);
                spec[i]->SetMarkerColor(colors[i]);
                spec[i]->SetMarkerStyle(markers[i]);

                specReb[i] = new TH1F(Form("specReb%d",i),"specReb",ptbinsJetMeasN,ptbinsJetMeasA);
                for(int j=1;j<specReb[i]->GetNbinsX()+1;j++){
                    double pt = specReb[i]->GetBinCenter(j);
                    int bin = spec[i]->GetXaxis()->FindBin(pt);
                    double value = spec[i]->GetBinContent(bin);
                    double error = spec[i]->GetBinError(bin);
                    specReb[i]->SetBinContent(j,value);
                    specReb[i]->SetBinError(j,error);
                }

                specReb[i]->SetTitle();
                specReb[i]->SetLineColor(colors[i]);
                specReb[i]->SetMarkerColor(colors[i]);
                specReb[i]->SetMarkerStyle(markers[i]);

                spec[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
                if(!i) spec[i]->Draw();
                else spec[i]->Draw("same");
                leg->AddEntry(spec[i],desc[i].Data());
            }
            leg->Draw("same");

          cspec->SaveAs(Form("%s/UnfoldingRangesComparison_ptMeas_%d_%d_ptTrue_%d_%d_reg%d.pdf",out.Data(),measmin,measmax,truemin,truemax,reg));
          cspec->SaveAs(Form("%s/UnfoldingRangesComparison_ptMeas_%d_%d_ptTrue_%d_%d_reg%d.png",out.Data(),measmin,measmax,truemin,truemax,reg));
          TLegend *leg2 = new TLegend(0.55,0.55,0.9,0.85);

            leg2->SetBorderSize(0);
            TCanvas *cspec2 = new TCanvas("cspec2","cspec2",800,400);
            TH1F *hratio[nFiles-1];
            for(int i=0; i<nFiles-1; i++){
                hratio[i] = (TH1F*)specReb[i+1]->Clone( Form("hratio_%d",i));
                hratio[i]->Divide(specReb[0]);
                hratio[i]->SetLineStyle(linestyle[i]);
                  hratio[i]->SetLineWidth(2);
                hratio[i]->GetXaxis()->SetRangeUser(ptbinsJetMeasA[0],ptbinsJetMeasA[ptbinsJetMeasN]);
                hratio[i]->GetYaxis()->SetRangeUser(0.95,1.4);
                hratio[i]->GetYaxis()->SetTitle(Form("ratio to central (%s)",desc[0].Data()));
                if(!i) hratio[i]->Draw("hist");
                else hratio[i]->Draw("samehist");
                leg2->AddEntry(hratio[i],desc[i+1].Data());
            }
            leg2->Draw("same");

            TLine *line = new TLine(ptbinsJetMeasA[0],1,ptbinsJetMeasA[ptbinsJetMeasN],1);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw("same");


            cspec2->SaveAs(Form("%s/UnfoldingRangesComparison_ptMeas_%d_%d_ptTrue_%d_%d_reg%d_ratio.pdf",out.Data(),measmin,measmax,truemin,truemax,reg));
            cspec2->SaveAs(Form("%s/UnfoldingRangesComparison_ptMeas_%d_%d_ptTrue_%d_%d_reg%d_ratio.png",out.Data(),measmin,measmax,truemin,truemax,reg));

            return;

}

void compareReg(TString inDir = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50")
{

            gStyle->SetOptStat(0000); //Mean and RMS shown
            gSystem->Exec(Form("mkdir %s/systematics",inDir.Data()));

            const int nFiles = 2;
            TString dirName[nFiles];
            int regList[nFiles] = {3,4};

            for (int i=0; i<nFiles; i++){
                dirName[i] = inDir;
                dirName[i] += "/unfolding_Bayes_";
                dirName[i] += regList[i];
            }

            TString desc[nFiles];
            for(int i=0; i<nFiles; i++){
              desc[i] = "reg=";
              desc[i] += regList[i];

            }

            double plotmin = ptJetbins[0], plotmax = ptJetbins[nJetBins];

            TFile *fproj[nFiles];
            for(int i=0; i<nFiles; i++) fproj[i] = new TFile(Form("%s/unfoldedSpectrum_unfoldedJetSpectrum.root",dirName[i].Data()),"READ");

            TCanvas *cspec = new TCanvas("cspec","cspec",800,600);
            cspec->SetLogy();

            TLegend *leg = new TLegend(0.5,0.6,0.85,0.8);
            leg->SetBorderSize(0);

            TH1F *spec[nFiles];
            for(int i=0; i<nFiles; i++) {
                spec[i] = (TH1F*)fproj[i]->Get("unfoldedSpectrum");
                spec[i]->Sumw2();
                spec[i] -> Scale(1,"width");
                spec[i]->SetTitle();
                spec[i]->SetLineColor(colors[i]);
                spec[i]->SetMarkerColor(colors[i]);
                spec[i]->SetMarkerStyle(markers[i]);
                spec[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
                if(!i) spec[i]->Draw();
                else spec[i]->Draw("same");
                leg->AddEntry(spec[i],desc[i].Data());
            }
            leg->Draw("same");

          cspec->SaveAs(Form("%s/systematics/RegComparison_reg.pdf",inDir.Data()));
          cspec->SaveAs(Form("%s/systematics/RegComparison_reg.png",inDir.Data()));
          TLegend *leg2 = new TLegend(0.55,0.75,0.9,0.85);

            leg2->SetBorderSize(0);
            TCanvas *cspec2 = new TCanvas("cspec2","cspec2",800,400);
            TH1F *hratio[nFiles-1];
            for(int i=0; i<nFiles-1; i++){
                hratio[i] = (TH1F*)spec[i+1]->Clone( Form("hratio_%d",i));
                hratio[i]->Divide(spec[0]);
                hratio[i]->SetLineStyle(linestyle[i]);
                hratio[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
                hratio[i]->GetYaxis()->SetRangeUser(0.8,1.15);
                hratio[i]->GetYaxis()->SetTitle(Form("ratio to central (%s)",desc[0].Data()));
                if(!i) hratio[i]->Draw("hist");
                else hratio[i]->Draw("histsame");
                leg2->AddEntry(hratio[i],desc[i+1].Data());
            }
            leg2->Draw("same");

            TLine *line = new TLine(plotmin,1,plotmax,1);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw("same");


            cspec2->SaveAs(Form("%s/systematics/RegComparison_reg_ratio.pdf",inDir.Data()));
            cspec2->SaveAs(Form("%s/systematics/RegComparison_reg_ratio.png",inDir.Data()));

            return;

}


void comparePriors(int reg = 4 , TString inDir = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50")
{

        gStyle->SetOptStat(0000); //Mean and RMS shown
        gSystem->Exec(Form("mkdir %s/systematics",inDir.Data()));

        const int nFiles = 10;
        TString dirName[nFiles];
        dirName[0] = inDir;
        dirName[0] += "/unfolding_Bayes_";
        dirName[0] += reg;
        for (int i=1; i<nFiles; i++){
            dirName[i] = inDir;
            dirName[i] += "/unfolding_Bayes_";
            dirName[i] += reg;
            dirName[i] += "_priorType";
            dirName[i] += i-1;

        }

        TString desc[nFiles];
        desc[0] = "Bayes reg=";
        desc[0] += reg;
        desc[0] += " prior true";
        desc[nFiles-1] = "prior meas fit";
        for(int i=1; i< nFiles-1; i++){
          desc[i] = "prior=";
          desc[i] += i-1;

        }

        double plotmin = ptJetbins[0], plotmax = ptJetbins[nJetBins];

        TFile *fproj[nFiles];
        for(int i=0; i<nFiles; i++) fproj[i] = new TFile(Form("%s/unfoldedSpectrum_unfoldedJetSpectrum.root",dirName[i].Data()),"READ");

        TCanvas *cspec = new TCanvas("cspec","cspec",800,600);
        cspec->SetLogy();

        TLegend *leg = new TLegend(0.5,0.6,0.85,0.8);
        leg->SetBorderSize(0);

        TH1F *spec[nFiles];
        for(int i=0; i<nFiles; i++) {
            spec[i] = (TH1F*)fproj[i]->Get("unfoldedSpectrum");
            spec[i]->Sumw2();
            spec[i] -> Scale(1,"width");
            spec[i]->SetTitle();
            spec[i]->SetLineColor(colors[i]);
            spec[i]->SetMarkerColor(colors[i]);
            spec[i]->SetMarkerStyle(markers[i]);
            spec[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
            if(!i) spec[i]->Draw();
            else spec[i]->Draw("same");
            leg->AddEntry(spec[i],desc[i].Data());
        }
        leg->Draw("same");

      cspec->SaveAs(Form("%s/systematics/PriorComparison_reg%d.pdf",inDir.Data(),reg));
      cspec->SaveAs(Form("%s/systematics/PriorComparison_reg%d.png",inDir.Data(),reg));

        TLegend *leg2 = new TLegend(0.55,0.5,0.9,0.88);
        leg2->SetBorderSize(0);
        TCanvas *cspec2 = new TCanvas("cspec2","cspec2",800,400);
        TH1F *hratio[nFiles-1];
        for(int i=0; i<nFiles-1; i++){
            hratio[i] = (TH1F*)spec[i+1]->Clone( Form("hratio_%d",i));
            hratio[i]->Divide(spec[0]);
            hratio[i]->SetLineStyle(linestyle[i]);
            hratio[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
            hratio[i]->GetYaxis()->SetRangeUser(0.8,1.5);
            hratio[i]->GetYaxis()->SetTitle(Form("ratio to central (%s)",desc[0].Data()));
            if(!i) hratio[i]->Draw("hist");
            else hratio[i]->Draw("histsame");
            leg2->AddEntry(hratio[i],desc[i+1].Data());
        }
        leg2->Draw("same");

        TLine *line = new TLine(plotmin,1,plotmax,1);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");


        cspec2->SaveAs(Form("%s/systematics/PriorComparison_reg%d_ratio.pdf",inDir.Data(),reg));
        cspec2->SaveAs(Form("%s/systematics/PriorComparison_reg%d_ratio.png",inDir.Data(),reg));

        TH1F *hsys = new TH1F("hsys","syst. rms; p_{T,ch jet};  sys [%] (rms)",nJetBins,ptJetbins);
        TH1F *hmean = (TH1F*)hsys->Clone("hmean");
        getRMS(nFiles,hratio,hmean,hsys);

        hsys->GetYaxis()->SetRangeUser(0,10);
        hsys->SetLineColor(kViolet+2);
        TCanvas *cspecRMS = new TCanvas("cspecRMS","cspecRMS",800,400);
        hsys->Draw("hist");

        cspecRMS->SaveAs(Form("%s/systematics/PriorComparison_reg%d_rms.pdf",inDir.Data(),reg));
        cspecRMS->SaveAs(Form("%s/systematics/PriorComparison_reg%d_rms.png",inDir.Data(),reg));

        hmean->GetYaxis()->SetRangeUser(0.93,1.08);
        hmean->SetLineColor(kMagenta+1);
        TCanvas *cspecMean = new TCanvas("cspecMean","cspecMean",800,400);
        hmean->Draw("hist");
        line->Draw("same");

        cspecMean->SaveAs(Form("%s/systematics/PriorComparison_reg%d_mean.pdf",inDir.Data(),reg));
        cspecMean->SaveAs(Form("%s/systematics/PriorComparison_reg%d_mean.png",inDir.Data(),reg));

        return;

}

void FDsys(int reg = 3 , TString inDir = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50")
{

        gStyle->SetOptStat(0000); //Mean and RMS shown
        gSystem->Exec(Form("mkdir %s/systematics",inDir.Data()));

        const int nFiles = 3;

        TString dirName[nFiles];
        dirName[0] = inDir;
        dirName[0] += "/unfolding_Bayes_";
        dirName[0] += reg;

        dirName[1] = inDir;
        dirName[1] += "/unfolding_Bayes_";
        dirName[1] += reg;
        dirName[1] += "_FDsysUp";

        dirName[2] = inDir;
        dirName[2] += "/unfolding_Bayes_";
        dirName[2] += reg;
        dirName[2] += "_FDsysDown";

        TString desc[nFiles] = {"central","FD Up","FD Down"};

        double plotmin = ptJetbins[0], plotmax = ptJetbins[nJetBins];
      //    double plotmin = 5;

        TFile *fproj[nFiles];
        for(int i=0; i<nFiles; i++) fproj[i] = new TFile(Form("%s/unfoldedSpectrum_unfoldedJetSpectrum.root",dirName[i].Data()),"READ");

        TCanvas *cspec = new TCanvas("cspec","cspec",800,600);
        cspec->SetLogy();

        TLegend *leg = new TLegend(0.5,0.65,0.85,0.8);
        leg->SetBorderSize(0);

        TH1F *spec[nFiles];
        for(int i=0; i<nFiles; i++) {
            spec[i] = (TH1F*)fproj[i]->Get("unfoldedSpectrum");
            spec[i]->Sumw2();
            spec[i] -> Scale(1,"width");
            spec[i]->SetTitle();
            spec[i]->SetLineColor(colors[i]);
            spec[i]->SetMarkerColor(colors[i]);
            spec[i]->SetMarkerStyle(markers[i]);
            spec[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
            if(!i) spec[i]->Draw();
            else spec[i]->Draw("same");
            leg->AddEntry(spec[i],desc[i].Data());
        }
        leg->Draw("same");

      cspec->SaveAs(Form("%s/systematics/FD_reg%d.pdf",inDir.Data(),reg));
      cspec->SaveAs(Form("%s/systematics/FD_reg%d.png",inDir.Data(),reg));

        TLegend *leg2 = new TLegend(0.6,0.3,0.85,0.45);
        leg2->SetBorderSize(0);
        TCanvas *cspec2 = new TCanvas("cspec2","cspec2",800,400);
        TH1F *hratio[nFiles-1];
        TH1F *hratiof[nFiles-1];

        for(int i=0; i<nFiles-1; i++){
            hratio[i] = (TH1F*)spec[i+1]->Clone( Form("hratio_%d",i));
            hratio[i]->Divide(spec[0]);
            hratio[i]->SetLineStyle(linestyle[i]);
            hratio[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
            hratio[i]->GetYaxis()->SetRangeUser(0.8,1.2);
            hratio[i]->GetYaxis()->SetTitle(Form("ratio to central (%s)",desc[0].Data()));
            if(!i) hratio[i]->Draw("hist");
            else hratio[i]->Draw("histsame");
            leg2->AddEntry(hratio[i],desc[i+1].Data());
        }
        leg2->Draw("same");

        TLine *line = new TLine(plotmin,1,plotmax,1);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");


        cspec2->SaveAs(Form("%s/systematics/FD_reg%d_ratio.pdf",inDir.Data(),reg));
        cspec2->SaveAs(Form("%s/systematics/FD_reg%d_ratio.png",inDir.Data(),reg));


//-------------- FD uncertanties
TH1D *hFDUnc = (TH1D*)*spec[0]->Clone("hDFUnc");
hFDUnc->GetYaxis()->SetTitle("FD sys. unc [%]");
hFDUnc->SetLineColor(kViolet+1);
for(int j=1; j<=nJetBins; j++ ){
        double unc1 = spec[1]->GetBinContent(j) - spec[0]->GetBinContent(j);
        double unc2 = spec[0]->GetBinContent(j) - spec[2]->GetBinContent(j);
        double unc = 0;
        if(unc1>unc2) unc = unc1;
        else unc = unc2;
        unc /= spec[0]->GetBinContent(j);
        hFDUnc->SetBinContent(j,unc*100);
        hFDUnc->SetBinError(j,0);

        //cout << "unc: " << unc << endl;
}
hFDUnc->GetYaxis()->SetRangeUser(3,20);
TCanvas *cspecf2 = new TCanvas("cspecf2","cspecf2",800,400);
hFDUnc->Draw("hist");

cspecf2->SaveAs(Form("%s/systematics/FD_reg%d_unc.pdf",inDir.Data(),reg));
cspecf2->SaveAs(Form("%s/systematics/FD_reg%d_unc.png",inDir.Data(),reg));

        return;

}

void compareBkgM(int reg = 4 , TString inDir = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50")
{


    gStyle->SetOptStat(0000); //Mean and RMS shown
    gSystem->Exec(Form("mkdir %s/systematics",inDir.Data()));

    const int nFiles = 12;
    TString dirName[nFiles];
    dirName[0] = inDir;
    dirName[0] += "/unfolding_Bayes_";
    dirName[0 ]+= reg;
    for (int i=1 ;i<nFiles; i++){
        dirName[i] = inDir;
        dirName[i] += "_BkgM";
        dirName[i] += i;
        dirName[i] += "/unfolding_Bayes_";
        dirName[i] += reg;

    }

    TString desc[nFiles];
    desc[0] = "Bayes reg=";
    desc[0] += reg;
    desc[0] += " def BkgM";
    for(int i=1; i< nFiles; i++){
      desc[i] = "Bayes BkgM ";
      desc[i] += i;

    }

    double plotmin = ptJetbins[0], plotmax = ptJetbins[nJetBins];

    TFile *fproj[nFiles];
    for(int i=0; i<nFiles; i++) fproj[i] = new TFile(Form("%s/unfoldedSpectrum_unfoldedJetSpectrum.root",dirName[i].Data()),"READ");

    TCanvas *cspec = new TCanvas("cspec","cspec",800,600);
    cspec->SetLogy();

    TLegend *leg = new TLegend(0.5,0.45,0.85,0.85);
    leg->SetBorderSize(0);

    TH1F *spec[nFiles];
    for(int i=0; i<nFiles; i++) {
        spec[i] = (TH1F*)fproj[i]->Get("unfoldedSpectrum");
        spec[i]->Sumw2();
        spec[i] -> Scale(1,"width");
        spec[i]->SetTitle();
        spec[i]->SetLineColor(colors[i]);
        spec[i]->SetMarkerColor(colors[i]);
        spec[i]->SetMarkerStyle(markers[i]);
        spec[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
        if(!i) spec[i]->Draw();
        else spec[i]->Draw("same");
        leg->AddEntry(spec[i],desc[i].Data());
    }
    leg->Draw("same");

  cspec->SaveAs(Form("%s/systematics/BkgComparison_reg%d.pdf",inDir.Data(),reg));
  cspec->SaveAs(Form("%s/systematics/BkgComparison_reg%d.png",inDir.Data(),reg));

    TLegend *leg2 = new TLegend(0.55,0.5,0.9,0.88);
    leg2->SetBorderSize(0);
    TCanvas *cspec2 = new TCanvas("cspec2","cspec2",800,400);
    TH1F *hratio[nFiles-1];
    for(int i=0; i<nFiles-1; i++){
        hratio[i] = (TH1F*)spec[i+1]->Clone( Form("hratio_%d",i));
        hratio[i]->Divide(spec[0]);
        hratio[i]->SetLineStyle(linestyle[i]);
        hratio[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
        hratio[i]->GetYaxis()->SetRangeUser(0.89,1.2);
        hratio[i]->GetYaxis()->SetTitle(Form("ratio to central (%s)",desc[0].Data()));
        if(!i) hratio[i]->Draw("hist");
        else hratio[i]->Draw("histsame");
        leg2->AddEntry(hratio[i],desc[i+1].Data());
    }
    leg2->Draw("same");

    TLine *line = new TLine(plotmin,1,plotmax,1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");


    cspec2->SaveAs(Form("%s/systematics/BkgComparison_reg%d_ratio.pdf",inDir.Data(),reg));
    cspec2->SaveAs(Form("%s/systematics/BkgComparison_reg%d_ratio.png",inDir.Data(),reg));

    TH1F *hsys = new TH1F("hsys","syst. rms; p_{T,ch jet};  sys [%] (rms)",nJetBins,ptJetbins);
    TH1F *hmean = (TH1F*)hsys->Clone("hmean");
    getRMS(nFiles,hratio,hmean,hsys);

    hsys->GetYaxis()->SetRangeUser(0,6);
    hsys->SetLineColor(kViolet+2);
    TCanvas *cspecRMS = new TCanvas("cspecRMS","cspecRMS",800,400);
    hsys->Draw("hist");

    cspecRMS->SaveAs(Form("%s/systematics/BkgComparison_reg%d_rms.pdf",inDir.Data(),reg));
    cspecRMS->SaveAs(Form("%s/systematics/BkgComparison_reg%d_rms.png",inDir.Data(),reg));

    hmean->GetYaxis()->SetRangeUser(0.93,1.08);
    hmean->SetLineColor(kMagenta+1);
    TCanvas *cspecMean = new TCanvas("cspecMean","cspecMean",800,400);
    hmean->Draw("hist");
    line->Draw("same");

    cspecMean->SaveAs(Form("%s/systematics/BkgComparison_reg%d_mean.pdf",inDir.Data(),reg));
    cspecMean->SaveAs(Form("%s/systematics/BkgComparison_reg%d_mean.png",inDir.Data(),reg));

    return;

/*

    TH1D *hsys = new TH1D("hsys","syst. rms; p_{T,ch jet};  sys [%] (rms)",nJetBins,ptJetbins);
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

    TH1D *hmean = (TH1D*)hsys->Clone("hmean");
    hmean->GetYaxis()->SetTitle("mean");
    hmean->GetYaxis()->SetRangeUser(0.95,1.1);
    hmean->SetMarkerStyle(20);
    hmean->SetLineStyle(1);


    double rms[nJetBins];
    double mean[nJetBins];
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
        hsys->SetBinContent(i+1,rms[i]);


        hmean->SetBinContent(i+1,mean[i]);
        cout << "RMS pT " << (ptJetbins[i]+ptJetbins[i+1])/2. << " GeV/c:\t" << rms[i]*100 << endl;
    }



TCanvas *csys = new TCanvas("csys","csys",800,400);
hsys->Draw("h");

TCanvas *cmean = new TCanvas("cmean","cmean",800,400);
hmean->Draw();


cspec->SaveAs("unfoldingSys/jetSpectra.pdf");
cspec->SaveAs("unfoldingSys/jetSpectra.png");

cspec2->SaveAs("unfoldingSys2/jetSpectraRatio.pdf");
cspec2->SaveAs("unfoldingSys2/jetSpectraRatio.png");

csys->SaveAs("unfoldingSys/jetSpectraSys.pdf");
csys->SaveAs("unfoldingSys/jetSpectraSys.png");

cmean->SaveAs("unfoldingSys/jetSpectraMean.pdf");
cmean->SaveAs("unfoldingSys/jetSpectraMean.png");

*/
return;

/*
if(ratio){
    cspec = new TCanvas("cspec","cspec",800,500);
    hratio->Draw();

}
else{
    cspec = new TCanvas("cspec","cspec",800,1000);
    cspec->Divide(1,2,0);
    cspec->cd(1);
    gPad->SetLogy();
    spec1->Draw();
    spec2->Draw("same");
    cspec->cd(2);
    hratio->Draw();
}
//fit->Draw("same");

*/

  //  pvProd->Draw("same");


//cspec->SaveAs(Form("%s/unfoldedSpectra_%s.pdf",outDir.Data(),out.Data()));
//cspec->SaveAs(Form("%s/unfoldedSpectra_%s.png",outDir.Data(),out.Data()));

// cProj2->SaveAs(Form("%s/DetMatrixResProjectionsComparison.pdf",outDir.Data()));



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
