
#include <string>
#include <sstream>
#include <iostream>

#include "sys.h"
#include "config.h"

Int_t colors2[] = {1,2,kGreen+3,kMagenta+2,4,6,kCyan+1,kYellow+2,kOrange-1,8,kGray+1,kViolet+5};
Int_t markers2[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};
Int_t linestyle2[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

const int nFiles = 7;
TString inDir[nFiles] = {
  "349",
  "359",
  "348",
  "347",
  "358",
  "249",
  "259"
};

TString desc[nFiles] = {
  "S:3#sigma, B:4-9#sigma (Def)",
  "S:3#sigma, B:5-9#sigma",
  "S:3#sigma, B:4-8#sigma",
  "S:3#sigma, B:4-7#sigma",
  "S:3#sigma, B:5-8#sigma",
  "S:2#sigma, B:4-9#sigma",
  "S:2#sigma, B:5-9#sigma"
};

//const int nFilesCut = 6;
//TString inDirCut[nFilesCut] = {
//  "Default_jetMeas3_50_jetTrue3_50_PbPbbinning",
//  "SB59_jetMeas3_50_jetTrue3_50_PbPbbinning",
//  "SB48_jetMeas3_50_jetTrue3_50_PbPbbinning",
//  "SB47_jetMeas3_50_jetTrue3_50_PbPbbinning",
//  "S2sigma_SB49_jetMeas3_50_jetTrue3_50_PbPbbinning",
//  "S2sigma_SB59_jetMeas3_50_jetTrue3_50_PbPbbinning"
//};
//
//TString descCut[nFilesCut] = {
//  "S:3#sigma, B:4-9#sigma (Def)",
//  "S:3#sigma, B:5-9#sigma",
//  "S:3#sigma, B:4-8#sigma",
//  "S:3#sigma, B:4-7#sigma",
//  "S:2#sigma, B:4-9#sigma",
//  "S:2#sigma, B:5-9#sigma"
//};

double plotmin = 5, plotmax = 50;
const int ptbinsN = 6;
//double ptbinsA[ptbinsN+1] = { 5,6,8,10,14,20,30,50 };
double ptbinsA[ptbinsN+1] = { 5,10,15,20,25,35,50 };

int nJetBins2 = 6;
//double ptJetbins2[] = {5,6,8,10,14,20,30,50};
double ptJetbins2[] = { 5,10,15,20,25,35,50 };

void BkgSRangesComparison2(int reg=3,  TString inDirBase = "/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet", TString inName = "systematics")
{

  TString inputDir = inDirBase;
//  inDir += "/";
//  inDir += input;
  gSystem->Exec(Form("mkdir %s/%s",inputDir.Data(),inName.Data()));

 unfSpectra(inName,inDirBase);
// cutSys(inName,inDirBase);

return;
}

void unfSpectra(TString inName, TString inDirBase)
{

            TString out = inDirBase;
            out += "/";
            out += inName;
            gStyle->SetOptStat(0000); //Mean and RMS shown
            gSystem->Exec(Form("mkdir %s",out.Data()));

            TString dirName[nFiles];

            for (int i=0; i<nFiles; i++){
                dirName[i] = inDirBase;
                dirName[i] += "/";
                dirName[i] += "results_";
                dirName[i] += inDir[i];
            }

            TFile *outFile = new TFile(Form("%s/SBRangesComparison.root",out.Data()),"RECREATE");

            TFile *fproj[nFiles];
            for(int i=0; i<nFiles; i++) {
              fproj[i] = new TFile(Form("%s/DzeroR03_pPbCuts/Default/unfolding_Bayes_4/unfoldedSpectrum_unfoldedJetSpectrum.root",dirName[i].Data()),"READ");
              //fproj[i] = new TFile(Form("%s/FDsubtraction/JetPtSpectrum_FDsub.root",dirName[i].Data()),"READ");

              if(!fproj[i]) { cout << "NO FILE !!!!!" << endl; return; }
            }

            TCanvas *cspec = new TCanvas("cspec","cspec",1200,800);
            cspec->SetLogy();

            TLegend *leg = new TLegend(0.55,0.5,0.85,0.85);
            leg->SetBorderSize(0);

            TH1F *spec[nFiles];
            TH1F *specUnc[nFiles];
            TH1F *specUnc2[nFiles];
            TH1F *specUnc3[nFiles];
            TH1F *specReb[nFiles];

            for(int i=0; i<nFiles; i++) {
                spec[i] = (TH1F*)fproj[i]->Get("unfoldedSpectrum");
                //spec[i] = (TH1F*)fproj[i]->Get("hData_binned_sub");
                spec[i]->Sumw2();
                spec[i] -> Scale(1,"width");
                spec[i]->SetTitle();
                spec[i]->SetLineColor(colors2[i]);
                spec[i]->SetMarkerColor(colors2[i]);
                spec[i]->SetMarkerStyle(markers2[i]);
                specUnc[i] = (TH1F*)spec[i]->Clone(Form("specUnc_%d",i));
                specUnc2[i] = (TH1F*)spec[i]->Clone(Form("specUnc2_%d",i));
                specUnc3[i] = (TH1F*)spec[i]->Clone(Form("specUnc3_%d",i));
                for(int j=1;j<spec[i]->GetNbinsX()+1;j++){
                    double err;
    								if(spec[i]->GetBinContent(j)) err = spec[i]->GetBinError(j)/spec[i]->GetBinContent(j);
    								else err = 0;
                    specUnc[i]->SetBinContent(j,err);
                    specUnc[i]->SetBinError(j,0);
                    specUnc2[i]->SetBinContent(j,1+err);
                    specUnc2[i]->SetBinError(j,0);
                    specUnc3[i]->SetBinContent(j,1-err);
                    specUnc3[i]->SetBinError(j,0);
                }
                specUnc[i]->GetYaxis()->SetRangeUser(0,0.40);
                specUnc[i]->GetYaxis()->SetTitle("rel.sta.unc.");

                specReb[i] = new TH1F(Form("specReb%d",i),"specReb",ptbinsN,ptbinsA);
                for(int j=1;j<specReb[i]->GetNbinsX()+1;j++){
                    double pt = specReb[i]->GetBinCenter(j);
                    int bin = spec[i]->GetXaxis()->FindBin(pt);
                    double value = spec[i]->GetBinContent(bin);
                    double error = spec[i]->GetBinError(bin);
                    specReb[i]->SetBinContent(j,value);
                    specReb[i]->SetBinError(j,error);
                }

                specReb[i]->SetTitle();
                specReb[i]->SetLineColor(colors2[i]);
                specReb[i]->SetMarkerColor(colors2[i]);
                specReb[i]->SetMarkerStyle(markers2[i]);

                spec[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
                if(!i) spec[i]->Draw();
                else spec[i]->Draw("same");
                leg->AddEntry(spec[i],desc[i].Data());
            }

            leg->Draw("same");

          cspec->SaveAs(Form("%s/SBRangesComparison.pdf",out.Data()));
          cspec->SaveAs(Form("%s/SBRangesComparison.png",out.Data()));

          TCanvas *cspecUnc = new TCanvas("cspecUnc","cspecUnc",1200,800);
          TLegend *leg3 = new TLegend(0.15,0.45,0.45,0.85);
          leg3->SetBorderSize(0);

          for(int i=0; i<nFiles; i++) {

              if(!i) specUnc[i]->Draw();
              else specUnc[i]->Draw("same");
              leg3->AddEntry(specUnc[i],desc[i].Data());
          }
          leg3->Draw("same");

        cspecUnc->SaveAs(Form("%s/SBRangesUncComparison.pdf",out.Data()));
        cspecUnc->SaveAs(Form("%s/SBRangesUncComparison.png",out.Data()));


            TLegend *leg2 = new TLegend(0.15,0.15,0.5,0.4);
            leg2->SetBorderSize(0);
            TLegend *leg4 = new TLegend(0.45,0.2,0.65,0.25);
            leg4->SetBorderSize(0);
            TLegend *leg5 = new TLegend(0.15,0.75,0.65,0.85,"Signal extraction with different S and SB ranges");
            leg5->SetBorderSize(0);
            TCanvas *cspec2 = new TCanvas("cspec2","cspec2",1200,800);
            TH1F *hratio[nFiles-1];
            for(int i=0; i<nFiles-1; i++){
                //hratio[i] = (TH1F*)specReb[i+1]->Clone( Form("hratio_%d",i));
                //hratio[i]->Divide(specReb[0]);
                hratio[i] = (TH1F*)spec[i+1]->Clone( Form("hratio_%d",i));
                hratio[i]->Divide(spec[0]);
                hratio[i]->SetLineStyle(linestyle2[i]);
                hratio[i]->SetLineWidth(2);
                hratio[i]->GetXaxis()->SetRangeUser(ptbinsA[0],ptbinsA[ptbinsN]);
                hratio[i]->GetYaxis()->SetRangeUser(0.6,1.4);
                hratio[i]->GetYaxis()->SetTitle(Form("ratio to central (%s)",desc[0].Data()));
                if(!i) hratio[i]->Draw("hist");
                else hratio[i]->Draw("samehist");
                leg2->AddEntry(hratio[i],desc[i+1].Data(),"l");
            }
            leg2->Draw("same");
            specUnc2[0]->SetLineColor(1);
            specUnc2[0]->SetLineStyle(2);
            specUnc2[0]->Draw("same");
            specUnc3[0]->SetLineColor(1);
            specUnc3[0]->SetLineStyle(2);
            specUnc3[0]->Draw("same");
            leg4->AddEntry(specUnc2[0],"stat.unc.","l");
            leg4->Draw("same");
            leg5->Draw("same");

            TLine *line = new TLine(ptbinsA[0],1,ptbinsA[ptbinsN],1);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw("same");

            cspec2->SaveAs(Form("%s/SBRangesComparison_ratio.pdf",out.Data()));
            cspec2->SaveAs(Form("%s/SBRangesComparison_ratio.png",out.Data()));

            TH1F *hsys = new TH1F("hsys","syst. rms; p_{T,ch jet};  RMS [%]",ptbinsN,ptbinsA);
            TH1F *hmean = (TH1F*)hsys->Clone("hmean");
            getRMS(nFiles,hratio,hmean,hsys);

            hsys->GetYaxis()->SetRangeUser(0,15);
            hsys->SetLineColor(kViolet+2);
            TCanvas *cspecRMS = new TCanvas("cspecRMS","cspecRMS",800,400);
            hsys->Draw("hist");
            leg5->Draw("same");

            cspecRMS->SaveAs(Form("%s/SBRangesComparison_rms.pdf",out.Data()));
            cspecRMS->SaveAs(Form("%s/SBRangesComparison_rms.png",out.Data()));

            hmean->GetYaxis()->SetRangeUser(0.93,1.08);
            hmean->SetLineColor(kMagenta+1);
            TCanvas *cspecMean = new TCanvas("cspecMean","cspecMean",800,400);
            hmean->Draw("hist");
            line->Draw("same");

            cspecMean->SaveAs(Form("%s/SBRangesComparison_mean.pdf",out.Data()));
            cspecMean->SaveAs(Form("%s/SBRangesComparison_mean.png",out.Data()));

            outFile->cd();
            cspec->Write();
            cspec2->Write();
            cspecRMS->Write();
            cspecMean->Write();
            hmean->Write();
            hsys->Write();
            outFile->Close();

            return;

}


void cutSys(TString inName, TString inDirBase)
{

            TString out = inDirBase;
            out += "/";
            out += inName;
            gStyle->SetOptStat(0000); //Mean and RMS shown
            gSystem->Exec(Form("mkdir %s",out.Data()));

            TFile *outFile = new TFile(Form("%s/SBRangesComparisonCutSys.root",out.Data()),"RECREATE");

            TString dirName[nFilesCut];
            for (int i=0; i<nFilesCut; i++){
                dirName[i] = inDirBase;
                dirName[i] += "/";
                dirName[i] += inDirCut[i];
            }

            TFile *fproj[nFilesCut];
            for(int i=0; i<nFilesCut; i++) {
              fproj[i] = new TFile(Form("%s/systematics/cutSystematics.root",dirName[i].Data()),"READ");
              if(!fproj[i]) { cout << "NO FILE !!!!!" << endl; return; }
            }

            TCanvas *cspec = new TCanvas("cspec","cspec",1200,800);
          //  cspec->SetLogy();

            TLegend *leg = new TLegend(0.15,0.55,0.55,0.85);
            leg->SetBorderSize(0);

            TH1F *spec[nFilesCut];
            TH1F *specReb[nFilesCut];

            for(int i=0; i<nFilesCut; i++) {
                spec[i] = (TH1F*)fproj[i]->Get("cutSysRMS");
                spec[i]->Sumw2();
                //spec[i] -> Scale(1,"width");
                spec[i]->SetTitle();
                spec[i]->SetLineColor(colors2[i]);
                spec[i]->SetMarkerColor(colors2[i]);
                spec[i]->SetMarkerStyle(markers2[i]);
                spec[i]->SetLineStyle(linestyle2[i]);
                spec[i]->GetYaxis()->SetRangeUser(0.,50);
                //spec[i]->GetYaxis()->SetRangeUser(0.,0.5);

                specReb[i] = new TH1F(Form("specReb%d",i),"specReb",ptbinsN,ptbinsA);
                for(int j=1;j<specReb[i]->GetNbinsX()+1;j++){
                    double pt = specReb[i]->GetBinCenter(j);
                    int bin = spec[i]->GetXaxis()->FindBin(pt);
                    double value = spec[i]->GetBinContent(bin);
                    double error = spec[i]->GetBinError(bin);
                    specReb[i]->SetBinContent(j,value);
                    specReb[i]->SetBinError(j,error);
                }

                specReb[i]->SetTitle();
                specReb[i]->SetLineColor(colors2[i]);
                specReb[i]->SetMarkerColor(colors2[i]);
                specReb[i]->SetMarkerStyle(markers2[i]);
                specReb[i]->Scale(0.01);

                spec[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
                if(!i) spec[i]->Draw("hist");
                else spec[i]->Draw("histsame");
                leg->AddEntry(spec[i],descCut[i].Data(),"l");


            }
            leg->Draw("same");

          cspec->SaveAs(Form("%s/SBRangesComparisonCutSys.pdf",out.Data()));
          cspec->SaveAs(Form("%s/SBRangesComparisonCutSys.png",out.Data()));

            TLegend *leg2 = new TLegend(0.55,0.55,0.89,0.85);
            leg2->SetBorderSize(0);
            TCanvas *cspec2 = new TCanvas("cspec2","cspec2",1200,800);
            TH1F *hratio[nFilesCut-1];
            for(int i=0; i<nFilesCut-1; i++){
                //hratio[i] = (TH1F*)specReb[i+1]->Clone( Form("hratio_%d",i));
                //hratio[i]->Divide(specReb[0]);
                hratio[i] = (TH1F*)spec[i+1]->Clone( Form("hratio_%d",i));
                hratio[i]->Divide(spec[0]);
                hratio[i]->SetLineStyle(linestyle2[i]);
                hratio[i]->SetLineWidth(2);
                hratio[i]->GetXaxis()->SetRangeUser(ptbinsA[0],ptbinsA[ptbinsN]);
                hratio[i]->GetYaxis()->SetRangeUser(0.5,1.5);
                hratio[i]->GetYaxis()->SetTitle(Form("ratio to central (%s)",desc[0].Data()));
                if(!i) hratio[i]->Draw("hist");
                else hratio[i]->Draw("samehist");
                leg2->AddEntry(hratio[i],descCut[i+1].Data(),"l");
            }
            leg2->Draw("same");

            TLine *line = new TLine(ptbinsA[0],1,ptbinsA[ptbinsN],1);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw("same");


            cspec2->SaveAs(Form("%s/SBRangesComparisonCutSys_ratio.pdf",out.Data(),));
            cspec2->SaveAs(Form("%s/SBRangesComparisonCutSys_ratio.png",out.Data()));

            TH1F *hsys = new TH1F("hsys","syst. rms; p_{T,ch jet};  RMS [%]",ptbinsN,ptbinsA);
            TH1F *hmean = (TH1F*)hsys->Clone("hmean");
            getRMS(nFilesCut+1,spec,hmean,hsys);

            hsys->GetYaxis()->SetRangeUser(0,15);
            hsys->SetLineColor(kViolet+2);
            TCanvas *cspecRMS = new TCanvas("cspecRMS","cspecRMS",800,400);
            hsys->Draw("hist");

            cspecRMS->SaveAs(Form("%s/SBRangesComparisonCutSys_rms.pdf",out.Data()));
            cspecRMS->SaveAs(Form("%s/SBRangesComparisonCutSys_rms.png",out.Data()));

            // /hmean->Scale(100);
            hmean->GetYaxis()->SetRangeUser(0.,30);
            hmean->SetLineColor(kMagenta+1);
            TCanvas *cspecMean = new TCanvas("cspecMean","cspecMean",800,400);
            hmean->Draw("hist");
            //line->Draw("same");

            cspecMean->SaveAs(Form("%s/SBRangesComparisonCutSys_mean.pdf",out.Data()));
            cspecMean->SaveAs(Form("%s/SBRangesComparisonCutSys_mean.png",out.Data()));

            outFile->cd();
            cspec->Write();
            cspec2->Write();
            cspecRMS->Write();
            cspecMean->Write();
            hmean->Write();
            hsys->Write();
            outFile->Close();

            return;

}



void getRMS(const int nFiles, TH1F **hratio, TH1F *hmean, TH1F *hsys)
{

  //TH1D *hsys = new TH1D("hsys","syst. rms; p_{T,ch jet};  sys [%] (rms)",nJetBins2,ptJetbins2);
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

  double *rms = new double[nJetBins2];
  double *mean = new double[nJetBins2];
  for(int i=0; i<nJetBins2; i++){
      rms[i] = 0;
       mean[i] = 0;
       for (int j=0; j<nFiles-1; j++){
         mean[i] = mean[i]+ ( hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) );
      //double m = ( hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) );
      //rms[i] = rms[i]+ ( 1-hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) ) * ( 1-hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) ) ;
      }
   mean[i] = mean[i]/(double)(nFiles-1);

   for (int j=0; j<nFiles-1; j++){
      //mean[i] = mean[i]+ ( hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) );
      //double m = ( hratios[j]->GetBinContent(hratios[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) );
    //  rms[i] = rms[i]+ ( mean[i]-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) ) * ( mean[i]-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) ) ;
      rms[i] = rms[i]+ ( 1-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) ) * ( 1-hratio[j]->GetBinContent(hratio[j]->FindBin( (ptJetbins2[i]+ptJetbins2[i+1])/2. )) ) ;

  }
      rms[i] = sqrt(rms[i]/(double)(nFiles-1));
      hsys->SetBinContent(i+1,rms[i]*100);

      hmean->SetBinContent(i+1,mean[i]);
      cout << "RMS pT " << (ptJetbins2[i]+ptJetbins2[i+1])/2. << " GeV/c:\t" << rms[i]*100 << endl;
      cout << "Mean pT " << (ptJetbins2[i]+ptJetbins2[i+1])/2. << " GeV/c:\t" << mean[i] << endl;
  }

}
