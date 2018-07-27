
#include <string>
#include <sstream>
#include <iostream>


Int_t colors2[] = {kRed+1,kBlue+1,kRed-9,kBlue-9,kGreen+3,kGreen-9,kMagenta+2,4,6,kCyan+1,8,kOrange-1,kGray+1,kViolet+5,kYellow+2};
Int_t markers2[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};
Int_t linestyle2[] = {1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};


const int nFiles = 2;
TString inDirData[nFiles] = {
  "pp5TeV/D0jet/results/DzeroR03_pPbCuts/Default/unfolding_Bayes_4/finalSpectra/",
  //"pp5TeV/D0jet/results/DzeroR03_cuts2/Default/unfolding_Bayes_4/finalSpectra/",
  "pPb_run2/D0jet/PreliminaryOut/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50_ppbinning_MAIN/unfolding_Bayes_3/finalSpectra/"
}

TString dataFile[nFiles] = {
  "JetPtSpectrum_final.root",
  "JetPtSpectrum_final.root"
}

TString histName[nFiles] = {
  "hData_binned",
  "hData_binned"
}


TString desc[nFiles] = {
  "D^{0}-jet, pp 5@TeV (#times A)",
  "D^{0}-jet, p-Pb@5TeV"
};

const int     ptbinsN = 7;
double        ptbinsA[ptbinsN+1] = { 5,6,8,10,14,20,30,50 };
double        plotmin = 5, plotmax = 50;

double        sysUnc_pPb[ptbinsN] = {1.949626, 1.018353, 0.4938436, 0.1498547, 0.03191447, 0.005496932, 0.0007480437};
double        sysUncErr_pPb[ptbinsN] = {0.3571682, 0.1144968, 0.04982845, 0.01508283, 0.003925845, 0.001088327, 0.0001712551};
double        sysUnc_pp[ptbinsN];
double        sysUncErr_pp[ptbinsN];
double        sysUnc[ptbinsN];
double        sysUncErr[ptbinsN];

void RpPb(TString inDirBase = "$HOME/Work/alice/analysis/", TString outName = "$HOME/Work/alice/analysis/RpPb")
{

  gSystem->Exec(Form("mkdir %s",outName.Data()));
  plotmin = ptbinsA[0], plotmax = ptbinsA[ptbinsN];

  compareData(outName,inDirBase,"RpPb_pPbcuts");

return;

}

void compareData(TString inName, TString inDirBase, TString outHistName)
{

            TString out = inName; //= inDirBase;
            gStyle->SetOptStat(0000); //Mean and RMS shown
            gSystem->Exec(Form("mkdir %s",out.Data()));

            TString dirName[nFiles];
            for (int i=0; i<nFiles; i++){
                dirName[i] = inDirBase;
                dirName[i] += "/";
                dirName[i] += inDirData[i];
            }

            TFile *fproj[nFiles];
            for(int i=0; i<nFiles; i++) {
              fproj[i] = new TFile(Form("%s/%s",dirName[i].Data(),dataFile[i].Data()),"READ");
              if(!fproj[i]) { cout << "NO FILE !!!!!" << endl; return; }
            }

            TCanvas *cspec = new TCanvas("cspec","cspec",900,800);
            cspec->SetLogy();

            TLegend *leg = new TLegend(0.5,0.55,0.75,0.65);
            leg->SetBorderSize(0);

            TH1F *spec[nFiles];
            TH1F *specReb[nFiles];

            TH1F *specSim[nFiles];
            TH1F *specSimReb[nFiles];

            for(int i=0; i<nFiles; i++) {
                spec[i] = (TH1F*)fproj[i]->Get(Form("%s",histName[i].Data()));
                spec[i]->Sumw2();
                spec[i]->SetTitle();
                spec[i]->SetLineColor(colors2[i]);
                spec[i]->SetMarkerColor(colors2[i]);
                spec[i]->SetMarkerStyle(markers2[i]);
                spec[i]->SetMinimum(0.001);

                specReb[i] = (TH1F*)spec[i]->Rebin(ptbinsN,Form("specReb_%d",i),ptbinsA);
                specReb[i]->SetTitle();
                specReb[i]->SetLineColor(colors2[i]);
                specReb[i]->SetMarkerColor(colors2[i]);
                specReb[i]->SetMarkerStyle(markers2[i]);
                if(!i) specReb[i]->Scale(208);
                specReb[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
                Double_t ptval[ptbinsN];
                Double_t ptvalunc[ptbinsN];
                for(int j=0; j<specReb[i]->GetNbinsX(); j++){
                  ptval[j] = (ptbinsA[j]+ptbinsA[j+1]) / 2.;
                  ptvalunc[j] = (ptbinsA[j+1]-ptbinsA[j]) / 2.;

                  double error = specReb[i]->GetBinError(j+1);
                  double value = specReb[i]->GetBinContent(j+1);
                  double relError = error/value;
                  //cout << "rel error " << relError << endl;
                  double totalError = TMath::Sqrt(relError*relError);
                  specReb[i]->SetBinError(j+1,totalError*value);
                  if(!i) {
                    sysUnc_pp[j] = specReb[i]->GetBinContent(j+1);
                    sysUncErr_pp[j] = specReb[i]->GetBinContent(j+1)*0.15;
                  }

                }
                leg->AddEntry(specReb[i],desc[i].Data());
            }

            TGraphAsymmErrors *graSys_pp  = new TGraphAsymmErrors(ptbinsN,ptval,sysUnc_pp,ptvalunc,ptvalunc,sysUncErr_pp,sysUncErr_pp);
            graSys_pp->SetFillColor(colors2[2]);
            graSys_pp->SetLineColor(colors2[2]);

            TGraphAsymmErrors *graSys_pPb  = new TGraphAsymmErrors(ptbinsN,ptval,sysUnc_pPb,ptvalunc,ptvalunc,sysUncErr_pPb,sysUncErr_pPb);
            graSys_pPb->SetFillColor(colors2[3]);
            graSys_pPb->SetLineColor(colors2[3]);


            TH1F *Graph_central_syst_unc1 = new TH1F("Graph_central_syst_unc1","",100,2.5,52.5);
            Graph_central_syst_unc1->SetMinimum(2.e-04);
            Graph_central_syst_unc1->SetMaximum(10);
            Graph_central_syst_unc1->SetDirectory(0);
            Graph_central_syst_unc1->SetStats(0);
            Graph_central_syst_unc1->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
            Graph_central_syst_unc1->GetXaxis()->SetLabelFont(42);
            Graph_central_syst_unc1->GetXaxis()->SetLabelSize(0.035);
            Graph_central_syst_unc1->GetXaxis()->SetTitleSize(0.035);
            Graph_central_syst_unc1->GetXaxis()->SetTitleFont(42);
            Graph_central_syst_unc1->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#it{#eta}} [mb (GeV/#it{c})^{-1}]");
            Graph_central_syst_unc1->GetYaxis()->SetLabelFont(42);
            Graph_central_syst_unc1->GetYaxis()->SetLabelSize(0.035);
            Graph_central_syst_unc1->GetYaxis()->SetTitleSize(0.035);
            Graph_central_syst_unc1->GetYaxis()->SetTitleFont(42);
            Graph_central_syst_unc1->GetZaxis()->SetLabelFont(42);
            Graph_central_syst_unc1->GetZaxis()->SetLabelSize(0.035);
            Graph_central_syst_unc1->GetZaxis()->SetTitleSize(0.035);
            Graph_central_syst_unc1->GetZaxis()->SetTitleFont(42);
            graSys_pPb->SetHistogram(Graph_central_syst_unc1);

            graSys_pPb->Draw("2A");
            graSys_pp->Draw("2");
            specReb[1]->Draw("same");
            specReb[0]->Draw("same");
            leg->Draw("same");

            TPaveText *pt = new TPaveText(0.4,0.7,0.85,0.9,"NB NDC");
            pt->SetBorderSize(0);
            pt->SetFillStyle(0);
            pt->SetTextAlign(13);
            pt->SetTextFont(43);
            pt->SetTextSize(22);
            TText *text = new TText;
            text = pt->AddText("pp, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
            text = pt->AddText(Form("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.%d, |#it{#eta}_{lab}^{jet}| < 0.%d",3,6));
            text = pt->AddText(Form ("with D^{0}, %d < #it{p}_{T,D} < %d GeV/#it{c}",3,36));
            pt->Draw();

          cspec->SaveAs(Form("%s/%s.pdf",out.Data(),outHistName.Data()));
          cspec->SaveAs(Form("%s/%s.png",out.Data(),outHistName.Data()));


            TCanvas *cspec2 = new TCanvas("cspec2","cspec2",1000,600);

                TH1F *hratio = (TH1F*)specReb[1]->Clone( Form("hratio_%d",i));
                hratio->Divide(specReb[0]);
                hratio->SetLineStyle(linestyle2[1]);
                hratio->SetLineColor(colors2[4]);
                hratio->SetMarkerColor(colors2[4]);
                hratio->SetLineWidth(2);
                hratio->GetYaxis()->SetRangeUser(0.,2);
                hratio->GetYaxis()->SetTitle("R_{pPb}");

                for(int j=0; j<hratio->GetNbinsX(); j++){
                    double relUncpp = sysUncErr_pp[j] / sysUnc_pp[j];
                    double relUncpPb = sysUncErr_pPb[j] / sysUnc_pPb[j];

                    sysUnc[j] = hratio->GetBinContent(j+1);
                    sysUncErr[j] = sqrt(relUncpp*relUncpp+relUncpPb*relUncpPb) * sysUnc[j];

                }

                TGraphAsymmErrors *graSys  = new TGraphAsymmErrors(ptbinsN,ptval,sysUnc,ptvalunc,ptvalunc,sysUncErr,sysUncErr);
                graSys->SetFillColor(colors2[5]);
                graSys->SetLineColor(colors2[5]);
                TH1F *Graph_central_syst_unc = new TH1F("Graph_central_syst_unc","",100,2.5,52.5);
                Graph_central_syst_unc->SetMinimum(0);
                Graph_central_syst_unc->SetMaximum(2.6);
                Graph_central_syst_unc->SetDirectory(0);
                Graph_central_syst_unc->SetStats(0);
                Graph_central_syst_unc->GetXaxis()->SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})");
                Graph_central_syst_unc->GetXaxis()->SetLabelFont(42);
                Graph_central_syst_unc->GetXaxis()->SetLabelSize(0.035);
                Graph_central_syst_unc->GetXaxis()->SetTitleSize(0.035);
                Graph_central_syst_unc->GetXaxis()->SetTitleFont(42);
                Graph_central_syst_unc->GetYaxis()->SetTitle("R_{pPb}");
                Graph_central_syst_unc->GetYaxis()->SetLabelFont(42);
                Graph_central_syst_unc->GetYaxis()->SetLabelSize(0.035);
                Graph_central_syst_unc->GetYaxis()->SetTitleSize(0.035);
                Graph_central_syst_unc->GetYaxis()->SetTitleFont(42);
                Graph_central_syst_unc->GetZaxis()->SetLabelFont(42);
                Graph_central_syst_unc->GetZaxis()->SetLabelSize(0.035);
                Graph_central_syst_unc->GetZaxis()->SetTitleSize(0.035);
                Graph_central_syst_unc->GetZaxis()->SetTitleFont(42);
                graSys->SetHistogram(Graph_central_syst_unc);

                graSys->Draw("2A");
                hratio->Draw("same");
                pt->Draw();

            TLine *line = new TLine(ptbinsA[0],1,ptbinsA[ptbinsN],1);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw("same");


            cspec2->SaveAs(Form("%s/%s_ratio.pdf",out.Data(),outHistName.Data()));
            cspec2->SaveAs(Form("%s/%s_ratio.png",out.Data(),outHistName.Data()));

            return;

}
