#include "style.C"

TString out = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50_ppbinning/systematics";
TString jetDataFile = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50_ppbinning/signalExtraction/JetPtSpectra_SB_eff.root";
const int ptbinsJetN = 9;
double ptbinsJet[ptbinsJetN+1] = { 3,4,5,6,8,10,14,20,30,50 };

// Pb-Pb binning
TString out = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50_PbPbbinning/systematics";
TString jetDataFile = "/home/basia/Work/alice/analysis/pPb_run2/DzeroR03_RefDPt3PythiaEff_BaseCuts/Default_jetMeas3_50_jetTrue3_50_PbPbbinning/signalExtraction/JetPtSpectra_SB_eff.root";
const int ptbinsJetN = 7;
double ptbinsJet[ptbinsJetN+1] = {3,5,10,15,20,25,35,50 };


void compareVariationAnalysis(char *outDir = "JetPtComparison" )
{

    style();

    double plotmin = 3;
    double plotmax = 50;

    TString jetSimFile = "FinalRawYieldCentralPlusSystUncertainty_Dzero.root";
    TFile *outFile = new TFile(Form("%s/YieldExtraction.root",out.Data()),"RECREATE");

	  TFile *jetPtFile2 = new TFile(jetDataFile,"read");
	  TFile *jetPtFile4 = new TFile(jetSimFile,"read");

    TH1F *hjetDataFile_tmp = (TH1F*)jetPtFile2->Get("hjetptspectrum");
    TH1F *hjetSimFile_tmp = (TH1F*)jetPtFile4->Get("JetRawYieldCentral");
    TH1F *hjetDataFile = (TH1F*)hjetDataFile_tmp->Rebin(ptbinsJetN,"hjetDataFile",ptbinsJet);
    TH1F *hjetSimFile = (TH1F*)hjetSimFile_tmp->Rebin(ptbinsJetN,"hjetSimFile",ptbinsJet);

    hjetDataFile->SetMarkerColor(2);
    hjetDataFile->SetMarkerStyle(20);;
    hjetDataFile->SetLineColor(2);
    hjetDataFile->SetLineWidth(2);
    hjetDataFile->SetMarkerSize(1.1);

    hjetSimFile->SetMarkerColor(2);
    hjetSimFile->SetMarkerStyle(24);;
    hjetSimFile->SetLineColor(2);
    hjetSimFile->SetLineWidth(2);
    hjetSimFile->SetMarkerSize(1.1);
     hjetSimFile->SetFillColor(kRed-2);
     hjetSimFile->SetFillStyle(3005);
     hjetSimFile->SetTitle();

    //hjetDataFile->SetMinimum(50);
    hjetDataFile->GetXaxis()->SetRangeUser(3,30);
    hjetDataFile->Scale(1,"width");
    hjetSimFile->Scale(1,"width");
    hjetSimFile->GetYaxis()->SetTitle("eff. corr yield");
    hjetSimFile->GetXaxis()->SetTitle("p_{T}^{ch,jet} (GeV/c)");
    hjetSimFile->GetXaxis()->SetTitleSize(0.05);
    hjetSimFile->GetXaxis()->SetRangeUser(plotmin, plotmax);
    hjetDataFile->GetXaxis()->SetTitle("p_{T}^{ch,jet} (GeV/c)");
    hjetDataFile->GetXaxis()->SetRangeUser(plotmin, plotmax);
    hjetDataFile->GetXaxis()->SetTitleSize(0.05);

    TCanvas *cJetPt = new TCanvas("cJetPt","cJetPt",800,600);
    cJetPt->SetLogy();
    hjetDataFile->Draw();
    hjetSimFile->Draw("e2 same");

    TLegend *leg = new TLegend(0.6,0.65,0.85,0.85);
    leg->AddEntry(hjetDataFile,"SB, analysis","p");
    leg->AddEntry(hjetSimFile,"SB, variation","p");
    leg->Draw("same");

    cJetPt->SaveAs(Form("%s/jetPtComparison_DataMVariation.pdf",out.Data()));
    cJetPt->SaveAs(Form("%s/jetPtComparison_DataMVariation.png",out.Data()));

    TCanvas *cjetDataFile2 = new TCanvas("cjetDataFile2","cjetDataFile2",800,600);
    cjetDataFile2->SetLogy();
    hjetDataFile->Draw();
    hjetSimFile->Draw("e2 same");

    TLegend *leg11 = new TLegend(0.6,0.65,0.85,0.85,"SB");
    leg11->AddEntry(hjetDataFile,"data","p");
    leg11->AddEntry(hjetSimFile,"yield extr.","p");
    leg11->Draw("same");

    cjetDataFile2->SaveAs(Form("%s/jetPtComparison_DataMVariation_SB.pdf",out.Data()));
    cjetDataFile2->SaveAs(Form("%s/jetPtComparison_DataMVariation_SB.png",out.Data()));

     hjetSimFile->SetTitle("Jet spectrum central values+syst. yield unc.");

    TLine *l = new TLine(4,1,40,1);
    l->SetLineStyle(2);

    TH1F *relUncDataSB = (TH1F*)hjetDataFile->Clone("relUncDataSB");
    for(int i=1; i<=hjetDataFile->GetNbinsX(); i++){
            double bin = hjetDataFile->GetBinContent(i);
            double err = hjetDataFile->GetBinError(i);
            if(err) {
                relUncDataSB->SetBinContent(i,err/bin*100);
		             cout << "bin" << i << "\tunc: " << err / bin << endl;
                relUncDataSB->SetBinError(i,0);
            }
            else {
                relUncDataSB->SetBinContent(i,0);
                 relUncDataSB->SetBinError(i,0);
            }
    }

    relUncDataSB->GetYaxis()->SetTitle("Data rel. stat. unc.");

    TCanvas *cJetDataUn = new TCanvas("cJetDataUn","cJetDataUn",800,600);
    relUncDataSB->SetName("dataStatUnc");
    relUncDataSB->Draw("h");
    //leg2->Draw("same");

     cJetDataUn->SaveAs(Form("%s/RawYieldDataStatUnc.pdf",out.Data()));
     cJetDataUn->SaveAs(Form("%s/RawYieldDataStatUnc.png",out.Data()));

	cout << "SB systematic unc: " << endl;
    TH1F *relUncVarSB = (TH1F*)hjetSimFile->Clone("relUncVarSB");
    for(int i=1; i<=hjetSimFile->GetNbinsX(); i++){
            double bin = hjetSimFile->GetBinContent(i);
            double err = hjetSimFile->GetBinError(i);
            if(err) {
                relUncVarSB->SetBinContent(i,err/bin*100);
		cout << "bin" << i << "\tunc: " << err / bin << endl;
                relUncVarSB->SetBinError(i,0);
            }
            else {
                relUncVarSB->SetBinContent(i,0);
                 relUncVarSB->SetBinError(i,0);
            }
    }

     relUncVarSB->GetYaxis()->SetTitleOffset(1.5);
     relUncVarSB->GetYaxis()->SetRangeUser(0,10);
    relUncVarSB->GetYaxis()->SetTitle("Yield ext. syst. unc.");
    relUncVarSB->SetTitle();



    TCanvas *cJetVarUn = new TCanvas("cJetVarUn","cJetVarUn",800,600);
    relUncVarSB->Draw("h");
    //relUncVarEff->Draw("hsame");

    TLegend *leg32 = new TLegend(0.15,0.7,0.45,0.9,"Yield extr. method");
    leg32->AddEntry(hjetDataFile,"SB","p");
    leg32->Draw("same");


     cJetVarUn->SaveAs(Form("%s/YieldExtSysUnc.pdf",out.Data()));
     cJetVarUn->SaveAs(Form("%s/YieldExtSysUnc.png",out.Data()));

     TCanvas *cJetVarUn_sb = new TCanvas("cJetVarUn_sb","cJetVarUn_sb",800,600);
     relUncVarSB->SetName("sysUnc");
    relUncVarSB->Draw("h");

     cJetVarUn_sb->SaveAs(Form("%s/YieldExtSysUnc_SB.pdf",out.Data()));
     cJetVarUn_sb->SaveAs(Form("%s/YieldExtSysUnc_SB.png",out.Data()));

     TH1F *hVarRatioSB = (TH1F*)hjetDataFile->Clone("hVarRatioSB");
    hVarRatioSB->Sumw2();
    hVarRatioSB->Divide(hjetSimFile);
    hVarRatioSB->SetLineColor(kMagenta+2);
    hVarRatioSB->SetMarkerColor(kMagenta+2);
    hVarRatioSB->SetMarkerStyle(24);
     hVarRatioSB->GetYaxis()->SetTitleSize(0.06);
    hVarRatioSB->GetYaxis()->SetTitleOffset(0.8);
    hVarRatioSB->GetYaxis()->SetTitle("SB data/yield. ext.");
    hVarRatioSB->GetYaxis()->SetRangeUser(0.8,1.2);

    TCanvas *cratio4 = new TCanvas("cratio4","cratio4",800,300);
    hVarRatioSB->Draw();
    l->Draw("same");

    cratio4->SaveAs(Form("%s/YieldExtSBRatio.pdf",out.Data()));
    cratio4->SaveAs(Form("%s/YieldExtSBRatio.png",out.Data()));

    outFile->cd();
    relUncVarSB->Write();
    relUncDataSB->Write();
    outFile->Close();

    return;



}
