//Macro to extra jetpt spectra of different output files and compare them including ratios
//Author: A Mohanty (auro.mohanty@cern.ch)

//standard macro,not used now.. 
//this comparison: Simulations comparison
//output files are expected to lie in "results" folder for respective variation.

#include <string>
#include <sstream>
#include <iostream>

void compaYieldFDsim(TString prod="compare", TString outDir = "/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results/simCompare"){

double ratiomin=0.5, ratiomax=1.8;
double comin=1.0, comax=200000.0;
const int SE = 1; const int DC = 2; //
TString dir[DC]={"", "_effScaled"};
TString dirL[DC]={"Not scaled", "Eff scaled"};
TString histoName = "hPt";
TString simnum = "1504197460";TString outName = "pdf 21200";
//TString simnum = "1504197966";TString outName = "mb = 4.5";
//TString simnum = "1504199953";TString outName = "pdf 10800";
//TString simnum = "1504202511";TString outName = "mb = 5";
//TString simnum = "1504212024";TString outName = "F=0.5, R=0.5";
//TString simnum = "1504259653";TString outName = "F=1.0, R=0.5";
//TString simnum = "1504284947";TString outName = "F=0.5, R=1.0";
//TString simnum = "1504296768";TString outName = "F=2.0, R=1.0";
//TString simnum = "1504318569";TString outName = "F=2.0, R=2.0";
//TString simnum = "1506803374";TString outName = "F=1.0, R=2.0";
TString simnum = "150593961473";TString outName = "central";

TString title = Form("Yo");
TString rootfile;
//TH1F* hhsUnc[SE][DC];
TFile* File[SE][DC];
TH1F* hh[SE][DC]; TH1F* hhs[SE][DC]; // 
TH1F* hratio[DC]; // number of ratio histos
TH1F* hratioC[DC]; // number of ratio histos
//TH1F* hr[DC], *hh1[DC];

//Getting the histograms from the root files.
//-------------------------------------------------------------------------------
for (int j=0; j<SE;j++){        //j=0-> Sideband, j=1-> EffScaled
        for (int i=0; i<DC; i++){
                rootfile=Form("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results/DzeroR03_pPbCuts_noEff/Simulations/BFeedDown/JetPt_AnalysisResults_FastSim_powheg+pythia6_beauty_%s_Dpt3_36%s_Dzero.root",simnum.Data(),dir[i].Data());
///home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results/DzeroR03_pPbCuts/Simulations/BFeedDown/JetPt_AnalysisResults_FastSim_powheg+pythia6_beauty_1504197460_Dpt3_36_effScaled_Dzero.root
                File[j][i] = new TFile(rootfile,"read");
                hh[j][i] = (TH1F*)File[j][i]->Get(histoName.Data());
                hhs[j][i] = (TH1F*)hh[j][i]->Clone();
        }
}

//Dividing the histos from one variation by another

TLegend *le2 = new TLegend(0.65,0.5,0.9,0.8);
Color_t colors[] = {kBlue+2, kRed+2, kGreen+2, kYellow+2, kOrange+2, kBlack, kBlue+10,kBlue-2, kRed-2, kGreen-2};
TCanvas* c1 = new TCanvas("c1","c1",800,600);
TPad* pad1 = new TPad("pad1","pad1",0,0.3,1,1.0);
pad1->SetBottomMargin(0);
pad1->Draw();
pad1->SetLogy();
pad1->cd();

c1->cd();
TPad* pad2 = new TPad("pad2","pad2",0,0.03,1,0.3);
pad2->SetTopMargin(0);
pad2->SetBottomMargin(0.2);
pad2->Draw();
pad2->cd();

for (int i=0; i< DC; i++){
	hratio[i] = (TH1F*) hh[0][i]->Clone();
	hratio[i]->Divide(hh[0][0]);
//	hratio[i]->GetYaxis()->SetRangeUser(ratiomin,ratiomax);
	hratio[0]->SetTitle("");

	hratio[i]->SetMarkerColor(colors[i]);
        hratio[i]->SetMarkerStyle(20);
        hratio[i]->SetLineColor(colors[i]);

	c1->cd();

	pad1->cd();
	hhs[0][i]->GetYaxis()->SetTitle("#frac{dN}{dp_{T}}");
	hhs[0][i]->GetYaxis()->SetTitleOffset(0.7);
	hhs[0][i]->SetTitle(title.Data());

	hhs[0][i]->SetMarkerColor(colors[i]);
        hhs[0][i]->SetMarkerStyle(20);
        hhs[0][i]->SetLineColor(colors[i]);

	le2->AddEntry(hhs[0][i],Form("#it{sim%s} ",dir[i].Data()),"pe");

//	hhs[0][i]->Scale(1,"width");
//	hhs[0][i]->GetYaxis()->SetRangeUser(comin,comax);

	if (i==0) hhs[0][i]->Draw();
	else hhs[0][i]->Draw("same");

	c1->cd();

	pad2->cd();
	if (!i){
	hratio[i]->GetYaxis()->SetTitle("ratio");
	hratio[i]->GetYaxis()->SetTitleSize(20);
	hratio[i]->GetYaxis()->SetTitleFont(43);
	hratio[i]->GetYaxis()->SetLabelFont(83); // Absolute font size in pixel (precision 3)
	hratio[i]->GetYaxis()->SetLabelSize(15);
	hratio[i]->GetXaxis()->SetLabelFont(83); // Absolute font size in pixel (precision 3)
	hratio[i]->GetXaxis()->SetLabelSize(15);
	hratio[i]->GetXaxis()->SetTitle("p_{T,ch jet} (GeV/c)");
	hratio[i]->GetXaxis()->SetTitleSize(18);
	hratio[i]->GetXaxis()->SetTitleFont(43);
	hratio[i]->GetXaxis()->SetTitleOffset(2.5);
	}
	hratio[i]->Draw("same");
//	hratio[i]->GetXaxis()->SetTitle("p_{T,ch jet}");
}

c1->cd();pad1->cd();
le2->Draw("same");
c1->SaveAs(Form("%s/plots/%s.pdf",outDir.Data(),outName.Data()));
c1->SaveAs(Form("%s/plots/%s.png",outDir.Data(),outName.Data()));

TLegend *le22 = new TLegend(0.45,0.7,0.7,0.9);
TCanvas* c2 = new TCanvas("c2","c2",800,600);
for (int i=0; i<DC; i++){
	hratioC[i] = (TH1F*) hratio[i]->Clone();
	hratioC[i]->GetYaxis()->SetTitle("ratio");
	hratioC[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	if (!i){
	hratioC[i]->GetYaxis()->SetTitleSize(20);
	hratioC[i]->GetYaxis()->SetTitleFont(43);
	hratioC[i]->GetYaxis()->SetLabelFont(83); // Absolute font size in pixel (precision 3)
	hratioC[i]->GetYaxis()->SetLabelSize(15);
	hratioC[i]->GetXaxis()->SetLabelFont(83); // Absolute font size in pixel (precision 3)
	hratioC[i]->GetXaxis()->SetLabelSize(15);
	hratioC[i]->GetXaxis()->SetTitle("p_{T,ch jet} (GeV/c)");
	hratioC[i]->GetXaxis()->SetTitleSize(18);
	hratioC[i]->GetXaxis()->SetTitleFont(43);
	hratioC[i]->GetXaxis()->SetTitleOffset(2.5);
	}
	if (i){
	le22->AddEntry(hratioC[i],Form("#it{%s} ",dirL[i].Data()),"pe");
	hratioC[i]->Draw("same");
	}
}

le22->Draw("same");
c2->SaveAs(Form("%s/plots/%s_ratio.pdf",outDir.Data(),outName.Data()));
c2->SaveAs(Form("%s/plots/%s_ratio.png",outDir.Data(),outName.Data()));

}
