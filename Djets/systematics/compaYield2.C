//Macro to extra jetpt spectra of different output files and compare them including ratios
//Author: A Mohanty (auro.mohanty@cern.ch)

//standard macro,not used now.. 
//this comparison: Signal and SB ranges variation and comparison.
//output files are expected to lie in "results" folder for respective variation.

#include <string>
#include <sstream>
#include <iostream>

void compaYield2(TString prod="compare", TString outdir = "compare"){

double ratiomin=0.5, ratiomax=1.8;
double comin=1.0, comax=200000.0;
const int SE = 1; const int DC = 7; //
TString dir[DC]={"349", "359", "348", "347", "358", "249", "259"};
TString dirL[DC]={"S:3#sigma,B:4-9#sigma", "S:3#sigma,B:5-9#sigma", "S:3#sigma,B:4-8#sigma", "S:3#sigma,B:4-7#sigma", "S:3#sigma,B:5-8#sigma", "S:2#sigma,B:4-9#sigma", "S:2#sigma,B:5-9#sigma"};
TString histoName = "unfoldedSpectrum";
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
                rootfile=Form("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results/DzeroR03_pPbCuts/Default_%s/unfolding_Bayes_4/unfoldedSpectrum_unfoldedJetSpectrum.root",dir[i].Data());
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

	le2->AddEntry(hhs[0][i],Form("#it{D^{0}%s} ",dir[i].Data()),"pe");

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

TLegend *le22 = new TLegend(0.65,0.5,0.9,0.8);
TCanvas* c2 = new TCanvas("c2","c2",800,600);
for (int i=0; i<DC; i++){
	hratioC[i] = (TH1F*) hratio[i]->Clone();
	hratioC[i]->GetYaxis()->SetTitle("ratio to central (S:3#sigma,B:4-9#sigma)");
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

/*
//efficiencies. comment me
TCanvas* c3 = new TCanvas("c3","c3",800,600);
TLegend *le3 = new TLegend(0.65,0.5,0.9,0.8);
TCanvas* c4 = new TCanvas("c4","c4",800,600);
TLegend *le4 = new TLegend(0.65,0.5,0.9,0.8);
for (int i=0; i< DC; i++){
	c3->cd();

//	heff[0][i]->GetYaxis()->SetTitle("prompt eff");
	heff[0][i]->GetYaxis()->SetTitleOffset(0.7);
	heff[0][i]->SetTitle("prompt eff");

	heff[0][i]->SetMarkerColor(colors[i]);
        heff[0][i]->SetMarkerStyle(20);
        heff[0][i]->SetLineColor(colors[i]);

	le3->AddEntry(heff[0][i],Form("#it{D^{0}%s cuts} ",dir[i].Data()),"pe");

	if (i==0) heff[0][i]->Draw();
	else heff[0][i]->Draw("same");
}
le3->Draw("same");
for (int i=0; i< DC; i++){
	c4->cd();

//	heffnon[0][i]->GetYaxis()->SetTitle("non prompt eff");
	heffnon[0][i]->GetYaxis()->SetTitleOffset(0.7);
	heffnon[0][i]->SetTitle("non prompt eff");

	heffnon[0][i]->SetMarkerColor(colors[i]);
        heffnon[0][i]->SetMarkerStyle(20);
        heffnon[0][i]->SetLineColor(colors[i]);

	le4->AddEntry(heffnon[0][i],Form("#it{D^{0}%s cuts} ",dir[i].Data()),"pe");

	if (i==0) heffnon[0][i]->Draw();
	else heffnon[0][i]->Draw("same");

}

le4->Draw("same");
*/
}
