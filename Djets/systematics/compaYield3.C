//Macro to extra jetpt spectra of different output files and compare them icnluding ratios
//Author: A Mohanty (auro.mohanty@cern.ch)

//standard macro,not used now.. 
//this comparison: Signal and SB ranges variation and comparison.
//output files are expected to lie in "results" folder for respective variation.

#include <string>
#include <sstream>
#include <iostream>

void compaYield3(TString prod="compare", TString outdir = "compare"){

int colors[] = {1,2,kGreen+3,kMagenta+2,4,6,kCyan+1,kYellow+2,kOrange-1,8,kGray+1,kViolet+5};



double ratiomin=0.5, ratiomax=1.8;
double comin=1.0, comax=200000.0;
const int DC = 2;
int Dbins = 10;
double DbinEdges[] = {3,4,5,6,7,8,10,12,16,24,36};
TH1D *h = new TH1D("","",Dbins,DbinEdges);
TH1D *hh[DC];

double mc[]  = {0.009975,0.01086,0.01171,0.01262,0.01299,0.01371,0.01461,0.01612,0.01745,0.0185367 }; // set up sigma of the D signal from MC
double mcE[] = {0.000077,0.00009,0.00010,0.00014,0.00018,0.00016,0.00024,0.00025,0.00044,0.0013};
//double da[]  = {4,5,6,7,8,10,12,16,24,36,3};
//double daE[] = {4,5,6,7,8,10,12,16,24,36,3};

TFile *fData = new TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results/DzeroR03_pPbCuts/Default_249/signalExtraction/JetPtSpectra_SB_eff.root","READ");
TH1D *hsigma = (TH1D*)fData->Get("hsigma");
hh[1] = (TH1D*)hsigma->Clone();
hh[1]->GetYaxis()->SetRangeUser(9,28);
//-------------------------------------------------------------------------------

TCanvas *c2 = new TCanvas("c2","c2",1200,800);
TLegend *leg2 = new TLegend(0.1,0.8,0.4,0.9);
/*
for (int j=0; j<DC; j++){
	hh[j]= (TH1D*)h->Clone();
	for (int i=0; i<Dbins; i++){
		double value;double err;
		if(!j) {value = mc[i]; err = mcE[i];}
		else {value = da[i]; err = daE[i];}
cout<<value<<"======="<<endl;
//return;
		hh[j]->SetBinContent(i+1,value);
        	hh[j]->SetBinError(i+1,err);
		hh[j]->SetLineColor(colors[j]);
		hh[j]->SetMarkerColor(colors[j]);
	}
	if(!j)hh[j]->Draw();
	else hh[j]->Draw("same");
}
*/
	hh[0]= (TH1D*)h->Clone();
	for (int i=0; i<Dbins; i++){
		double value;double err;
		value = mc[i]*1000; err = mcE[i]*1000;
		hh[0]->SetBinContent(i+1,value);
        	hh[0]->SetBinError(i+1,err);
		hh[0]->SetLineColor(colors[3]);
		hh[0]->SetMarkerColor(colors[3]);
                hh[0]->SetMarkerStyle(20);
		hh[0]->SetLineWidth(2);

	}
leg2->AddEntry(hh[0],"MC");
leg2->AddEntry(hh[1],"Data");

hh[1]->Draw();
hh[0]->Draw("same");
leg2->Draw("same");

return;
}
