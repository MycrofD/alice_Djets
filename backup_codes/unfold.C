//Charm-jet pT unfolding for p-Pb
//Barbara Trzeciak

#if !defined( __CINT__) || defined(__MAKECINT__)

#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TProfile.h"
#include <iostream>
#include "TPaveText.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "THn.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBinByBin.h"

#endif


//====================== global =========================
TH2D* fMatrixPP;
TH1D* fTrueSpectrum;
TH1D* fMeasSpectrum;
TH2D* fMatrixDeltaPt;
TH1D* fRawSpectrum;

	/***********************************
	############# define your bins #####################
	************************************/

	// define your bins for the x-axis, the measured axis
	const int nbinsx=9;
	double bins_meas_fine2[nbinsx+1] = {3,4,5,6,8,10,14,20,30,50};

	// define your bins for the y-axis, the true axis, shoulde be less bins. Remember to keep over and underflow (should be taken care in this code) and to include the true overflow option in the unfolding
	const int nbinsy=9;
	double bins_true_fine2[nbinsy+1] = {3,4,5,6,8,10,14,20,30,50}; //illegal pointer :371 for nbinsy=7, from 5;fRawRebin:181

	// Check with the direct_extraction macro, why the following doesn't work 	// nbinsx = (int)(sizeof(bins_meas_fine2)/sizeof(*bins_meas_fine2)-1);//------- x: measured bins 	// nbinsy = (int)(sizeof(bins_true_fine2)/sizeof(*bins_true_fine2)-1);//------- y: true bins
//	int colortable[] = {1,2,4,kGreen+2,kGray+2, kMagenta+3,kOrange-3,kViolet+7,kCyan+3,kMagenta,kGray,3,5,kGreen+1,kGray+1, kMagenta+1,kOrange-1, kRed+2, kBlue+2, kGray-2, kGray-1,kMagenta-2};
	int colortable[] = {kMagenta, kViolet, kBlue, kCyan+2, kGreen+4, kGreen+1, kYellow+1, kOrange+1, kRed, kRed+2};
	//int markertable[] = { 20,21,22,23,29,24,25,26,22,23,33,29,20,21,22,23,33,29,24,25,26,22,23,33,29,20,21};
	//int markertable[] = { 20,21,22,23,29,47,45,41,39,33};
	int markertable[] = { 22,22,22,22,22,22,22,22,22,22};

/***********************************
############# begining of the macro ##################
************************************/
void unfold
(
TString outDir = "t_plots", // output directory
bool useDeltaPt = 1,  // if to use a separate bkg. fluctuation matrix
int bkgType = 4,  // type of bkg: Djet5Excl, .... (see below)
bool fDoWeighting = 1,
bool fdivide = 1,
bool overflow = 1,  // if to use overflow in the unfolding
bool bayesUnfolding = 0, // bayes unfolding
bool svdUnfolding = 1, //  svd unfolding
const int NTrials = nbinsy+1,//10,  //number of total trials
const int regBayes = 5,  // default reg. parameter for the bayes unfolding
bool fPrior = 0,  // if to use prior different than the true spectrum from the sim
int priorType = 1,   // if fPrior == 1, choose type of the prior
TString priorName = "priorPowerLaw6", // and the prior name to be included in the output file name
bool debug = 0
)
{	printf("%d\n",nbinsx);

  gSystem->Load("$HOME/ALICE_HeavyFlavour/RooUnfold-1.1.1/libRooUnfold");
  gStyle->SetOptStat(0000); //Mean and RMS shown
  //gStyle->SetPalette(1); //for 3D histograms

	gSystem->Exec(Form("mkdir %s",outDir.Data()));

	double plotmin = bins_true_fine2[0] ;
	double plotmax = bins_true_fine2[nbinsy+1];

	TString bkgName="";
	switch(bkgType){
			case 1: bkgName = "Djet5Excl";//D-tagged jet > 5 excluded
					break;
			case 2: bkgName = "Incljet5Excl";//no D tagging needed. jets > 5 excl
					break;
			case 3: bkgName = "Djet10Excl";
					break;
			case 4: bkgName = "Incljet10Excl";
					break;
	}

	TString outName;
	outName = "PythiaRM_";    //outName += "_";
	outName += bkgName; outName += "_";
	if(bayesUnfolding) { outName += "bayes"; outName += regBayes; }//svdUnfolding
	else if(svdUnfolding) { outName += "svd"; outName += regBayes; }//svdUnfolding
	else {cout << "!!!!!!!!!! other options for the unfolding to be included !!!!" << endl; return;}
	if(fDoWeighting) outName += "_weight";
	if(fPrior) { outName += "_"; outName += priorName; }
	outName += "_"; outName += bins_meas_fine2[0];	outName += "_"; outName += bins_meas_fine2[nbinsx];	outName += "_"; outName += bins_true_fine2[0];	outName += "_"; outName += bins_true_fine2[nbinsy];//extra line of code
/***********************************
############# load raw spectrum and response matrices ##################
************************************/

// load raw jet pT spectrum (after the feed-down subtraction), the first argument is the output file name, and the second name of the histogram; take histogram without the bin width scaling, this is done at the end after the unfolding
LoadRawSpectrum("$HOME/ALICE_HeavyFlavour/work/jets1/p3/simulation/Bfeed22_var.root","hDataCorr_unsc",rebin=0);//LoadRawSpectrum("PATH_TO_RAW_SPECTRUM/JetPtSpectrum_FDsub.root","YOUR_HIST_NAME",0);
// background fluctuation matrix; the first argument is the output file name, and the second name of the histogram;
LoadBackgroundMatrix(Form("$HOME/ALICE_HeavyFlavour/work/jets1/p3/unfolding/matrices/RandCones_BkgM_%s.root",bkgName.Data()),"hBkgM");//LoadBackgroundMatrix("PATH_TO_THE_FILE/RandCones_BkgM_Djet5Excl.root","hBkgM");
// detector reponse matrix for prompt D, from Pythia; the first argument is the output file name, and the second name of the response matrix, third: D-jet spectrum at the generator level, fourth: D-jet spectrum at the reconstruction level;
//LoadDetectorMatrix("$HOME/ALICE_HeavyFlavour/work/jets1/p3/unfolding/matrices/DetMatrix_Dpt0_100.root","hPtJet2d","hPtJetGen","hPtJetRec",0);//LoadDetectorMatrix("PATH_TO_THE_FILE/DetMatrix_Dpt0_100.root","hPtJet2d","hPtJetGen","hPtJetRec",0);
LoadDetectorMatrix("$HOME/ALICE_HeavyFlavour/work/jets1/p3/unfolding/matrices/DetMatrix_Dpt3_36.root","hPtJet2d","hPtJetGen","hPtJetRec",0);//Remember to change the output directory

	if (!fRawSpectrum)  { Error("Unfold", "No raw spectrum!"); return 0;	}
	if (!fTrueSpectrum) { Error("Unfold", "No true spectrum!");	return 0; }
	if (!fMeasSpectrum) { Error("Unfold", "No reconstructed spectrum!"); return 0; }

	TH1D* hDeltaPtFlat = (TH1D*)fMatrixDeltaPt->ProjectionY(); //fMatrixDeltaPt  is an (!un)initialized 2D histogram.
	TH1D* truespec = (TH1D*)fTrueSpectrum->Clone("truespec"); //fTrueSpectrum is an (!un)initialized 1D histogram. 		//They are initialised by the LoadR/B/D....() functions above.
	truespec->Scale(1./truespec->Integral());
	fMatrixDeltaPt->Sumw2();
	//WeightMatrixY(fMatrixDeltaPt,fTrueSpectrum,0);

	TH2D* MatrixDeltaReb = Rebin2D("MatrixDeltaReb", fMatrixDeltaPt, nbinsx, bins_meas_fine2, nbinsy, bins_true_fine2,0);
	TH2D* MatrixComb = ProductMatrix(fMatrixDeltaPt,fMatrixPP); //product of two matrices

    TH2D* fMatrixProd;
	if(useDeltaPt) fMatrixProd = (TH2D*)MatrixComb->Clone("fMatrixProd");
	else fMatrixProd = (TH2D*)fMatrixPP->Clone("fMatrixProd");
	if (!fMatrixProd) { Error("Unfold", "Error getting product matrix!"); return 0;	}

	TH1D* priorhisto = (TH1D*) fTrueSpectrum->Clone("priorhisto");

        TF1* fPriorFunction ;
        if(fPrior) {
            fPriorFunction = getPriorFunction(fPrior, priorhisto,priorType);
TCanvas* cPrior = new TCanvas("cPrior0", "cPrior0", 800, 600);
cPrior->SetLogy();
TH1D* histoPrior=(TH1D*)priorhisto->Clone();
histoPrior->Draw();
fPriorFunction->Draw("same");
}

    TH1D* hNormY;
    // weighting the matrix
	if (fDoWeighting) {
        cout << "==== weighting ==== " << endl;
		hNormY=(TH1D*)fMatrixProd->ProjectionY("hNormY");
		if (fPrior){
	        cout << "=== using prior function ====" << endl;
			if (! hNormY->Divide(fPriorFunction) ) { cout << "\"divide\" failed "; return; }
        	}
		else{
	        cout << "==== dividing ==== " << endl;
			hNormY->Divide(priorhisto);
        }
		WeightMatrixY(fMatrixProd,hNormY,fdivide);
	}

    TH1D* hProjYeff=(TH1D*)fMatrixProd->ProjectionY("hProjYeff");
    TH1D* hProjXeff=(TH1D*)fMatrixProd->ProjectionX("hProjXeff");

    TH2D* Matrix = Rebin2D("Matrix", fMatrixProd, nbinsx, bins_meas_fine2, nbinsy, bins_true_fine2,0);
    TH1D* hProjXeffRebin=(TH1D*)hProjXeff->Rebin(nbinsx, "hProjXeffRebin", bins_meas_fine2);
    TH1D* hProjYeffRebin=(TH1D*)hProjYeff->Rebin(nbinsy, "hProjYeffRebin", bins_true_fine2);
    TH1D* fRawRebin=(TH1D*)fRawSpectrum->Rebin(nbinsx, "fRawRebin", bins_meas_fine2);

/**************************************************
############# unfolding settings ##################
**************************************************/
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    RooUnfoldResponse response(hProjXeffRebin,hProjYeffRebin, Matrix, "response","response");
    response.UseOverflow(overflow);

    TH1D* fUnfoldedBayes[NTrials];
    TH1D* folded[NTrials];
    TH2D* fPearsonCoeffs[NTrials];
    TH1D* hRatioSpectrum[NTrials];
    TH1D* hRatio[NTrials];

    TCanvas* cUnfolded = new TCanvas("cUnfolded","cUnfolded",800,600);
    cUnfolded->SetLogy();
    TLegend* leg =  new TLegend(0.6,0.4,0.85,0.85);
    leg->SetBorderSize(0);

//------------ do unfolding NTrials times ------------
	for(Int_t ivar=0; ivar<NTrials; ivar++){//changes
/***********************************
############# unfolding ##################
************************************/
//	if (bayesUnfolding){
//	RooUnfoldBayes unfoldB (&response, fRawRebin, ivar+1);
//		fUnfoldedBayes[ivar] = (TH1D*)unfoldB.Hreco();
//	}
//	else if (svdUnfolding){
	RooUnfoldSvd unfold (&response, fRawRebin, ivar+1);
//		fUnfoldedBayes[ivar] = (TH1D*)unfoldS.Hreco();
//	}
//	if 
		fUnfoldedBayes[ivar] = (TH1D*)unfold.Hreco();
		folded[ivar] = (TH1D*)response.ApplyToTruth(fUnfoldedBayes[ivar]);

		// ------------ Get Person coefficient ------------
		fPearsonCoeffs[ivar] = getPearsonCoeffs( unfold.Ereco(RooUnfold::kCovariance) );
		fPearsonCoeffs[ivar]->SetName(Form("PearsonCoeffs%d",ivar));

		fUnfoldedBayes[ivar]->SetLineColor(colortable[ivar]);
		fUnfoldedBayes[ivar]->SetMarkerColor(colortable[ivar]);
		fUnfoldedBayes[ivar]->SetMarkerStyle(20);
		fUnfoldedBayes[ivar]->SetLineWidth(2);
		fUnfoldedBayes[ivar]->GetXaxis()->SetTitle("p_{T,jet}^{ch} (GeV/#it{c})");
		fUnfoldedBayes[ivar]->GetYaxis()->SetTitle("dN/dp_{T}");
		if(bayesUnfolding) fUnfoldedBayes[ivar]->SetTitle("Bayes unfolding");
		else if(svdUnfolding) fUnfoldedBayes[ivar]->SetTitle("Svd unfolding");
		fUnfoldedBayes[ivar]->GetXaxis()->SetRangeUser(plotmin,plotmax);

		if(ivar == 0) { fUnfoldedBayes[ivar]->Draw(); }
		fUnfoldedBayes[ivar]->Draw("same");
		leg->AddEntry(fUnfoldedBayes[ivar],Form("Reg=%d",ivar+1),"p");
	
	}
	leg->Draw("same");

int regCount = 0;
double old = 1; double new = 0;

for (int i =0; i<NTrials-1; i++){ //for every unfolded histo to the penultimate histo
	int secCount = 0;
	for (int ll =0; ll<nbinsy; ll++){
		old=fUnfoldedBayes[i]->GetBinContent(ll);
		new=fUnfoldedBayes[i+1]->GetBinContent(ll);
		
		if(fabs((old - new) / old)>0.1) break;
		else secCount += 1;
	}
	if (secCount == nbinsy) {regCount = i;break;}
}
cout<<"-------------regCount is: "<<regCount<<endl;


	//cUnfolded->SaveAs(Form("%s/%s_unfSpectra.png",outDir.Data(),outName.Data()));
	cUnfolded->SaveAs(Form("%s/%s_unfSpectra.pdf",outDir.Data(),outName.Data()));

	TCanvas *cPearson = new TCanvas("cPearson","cPearson",1800,1800);
	cPearson->Divide(3,4);
	for(Int_t ivar=1; ivar<NTrials+1; ivar++)
	{
		cPearson->cd(ivar);
		fPearsonCoeffs[ivar-1]->SetTitle(Form("k=%d",ivar));
		fPearsonCoeffs[ivar-1]->SetMaximum(1);
		fPearsonCoeffs[ivar-1]->SetMinimum(-1);
		fPearsonCoeffs[ivar-1]->Draw("colz");
	}

	cPearson->SaveAs(Form("%s/%s_Pearson.pdf",outDir.Data(),outName.Data()));

    TCanvas *cRatio2 = new TCanvas("cRatio2","cRatio2",800,600);
    TLegend *legRatio2 =  new TLegend(0.6,0.5,0.7,0.85);//0.6,0.5,0.65,0.85
    legRatio2->SetBorderSize(0);
	TH1D* hBaseMeasure = (TH1D*)fRawRebin->Clone("hBaseSpectrum");

    for(Int_t ivar=0; ivar<NTrials; ivar++){
        hRatio[ivar] = (TH1D*) folded[ivar]->Clone(Form("hRatio%d",ivar));
        hRatio[ivar]->Divide(hBaseMeasure);
        hRatio[ivar]->GetYaxis()->SetTitle("folded/measured");
        //hRatio[ivar]->GetYaxis()->SetRangeUser(0,3.5);
        hRatio[ivar]->SetTitle();

        hRatio[ivar]->SetMarkerStyle(markertable[ivar]);
        hRatio[ivar]->SetMarkerColor(colortable[ivar]);//(colortable[ivar+1]);
        hRatio[ivar]->SetMarkerSize(1.2);
        hRatio[ivar]->SetLineColor(colortable[ivar]);//(colortable[ivar+1]);

		int counter = 0;
        for(int kk=1; kk<nbinsx; kk++){
                double val = fabs(1-hRatio[ivar]->GetBinContent(kk));
                //cout << "==val: " << val << endl;
                if(val<0.1) counter++;
        }
       //if(counter == nbinsx-2) cout << "!!!!!! bin found: " << ivar+1 << endl;
        hRatio[ivar]->GetXaxis()->SetRangeUser(plotmin,plotmax);

		if(ivar==0){ //hRatio[ivar]->GetXaxis()->SetRangeUser(plotmin,plotmax);
			hRatio[ivar]->Draw();
		}
		else hRatio[ivar]->Draw("same");
		legRatio2->AddEntry(hRatio[ivar],Form("Reg=%d",ivar+1),"p");

    }
    legRatio2->Draw("same");
    TLine *line = new TLine(4,1,40,1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");

    //cRatio2->SaveAs(Form("%s/%s_foldedRatio.png",outDir.Data(),outName.Data()));
    cRatio2->SaveAs(Form("%s/%s_foldedRatio.pdf",outDir.Data(),outName.Data()));

	TCanvas *cRatio = new TCanvas("cRatio","cRatio",800,600);
    TLegend *legRatio =  new TLegend(0.6,0.5,0.65,0.85);
    legRatio->SetBorderSize(0);
    TH1D *hBaseSpectrum = (TH1D*)fUnfoldedBayes[regBayes-1]->Clone("hBaseSpectrum");
//    int iter = 0;
    for(Int_t ivar=0; ivar<NTrials; ivar++){
        hRatioSpectrum[ivar] = (TH1D*) fUnfoldedBayes[ivar]->Clone(Form("hRatioSpectrum%d",ivar));
        hRatioSpectrum[ivar]->Divide(hBaseSpectrum);
        hRatioSpectrum[ivar]->GetYaxis()->SetTitle(Form("RegX/Reg%d",regBayes));
        hRatioSpectrum[ivar]->SetMarkerStyle(markertable[ivar]);//[iter]);
        hRatioSpectrum[ivar]->SetMarkerColor(colortable[ivar]);//[iter]);//(colortable[iter+1]);
        hRatioSpectrum[ivar]->SetMarkerSize(1.2);
        hRatioSpectrum[ivar]->SetLineColor(colortable[ivar]);//[iter]);//(colortable[iter+1]);

        if((ivar == 0) ||  (ivar == 1) || (ivar == regBayes-1)) continue;
        //if((ivar == baseReg-1)) continue;

        if(ivar==2) hRatioSpectrum[ivar]->Draw();
        else hRatioSpectrum[ivar]->Draw("same");
        legRatio->AddEntry(hRatioSpectrum[ivar],Form("Reg=%d",ivar+1),"p");

  //      iter++;

    }
	legRatio->Draw("same");
    line->Draw("same");

	//cRatio->SaveAs(Form("%s/%s_unfRatio.png",outDir.Data(),outName.Data()));
	cRatio->SaveAs(Form("%s/%s_unfRatio.pdf",outDir.Data(),outName.Data()));

    fRawRebin->SetLineColor(kBlue+1);
    fRawRebin->SetMarkerColor(kBlue+1);
    fRawSpectrum->SetLineColor(kBlue+1);
    fRawSpectrum->SetMarkerColor(kBlue+1);
    fUnfoldedBayes[regBayes-1]->SetLineColor(kRed+1);
    fUnfoldedBayes[regBayes-1]->SetMarkerColor(kRed+1);
    fUnfoldedBayes[regBayes-1]->SetMarkerStyle(21);
    folded[regBayes-1]->SetLineColor(kGreen+1);//(kRed+1);
    folded[regBayes-1]->SetMarkerColor(kGreen+1);//(kRed+1);

    fUnfoldedBayes[regBayes-1]->SetName("unfoldedSpectrum");
    folded[regBayes-1]->SetName("foldedSpectrum");
    fRawRebin->SetName("fRawRebin");

    TFile *outSpectra = new TFile(Form("%s/%s_unfoldedJetSpectrum.root",outDir.Data(),outName.Data()),"recreate");
    hProjYeff ->Write();
    hProjXeff ->Write();
    hProjYeffRebin ->Write();
    hProjXeffRebin ->Write();
    fMatrixProd->Write();
    Matrix->Write();


/* ==================
 * Scale with the bin width
 ================== */
    fRawRebin->Scale(1,"width");
    for(Int_t ivar=0; ivar<NTrials; ivar++){fUnfoldedBayes[ivar]->Scale(1,"width");}
    folded[regBayes-1]->Scale(1,"width");

    fRawRebin->SetTitle();
//    fUnfoldedBayes[regBayes-1]->SetLineColor(kRed+1);
//    fUnfoldedBayes[regBayes-1]->SetMarkerColor(kRed+1);

    fRawRebin ->Write();
    fUnfoldedBayes[regBayes-1]->Write();
    folded[regBayes-1]->Write();

    outSpectra->Close();

/* ==============================
 * Drawing, etc., below
 ================================*/

    TCanvas* cSpectra = new TCanvas("cSpectra","cSpectra",800,600);
    cSpectra->SetLogy();
    //fRawSpectrum ->Draw();
    fUnfoldedBayes[regBayes-1]->Draw();
    folded[regBayes-1]->Draw("same");
    fRawRebin ->Draw("same");

    TLegend* ls = new TLegend(0.55,0.6,0.85,0.85);
    ls->SetBorderSize(0);
    ls->AddEntry(fRawRebin,"measured","p");
    ls->AddEntry(fUnfoldedBayes[regBayes-1],"unfolded","p");
    ls->AddEntry(folded[regBayes-1],"folded","l");
    ls->Draw("same");

    cSpectra->SaveAs(Form("%s/%s_UnfSpectrum.pdf",outDir.Data(),outName.Data()));


    TH1F* hratio = (TH1F*)fUnfoldedBayes[regBayes-1]->Rebin(nbinsy,"hratio",bins_true_fine2);
	if (nbinsy != nbinsx){
	    TH1D* fRawRebinClone = new TH1D("fRawRebinClone","Measured hist, True rebinned",nbinsy,bins_true_fine2);
		int istart = 0;
		while (bins_true_fine2[0] != bins_meas_fine2[istart]){
			istart++;
		}
        for(int j=0; j<=fRawRebinClone->GetNbinsX()+1;j++){
        	double cont = fRawRebin->GetBinContent(j+istart);
        	double err = fRawRebin->GetBinError(j+istart);
        	fRawRebinClone->SetBinContent(j, cont);
        	fRawRebinClone->SetBinError(j, err);
//			printf("#_#_#_#_#_#_#_# %d-",j);
//			printf("%f\n",cont);
        }
	}
	else    TH1D* fRawRebinClone = (TH1D*)fRawRebin->Clone("fRawRebinClone");
//    TH1D* fRawRebinClone = (TH1D*)hBaseMeasure->Clone("fRawRebinClone");
	hratio->Divide(fRawRebinClone);
	hratio->SetLineColor(kMagenta+1);
	hratio->SetMarkerColor(kMagenta+1);

    TCanvas* cr= new TCanvas("cr","cr",800,600);
	hratio->Draw();hratio->SetTitle(Form("Unfolded(Reg=%d)/Measured",regBayes));
    cr->SaveAs(Form("%s/%s_UnfMeasRatio.png",outDir.Data(),outName.Data()));

    hProjXeff->GetXaxis()->SetRangeUser(plotmin,plotmax);
    hProjYeff->GetXaxis()->SetRangeUser(plotmin,plotmax);

    TCanvas* cProjMatrix= new TCanvas("cProjMatrix","cProjMatrix",800,600);
    cProjMatrix->SetLogz();
    fMatrixProd->Draw("colz");
    cProjMatrix->SaveAs(Form("%s/%s_MatrixProd.png",outDir.Data(),outName.Data()));

    TCanvas* cMatrix= new TCanvas("cMatrix","cMatrix",800,600);
    cMatrix->SetLogz();
    Matrix->Draw("colz");
    cMatrix->SaveAs(Form("%s/%s_Matrix.png",outDir.Data(),outName.Data()));

    TCanvas* cProjMReb= new TCanvas("cProjMReb","cProjMReb",800,600);
    cProjMReb->SetLogy();
    hProjYeffRebin->GetXaxis()->SetRangeUser(plotmin,plotmax);
    hProjXeffRebin->GetXaxis()->SetRangeUser(plotmin,plotmax);

    hProjYeffRebin->SetLineColor(2);
    hProjYeffRebin->Draw();
    hProjXeffRebin->SetLineColor(4);
    hProjXeffRebin->Draw("same");

    cProjMReb->SaveAs(Form("%s/%s_MatrixProdProjReb.png",outDir.Data(),outName.Data()));
}


int LoadDetectorMatrix(TString fn, TString mxname, TString tsname, TString msname, bool norm = 1, TString spostfix="") {
	TFile *f  = TFile::Open(fn);
	if (!f) { Error("LoadDetectorMatrix","Detector matrix file %s not found.",fn.Data()); return 0; }

    if(norm){
	if (mxname == "") mxname = "Detectormatrix";
        //fMatrixPP = (TH2D*)f->Get(mxname);
        TH2D *matrix = (TH2D*)f->Get(mxname);
        //if (!fMatrixPP) {
        if (!matrix) {
            Error("LoadDetectorMatrix","Detector matrix %s could not be gotten from file.",mxname.Data());
            return 0;
        }
        fMatrixPP = NormMatrixY(mxname,matrix);
    }
    else {
        fMatrixPP = (TH2D*)f->Get(mxname);
        if (!fMatrixPP) {
            Error("LoadDetectorMatrix","Detector matrix %s could not be gotten from file.",mxname.Data());
            return 0;
        }

    }

    for(int i=0; i<=fMatrixPP->GetNbinsX()+1;i++){
        for(int j=0; j<=fMatrixPP->GetNbinsX()+1;j++){

            double cont = fMatrixPP->GetBinContent(i,j);
            if(i==0 && j==0)fMatrixPP->SetBinContent(i,j,0);
            else if(i==fMatrixPP->GetNbinsX()+1 && j==fMatrixPP->GetNbinsY()+1)fMatrixPP->SetBinContent(i,j,0);
            else fMatrixPP->SetBinContent(i,j,cont);

        }
    }

	if (tsname == "") tsname = "hTrue";
	fTrueSpectrum = (TH1D*) f->Get(tsname);
	//fTrueSpectrum = (TH1D*) fMatrixPP->ProjectionY("fTrueSpectrum");
	tsname+=spostfix;
	if (!fTrueSpectrum) {
		Error("LoadDetectorMatrix","True spectrum %s could not be gotten from file.",tsname.Data());
		return 0;
	}

	if (msname == "") msname = "hReco";
	msname+=spostfix;
	fMeasSpectrum = (TH1D*) f->Get(msname);
	//fMeasSpectrum = (TH1D*) fMatrixPP->ProjectionX("fMeasSpectrum");
    if (!fMeasSpectrum) {
		Error("LoadDetectorMatrix","True spectrum %s could not be gotten from file.",msname.Data());
		return 0;
	}

	Info("LoadDetectorMatrix", "mtx=%s, ts=%s, ms=%s loaded.",
		mxname.Data(), tsname.Data(), msname.Data());


	return 1;
}

/// load backround matrix
int LoadBackgroundMatrix(TString fn, TString mxname) {
	TFile *f  = TFile::Open(fn);
	if (!f) { Error("LoadBackgroundMatrix","Background matrix file %s not found.",fn.Data()); return 0; }

	fMatrixDeltaPt = (TH2D*)f->Get(mxname);
	if (!fMatrixDeltaPt) {
		Error("LoadBackgroundMatrix","Background matrix %s could not be gotten from file.",mxname.Data());
		return 0;
	}
	Info("LoadBackgroundMatrix", "%s loaded.", mxname.Data());

     for(int i=0; i<=fMatrixDeltaPt->GetNbinsX()+1;i++){
        for(int j=0; j<=fMatrixDeltaPt->GetNbinsX()+1;j++){

            double cont = fMatrixDeltaPt->GetBinContent(i,j);
            if(i==0 && j==0)fMatrixDeltaPt->SetBinContent(i,j,0);
            else if(i==fMatrixDeltaPt->GetNbinsX()+1 && j==fMatrixDeltaPt->GetNbinsY()+1)fMatrixDeltaPt->SetBinContent(i,j,0);
            else fMatrixDeltaPt->SetBinContent(i,j,cont);

        }
    }


	return 1;

}


/// load raw spectrum (that is to be unfolded)
/// if sname is not specified, try to load according to default naming conventions
int LoadRawSpectrum(TString fn, TString sname, int rebin = 0, TString spostfix="") {
	TFile *f  = TFile::Open(fn);
	if (!f) { Error("LoadRawSpectrum","Raw spectrum file %s not found.",fn.Data());	return 0; }

	if (sname == "") sname = "hReco";
	sname += spostfix;


	TH1D *spectrum = (TH1D*)f->Get(sname);
	if (!spectrum) {
		Error("LoadRawSpectrum","Raw spectrum %s could not be gotten from file.",sname.Data());
		return 0;
	}
	Info("LoadRawSpectrum", "%s loaded.", sname.Data());

	if(rebin)  fRawSpectrum = (TH1D*)spectrum->Rebin(nbinsx,"fRawSpectrum",bins_meas_fine2);
	else fRawSpectrum = (TH1D*)spectrum->Clone("fRawSpectrum");


    fRawSpectrum->Sumw2();

	return 1;
}


/// get response matrix. If no background, use just detector matrix, otherwise get product. Normalize if requested.
TH2D * getResponseMatrix(bool useDeltaPt) {
	TH2D *mtx;
	TH2D *mtx2;
	TH2D *fMatrixPP2;
	TH2D *fMatrixDeltaPt2;

	if (!useDeltaPt) {
		if (!fMatrixPP) { Error("getResponseMatrix","No unfolding matrix present."); return 0; }
 		mtx = (TH2D*) fMatrixPP->Clone();
 		//mtx = (TH2D*) fMatrixDeltaPt->Clone();
		mtx->Sumw2();
	}
	else {

		mtx = ProductMatrix(fMatrixPP, fMatrixDeltaPt);

	}

    TFile *outFileM = new TFile("outMatrix.root","recreate");
	outFileM->cd();
	fMatrixPP->Write();
	mtx->Write();
	outFileM->Close();
	delete outFileM;

	return mtx;
}



/// get product of two matrices
TH2D * ProductMatrix(TH2D * MtxA, TH2D * MtxB) {
	// make sure the matrices exist
	if (!MtxA) {
		cerr << "Error in <AliHeavyUnfoldTools::ProductMatrix> : MtxA==0." << endl;
		return 0;
	}
	if (!MtxB) {
		cerr << "Error in <AliHeavyUnfoldTools::ProductMatrix> : MtxB==0." << endl;
		return 0;
	}

	Int_t binx_a=MtxA->GetNbinsX();
	Int_t biny_a=MtxA->GetNbinsY();
	Int_t binx_b=MtxB->GetNbinsX();
	Int_t biny_b=MtxB->GetNbinsY();

	// make sure the matrices are of the same size
	if (binx_b!=binx_a) {
		cerr << "Error in <AliHeavyUnfoldTools::ProductMatrix> : MtxA--MtxB dimension mismatch." << endl;
		return 0;
	}

	Double_t x_low = MtxA->GetXaxis()->GetBinLowEdge(1);
	Double_t x_up = MtxA->GetXaxis()->GetBinUpEdge(binx_a);
	Double_t y_low = MtxB->GetYaxis()->GetBinLowEdge(1);
	Double_t y_up = MtxB->GetYaxis()->GetBinUpEdge(biny_a);

	TH2D * MtxC=new TH2D("Matrix_prod","Product Matrix", binx_a, x_low, x_up, biny_b, y_low, y_up);
	Double_t c=0;
	for(Int_t k=0; k<=binx_a+1; k++){
		for(Int_t i=0; i<=biny_b+1; i++){
			c=0;
			for (Int_t j=0; j<=biny_a+1; j++){
				Double_t a = MtxB->GetBinContent(k,j);
				Double_t b = MtxA->GetBinContent(j,i);
				c+=a*b;
			}
			//MtxC->Fill(Float_t(k-0.5),Float_t(i-0.5),c);
			MtxC->SetBinContent(k,i,c);
		}
	}

	return MtxC;
}



/// Weight matrix along y axis by histo values
void WeightMatrixY(TH2D * Mtx, TH1D * h, bool divide) {
	if (!Mtx) {
		cerr << "Warning in <AliHeavyUnfoldTools::WeightMatrixY> : Mtx==0." << endl;
		return;
	}
	if (!h) {
		cerr << "Warning in <AliHeavyUnfoldTools::WeightMatrixY> : h==0." << endl;
		return;
	}

	for(int j=1; j<=Mtx->GetNbinsY()+1; j++) {

		double value = Mtx->GetYaxis()->GetBinCenter(j);
		double c = h->GetBinContent(h->GetXaxis()->FindBin(value));
		if (divide && c)  c = 1./c;
		//else c = 1.;
		for(int i=1; i<=Mtx->GetNbinsX()+1; i++) {
			Mtx->SetBinContent(i, j, Mtx->GetBinContent(i,j)*c);
			//Mtx->SetBinError(i, j, Mtx->GetBinError(i,j)*c);
		}
	}
}

/// rebin in 2d variable size - no such routine in Root
TH2D * Rebin2D(const char* name, TH2D *h, int nx, const double *binx, int ny, const double *biny, bool crop) {
	if (!h) {
		cerr << "Warning in <AliHeavyUnfoldTools::Rebin2D> : h==0." << endl;
		return 0;
	}

	TAxis *xaxis = h->GetXaxis();
    TAxis *yaxis = h->GetYaxis();

	TH2D * hre = new TH2D(name,name,nx,binx,ny,biny);
    hre->Sumw2();
    for (int i=1; i<=xaxis->GetNbins();i++) {
        for (int j=1; j<=yaxis->GetNbins();j++) {
            hre->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),h->GetBinContent(i,j));
        }
    }

	/*
	const double epsilon = 0.00001;
	int ixmin = crop ? 0  : -1;
	int ixmax = crop ? nx : nx+1;
	int iymin = crop ? 0  : -1;
	int iymax = crop ? ny : ny+1;

	for(int ix=ixmin; ix<ixmax; ix++) {
		double xlo = (ix!=-1) ? binx[ix]   : h->GetXaxis()->GetBinCenter(0)-epsilon;
		double xhi = (ix!=nx) ? binx[ix+1] : h->GetXaxis()->GetBinCenter(h->GetNbinsX()+1)+epsilon;
		for(int iy=iymin; iy<iymax; iy++) {
			double ylo = (iy!=-1) ? binx[iy]   : h->GetYaxis()->GetBinCenter(0)-epsilon;
			double yhi = (iy!=ny) ? binx[iy+1] : h->GetYaxis()->GetBinCenter(h->GetNbinsY()+1)+epsilon;
			int k = 0;
			double hits = 0;
			double ersq = 0;
			for(int jx=0; jx<=h->GetNbinsX()+1; jx++) {
				if (h->GetXaxis()->GetBinCenter(jx) <  xlo) continue;
				if (h->GetXaxis()->GetBinCenter(jx) >= xhi) continue;
				for(int jy=0; jy<=h->GetNbinsY()+1; jy++) {
					if (h->GetYaxis()->GetBinCenter(jy) <  ylo) continue;
					if (h->GetYaxis()->GetBinCenter(jy) >= yhi) continue;

					hits += h->GetBinContent(jx,jy);
					ersq += pow(h->GetBinError(jx,jy)/h->GetBinContent(jx,jy),2);
					k++;
				}
			}

			hre->SetBinContent(ix+1, iy+1, hits);
			hre->SetBinError(ix+1, iy+1, hits*sqrt(ersq));
			//hre->SetBinError(ix+1, iy+1, 0);
		}
	}
*/

    //for(int i=0;i<=hre->GetNbinsX();i++){
        for(int j=0;j<=hre->GetNbinsY()+1;j++){
            hre->SetBinContent(0,j,0);
            hre->SetBinError(0,j,0);

			hre->SetBinContent(hre->GetNbinsX()+1,j,0);
            hre->SetBinError(hre->GetNbinsX()+1,j,0);

    }
    //}

	return hre;
}



/// Prior for unfolding: try different prior functions wich best describe the raw specrum.
TF1* getPriorFunction(int prior, TH1D* spect, int priorType = 0) {

	if (!spect) {
		Error("getPriorFunction","Required spectrum does not exist.");
		return 0;
	}


	double fitlo = 3; // fFitPtMin;
	double fithi = 30; // fFitPtMax;

	TF1 *PriorFunction = 0;

	TString fitopt = "NR";
	//if (fLHfit) fitopt += "L";


	//PriorFunction = new TF1("PriorFunction","[0]*pow(1+(sqrt([1]*[1]+x*x)-[1])/([2]*[3]),-1*[2])",0,50);
	//PriorFunction->SetParameters(10,0.14,6.6,0.145);/
    PriorFunction = new TF1("PriorFunction","[0]* pow(x,-[1]) * exp(-[1]*[2]/x)",fitlo,fithi);
    //PriorFunction->SetParameters(10,3,3);
    PriorFunction->FixParameter(2,6);
    if(priorType == 0) PriorFunction->FixParameter(1,3);
    else if(priorType == 1) PriorFunction->FixParameter(1,6);
    //PriorFunction->FixParameter(2,2);
    spect->Fit(PriorFunction, fitopt,"",fitlo,fithi);


	//PriorFunction = new TF1("PriorFunction","[0]*x+[1]",0,50);
	//PriorFunction->SetParameters(10,20);
	//spect->Fit(PriorFunction, fitopt,"",fitlo,fithi);


	PriorFunction->SetTitle("Chosen prior function");
	return PriorFunction;
}

/// Create a new, y-normalized matrix
TH2D* NormMatrixY(const char* name, TH2D* Mtx) {
	if (!Mtx) {
		cerr << "Warning in <AliHeavyUnfoldTools::NormMatrixY> : Mtx==0." << endl;
		return 0;
	}

	TH2D * Mre = (TH2D*)Mtx->Clone(name);
	NormMatrixY(Mre);

	return Mre;
}


/// Normalize matrix along y axis projection
void NormMatrixY(TH2D * Mtx) {
	if (!Mtx) {
		cerr << "Warning in <AliHeavyUnfoldTools::NormMatrixY> : Mtx==0." << endl;
		return;
	}

	TH1D* h = Mtx->ProjectionY();
	WeightMatrixY(Mtx, h, true);
}


/// get pearson coeffs from covariance matrix
TH2D * getPearsonCoeffs(const TMatrixD &covMatrix) {

	Int_t nrows = covMatrix.GetNrows();
	Int_t ncols = covMatrix.GetNcols();

	TH2D* PearsonCoeffs = new TH2D("PearsonCoeffs","Pearson Coefficients", nrows, 0, nrows, ncols, 0, ncols);
	for(Int_t row = 0; row<nrows; row++) {
		for(Int_t col = 0; col<ncols; col++) {
			Double_t pearson = 0.;
			if(covMatrix(row,row)!=0. && covMatrix(col,col)!=0.)
				pearson = covMatrix(row,col)/TMath::Sqrt(covMatrix(row,row)*covMatrix(col,col));
			PearsonCoeffs->SetBinContent(row+1,col+1, pearson);
		}
	}

	return PearsonCoeffs;
}
