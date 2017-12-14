//Charm-jet pT combined detector and bkg fluctuations matrxi for p-Pb
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
#include "TH3F.h"
#include "THn.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h" 
#include "RooUnfoldSvd.h" 
#include "RooUnfoldBinByBin.h" 

#endif

#include "style.C"

//====================== global =========================
TH2D *fMatrixPP;
TH1D *fTrueSpectrum;
TH1D *fMeasSpectrum;
TH2D *fMatrixDeltaPt;
TH1D *fRawSpectrum;

	
	/***********************************
	############# define your bins #####################
	************************************/
	
	const Int_t nbinsx = 9;
	double bins_meas_fine2[nbinsx+1] = { 3,4,5,6,8,10,14,20,30,50 };
	const Int_t nbinsy = 9;
	double bins_true_fine2[nbinsy+1] = { 3,4,5,6,8,10,14,20,30,50 };
	

void combineRM (
int bkgType = 1,  // bkg type: Djet5Excl, ... (see below)
bool isFD = 0,	// if prompt of non-prompt detector response matrix taken
bool useDeltaPt = 1, 
bool fDoWeighting = 1, 
bool fdivide = 1, 
bool fPrior = 0, int priorType = 1,  TString priorName = "priorPowerLaw6",
bool debug = 0)
{
  
    gSystem->Load("PATH_TO_YOUR_ROOUNFOLD_PACKAGE/libRooUnfold.so");
    gStyle->SetOptStat(0000); //Mean and RMS shown
    style();
	
	TString outDir;
	if(isFD) outDir = "combinedMatrix_FD";
	else outDir = "combinedMatrix";
	gSystem->Exec(Form("mkdir %s",outDir.Data()));
	
	TString bkgName="";
	switch(bkgType){
			case 1: bkgName = "Djet5Excl";
					break;
			case 2: bkgName = "Incljet5Excl";
					break;
			case 3: bkgName = "Djet10Excl";
					break;
			case 4: bkgName = "Incljet10Excl";
					break;
	}
 
	
    TString outName = "";
	if(isFD) outName += "FD_";
    outName += "PythiaRM_";
    outName += bkgName;
    if(fDoWeighting) outName += "_weight"; 
    if(fPrior) { outName += "_"; outName += priorName; }
	

LoadBackgroundMatrix(Form("bkgRM/RandCones_BkgM_%s.root",bkgName.Data()),"hBkgM");

if(isFD) LoadDetectorMatrix("outResDetMatrixPythia_FD/DetMatrix_Dpt0_100.root","hPtJet2d","hPtJetGen","hPtJetRec",0);
else LoadDetectorMatrix("outResDetMatrixPythia/DetMatrix_Dpt0_100.root","hPtJet2d","hPtJetGen","hPtJetRec",0);

	if (!fTrueSpectrum) { Error("Unfold", "No true spectrum!");	return 0; }
	if (!fMeasSpectrum) { Error("Unfold", "No reconstructed spectrum!"); return 0; }

    
	TH1D *hDeltaPtFlat = (TH1D*)fMatrixDeltaPt->ProjectionY();
	TH1D *truespec = (TH1D*)fTrueSpectrum->Clone("truespec");
	truespec->Scale(1./truespec->Integral());
	
	TCanvas *ctrue = new TCanvas;
	fTrueSpectrum->Draw();
	
	fMatrixDeltaPt->Sumw2();
	//WeightMatrixY(fMatrixDeltaPt,fTrueSpectrum,0);
	
	fMatrixDeltaPt->GetXaxis()->SetTitle("p_{T,ch jet}^{rec.} (GeV/#it{c})");
    fMatrixDeltaPt->GetYaxis()->SetTitle("p_{T,ch jet}^{gen.} (GeV/#it{c})");
	fMatrixDeltaPt->SetTitle("Bkg.Fluc. Matrix");
	
	fMatrixPP->GetXaxis()->SetTitle("p_{T,ch jet}^{rec.} (GeV/#it{c})");
    fMatrixPP->GetYaxis()->SetTitle("p_{T,ch jet}^{gen.} (GeV/#it{c})");
	fMatrixPP->SetTitle("Det.Res. Matrix");
	
	TCanvas *cMatrix = new TCanvas("cMatrix","cMatrix",1200,800);
	cMatrix->Divide(2,1);
	cMatrix->cd(1);
	gPad->SetLogz();
	fMatrixDeltaPt->Draw("colz");
	cMatrix->cd(2);
	gPad->SetLogz();
	fMatrixPP->Draw("colz");
	
	
	TCanvas *cMatrixProj = new TCanvas;
	
	TH1D *hDeltaProj = (TH1D*)fMatrixDeltaPt->ProjectionY();
	TCanvas *cMProj = new TCanvas;
	hDeltaProj->Draw();
	
	TH2D *MatrixComb = ProductMatrix(fMatrixDeltaPt,fMatrixPP);
	TH1D *hProj = (TH1D*)MatrixComb->ProjectionY();
	hProj->Draw();
	
	
	TCanvas *cMatrix2 = new TCanvas;
	cMatrix2->SetLogz();
    TH2D *fMatrixProd;
	
	if(useDeltaPt) fMatrixProd = (TH2D*)MatrixComb->Clone("fMatrixProd");
	else fMatrixProd = (TH2D*)fMatrixPP->Clone("fMatrixProd");
		
	if (!fMatrixProd) { Error("Unfold", "Error getting product matrix!"); return 0;	}

    
	TH1D* priorhisto = (TH1D*) fTrueSpectrum->Clone("priorhisto");
    TF1 *fPriorFunction ;
    if(fPrior) { 
		TCanvas *cfit = new TCanvas;
        fPriorFunction = getPriorFunction(fPrior, priorhisto,priorType);
        priorhisto->Draw();
        fPriorFunction->Draw("same");
    }
    
    TH1D* hNormY;
    // weighting the matrix
	if (fDoWeighting) {
        cout << "==== weighting ==== " << endl;
		hNormY=(TH1D*)fMatrixProd->ProjectionY("hNormY");
	
		if (fPrior){
            cout << "=== using prior function ====" << endl;
			if (! hNormY->Divide(fPriorFunction) ) { cout << "divide faild "; return; } 
        }
		else{
            cout << "==== dividing ==== " << endl;
			hNormY->Divide(priorhisto);

        }
            
		WeightMatrixY(fMatrixProd,hNormY,fdivide);
		
	}
	
	fMatrixProd->GetXaxis()->SetTitle("p_{T,ch jet}^{rec.} (GeV/#it{c})");
    fMatrixProd->GetYaxis()->SetTitle("p_{T,ch jet}^{gen.} (GeV/#it{c})");
	fMatrixProd->SetTitle("Combined Matrix");
	
	TCanvas *cMatrixProd = new TCanvas();
	cMatrixProd->SetLogz();
	fMatrixProd->Draw("colz");

    TH1D* hProjYeff=(TH1D*)fMatrixProd->ProjectionY("hProjYeff");
	TH1D* hProjXeff=(TH1D*)fMatrixProd->ProjectionX("hProjXeff");
   
	TH2D *Matrix = Rebin2D("Matrix", fMatrixProd, nbinsx, bins_meas_fine2, nbinsy, bins_true_fine2,0);
	Matrix->GetXaxis()->SetTitle("p_{T,ch jet}^{rec.} (GeV/#it{c})");
    Matrix->GetYaxis()->SetTitle("p_{T,ch jet}^{gen.} (GeV/#it{c})");
	Matrix->SetTitle("Combined Matrix");
	TCanvas *cMatrixProdReb = new TCanvas();
	cMatrixProdReb->SetLogz();
	Matrix->Draw("colz");
  
    TFile *outMatrix = new TFile(Form("%s/combineMatrix%s.root",outDir.Data(),outName.Data()),"recreate");
	fMatrixProd->Write();
    Matrix->Write();
    outMatrix->Close();

	cMatrix->SaveAs(Form("%s/Matrices%s.png",outDir.Data(),outName.Data()));
	cMatrixProd->SaveAs(Form("%s/ProdMatrix%s.png",outDir.Data(),outName.Data()));
	cMatrixProdReb->SaveAs(Form("%s/ProdMatrixRebin%s.png",outDir.Data(),outName.Data()));

	return;

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
    
    
    /*
    TH1D *proj[50];
  fMatrixPP = new TH2D("fMatrixPP","fMatrixPP",50,0,50,50,0,50);
   for(int i=0; i<=50;i++){
       proj[i] = (TH1D*)matrix->ProjectionX(Form("proj%d",i), i, i+1 );
       proj[i]->SetName(Form("proj%d",i));
       proj[i]->Scale(1./proj[i]->Integral());
      
       for(int j=1;j<=50;j++){
            fMatrixPP->SetBinContent(j,i+1,proj[i]->GetBinContent(j));
           
        }
    }
    */
    
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
	//fMatrixPP2->Write();
	//fMatrixDeltaPt->Write();
	//fMatrixDeltaPt2->Write();
	//mtx2->Write();
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
/// these shapes are from Gyulnara
TF1 * getPriorFunction(int prior, TH1D * spect, int priorType = 0) {

	if (!spect) {
		Error("getPriorFunction","Required spectrum does not exist.");	
		return 0;
	}


	double fitlo = 3; // fFitPtMin;
	double fithi = 50; // fFitPtMax;

	TF1 *PriorFunction = 0;

	TString fitopt = "NR";
	//if (fLHfit) fitopt += "L";


   // PriorFunction = new TF1("PriorFunction","[0]*pow(1+(sqrt([1]*[1]+x*x)-[1])/([2]*[3]),-1*[2])",0,50);
   //PriorFunction->SetParameters(10,0.14,6.6,0.145);/
    PriorFunction = new TF1("PriorFunction","[0]* pow(x,-[1]) * exp(-[1]*[2]/x)",3,30);
    //PriorFunction->SetParameters(10,3,3);
    PriorFunction->FixParameter(2,6);
    if(priorType == 0) PriorFunction->FixParameter(1,3);
    else if(priorType == 1) PriorFunction->FixParameter(1,6);
    //PriorFunction->FixParameter(2,2);
    spect->Fit(PriorFunction, fitopt,"",fitlo,fithi);
   
    
   // PriorFunction = new TF1("PriorFunction","[0]*x+[1]",0,50);
   // PriorFunction->SetParameters(10,20);
    //spect->Fit(PriorFunction, fitopt,"",fitlo,fithi);

	
	PriorFunction->SetTitle("Chosen prior function");
	return PriorFunction;
}

/// Create a new, y-normalized matrix
TH2D * NormMatrixY(const char* name, TH2D * Mtx) {
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
	
	TH1D * h = Mtx->ProjectionY();
	WeightMatrixY(Mtx, h, true);
}


/// get pearson coeffs from covariance matrix
TH2D * getPearsonCoeffs(const TMatrixD &covMatrix) {

	Int_t nrows = covMatrix.GetNrows();
	Int_t ncols = covMatrix.GetNcols();

	TH2D * PearsonCoeffs = new TH2D("PearsonCoeffs","Pearson Coefficients", nrows, 0, nrows, ncols, 0, ncols);
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
