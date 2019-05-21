//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//  [Modified] A.Mohanty
//  Utrecht University
//  auro.mohanty@cern.ch
//-----------------------------------------------------------------------
//https://arxiv.org/pdf/1105.1160.pdf
//http://hep.physics.utoronto.ca/~orr/wwwroot/Unfolding/d-agostini.pdf

#include "../DsignalExtraction/configDzero_ppz.h"

//====================== global =========================
double plotmin = fptbinsZMeasA[0], plotmax =  fptbinsZMeasA[fptbinsZMeasN];

const int fptbinsZGenN=10, fJetptbinsGenN=9;
double fptbinsZGenA[fptbinsZGenN+1]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.02};
double fJetptbinsGenA[fJetptbinsGenN+1]={0,2,3,4,5,7,10,15,50,60};
const int nDim    = 4;//the four dimensions
double bins[nDim] = {fptbinsZGenN,fptbinsZGenN};//for creating 4D ResponseMatrix
int dim[nDim]     = {0,1,5,6};//for extacting 4D info from THnSparse
TH1D *fRawSpec2Dproj[fJetptbinsN];

TH1D *LoadRawSpec(TString fn, TString sname, TString spostfix="");
TH2D *Rebin2D(const char* name, TH2D *h, int nx, const double *binx, int ny, const double *biny, bool crop);
////===================== for 2D unfolding
/***********************************
############# define your bins #####################
************************************/
int colortable[] = {kMagenta, kViolet, kBlue, kCyan+2, kGreen+4, kGreen+1, kYellow+1, kOrange+1, kRed, kRed+2};
int linesytle[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
/***********************************
############# begining of the macro ##################
************************************/
void unfold_Bayeszjet(
  //TString roounfoldpwd = "",
  TString listName = "FD",
  bool isPrompt=1,
  bool postfix=0,
  bool isprefix=0,
  TString effFile = "",
  TString datafile = "file.root",
  TString detRMfile = "detRM.root",
  TString bkgRMfile = "bkgRM.root",
  TString outDir = "out", // output directory
  const int regBayes = 5,  // default reg. parameter for the bayes unfolding
  bool isPrior = 0,  // if to use prior different than the true spectrum from the sim
  int priorType = 1,   // if isPrior == 1, choose type of the prior
  bool useDeltaPt = 1,  // if to use a separate bkg. fluctuation matrix
  bool isFDUpSpec = 0,
  bool isFDDownSpec = 0,
  bool fDoWeighting = 1,
  bool fdivide = 1,
  bool overflow = 1,  // if to use overflow in the unfolding
  const int NTrials = 10,//10,  //number of total trials
  bool debug = 0
)
{
    bool tempo_data_plot = 0;
    //-----------------------------
    gStyle ->SetOptStat(0000);
    gSystem->Exec(Form("mkdir  %s",outDir.Data()));
    gSystem->Exec(Form("mkdir  %s/Bayes",outDir.Data()));
    outDir+="/Bayes";
    gSystem->Exec(Form("mkdir  %s",outDir.Data()));
    gSystem->Exec(Form("mkdir  %s/plots",outDir.Data()));
    gSystem->Exec(Form("mkdir  %s/alljetz2D",outDir.Data()));
    //-----------------------------
    // Extracting data 2D histogram in z and jetpt. Before it was just in z...  for unfolding.
    TH2* hData2D = new TH2D("hFD_zjet", "hFD_zjet", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    TH2* hData2DS = new TH2D("hFD_zjet", "hFD_zjet", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    TString data2D[fJetptbinsN];
    for (int binjet=0; binjet < fJetptbinsN; binjet++){
        //reading data from each 1D
        data2D[binjet]=Form("%s",datafile.Data());
        data2D[binjet]+=Form("/Jetbin_%d_%d",(int)fJetptbinsA[binjet], (int)fJetptbinsA[binjet+1]);
        data2D[binjet]+="/plots/JetPtSpectrum_FDsub.root";
        //reading 2D DetRM: need to create 2D RM file in ResponseMatrix folder
    }
    for (int binjet=0; binjet < fJetptbinsN; binjet++){
        //cout<<fJetptbinsA[binjet]<<endl;
        fRawSpec2Dproj[binjet] = (TH1D*)LoadRawSpec(data2D[binjet].Data(),"hData_binned_sub");
        TH1D *fRawSpec2DprojScale=(TH1D*)fRawSpec2Dproj[binjet]->Clone();
        fRawSpec2DprojScale->Scale(1,"width");
        for (Int_t binz=0; binz < fptbinsZMeasN+1; binz++){
            double cont  = fRawSpec2Dproj[binjet]->GetBinContent(binz);
            double contS = fRawSpec2DprojScale->GetBinContent(binz);
        	  double contErr = fRawSpec2Dproj[binjet]->GetBinError(binz);
        	  double contErrS = fRawSpec2DprojScale->GetBinError(binz);
        	  hData2D->SetBinContent(binz,binjet+1,cont);
        	  hData2DS->SetBinContent(binz,binjet+1,contS);
        	  hData2D->SetBinError(binz,binjet+1,contErr);
        	  hData2DS->SetBinError(binz,binjet+1,contErrS);
        }
    }
    // Saving the 2D data file.
    TCanvas *cFD_2D = new TCanvas("Measured2D", "Measured2D", 800, 600);
    cFD_2D  ->SetLogz();
    hData2DS->SetTitle("z-jet p_{T} spectrum (before unfolding)");
    hData2DS->GetYaxis()->SetTitle("jet p_{T}");
    hData2DS->GetXaxis()->SetTitle("z_{||}^{ch}");
    hData2DS->Draw("colz");
    hData2DS->Draw("TEXT SAME");
    cFD_2D  ->SaveAs(Form("%s/alljetz2D/FDdata2D.pdf",outDir.Data()));
    cFD_2D  ->SaveAs(Form("%s/alljetz2D/FDdata2D.png",outDir.Data()));
    cFD_2D  ->SaveAs(Form("%s/alljetz2D/FDdata2D.svg",outDir.Data()));
    TFile *outFile = new TFile(Form("%s/alljetz2D/outFD.root",outDir.Data()),"recreate");
    hData2D ->Write();
    hData2DS->Write();
    outFile ->Close();
    //-----------------------------
    // Time for Response Matrix in 2D
    TFile          *File = new TFile(effFile,"read");
    TDirectoryFile *dir  = (TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
    TString histName;
    if(!isprefix){
        if(fDmesonSpecie) histName = "histosDStarMBN";
        else histName              = "histosD0MBN";}
    else{
        if(fDmesonSpecie) histName = "histosDStarMBN";
        else histName              = "histosD0";}
    TPaveText *pv2 = new TPaveText(0.15,0.8,0.25,0.9,"brNDC");
    pv2->SetFillStyle(0);
    pv2->SetBorderSize(0);
    pv2->AddText(Form("R=0.%d",Rpar));

    THnSparseD *hZjetRecGen;THnSparseD *hZjetRG[NDMC];
    THnD *hZ4d; //TH2D *;

    TList *histList[NDMC];THnSparseD *sparseMC[NDMC];THnSparseD *sparsereco[NDMC];

    //Empty 2D histograms to define binning of response matrix at reco and gen level
    TH2D* hZjetRRebin = new TH2D("hRecRebin","hRecRebin", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    hZjetRRebin->Sumw2();
    TH2D* hZjetGRebin = new TH2D("hGenRebin","hGenRebin", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsN, fJetptbinsA);
    hZjetGRebin->Sumw2();

    for(int i=0; i<NDMC; i++){
        if(!isprefix){
            if(postfix) {
                histList[i] =  (TList*)dir->Get(Form("%s%d%sMCrec",histName.Data(),i,listName.Data())); }
            else {
                if(isPrompt) histList[i] =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
                else histList[i] =  (TList*)dir->Get(Form("%s%dFDMCrec",histName.Data(),i));
            }
        }
        else{
            if(postfix) {
                if(isPrompt){ histList[i] =  (TList*)dir->Get(Form("%s%sMBN%dMCrec",histName.Data(),listName.Data(),i)); }
                else{    histList[i] =  (TList*)dir->Get(Form("%s%sMBN%dFDMCrec",histName.Data(),listName.Data(),i)); }
            }
            else { cout<<"-----postfix has to be true if prefix is true!! check again----------------"<<endl; return;       }
        }

        sparseMC[i] = (THnSparseD*)histList[i]->FindObject("ResponseMatrix");
        //sparseMC[i]->GetAxis(5)->SetRangeUser(fptbinsZGenA[0],fptbinsZGenA[fptbinsZGenN]);
        sparseMC[i]->GetAxis(1)->SetRangeUser(fJetptbinsGenA[0],fJetptbinsGenA[fJetptbinsGenN]);

        //----- getting min and max of each dimension
        for (int idim=0; idim< sparseMC[i]->GetNdimensions(); idim++){
        	auto axis = sparseMC[i]->GetAxis(i);
        	int min = 0; int max = axis->GetNbins()+1;
        	axis->SetRange(min,max);
        }
        hZjetRG[i] = (THnSparseD*)sparseMC[i]->Projection(nDim,dim,"E");
        //-----
        if (!i){
                hZjetRecGen = (THnSparseD*)hZjetRG[0]->Clone("hZjetRecGen");
        }
        else {
                hZjetRecGen->Add(hZjetRG[i]);
        }
    }//end of NDMC for loop //hZjetRecGen->SaveAs("proj.root");

    /**************************
    #### normalizing
    **************************/
    TH2D *hZjetGen = (TH2D*)hZjetRecGen->Projection(3,2,"E");
    //TH2D *hZjetGen = (TH2D*)hZjetRecGen->Projection(1,0,"E");
    TH2D *hZjetGenRebin = Rebin2D("hZjetGenRebin", hZjetGen, fptbinsZGenN, fptbinsZGenA, fJetptbinsGenN, fJetptbinsGenA, 0);
    TCanvas* cZGen = new TCanvas();
    cZGen->SetLogz();
    hZjetGenRebin->Draw("colz");
    hZjetGenRebin->Draw("TEXT SAME");
    //WeightMatrixY1Y2(hZjetRecGen,hZjetGen);

    /***************************
    #### unfolding settings ####
    ***************************/
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    RooUnfoldResponse response (hZjetRRebin, hZjetGRebin);
    //----------- fill 4D histo response matrix
    //int eventcount = 0;
    for (int z = 0; z< hZjetRecGen->GetNbins();z++) {
       int coord[4]={0,0,0,0};
       double content = hZjetRecGen->GetBinContent(z,coord);
       int i = coord[0], j = coord[1], k = coord[2], m = coord[3];
       double weight = content;
       double i_center = hZjetRecGen->GetAxis(0)->GetBinCenter(i);
       double j_center = hZjetRecGen->GetAxis(1)->GetBinCenter(j);
       double k_center = hZjetRecGen->GetAxis(2)->GetBinCenter(k);
       double m_center = hZjetRecGen->GetAxis(3)->GetBinCenter(m);

       bool measurement_ok = kTRUE;
       //if (i==0 || j ==0 || k ==0 || m == 0){
       if (i_center<0 || j_center<0){
           measurement_ok = kFALSE;
           cout<<i<<"--"<<j<<endl;//"--"<<k<<"--"<<m<<endl;
           //eventcount+=1;
       }
       if (measurement_ok){
           //double newi = hZjetRRebin->GetXaxis()->GetBinCenter(hZjetRRebin->GetXaxis()->FindBin(i_center));
           //double newj = hZjetRRebin->GetYaxis()->GetBinCenter(hZjetRRebin->GetYaxis()->FindBin(j_center));
           //double newk = hZjetGRebin->GetXaxis()->GetBinCenter(hZjetGRebin->GetXaxis()->FindBin(k_center));
           //double newm = hZjetGRebin->GetYaxis()->GetBinCenter(hZjetGRebin->GetYaxis()->FindBin(m_center));
           //response.Fill(newi,newj,newk,newm,weight);
           response.Fill(i_center,j_center,k_center,m_center,weight);
       }
       else{
           //double newk = hZjetGRebin->GetXaxis()->GetBinCenter(hZjetGRebin->GetXaxis()->FindBin(k_center));
           //double newm = hZjetGRebin->GetYaxis()->GetBinCenter(hZjetGRebin->GetYaxis()->FindBin(m_center));
           //response.Miss(newk,newm,weight);
           response.Miss(k_center,m_center,weight);
       }//cout<<i<<"-"<<j<<"-"<<k<<"-"<<m<<endl;

    }
    //cout<<eventcount<<endl;

    //response.UseOverflow(overflow);

    //THnD *responseplot = (THnD*)response.Hresponse();
    //TH2D *responseplot1 = (TH2D*)responseplot->Projection(0,1);
    TH2D *fUnfoldedBayes[NTrials];
    TH2D *folded[NTrials];
    TString outName = "unfoldedSpectrum";
    TCanvas* cUnfolded = new TCanvas("cUnfolded","cUnfolded",800,600);
    cUnfolded->SetLogz();
    TLegend* leg =  new TLegend(0.15,0.5,0.30,0.85);
    leg->SetBorderSize(0);
    //------------ do unfolding NTrials times ------------
    for(Int_t ivar=0; ivar<NTrials; ivar++){//changes
    /***********************************
    ############# unfolding ##################
    ************************************/
		    RooUnfoldBayes unfold (&response, hData2D, ivar+1);
        fUnfoldedBayes[ivar] = (TH2D*)unfold.Hreco();
        folded[ivar] = (TH2D*)response.ApplyToTruth(fUnfoldedBayes[ivar]);

        fUnfoldedBayes[ivar]->SetTitle("unfolded z-jetpt spectrum");
        if(ivar == regBayes-1){fUnfoldedBayes[ivar]->Draw("colz");}
        if(ivar == regBayes-1){fUnfoldedBayes[ivar]->Draw("TEXT SAME");}
        //else{fUnfoldedBayes[ivar]->Draw("TEXT SAME");}
        //leg->AddEntry(fUnfoldedBayes[ivar],Form("Reg=%d",ivar+1),"p");
    }
    //leg->Draw("same");
	  cUnfolded->SaveAs(Form("%s/plots/%s_unfSpectra.pdf",outDir.Data(),outName.Data()));
	  cUnfolded->SaveAs(Form("%s/plots/%s_unfSpectra.png",outDir.Data(),outName.Data()));
	  cUnfolded->SaveAs(Form("%s/plots/%s_unfSpectra.svg",outDir.Data(),outName.Data()));
    //TCanvas* cUnfolded_proj = new TCanvas("cUnfolded_proj","cUnfolded_proj",800,600);
    //cUnfolded_proj->SetLogy();
    //TH1D *proj_hist = (TH1D*)fUnfoldedBayes[regBayes-1]->ProjectionX("projunfold",1,1,"E");
    //responseplot1->Draw();
    TH2D *fUnfoldedBayesScale = (TH2D*)fUnfoldedBayes[regBayes-1]->Clone();
    TH1D *UnfProjectX[4];
    TH1D *hDataProjectX[4];
    TH1D *UnfProjectXScale[4];
    TFile *unfold2DoutFile = new TFile(Form("%s/alljetz2D/unfold2DoutFile.root",outDir.Data()),"recreate");
    for (int binjet=1; binjet<= fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
      UnfProjectX[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes-1]->ProjectionX(Form("UnfProjectX_%d",binjet), binjet, binjet, "E");
      hDataProjectX[binjet-1] = (TH1D*)hData2D->ProjectionX(Form("hDatarojectX_%d",binjet), binjet, binjet, "E");
      UnfProjectXScale[binjet-1] = (TH1D*)UnfProjectX[binjet-1]->Clone();
      UnfProjectXScale[binjet-1]->Scale(1,"width");
      UnfProjectX[binjet-1]     ->Write(Form("UnfProjectX_%d",binjet-1));
      hDataProjectX[binjet-1]   ->Write(Form("hDataProjectX_%d",binjet-1));
      //UnfProjectXScale[binjet-1]->Write(Form("UnfProjectXScale_%d",binjet-1));
      for(int binz=1;binz<=fUnfoldedBayesScale->GetNbinsX();binz++){
        fUnfoldedBayesScale->SetBinContent(binz,binjet,UnfProjectXScale[binjet-1]->GetBinContent(binz));
        fUnfoldedBayesScale->SetBinError(binz,binjet,UnfProjectXScale[binjet-1]->GetBinError(binz));
      }
    }
    unfold2DoutFile->Close();
    TCanvas *cUnfoldedScale = new TCanvas("cUnfoldedScale", "cUnfoldedScale", 800, 600);
    cUnfoldedScale->SetLogz();
    fUnfoldedBayesScale->SetTitle("z-jet p_{T} spectrum (after unfolding)");
    fUnfoldedBayesScale->GetYaxis()->SetTitle("jet p_{T}");
    fUnfoldedBayesScale->GetXaxis()->SetTitle("z_{||}^{ch}");
    fUnfoldedBayesScale->Draw("colz");
    fUnfoldedBayesScale->Draw("TEXT SAME");
cout<<hData2D->GetNbinsY()<<endl;
cout<<hData2D->GetNbinsX()<<endl;

return;
}


// All functions-----------------------------------------------------------------------------------------------------------
/// Rebin 2d variable size - no such routine in Root
TH2D *Rebin2D(const char* name, TH2D *h, int nx, const double *binx, int ny, const double *biny, bool crop) {
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
    for(int j=0;j<=hre->GetNbinsY()+1;j++){
            hre->SetBinContent(0,j,0);
            hre->SetBinError(0,j,0);

						hre->SetBinContent(hre->GetNbinsX()+1,j,0);
            hre->SetBinError(hre->GetNbinsX()+1,j,0);

    }
	return hre;
}
/// load 1D raw spectrum (that is to be unfolded)
TH1D *LoadRawSpec(TString fn, TString sname, TString spostfix="") {
	TFile *f  = TFile::Open(fn);
	if (!f) { Error("LoadRawSpectrum","Raw spectrum file %s not found.",fn.Data());	return 0; }
	TH1D *spectrum = (TH1D*)f->Get(sname);
	if (!spectrum) {
		Error("LoadRawSpectrum","Raw spectrum %s could not be gotten from file.",sname.Data());
		return 0;
	}
	Info("LoadRawSpectrum", "%s loaded.", sname.Data());
	TH1D *fRawSpectrum = (TH1D*)spectrum->Clone("fRawSpectrum");
  fRawSpectrum->Sumw2();
  return fRawSpectrum;
}
/// Weigh matrix along y1(Gen Z) and y2(Gen jetpt)
//THnSparseD *WeightMatrixY1Y2(THnSparseD *hZjetRecGen, TH2D *hZjetGen){
//  for (int j = 1; j <= )
///}
//
/*
THnSparseD *ConvertToResponse(THnSparseD* rawTHnS, TH2D *BinnedRec, TH2D *BinnedGen){
  THnSparseD *BinnedRecGen;
  for(int z = 0; z < rawTHnS->GetNbins(); z++){
    int coord[4]={};
  }
  return;
}
*/
