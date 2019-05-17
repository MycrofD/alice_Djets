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
TH1D* fRawSpec2Dproj[fJetptbinsN];

TH1D *LoadRawSpec2D(TString fn, TString sname, TString spostfix="");
////===================== for 2D unfolding
////=====================
/***********************************
############# define your bins #####################
************************************/
int colortable[] = {kMagenta, kViolet, kBlue, kCyan+2, kGreen+4, kGreen+1, kYellow+1, kOrange+1, kRed, kRed+2};
int linesytle[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
/***********************************
############# begining of the macro ##################
************************************/
void unfold_Bayeszjet
(
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
){
    bool tempo_data_plot = 0;
    //-----------------------------
    gStyle->SetOptStat(0000);
    gSystem->Load(Form("%s",roounfoldpwd.Data()));
    gSystem->Exec(Form("mkdir  %s",outDir.Data()));
    gSystem->Exec(Form("mkdir  %s/Bayes",outDir.Data()));
    outDir+="/Bayes";
    gSystem->Exec(Form("mkdir  %s",outDir.Data()));
    gSystem->Exec(Form("mkdir  %s/plots",outDir.Data()));
    gSystem->Exec(Form("mkdir  %s/alljetz2D",outDir.Data()));
    //-----------------------------
    // Extracting data 2D histogram in z and jetpt. Before it was just in z...  for unfolding.
    TH2* hData2D = new TH2D("hFD_zjet", "hFD_zjet", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
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
        fRawSpec2Dproj[binjet] = (TH1D*)LoadRawSpec2D(data2D[binjet].Data(),"hData_binned_sub");
        fRawSpec2Dproj[binjet]->Scale(1,"width");
        for (Int_t binz=0; binz < fptbinsZMeasN+1; binz++){
        	double cont = fRawSpec2Dproj[binjet]->GetBinContent(binz);
        	double contErr = fRawSpec2Dproj[binjet]->GetBinError(binz);
        	hData2D->SetBinContent(binz,binjet+1,cont);
        	hData2D->SetBinError(binz,binjet+1,contErr);
        }
    }
    // Saving the 2D data file.
    TCanvas* cFD_2D = new TCanvas();
    cFD_2D->SetLogz();
    hData2D->Draw("colz");
    hData2D->Draw("TEXT SAME");
    cFD_2D->SaveAs(Form("%s/alljetz2D/FDdata2D.pdf",outDir.Data()));
    cFD_2D->SaveAs(Form("%s/alljetz2D/FDdata2D.png",outDir.Data()));
    cFD_2D->SaveAs(Form("%s/alljetz2D/FDdata2D.svg",outDir.Data()));
    TFile *outFile = new TFile(Form("%s/alljetz2D/outFD.root",outDir.Data()),"recreate");
    hData2D->Write();
    outFile->Close();
    //-----------------------------
    // Time for Response Matrix in 2D
    TFile *File = new TFile(effFile,"read");
    TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
    TString histName;
    if(!isprefix){
        if(fDmesonSpecie) histName = "histosDStarMBN";
        else histName = "histosD0MBN";}
    else{
        if(fDmesonSpecie) histName = "histosDStarMBN";
        else histName = "histosD0";}
    double jetmin = 0, jetmax = 60;
    TPaveText *pv2 = new TPaveText(0.15,0.8,0.25,0.9,"brNDC");
    pv2->SetFillStyle(0);
    pv2->SetBorderSize(0);
    pv2->AddText(Form("R=0.%d",Rpar));

    TH2D* hZjetG[NDMC];TH2D* hZjetR[NDMC];
    TH2D* hZjetGen;TH2D* hZjetRec;
    THnSparseD* hZjetRecGen;THnSparseD* hZjetRG[NDMC];
    TList *histList[NDMC];THnSparseF *sparseMC[NDMC];THnSparseF *sparsereco[NDMC];

    TH2D* hZjetGRebin = new TH2D("hGenRebin","hGenRebin", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    hZjetGRebin->Sumw2();
    TH2D* hZjetRRebin = new TH2D("hRecRebin","hRecRebin", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsN, fJetptbinsA);
    hZjetRRebin->Sumw2();

    const int nDim = 4;
    int dim[nDim] = {0,1,5,6};

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

        sparseMC[i] = (THnSparseF*)histList[i]->FindObject("ResponseMatrix");
        //sparseMC[i]->GetAxis(6)->SetRangeUser(jetmin,jetmax);
        hZjetG[i] = (TH2D*)sparseMC[i]->Projection(6,5,"E");
        hZjetR[i] = (TH2D*)sparseMC[i]->Projection(1,0,"E");

        hZjetR[i]->Sumw2();
        hZjetG[i]->Sumw2();
        hZjetG[i]->SetName(Form("hZjetG_%d",i));
        hZjetR[i]->SetName(Form("hZjetR_%d",i));

        //----- getting min and max of each dimension
        for (int idim=0; idim< sparseMC[i]->GetNdimensions(); idim++){
        	auto axis = sparseMC[i]->GetAxis(i);
        	int min = 0; int max = axis->GetNbins()+1;
        	axis->SetRange(min,max);
        }
        hZjetRG[i] = (THnSparseD*)sparseMC[i]->Projection(nDim,dim,"E");
        //-----
        if (!i){
                hZjetGen    = (TH2D*)hZjetG[0] ->Clone("hZjetGen");
                hZjetRec    = (TH2D*)hZjetR[0] ->Clone("hZjetRec");
                hZjetRecGen = (THnSparseD*)hZjetRG[0]->Clone("hZjetRecGen");
        }
        else {
                hZjetGen   ->Add(hZjetG[i]);
                hZjetRec   ->Add(hZjetR[i]);
                hZjetRecGen->Add(hZjetRG[i]);
        }
    }//end of NDMC for loop //hZjetRecGen->SaveAs("proj.root");

    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    RooUnfoldResponse response (hZjetRRebin, hZjetGRebin);
    //----------- fill 4D histo response matrix
    for (int z = 0; z< hZjetRecGen->GetNbins();z++) {
       int coord[4]={0,0,0,0};
       double content = hZjetRecGen->GetBinContent(z,coord);
       int i = coord[0];
       int j = coord[1];
       int k = coord[2];
       int m = coord[3];

       double i_center = hZjetRecGen->GetAxis(0)->GetBinCenter(i);
       double j_center = hZjetRecGen->GetAxis(1)->GetBinCenter(j);
       double k_center = hZjetRecGen->GetAxis(2)->GetBinCenter(k);
       double m_center = hZjetRecGen->GetAxis(3)->GetBinCenter(m);

       bool measurement_ok = kTRUE;
       if (i==0 || j ==0 || k ==0 || m == 0)
           measurement_ok = kFALSE;

       if (measurement_ok){
           response.Fill(i_center,j_center,k_center,m_center,content);
       }
       else{
           response.Miss(k_center,m_center,content);
       }//cout<<i<<"-"<<j<<"-"<<k<<"-"<<m<<endl;
    }
    response.UseOverflow(overflow);
return;
}


// All functions----------------------------------------------------------------
/// load raw spectrum (that is to be unfolded)
TH1D *LoadRawSpec2D(TString fn, TString sname, TString spostfix="") {
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
