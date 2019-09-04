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
// to fill the response matrix
const int fJetptbinsGenN=9; TString truebinnum="";
double fJetptbinsGenA[fJetptbinsGenN+1]={0,2,3,4,5,7,10,15,50,100};
//const int fJetptbinsGenN=6; TString truebinnum="550";
//double fJetptbinsGenA[fJetptbinsGenN+1]={5,7,10,15,50,50};
//const int fJetptbinsGenN=8; TString truebinnum="360";
//double fJetptbinsGenA[fJetptbinsGenN+1]={3,4,5,7,10,15,50,60};
//const int fJetptbinsGenN=7; TString truebinnum="350";
//double fJetptbinsGenA[fJetptbinsGenN+1]={3,4,5,7,10,15,50};
const int nDim    = 4;//the four dimensions
//int dim[nDim]     = {0,1,5,6};//for extacting 4D info from THnSparse
const int DnDim   = 6;//the six dimensions
//int Ddim[DnDim]   = {0,1,5,6,2,7};//for extacting 6D info from THnSparse
int zRec = 0, jetRec = 1, DRec = 2;
int zGen = 5, jetGen = 6, DGen = 7; 
int Ddim[DnDim]   = {zRec,jetRec,zGen,jetGen,DRec,DGen};//for extacting 6D info from THnSparse
int dim[nDim]   = {zRec,jetRec,zGen,jetGen};//for extacting 4D info from THnSparse
TH1D *fRawSpec2Dproj[fJetptbinsN];
TH1D *fGenSpec2Dproj[fJetptbinsN];
TH1D *fRawSpec2DprojUp[fJetptbinsN];//FD systematics
TH1D *fRawSpec2DprojDo[fJetptbinsN];//FD systematics

TH1D *LoadRawSpec(TString fn, TString sname, TString spostfix="");
TH2D *Rebin2D(const char* name, TH2D *h, int nx, const double *binx, int ny, const double *biny, bool crop);
////===================== for 2D unfolding
int colortable[] = {kMagenta, kViolet, kBlue, kCyan+2, kGreen+4, kGreen+1, kYellow+1, kOrange+1, kRed, kRed+2};
int linesytle[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
/*****************************************************
 ############ begining of the macro ##################
 *****************************************************/
void closureMCadv(
  //TString roounfoldpwd = "",
  TString listName = "FD",
  bool isPrompt=1,
  bool postfix=0,
  bool isprefix=0,
  TString effFile = "/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_642_pp5TeV_z.root",
  TString datafile = "file.root",
  TString detRMfile = "detRM.root",
  TString bkgRMfile = "bkgRM.root",
  TString outDir = "MCclosure", // output directory
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
    bool DptcutTrue=1;
    //-----------------------------
    gStyle ->SetOptStat(0000);
//    gSystem->Exec(Form("mkdir  %s",outDir.Data()));
//    gSystem->Exec(Form("mkdir  %s/BayesClosure",outDir.Data()));
//    outDir+="/BayesClosure";
//    gSystem->Exec(Form("mkdir  %s",outDir.Data()));
//    gSystem->Exec(Form("mkdir  %s/plots",outDir.Data()));
//    gSystem->Exec(Form("mkdir  %s/alljetz2D",outDir.Data()));
    //===============================
    // Time for Response Matrix in 2D
    //===============================
    //------ reading the MC (effFile)-------------------------
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

    //--- Getting the THnSparse --------------------
    THnSparseD *hZjetRecGen;THnSparseD *hZjetRG[NDMC]; // RecGen and RG[NDMC]
    //THnSparseD *hZjetRecGen;THnSparseD *hZjetRecGen_D2;THnSparseD *hZjetRecGen_D3;THnSparseD *hZjetRecGen_D5;
    //THnSparseD *hZjetRG_D2[NDMC];THnSparseD *hZjetRG_D3[NDMC];THnSparseD *hZjetRG_D5[NDMC];
    THnSparseD *hZjetRecGenD;
    THnSparseD *hZjetRGD[NDMC];
    //---- Storing the response matrix for drawing purposes ------------

    TList *histList[NDMC];THnSparseD *sparseMC[NDMC];THnSparseD *sparsereco[NDMC];
    //--- Empty 2D histograms to define binning of response matrix at reco and gen level --------------------
    TH2D* hZjetRRebin = new TH2D("hRecRebin","hRecRebin", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    //TH2D* hZjetRRebin = new TH2D("hRecRebin","hRecRebin", fptbinsZTrueNN, fptbinsZTrueAA, fJetptbinsN, fJetptbinsA);
    hZjetRRebin->Sumw2();
    //TH2D* hZjetGRebin = new TH2D("hGenRebin","hGenRebin", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsN, fJetptbinsA);
    //TH2D* hZjetGRebin = new TH2D("hGenRebin","hGenRebin", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    //TH2D* hZjetGRebin = new TH2D("hGenRebin","hGenRebin", fptbinsZTrueN, fptbinsZTrueA, fJetptbinsN, fJetptbinsA);
    TH2D* hZjetGRebin = new TH2D("hGenRebin","hGenRebin", fptbinsZTrueNN, fptbinsZTrueAA, fJetptbinsGenN, fJetptbinsGenA);
    hZjetGRebin->Sumw2();
    
    //Kinematic efficiency
    //--------------------
    TH2D* kineEffCut  = new TH2D("hkineEffCut","hkineEffCut", fptbinsZTrueNN, fptbinsZTrueAA, fJetptbinsGenN, fJetptbinsGenA);
    kineEffCut->Sumw2();
    TH2D* kineEffFull  = new TH2D("hkineEffFull","hkineEffFull", fptbinsZTrueNN, fptbinsZTrueAA, fJetptbinsGenN, fJetptbinsGenA);
    kineEffFull->Sumw2();

    //---- Reading from the ResponseMatrix of each D meson in a loop
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
        sparseMC[i]->GetAxis(1)->SetRangeUser(fJetptbinsGenA[0],fJetptbinsGenA[fJetptbinsGenN]);
        sparseMC[i]->GetAxis(4)->SetRangeUser(-0.9+fRpar,0.9-fRpar);//pseudo-rapidity cuts
        sparseMC[i]->GetAxis(9)->SetRangeUser(-0.9+fRpar,0.9-fRpar);//pseudo-rapidity cuts

        //----- getting min and max of each dimension.... //is it necessary?
        for (int idim=0; idim< sparseMC[i]->GetNdimensions(); idim++){
        	auto axis = sparseMC[i]->GetAxis(i);
        	int min = 0; int max = axis->GetNbins()+1;
        	axis->SetRange(min,max);
        }
        hZjetRG[i] = (THnSparseD*)sparseMC[i]->Projection(nDim,dim,"E");// Just z and jet components
        hZjetRGD[i] = (THnSparseD*)sparseMC[i]->Projection(DnDim,Ddim,"E");// z, jet and D components
        //-----
        if (!i){
                hZjetRecGen = (THnSparseD*)hZjetRG[0]->Clone("hZjetRecGen");
                hZjetRecGenD = (THnSparseD*)hZjetRGD[0]->Clone("hZjetRecGenD");
        }
        else {
                hZjetRecGen->Add(hZjetRG[i]);
                hZjetRecGenD->Add(hZjetRGD[i]);
        }
    }//end of NDMC for loop //hZjetRecGen->SaveAs("proj.root");


    // Separating the MC for response matrix and reco data.
    // -----------------------------------------------------

    //int eventcount = 0;
    double binval = 0;
    double splitPC = 0.7;
    //----------- fill 4D histo response matrix
    int totalBins = hZjetRecGenD->GetNbins();
    int responseBins = (int)(totalBins*splitPC);
    int MCdataBins =  totalBins - responseBins;
    int sortsize = 0;

    //Getting an array of random numbers
    int responseBinsFinal = 0;

    mt19937 mt_resp(random_device{}());
    uniform_int_distribution<int> dist_resp(0,totalBins);
    vector<int> result_resp;
    result_resp.reserve(responseBins);
    vector<int> result_reco;
//    result_reco.reserve(MCdataBins);
    set<int> seen_resp;

    for (int i = 0; i< responseBins;){
        auto n = dist_resp(mt_resp);
        if(seen_resp.insert(n).second){
            result_resp.push_back(n);
            ++i;
        }
    }

    for (int i=1; i<totalBins ; i++){
        if(find(result_resp.begin(), result_resp.end(), i)!=result_resp.end())
            i+=1;
        else result_reco.insert (result_reco.end(),i);
    }
    cout<<result_resp.size()<<endl;
    cout<<result_reco.size()<<endl;
    cout<<result_reco.size() + result_resp.size()<<endl;
    cout<<hZjetRecGenD->GetNbins()<<endl;

    /***************************
    #### MC data settings ####
    ***************************/
    const int minDim = 6;
    int bins[minDim] = {200,110,220,200,110,220}; //z,jet,d rec and gen
    double xmin[minDim] = {-2,-10,-10,-2,-10,-10};
    double xmax[minDim] = {2,100,100,2,100,100};
    /*
//    THnSparseD mcdata00[5]("mcdata","mcdata",minDim,bins,xmin,xmax);
    THnSparseD mcdata[5];
    mcdata[0] = new THnSparseD("mcdata","mcdata",minDim,bins,xmin,xmax);

    for (int i=0; i<5; i++){
//        mcdata[i] = (THnSparseD)mcdata00.Clone();
    
    }
    */
    THnSparseD mcdata0("mcdata","mcdata",minDim,bins,xmin,xmax);
    THnSparseD mcdata1("mcdata","mcdata",minDim,bins,xmin,xmax);
    THnSparseD mcdata2("mcdata","mcdata",minDim,bins,xmin,xmax);
    THnSparseD mcdata3("mcdata","mcdata",minDim,bins,xmin,xmax);
    mcdata0.Sumw2();
    mcdata1.Sumw2();
    mcdata2.Sumw2();
    mcdata3.Sumw2();

    for(int z=1; z<=result_reco.size(); z++){
        int coord[DnDim]={0,0,0,0,0,0};
        double content = hZjetRecGenD->GetBinContent(result_reco[z-1],coord);
        int i = coord[0], j = coord[1], k = coord[2], m = coord[3], n = coord[4], p = coord[5];
        double weight = content;
        double zR_center = hZjetRecGenD->GetAxis(0)->GetBinCenter(i);
        double jR_center = hZjetRecGenD->GetAxis(1)->GetBinCenter(j);
        double zG_center = hZjetRecGenD->GetAxis(2)->GetBinCenter(k);
        double jG_center = hZjetRecGenD->GetAxis(3)->GetBinCenter(m);
        double DR_center = hZjetRecGenD->GetAxis(4)->GetBinCenter(n);
        double DG_center = hZjetRecGenD->GetAxis(5)->GetBinCenter(p);
        double xbin[] = {zR_center,jR_center,DR_center,zG_center,jG_center,DG_center};
        mcdata0.Fill(xbin,weight);
        mcdata1.Fill(xbin,weight);
        mcdata2.Fill(xbin,weight);
        mcdata3.Fill(xbin,weight);
    }

    //---------------------------------
    // Creating data histos from mcdata
    //----------------------------------
    //const int     fJetptbinsN = 4;  
    //double        fJetptbinsA[fJetptbinsN+1] = {3.0, 5.0, 7.0, 10.0, 15.0, 50.0};
    //
    TH1D* hZR[4]; TH1D* hZG[4]; // Z distributions for different jetpt intervals

    double dmin[] = {2.0,2.0,3.0,5.0,5.0};
    mcdata0.GetAxis(1)->SetRangeUser(fJetptbinsA[0],fJetptbinsA[1]);
    mcdata1.GetAxis(1)->SetRangeUser(fJetptbinsA[1],fJetptbinsA[2]);
    mcdata2.GetAxis(1)->SetRangeUser(fJetptbinsA[2],fJetptbinsA[3]);
    mcdata3.GetAxis(1)->SetRangeUser(fJetptbinsA[3],fJetptbinsA[4]);
    mcdata0.GetAxis(2)->SetRangeUser(dmin[0],fJetptbinsA[1]);
    mcdata1.GetAxis(2)->SetRangeUser(dmin[1],fJetptbinsA[2]);
    mcdata2.GetAxis(2)->SetRangeUser(dmin[2],fJetptbinsA[3]);
    mcdata3.GetAxis(2)->SetRangeUser(dmin[3],fJetptbinsA[4]);

    hZR[0]=(TH1D*)mcdata0.Projection(0);
    hZG[0]=(TH1D*)mcdata0.Projection(3);
    hZR[1]=(TH1D*)mcdata1.Projection(0);
    hZG[1]=(TH1D*)mcdata1.Projection(3);
    hZR[2]=(TH1D*)mcdata2.Projection(0);
    hZG[2]=(TH1D*)mcdata2.Projection(3);
    hZR[3]=(TH1D*)mcdata3.Projection(0);
    hZG[3]=(TH1D*)mcdata3.Projection(3);

    TH1D* hzRreb_tmp[4];
    TH1D* hzGreb_tmp[4];
    TH1D* hzRreb[4];
    TH1D* hzGreb[4];
    for (int i=0; i<4; i++){
        hzRreb_tmp[i] = (TH1D*)hZR[i]->Clone(Form("hzRreb_tmp%d",i));
        hzGreb_tmp[i] = (TH1D*)hZG[i]->Clone(Form("hzGreb_tmp%d",i));
        hzRreb[i] = (TH1D*)hzRreb_tmp[i]->Rebin(fptbinsZMeasN,"hzRreb",fptbinsZMeasA);
        hzGreb[i] = (TH1D*)hzGreb_tmp[i]->Rebin(fptbinsZMeasN,"hzGreb",fptbinsZMeasA);
    
    }
    
    //==========================================================
    // Extracting data 2D histogram in z and jetpt. Before it was just in z...  for unfolding.
    //==========================================================
    TH2* hData2D = new TH2D("hFD_zjet", "hFD_zjet", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    TH2* hData2DS = new TH2D("hFD_zjet", "hFD_zjet", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    TH2* hData2DUp = new TH2D("hFD_zjetUp", "hFD_zjetUp", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    TH2* hData2DDo = new TH2D("hFD_zjetDo", "hFD_zjetDo", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    TString data2D[fJetptbinsN];
    // ---------------------------------------
    // reading data from each 1D
    // ---------------------------------------
    for (int binjet=0; binjet < fJetptbinsN; binjet++){
        // ---------------------------------------
        // 1D projections of FD subtracted spectra
        // ---------------------------------------
        //cout<<fJetptbinsA[binjet]<<endl;
        fRawSpec2Dproj[binjet] = (TH1D*)hzRreb[binjet]->Clone();
        fGenSpec2Dproj[binjet] = (TH1D*)hzGreb[binjet]->Clone();
        TH1D *fRawSpec2DprojScale=(TH1D*)fRawSpec2Dproj[binjet]->Clone();
        fRawSpec2DprojScale->Scale(1,"width");
        // --------------------------------------------------
        // combining those 1D projections into a 2D histogram
        // --------------------------------------------------
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
    // ---------------------------------------
    // Saving the 2D data file.
    // ---------------------------------------
    TCanvas *cFD_2D = new TCanvas("Measured2D", "Measured2D", 800, 600);
    cFD_2D  ->SetLogz();
    hData2DS->SetTitle("z-jet p_{T} spectrum (before unfolding)");
    hData2DS->GetYaxis()->SetTitle("jet p_{T}");
    hData2DS->GetXaxis()->SetTitle("z_{||}^{ch}");
    hData2DS->Draw("colz");
    hData2DS->Draw("TEXT SAME");
    cFD_2D  ->SaveAs(Form("%s/FDdata2D.pdf",outDir.Data()));
    cFD_2D  ->SaveAs(Form("%s/FDdata2D.png",outDir.Data()));
    cFD_2D  ->SaveAs(Form("%s/FDdata2D.svg",outDir.Data()));
    TFile *outFile = new TFile(Form("%s/outFD.root",outDir.Data()),"recreate");
    hData2D ->Write();//raw 2D
    hData2DS->Write();//scaled 2D
    outFile ->Close();
    /***************************
    #### unfolding settings ####
    ***************************/
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    RooUnfoldResponse response (hZjetRRebin, hZjetGRebin);
    //response.UseOverflow(1);
    //----------- fill 4D histo response matrix

    if(DptcutTrue){
//        for (int z = 0; z < hZjetRecGenD->GetNbins();z++) {
        for (int z = 0; z < responseBins;z++) {
            int coord[DnDim]={0,0,0,0,0,0};
            //double content = hZjetRecGenD->GetBinContent(responseResult[z],coord);
            double content = hZjetRecGenD->GetBinContent(result_resp[z],coord);
            //int Ddim[DnDim]   = {zRec,jetRec,zGen,jetGen,DRec,DGen};//for extacting 6D info from THnSparse
            int i = coord[0], j = coord[1], k = coord[2], m = coord[3], n = coord[4], p = coord[5];
            double weight = content;
            double zR_center = hZjetRecGenD->GetAxis(0)->GetBinCenter(i);
            double jR_center = hZjetRecGenD->GetAxis(1)->GetBinCenter(j);
            double zG_center = hZjetRecGenD->GetAxis(2)->GetBinCenter(k);
            double jG_center = hZjetRecGenD->GetAxis(3)->GetBinCenter(m);
            double DR_center = hZjetRecGenD->GetAxis(4)->GetBinCenter(n);
            double DG_center = hZjetRecGenD->GetAxis(5)->GetBinCenter(p);

            bool measurement_ok = kTRUE;
            if( DR_center<2 && jR_center>=5 ){
                  measurement_ok = kFALSE;
            }
            if( DR_center<3 && jR_center>=7 ){
                  measurement_ok = kFALSE;
            }
            if(DR_center<5 && jR_center>=10 ){
                  measurement_ok = kFALSE;
            }
            if (measurement_ok){
                response.Fill(zR_center,jR_center,zG_center,jG_center,weight);
            }
            //else{
            //    response.Miss(zG_center,jG_center,weight);
            //}//cout<<i<<"-"<<j<<"-"<<k<<"-"<<m<<endl;
        }
    }

// checking to shorten the code
    else{
        for (int z = 0; z< hZjetRecGen->GetNbins();z++) {
            int coord[nDim]={0,0,0,0};
            double content = hZjetRecGen->GetBinContent(z,coord);
            int i = coord[0], j = coord[1], k = coord[2], m = coord[3];
            double weight = content;
            double zR_center = hZjetRecGen->GetAxis(0)->GetBinCenter(i);
            double jR_center = hZjetRecGen->GetAxis(1)->GetBinCenter(j);
            double zG_center = hZjetRecGen->GetAxis(2)->GetBinCenter(k);
            double jG_center = hZjetRecGen->GetAxis(3)->GetBinCenter(m);

            bool measurement_ok = kTRUE;
            if (measurement_ok){
                response.Fill(zR_center,jR_center,zG_center,jG_center,weight);
            }
        }
    }

		
    TH2D *fUnfoldedBayes[NTrials];
    TH2D *folded[NTrials];
    TString outName = "unfoldedSpectrum";
//    TCanvas* cUnfolded = new TCanvas("cUnfolded","cUnfolded",800,600);
//    cUnfolded->SetLogz();
//    TLegend* leg =  new TLegend(0.15,0.5,0.30,0.85);
//    leg->SetBorderSize(0);
    //------------ do unfolding NTrials times ------------
    for(Int_t ivar=0; ivar<NTrials; ivar++){//changes
        /***********************************
        ############# unfolding ##################
        ************************************/
        RooUnfoldBayes unfold (&response, hData2D, ivar+1);
        fUnfoldedBayes[ivar] = (TH2D*)unfold.Hreco();
        folded[ivar] = (TH2D*)response.ApplyToTruth(fUnfoldedBayes[ivar]);
    
        fUnfoldedBayes[ivar]->SetTitle("unfolded z-jetpt spectrum");
//        if(ivar == regBayes-1){fUnfoldedBayes[ivar]->Draw("colz");}
//        if(ivar == regBayes-1){fUnfoldedBayes[ivar]->Draw("TEXT SAME");}
            //else{fUnfoldedBayes[ivar]->Draw("TEXT SAME");}
            //leg->AddEntry(fUnfoldedBayes[ivar],Form("Reg=%d",ivar+1),"p");
    }
    //leg->Draw("same");
//    cUnfolded->SaveAs(Form("%s/plots/%s_unfSpectra.pdf",outDir.Data(),outName.Data()));
//    cUnfolded->SaveAs(Form("%s/plots/%s_unfSpectra.png",outDir.Data(),outName.Data()));
//    cUnfolded->SaveAs(Form("%s/plots/%s_unfSpectra.svg",outDir.Data(),outName.Data()));
//
//
//    /***********************************
//    // FD systematics begin
//    ************************************/
//    /***********************************
//    // FD systematics end
//    ************************************/
//
//    TH2D *fUnfoldedBayesScale = (TH2D*)fUnfoldedBayes[regBayes-1]->Clone();
    TH2D *fUnfoldedBayesScale = (TH2D*)hData2D->Clone();
    TH1D *UnfProjectX[fJetptbinsN];
    TH1D *hDataProjectX[fJetptbinsN];
    TH1D *UnfProjectXScale[fJetptbinsN];
    //************************************Unfolding systematics
    TH1D *UnfProjectXIterPre[fJetptbinsN];
    TH1D *UnfProjectXIterPost[fJetptbinsN];
    //************************************
    TH1D *UnfProjectXUp[fJetptbinsN];
    TH1D *hDataProjectXUp[fJetptbinsN];
    TH1D *UnfProjectXScaleUp[fJetptbinsN];
    TH1D *UnfProjectXDo[fJetptbinsN];
    TH1D *hDataProjectXDo[fJetptbinsN];
    TH1D *UnfProjectXScaleDo[fJetptbinsN];
    //************************************
    // Saving response matrix    *********
    //************************************

    TH2D* respobj = (TH2D*)response.Hresponse()->Clone();
    TFile *responseobject = new TFile(Form("%s/responseobject.root",outDir.Data()),"recreate");
    //response.SetName("responseNow");
    respobj->SetName("responseNow");
    response.Write();
    response.Hresponse()->Write();
 //   response.Hresponse()->SetName("resp");
    respobj->Write();
    responseobject->Close();
    TCanvas *cRespon = new TCanvas("","",800,600);
    TH2D* respobjP = (TH2D*)respobj->Clone();
    respobjP->Scale(1.0/respobj->Integral());
    cRespon->SetLogz();
    respobjP->Draw("");
    //respobjP->Draw("colz");
    cRespon->SaveAs(Form("%s/responseobj.pdf",outDir.Data()));


    //-----------Kinematic eff and canvas
    if(DptcutTrue){
        for (int z = 0; z < responseBins;z++) {
            int coord[DnDim]={0,0,0,0,0,0};
            //double content = hZjetRecGenD->GetBinContent(responseResult[z],coord);
            double content = hZjetRecGenD->GetBinContent(result_resp[z],coord);
            //int Ddim[DnDim]   = {zRec,jetRec,zGen,jetGen,DRec,DGen};//for extacting 6D info from THnSparse
            int i = coord[0], j = coord[1], k = coord[2], m = coord[3], n = coord[4], p = coord[5];
            double weight = content;
            double zR_center = hZjetRecGenD->GetAxis(0)->GetBinCenter(i);
            double jR_center = hZjetRecGenD->GetAxis(1)->GetBinCenter(j);
            double zG_center = hZjetRecGenD->GetAxis(2)->GetBinCenter(k);
            double jG_center = hZjetRecGenD->GetAxis(3)->GetBinCenter(m);
            double DR_center = hZjetRecGenD->GetAxis(4)->GetBinCenter(n);
            double DG_center = hZjetRecGenD->GetAxis(5)->GetBinCenter(p);

            bool measurement_ok = kTRUE;
            if( DR_center<2 && jR_center>=5 ){
                  measurement_ok = kFALSE;
            }
            if( DR_center<3 && jR_center>=7 ){
                  measurement_ok = kFALSE;
            }
            if( DR_center<5 && jR_center>=10 ){
                  measurement_ok = kFALSE;
            }
            if( measurement_ok ){
                kineEffFull->Fill(zG_center,jG_center,weight);
            }
            if( jR_center < 5  || jR_center > 50){
                measurement_ok = kFALSE;           
            }
            if( measurement_ok ){
                kineEffCut->Fill(zG_center,jG_center,weight);
            }
        }
    }
    TCanvas* cKine = new TCanvas("hKine","hKine",800,600);
    TH2D* hKine = (TH2D*) kineEffCut->Clone();
    hKine->Divide(kineEffFull);
    hKine->GetYaxis()->SetRangeUser(5,50);
    hKine->GetXaxis()->SetRangeUser(0.4,1.02);
    hKine->Draw("colz");
    hKine->Draw("TEXT same");

    //************************************
    //Finding the projected and also scaled histograms of unfolded ones
    //************************************
    int deltabin = 4;
    for (int binjet=1; binjet<= hData2D->GetNbinsY(); binjet++){
        UnfProjectX[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes-1]->ProjectionX(Form("UnfProjectX_%d",binjet-1), binjet+deltabin, binjet+deltabin, "E");
	    //UnfProjectX[binjet-1]->SetName(Form("UnfProjectX_%d",binjet-1));
        hDataProjectX[binjet-1] = (TH1D*)hData2D->ProjectionX(Form("hDataProjectX_%d",binjet-1), binjet+deltabin, binjet+deltabin, "E");
        UnfProjectXScale[binjet-1] = (TH1D*)UnfProjectX[binjet-1]->Clone();
        UnfProjectXScale[binjet-1]->Scale(1,"width");
        UnfProjectXIterPre[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes-2]->ProjectionX(Form("UnfProjectXIterPre_%d",binjet-1), binjet+deltabin, binjet+deltabin, "E");
        //************************************Unfolding systematics
        UnfProjectXIterPost[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes]->ProjectionX(Form("UnfProjectXIterPost_%d",binjet-1), binjet+deltabin, binjet+deltabin, "E");
        //************************************
        for(int binz=1;binz<=hData2D->GetNbinsX();binz++){
            fUnfoldedBayesScale->SetBinContent(binz,binjet,UnfProjectXScale[binjet-1]->GetBinContent(binz+4));
            fUnfoldedBayesScale->SetBinError(binz,binjet,UnfProjectXScale[binjet-1]->GetBinError(binz+4));
        }
    }
    //******************
    //Writing the above histograms separately in a file
    //******************
    TFile *unfold2DoutFile = new TFile(Form("%s/unfold2DoutFile%sAPW.root",outDir.Data(),truebinnum.Data()),"recreate");
    for (int binjet=1; binjet<= hData2D->GetNbinsY(); binjet++){
        UnfProjectX[binjet-1]     ->Write();
        hDataProjectX[binjet-1]   ->Write();
        UnfProjectXScale[binjet-1]->Write();
        UnfProjectXIterPre[binjet-1]     ->Write();
        //************************************Unfolding systematics
        UnfProjectXIterPost[binjet-1]     ->Write();
    }
    unfold2DoutFile->Close();
    // Plotting scaled 2D unfolded spectrum
    // ------------------------------------
    TCanvas *cUnfoldedScale = new TCanvas("cUnfoldedScale", "cUnfoldedScale", 800, 600);
    cUnfoldedScale->SetLogz();
    fUnfoldedBayesScale->SetTitle("z-jet p_{T} spectrum (after unfolding)");
    fUnfoldedBayesScale->GetYaxis()->SetTitle("jet p_{T}");
    fUnfoldedBayesScale->GetXaxis()->SetTitle("z_{||}^{ch}");
    fUnfoldedBayesScale->Draw("colz");
    fUnfoldedBayesScale->Draw("TEXT SAME");
    //fUnfoldedBayes[regBayes-1]->Draw("TEXT SAME");
    cUnfoldedScale->SaveAs(Form("%s/%s_unfSpectraScale.pdf",outDir.Data(),outName.Data()));
    cUnfoldedScale->SaveAs(Form("%s/%s_unfSpectraScale.png",outDir.Data(),outName.Data()));
    cUnfoldedScale->SaveAs(Form("%s/%s_unfSpectraScale.svg",outDir.Data(),outName.Data()));

    // Unfolded Reco v. Gen
    // --------------------
    TH1D* hUnfP[fJetptbinsN];
    TCanvas* cGenUnf = new TCanvas("MCClosure","MClosure",1400,900);
    cGenUnf->Divide(3,2);
    for (int binjet=1; binjet <= hData2D->GetNbinsY(); binjet++){
        cGenUnf->cd(binjet);
        cGenUnf->cd(binjet)->SetLogy();
        hUnfP[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes-1]->ProjectionX(Form("UnfProjectX_%d",binjet-1), binjet+deltabin, binjet+deltabin, "E");
        hUnfP[binjet-1]->Sumw2();
   //     hData2Dplot->SetLineWidth(2);
   //     hData2Dplot->SetMarkerStyle(20);
   //     hData2Dplot->SetMarkerSize(1);
        hUnfP[binjet-1]->Draw();
        fGenSpec2Dproj[binjet-1]->Draw("same");// = (TH1D*)LoadRawSpec(data2D[binjet].Data(),"hzGreb");
    }
    cout<<deltabin<<endl;
    return;
    TLegend* lGenUnf = new TLegend(0.25,0.25,0.85,0.85);
    cGenUnf->cd(fUnfoldedBayes[regBayes-1]->GetNbinsY()+1);
    lGenUnf->AddEntry(hUnfP[fUnfoldedBayes[regBayes-1]->GetNbinsY()-1],"Unfolded","l");
    lGenUnf->AddEntry(fGenSpec2Dproj[fUnfoldedBayes[regBayes-1]->GetNbinsY()-1],"True","l");
    lGenUnf->Draw();

//    TCanvas* ctrial = new TCanvas("a","a",800,800);
//    hUnfP[4]->Draw();
//    fGenSpec2Dproj[4]->Draw("same");
//
//    return;


    TString jetpttitle[5] = {"3-5","5-7","7-10","10-15","15-50"};

//--------- TRatioPlot
    TCanvas* cTRp = new TCanvas("Closure","Closure",1400,900);
    cTRp->Divide(3,2);
    for (int binjet=1; binjet <= fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
        cTRp->cd(binjet);
        //hUnfP[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes-1]->ProjectionX(Form("UnfProjectX_%d",binjet-1), binjet, binjet, "E");
        hUnfP[binjet-1]->SetTitle(Form("Jet pt: %s",jetpttitle[binjet-1].Data()));
        hUnfP[binjet-1]->SetLineColor(kRed+2);
        hUnfP[binjet-1]->SetMarkerColor(kRed+2);
        hUnfP[binjet-1]->SetMarkerStyle(20);
//        fGenSpec2Dproj[binjet-1]->SetMarkerColor(kBlue+2);
        fGenSpec2Dproj[binjet-1]->SetMarkerStyle(23);
        auto rp = new TRatioPlot(hUnfP[binjet-1], fGenSpec2Dproj[binjet-1]);//,"divsym");
        rp->Draw();
        rp->GetUpperPad()->cd();
        hUnfP[binjet-1]->Draw("same");
    cTRp->cd(binjet-1)->Update();
    }
    TLegend* lGenUnf2 = new TLegend(0.4,0.4,0.85,0.85);
    cTRp->cd(fUnfoldedBayes[regBayes-1]->GetNbinsY()+1);
    lGenUnf2->AddEntry(hUnfP[fUnfoldedBayes[regBayes-1]->GetNbinsY()-1],"Unfolded","l");
    lGenUnf2->AddEntry(fGenSpec2Dproj[fUnfoldedBayes[regBayes-1]->GetNbinsY()-1],"True","l");
    lGenUnf2->Draw();

return;
    // Folded/Measured spectra comparison
    // ----------------------------------
    TCanvas *cFolded = new TCanvas("Folded/Measured","", 1400, 900);
    cFolded->Divide(3,2);
    for (int binjet=1; binjet<= fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
        /***********************************
        ############# folding ##################
        ************************************/
        cFolded->cd(binjet);

        TH1D* hData2Dplot = (TH1D*)hData2D->ProjectionX("",binjet,binjet,"E")->Clone();
        //hData2Dplot->SetLineColor(kBlack);
        hData2Dplot->SetLineWidth(2);
        hData2Dplot->SetMarkerStyle(20);
        hData2Dplot->SetMarkerSize(1);
        hData2Dplot->Draw();
        for(int ivar=0; ivar<NTrials; ivar++){//changes
            folded[ivar]->SetLineColor(colortable[ivar]);
            folded[ivar]->ProjectionX("",binjet,binjet,"E")->Draw("same");
        }
    }
    TLegend *lFolded2 = new TLegend(0.25,0.25, 0.85, 0.85);
    cFolded->cd(fUnfoldedBayes[regBayes-1]->GetNbinsY()+1);
    for(int ivar=0; ivar<NTrials; ivar++){//changes
        TH1D* hRatio = (TH1D*)folded[ivar]->ProjectionX("",1,1,"E");
        hRatio->SetLineColor(colortable[ivar]);
        lFolded2->AddEntry(hRatio,Form("Bayes iter = %d",ivar+1),"l");
    }
    lFolded2->Draw();
 
    // Folded/Measured spectra comparison ratio
    // ----------------------------------------
    TLine *lm = new TLine(0.4, 1, 1.02, 1);
    lm->SetLineStyle(2);

    TCanvas *cFoldedR = new TCanvas("Folded/Measured Ratio","", 1400, 900);
    cFoldedR->Divide(3,2);
    TH1D* hRatio[fUnfoldedBayes[regBayes-1]->GetNbinsY()][NTrials];
    for (int binjet=1; binjet<= fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
        //TLegend *lFoldedR = new TLegend(0.65,0.65, 0.85, 0.85);
        /***********************************
        ############# folding ##################
        ************************************/
        cFoldedR->cd(binjet);
        TH1D* hDbase = (TH1D*)hData2D->ProjectionX("",binjet,binjet,"E");
        for(int ivar=0; ivar<NTrials; ivar++){//changes
            hRatio[binjet-1][ivar] = (TH1D*)folded[ivar]->ProjectionX("",binjet,binjet,"E");
            hRatio[binjet-1][ivar]->Divide(hDbase);
            hRatio[binjet-1][ivar]->SetLineColor(colortable[ivar]);
            for (int binN=1; binN <=hRatio[binjet-1][ivar]->GetNbinsX(); binN++){
                hRatio[binjet-1][ivar]->SetBinError(binN,0);
            }
            if(!ivar)hRatio[binjet-1][ivar]->Draw();
            else hRatio[binjet-1][ivar]->Draw("same");
            hRatio[binjet-1][ivar]->GetYaxis()->SetRangeUser(0.92,1.08);
            //lFoldedR->AddEntry(hRatio,Form("Bayes iter = %d",ivar+1),"l");
        }
        //lFoldedR->Draw("same");
        lm->Draw("same");
    }
    TLegend *lFoldedR2 = new TLegend(0.25,0.25, 0.85, 0.85);
    cFoldedR->cd(fUnfoldedBayes[regBayes-1]->GetNbinsY()+1);
    for(int ivar=0; ivar<NTrials; ivar++){//changes
        TH1D* hRatio = (TH1D*)folded[ivar]->ProjectionX("",1,1,"E");
        hRatio->SetLineColor(colortable[ivar]);
        lFoldedR2->AddEntry(hRatio,Form("Bayes iter = %d",ivar+1),"l");
    }
    lFoldedR2->Draw();
    //
    //}
    
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
