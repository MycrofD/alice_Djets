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
TString jetpttitle[5] = {"3-5","5-7", "7-10","10-15","15-50"};
// to fill the response matrix
bool BayesParameterSys=1;
const int fJetptbinsGenNN=7;
double fJetptbinsGenAA[fJetptbinsGenNN+1]={0,3,5,7,10,15,50,100};
         //const int fJetptbinsGenN=6; TString truebinnum="";
         //double fJetptbinsGenA[fJetptbinsGenN+1]={0,5,7,10,15,50,100};
         const int fJetptbinsGenN=5; TString truebinnum="250";
         double fJetptbinsGenA[fJetptbinsGenN+1]={2,5,7,10,15,50};
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
  int priorType = 1,   // if isPrior == 1, choose type of the prior prior 1 is steep, prior 2 is gentle
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
    bool FDsys = 1;
    bool DptcutTrue=1;
    //-----------------------------
    gStyle ->SetOptStat(0000);
    gSystem->Exec(Form("mkdir  %s",outDir.Data()));
    gSystem->Exec(Form("mkdir  %s/Bayes",outDir.Data()));
    outDir+="/Bayes";
    gSystem->Exec(Form("mkdir  %s",outDir.Data()));
    gSystem->Exec(Form("mkdir  %s/plots",outDir.Data()));
    gSystem->Exec(Form("mkdir  %s/alljetz2D",outDir.Data()));
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
        data2D[binjet]=Form("%s",datafile.Data());
        data2D[binjet]+=Form("/Jetbin_%d_%d",(int)fJetptbinsA[binjet], (int)fJetptbinsA[binjet+1]);
        data2D[binjet]+="/plots/JetPtSpectrum_FDsub.root";
        //reading 2D DetRM: need to create 2D RM file in ResponseMatrix folder
	//no, we are creating RM here in this macro itself
    }
    for (int binjet=0; binjet < fJetptbinsN; binjet++){
        // ---------------------------------------
        // 1D projections of FD subtracted spectra
        // ---------------------------------------
        //cout<<fJetptbinsA[binjet]<<endl;
        fRawSpec2Dproj[binjet] = (TH1D*)LoadRawSpec(data2D[binjet].Data(),"hData_binned_sub");
        TH1D *fRawSpec2DprojScale=(TH1D*)fRawSpec2Dproj[binjet]->Clone();
        fRawSpec2DprojScale->Scale(1,"width");
        if(FDsys){
            fRawSpec2DprojUp[binjet] = (TH1D*)LoadRawSpec(data2D[binjet].Data(),"hData_binned_sub_up");
            fRawSpec2DprojDo[binjet] = (TH1D*)LoadRawSpec(data2D[binjet].Data(),"hData_binned_sub_down");
            TH1D *fRawSpec2DprojUpScale=(TH1D*)fRawSpec2DprojUp[binjet]->Clone();
            TH1D *fRawSpec2DprojDoScale=(TH1D*)fRawSpec2DprojDo[binjet]->Clone();
            fRawSpec2DprojUpScale->Scale(1,"width");
            fRawSpec2DprojDoScale->Scale(1,"width");
        }
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
            if(FDsys){
                double contUp  = fRawSpec2DprojUp[binjet]->GetBinContent(binz);
                double contErrUp = fRawSpec2DprojUp[binjet]->GetBinError(binz);
                hData2DUp->SetBinContent(binz,binjet+1,contUp);
                hData2DUp->SetBinError(binz,binjet+1,contErrUp);
                double contDo  = fRawSpec2DprojDo[binjet]->GetBinContent(binz);
                double contErrDo = fRawSpec2DprojDo[binjet]->GetBinError(binz);
                hData2DDo->SetBinContent(binz,binjet+1,contDo);
                hData2DDo->SetBinError(binz,binjet+1,contErrDo);
            }
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
    cFD_2D  ->SaveAs(Form("%s/alljetz2D/FDdata2D.pdf",outDir.Data()));
    cFD_2D  ->SaveAs(Form("%s/alljetz2D/FDdata2D.png",outDir.Data()));
    cFD_2D  ->SaveAs(Form("%s/alljetz2D/FDdata2D.svg",outDir.Data()));
    TFile *outFile = new TFile(Form("%s/alljetz2D/outFD.root",outDir.Data()),"recreate");
    hData2D ->Write();//raw 2D
    hData2DS->Write();//scaled 2D
    outFile ->Close();
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
    hZjetRRebin->Sumw2();
    TH2D* hZjetGRebin = new TH2D("hGenRebin","hGenRebin", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsGenN, fJetptbinsGenA);
    hZjetGRebin->Sumw2();

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
        //sparseMC[i]->GetAxis(1)->SetRangeUser(fJetptbinsGenA[0],fJetptbinsGenA[fJetptbinsGenN]);
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

    cout<<hZjetRecGenD->GetAxis(0)->GetXmin()<<endl;
    cout<<hZjetRecGenD->GetAxis(1)->GetXmin()<<endl;
    cout<<hZjetRecGenD->GetAxis(2)->GetXmin()<<endl;
//    cout<<hZjetRecGenD->GetAxis(1)->GetXmax()<<endl;
//    return;

    /***************************
    #### unfolding settings ####
    ***************************/
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    RooUnfoldResponse response (hZjetRRebin, hZjetGRebin);
    //----------- fill 4D histo response matrix
    //int eventcount = 0;
    double binval = 0;
    //----------- fill 4D histo response matrix
    int display = 0;
    if(DptcutTrue){
        for (int z = 0; z < hZjetRecGenD->GetNbins();z++) {
            int coord[DnDim]={0,0,0,0,0,0};
            double content = hZjetRecGenD->GetBinContent(z,coord);
            //int Ddim[DnDim]   = {zRec,jetRec,zGen,jetGen,DRec,DGen};//for extacting 6D info from THnSparse
            int i = coord[0], j = coord[1], k = coord[2], m = coord[3], n = coord[4], p = coord[5];
            double weight = content;
            double zR_center = hZjetRecGenD->GetAxis(0)->GetBinCenter(i);
            double jR_center = hZjetRecGenD->GetAxis(1)->GetBinCenter(j);
            double zG_center = hZjetRecGenD->GetAxis(2)->GetBinCenter(k);
            double jG_center = hZjetRecGenD->GetAxis(3)->GetBinCenter(m);
            double DR_center = hZjetRecGenD->GetAxis(4)->GetBinCenter(n);
            double DG_center = hZjetRecGenD->GetAxis(5)->GetBinCenter(p);

            if(isPrior){
                if(priorType==1){
                weight = weight*zG_center*jG_center;}
                else if(priorType==2){
                weight = weight/(zG_center*jG_center);}
                else if(priorType==3){
                weight = weight/(zG_center);}
                else if(priorType==4){
                weight = weight*(zG_center);}
                else if(priorType==5){
                weight = weight/(jG_center);}
                else if(priorType==6){
                weight = weight*(jG_center);}
                else if(priorType==7){
                weight = weight*(zG_center/jG_center);}
                else if(priorType==8){
                weight = weight*(jG_center/zG_center);}
                else if(priorType==9){
                weight = 10*weight/(zG_center);}
                else if(priorType==10){
                weight = 0.1*weight*(zG_center);}
                else weight = weight;
            }
            bool measurement_ok = kTRUE;
            if(DR_center<2 || DG_center<2 || jR_center<2){
                  measurement_ok = kFALSE;
            }
            if((DG_center<2 && jG_center>=2 )||(DR_center<2 && jR_center>=2 )){
                  measurement_ok = kFALSE;
            }
            if((DG_center<2 && jG_center>=5 )||(DR_center<2 && jR_center>=5 )){
                  measurement_ok = kFALSE;
            }
            if((DG_center<3 && jG_center>=7 )||(DR_center<3 && jR_center>=7 )){
                  measurement_ok = kFALSE;
            }
            if((DG_center<5 && jG_center>=10 )||(DR_center<5 && jR_center>=10) ){
                  measurement_ok = kFALSE;
            }
            if((DR_center>36 || jR_center>50 || zR_center<0.4 )){
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
            //if (i==0 || j ==0 || k ==0 || m == 0){
            // if (i_center<0 || j_center<0){
            //     measurement_ok = kFALSE;
            //     cout<<i<<"--"<<j<<endl;//"--"<<k<<"--"<<m<<endl;
            //     //eventcount+=1;
            // }
            if (measurement_ok){
                response.Fill(zR_center,jR_center,zG_center,jG_center,weight);
            }
            // else{
            //     response.Miss(k_center,m_center,weight);
            // }//cout<<i<<"-"<<j<<"-"<<k<<"-"<<m<<endl;
        }
    }

		
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


    /***********************************
    // FD systematics begin
    ************************************/
    TH2D *fUnfoldedBayesUp[NTrials];
    TH2D *foldedUp[NTrials];
    TH2D *fUnfoldedBayesDo[NTrials];
    TH2D *foldedDo[NTrials];
    if(FDsys){
        for(int ivar=0; ivar<NTrials; ivar++){//changes
            RooUnfoldBayes unfoldUp (&response, hData2DUp, ivar+1);
            RooUnfoldBayes unfoldDo (&response, hData2DDo, ivar+1);
            fUnfoldedBayesUp[ivar] = (TH2D*)unfoldUp.Hreco();
            fUnfoldedBayesDo[ivar] = (TH2D*)unfoldDo.Hreco();
            foldedUp[ivar] = (TH2D*)response.ApplyToTruth(fUnfoldedBayesUp[ivar]);
            foldedDo[ivar] = (TH2D*)response.ApplyToTruth(fUnfoldedBayesDo[ivar]);
        }
    }
    /***********************************
    // FD systematics end
    ************************************/

    //*********************************
    // Kinematic efficiency begins
    //*********************************/
    TH2D* kineEffCut  = new TH2D("hkineEffCut","hkineEffCut", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsGenN, fJetptbinsGenA);
    kineEffCut->Sumw2();
    TH2D* kineEffFull  = new TH2D("hkineEffFull","hkineEffFull", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsGenN, fJetptbinsGenA);
    kineEffFull->Sumw2();
    TH2D* kineEffCutAll  = new TH2D("hkineEffCuta","hkineEffCuta", fptbinsZTrueNN, fptbinsZTrueAA, fJetptbinsGenNN, fJetptbinsGenAA);
    kineEffCutAll->Sumw2();
    TH2D* kineEffFullAll  = new TH2D("hkineEffFulla","hkineEffFulla", fptbinsZTrueNN, fptbinsZTrueAA, fJetptbinsGenNN, fJetptbinsGenAA);
    kineEffFullAll->Sumw2();
    if(DptcutTrue){
        for (int z = 0; z < hZjetRecGenD->GetNbins();z++) {
            int coord[DnDim]={0,0,0,0,0,0};
            double content = hZjetRecGenD->GetBinContent(z,coord);
            //int Ddim[DnDim]   = {zRec,jetRec,zGen,jetGen,DRec,DGen};//for extacting 6D info from THnSparse
            int i = coord[0], j = coord[1], k = coord[2], m = coord[3], n = coord[4], p = coord[5];
            double weight = content;
            double zR_center = hZjetRecGenD->GetAxis(0)->GetBinCenter(i);
            double jR_center = hZjetRecGenD->GetAxis(1)->GetBinCenter(j);
            double zG_center = hZjetRecGenD->GetAxis(2)->GetBinCenter(k);
            double jG_center = hZjetRecGenD->GetAxis(3)->GetBinCenter(m);
            double DR_center = hZjetRecGenD->GetAxis(4)->GetBinCenter(n);
            double DG_center = hZjetRecGenD->GetAxis(5)->GetBinCenter(p);

            if(isPrior){
                if(priorType==1){
                weight = weight*zG_center*jG_center;}
                else if(priorType==2){
                weight = weight/(zG_center*jG_center);}
                else if(priorType==3){
                weight = weight/(zG_center);}
                else if(priorType==4){
                weight = weight*(zG_center);}
                else if(priorType==5){
                weight = weight/(jG_center);}
                else if(priorType==6){
                weight = weight*(jG_center);}
                else if(priorType==7){
                weight = weight*(zG_center/jG_center);}
                else if(priorType==8){
                weight = weight*(jG_center/zG_center);}
                else weight = weight;
            }
            bool measurement_ok = kTRUE;
            if(DR_center<1 || DG_center<1 || jR_center<2){
            //if(DR_center<1){
            //if(DG_center<2){
                  measurement_ok = kFALSE;
            }
            //if((DG_center<2 && jG_center>=4.25 )||(DR_center<2 && jR_center>=4.25 )){
            if((DG_center<2 && jG_center>=2 )||(DR_center<2 && jR_center>=2 )){
            //if(( DG_center<3) && jG_center>=7 ){
                  measurement_ok = kFALSE;
            }
            if((DG_center<2 && jG_center>=5 )||(DR_center<2 && jR_center>=5 )){
            //if((DR_center<2 && jR_center>=4.25 )){
                  measurement_ok = kFALSE;
            }
            //if((DR_center<3 || DG_center<3) && jG_center>=7 ){
            if((DG_center<3 && jG_center>=7 )||(DR_center<3 && jR_center>=7 )){
            //if((DR_center<3 && jR_center>=7 )){
            //if(( DG_center<3) && jG_center>=7 ){
                  measurement_ok = kFALSE;
            }
            //if((DR_center<5 || DG_center<5) && jG_center>=10 ){
            if((DG_center<5 && jG_center>=10 )||(DR_center<5 && jR_center>=10) ){
            //if((DR_center<5 && jR_center>=10) ){
            //if((DG_center<5) && jG_center>=10 ){
                  measurement_ok = kFALSE;
            }
            if (measurement_ok){
                kineEffFull->Fill(zG_center,jG_center,weight);
                kineEffFullAll->Fill(zG_center,jG_center,weight);
            }
            if(jR_center<3 || jR_center>50 || zR_center<0.4){
            //if((DG_center<5) && jG_center>=10 ){
                  measurement_ok = kFALSE;
            }
            if (measurement_ok){
                kineEffCut->Fill(zG_center,jG_center,weight);
                kineEffCutAll->Fill(zG_center,jG_center,weight);
            }
            //else{
            //    response.Miss(zG_center,jG_center,weight);
            //}//cout<<i<<"-"<<j<<"-"<<k<<"-"<<m<<endl;
        }
    }
    TCanvas* cKine = new TCanvas("hKine","hKine",800,600);
    TH2D* hKine = (TH2D*) kineEffCut->Clone();
    hKine->Divide(kineEffCut,kineEffFull,1,1,"B");
    
    hKine->SetTitle("");
    hKine->Draw("colz");
    hKine->Draw("TEXT same");
//    hKine->SaveAs(Form("%s/plots/_hKine.pdf",outDir.Data()));
//    hKine->SaveAs(Form("%s/plots/_hKine.png",outDir.Data()));
//    hKine->SaveAs(Form("%s/plots/_hKine.svg",outDir.Data()));

    TCanvas* cKineA = new TCanvas("hKinea","hKinea",800,600);
    TH2D* hKineAll = (TH2D*) kineEffCutAll->Clone();
    hKineAll->Divide(kineEffFullAll);
    hKineAll->GetYaxis()->SetRangeUser(3,50);
    hKineAll->GetXaxis()->SetRangeUser(0.4,1.02);

    hKineAll->SetTitle("");
    hKineAll->Draw("colz");
    hKineAll->Draw("TEXT same");
//
    //*********************************
    // Kinematic efficiency end
    //*********************************/

    TH2D *fUnfoldedBayesScale = (TH2D*)fUnfoldedBayes[regBayes-1]->Clone();
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
    TH2D* respTruth = (TH2D*)response.Htruth()->Clone();
    TH2D* respMeas = (TH2D*)response.Hmeasured()->Clone();
    respTruth->Scale(1/response.Htruth()->Integral());
    respMeas->Scale(1/response.Hmeasured()->Integral());
    TFile *responseobject = new TFile(Form("%s/alljetz2D/responseobject%s.root",outDir.Data(),  isPrior?Form("%d",priorType):""  ),"recreate");
    //response.SetName("responseNow");
    respobj->SetName("responseNow");
    response.Write();
//    response.Hresponse()->Write();
 //   response.Hresponse()->SetName("resp");
    respobj->Write();
    respTruth->Write();
    respMeas->Write();
    responseobject->Close();
    //TCanvas *cRespon = new TCanvas("Resp","Resp",950,800);
    TCanvas *cRespon = new TCanvas("","",950,800);
    TH2D* respobjP = (TH2D*)respobj->Clone();
    //respobjP->Scale(1.0/respobj->Integral());
    respobjP->GetYaxis()->SetLabelSize(0.);
    respobjP->GetXaxis()->SetLabelSize(0.);
    cRespon->SetLogz();
    respobjP->SetTitle("");
    respobjP->GetYaxis()->SetTitle("z, jetpt, GEN");
    respobjP->GetYaxis()->SetTitleOffset(1.5);
    respobjP->GetXaxis()->SetTitle("z, jetpt, RECO");
    respobjP->GetXaxis()->SetTitleOffset(1.5);
    respobjP->Draw("colz");
    //setting axes and labels properly---------
    TGaxis *axisV1 = new TGaxis(0, 0, 0, 1, 0.0, 0.4, 1, "-L");
    axisV1->SetLabelFont(43);    axisV1->SetLabelSize(15);    axisV1->SetLabelOffset(0.03);    axisV1->Draw();
    TGaxis *axisV2 = new TGaxis(0, 1, 0, 2, 0.4, 0.6, 1, "-L");
    axisV2->SetLabelFont(43);    axisV2->SetLabelSize(15);    axisV2->SetLabelOffset(0.03);    axisV2->Draw();
    TGaxis *axisV3 = new TGaxis(0, 2, 0, 6, 0.6, 1.0, 4, "-L");
    axisV3->SetLabelFont(43);    axisV3->SetLabelSize(15);    axisV3->SetLabelOffset(0.03);    axisV3->Draw();
    TGaxis *axisV4 = new TGaxis(0, 8, 0, 12, 0.6, 1.0, 4, "-L");
    axisV4->SetLabelFont(43);    axisV4->SetLabelSize(15);    axisV4->SetLabelOffset(0.03);    axisV4->Draw();
    TGaxis *axisV5 = new TGaxis(0, 14, 0, 18, 0.6, 1.0, 4, "-L");
    axisV5->SetLabelFont(43);    axisV5->SetLabelSize(15);    axisV5->SetLabelOffset(0.03);    axisV5->Draw();
    TGaxis *axisV6 = new TGaxis(0, 20, 0, 24, 0.6, 1.0, 4, "-L");
    axisV6->SetLabelFont(43);    axisV6->SetLabelSize(15);    axisV6->SetLabelOffset(0.03);    axisV6->Draw();
    TGaxis *axisV7 = new TGaxis(0, 26, 0, 30, 0.6, 1.0, 4, "-L");
    axisV7->SetLabelFont(43);    axisV7->SetLabelSize(15);    axisV7->SetLabelOffset(0.03);    axisV7->Draw();
    TGaxis *axisV8 = new TGaxis(0, 32, 0, 36, 0.6, 1.0, 4, "-L");
    axisV8->SetLabelFont(43);    axisV8->SetLabelSize(15);    axisV8->SetLabelOffset(0.03);    axisV8->Draw();
    TGaxis *axisV9 = new TGaxis(0, 0, 0, 6, 0, 5, 1, "-L");
    axisV9->SetLabelFont(43);    axisV9->SetLabelSize(19);    axisV9->SetLabelOffset(0.05);    axisV9->Draw();
    TGaxis *axisV10 = new TGaxis(0, 12, 0, 18, 7, 10, 1, "-L");
    axisV10->SetLabelFont(43);    axisV10->SetLabelSize(19);    axisV10->SetLabelOffset(0.06);    axisV10->Draw();
    TGaxis *axisV11 = new TGaxis(0, 24, 0, 30, 15, 50, 1, "-L");
    axisV11->SetLabelFont(43);    axisV11->SetLabelSize(19);    axisV11->SetLabelOffset(0.06);    axisV11->Draw();
    TGaxis *axisV12 = new TGaxis(0, 30, 0, 36, 50, 100, 1, "-L");
    axisV12->SetLabelFont(43);    axisV12->SetLabelSize(19);    axisV12->SetLabelOffset(0.06);    axisV12->Draw();
    TGaxis *axisH1 = new TGaxis(0, 0, 6, 0, 0.4, 1.0, 3, "-+L");
    axisH1->SetLabelFont(43);    axisH1->SetLabelSize(15);    axisH1->SetLabelOffset(0.005);    axisH1->Draw();
    TGaxis *axisH2 = new TGaxis(6, 0, 12, 0, 0.4, 1.0, 3, "-+L");
    axisH2->SetLabelFont(43);    axisH2->SetLabelSize(15);    axisH2->SetLabelOffset(0.02);    axisH2->Draw();
    TGaxis *axisH3 = new TGaxis(12, 0, 18, 0, 0.4, 1.0, 3, "-+L");
    axisH3->SetLabelFont(43);    axisH3->SetLabelSize(15);    axisH3->SetLabelOffset(0.005);    axisH3->Draw();
    TGaxis *axisH4 = new TGaxis(18, 0, 24, 0, 0.4, 1.0, 3, "-+L");
    axisH4->SetLabelFont(43);    axisH4->SetLabelSize(15);    axisH4->SetLabelOffset(0.02);    axisH4->Draw();
    TGaxis *axisH5 = new TGaxis(24, 0, 30, 0, 0.4, 1.0, 3, "-+L");
    axisH5->SetLabelFont(43);    axisH5->SetLabelSize(15);    axisH5->SetLabelOffset(0.005);    axisH5->Draw();
    TGaxis *axisH6 = new TGaxis(30, 0, 36, 0, 0.4, 1.0, 3, "-+L");
    axisH6->SetLabelFont(43);    axisH6->SetLabelSize(15);    axisH6->SetLabelOffset(0.02);    axisH6->Draw();
    TGaxis *axisH7 = new TGaxis(0, 0, 6, 0, 3, 4.25, 1, "-+L");
    axisH7->SetLabelFont(43);    axisH7->SetLabelSize(19);    axisH7->SetLabelOffset(0.04);    axisH7->Draw();
    TGaxis *axisH8 = new TGaxis(12, 0, 18, 0, 5, 7, 1, "-+L");
    axisH8->SetLabelFont(43);    axisH8->SetLabelSize(19);    axisH8->SetLabelOffset(0.04);    axisH8->Draw();
    TGaxis *axisH9 = new TGaxis(24, 0, 30, 0, 10, 15, 1, "-+L");
    axisH9->SetLabelFont(43);    axisH9->SetLabelSize(19);    axisH9->SetLabelOffset(0.04);    axisH9->Draw();
    TGaxis *axisH10 = new TGaxis(30, 0, 36, 0, 15, 50, 1, "-+L");
    axisH10->SetLabelFont(43);    axisH10->SetLabelSize(19);    axisH10->SetLabelOffset(0.04);   axisH10->Draw();
    //respobjP->Draw();
    cRespon->SaveAs(Form("%s/alljetz2D/responseobj%s.pdf",outDir.Data(), isPrior?Form("%d",priorType):"" ));
    //cRespon->SaveAs(Form("%s/alljetz2D/responseobj%s.png",outDir.Data(), isPrior?Form("%d",priorType):"" ));
    cRespon->SaveAs(Form("%s/alljetz2D/responseobj%s.svg",outDir.Data(), isPrior?Form("%d",priorType):"" ));


    //************************************
    TFile *unfold2DoutFile = new TFile(Form("%s/alljetz2D/unfold2DoutFile%sAPW%s.root",outDir.Data(),truebinnum.Data(), isPrior?Form("%d",priorType):"" ),"recreate");
    for (int binjet=1; binjet<= fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
        UnfProjectX[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes-1]->ProjectionX(Form("UnfProjectX_%d",binjet-1), binjet, binjet, "E");
	    //UnfProjectX[binjet-1]->SetName(Form("UnfProjectX_%d",binjet-1));
        hDataProjectX[binjet-1] = (TH1D*)hData2D->ProjectionX(Form("hDataProjectX_%d",binjet-1), binjet, binjet, "E");
        UnfProjectXScale[binjet-1] = (TH1D*)UnfProjectX[binjet-1]->Clone();
        UnfProjectXScale[binjet-1]->Scale(1,"width");
        UnfProjectX[binjet-1]     ->Write();
        hDataProjectX[binjet-1]   ->Write();
        //UnfProjectXScale[binjet-1]->Write();
        UnfProjectXIterPre[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes-2]->ProjectionX(Form("UnfProjectXIterPre_%d",binjet-1), binjet, binjet, "E");
        UnfProjectXIterPre[binjet-1]     ->Write();
        //************************************Unfolding systematics
        UnfProjectXIterPost[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes]->ProjectionX(Form("UnfProjectXIterPost_%d",binjet-1), binjet, binjet, "E");
        UnfProjectXIterPost[binjet-1]     ->Write();
        //************************************
        if(FDsys){
            UnfProjectXUp[binjet-1] = (TH1D*)fUnfoldedBayesUp[regBayes-1]->ProjectionX(Form("UnfProjectXUp_%d",binjet-1), binjet, binjet, "E");
            UnfProjectXDo[binjet-1] = (TH1D*)fUnfoldedBayesDo[regBayes-1]->ProjectionX(Form("UnfProjectXDo_%d",binjet-1), binjet, binjet, "E");
            hDataProjectXUp[binjet-1] = (TH1D*)hData2DUp->ProjectionX(Form("hDataProjectXUp_%d",binjet-1), binjet, binjet, "E");
            hDataProjectXDo[binjet-1] = (TH1D*)hData2DDo->ProjectionX(Form("hDataProjectXDo_%d",binjet-1), binjet, binjet, "E");
            UnfProjectXScaleUp[binjet-1] = (TH1D*)UnfProjectXUp[binjet-1]->Clone();
            UnfProjectXScaleDo[binjet-1] = (TH1D*)UnfProjectXDo[binjet-1]->Clone();
            UnfProjectXScaleUp[binjet-1]->Scale(1,"width");
            UnfProjectXScaleDo[binjet-1]->Scale(1,"width");
            UnfProjectXUp[binjet-1]     ->Write();
            hDataProjectXUp[binjet-1]   ->Write();
            UnfProjectXDo[binjet-1]     ->Write();
            hDataProjectXDo[binjet-1]   ->Write();
        }
        for(int binz=1;binz<=fUnfoldedBayes[regBayes-1]->GetNbinsX();binz++){
            fUnfoldedBayesScale->SetBinContent(binz,binjet,UnfProjectXScale[binjet-1]->GetBinContent(binz));
            fUnfoldedBayesScale->SetBinError(binz,binjet,UnfProjectXScale[binjet-1]->GetBinError(binz));
            if(FDsys){
                TH2D *fUnfoldedBayesScaleUp = (TH2D*)fUnfoldedBayesUp[regBayes-1]->Clone();
                TH2D *fUnfoldedBayesScaleDo = (TH2D*)fUnfoldedBayesDo[regBayes-1]->Clone();
                fUnfoldedBayesScaleUp->SetBinContent(binz,binjet,UnfProjectXScaleUp[binjet-1]->GetBinContent(binz));
                fUnfoldedBayesScaleUp->SetBinError(binz,binjet,UnfProjectXScaleUp[binjet-1]->GetBinError(binz));
                fUnfoldedBayesScaleDo->SetBinContent(binz,binjet,UnfProjectXScaleDo[binjet-1]->GetBinContent(binz));
                fUnfoldedBayesScaleDo->SetBinError(binz,binjet,UnfProjectXScaleDo[binjet-1]->GetBinError(binz));
            }
        }
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
    cUnfoldedScale->SaveAs(Form("%s/plots/%s_unfSpectraScale.pdf",outDir.Data(),outName.Data()));
    cUnfoldedScale->SaveAs(Form("%s/plots/%s_unfSpectraScale.png",outDir.Data(),outName.Data()));
    cUnfoldedScale->SaveAs(Form("%s/plots/%s_unfSpectraScale.svg",outDir.Data(),outName.Data()));

    // Bayesian Parameter Systematics
    // ------------------------------
    if (BayesParameterSys){
      TFile *unfold2DoutFileRangeSys = new TFile(Form("%s/alljetz2D/unfold2DoutFileRegBayesSys%s.root",outDir.Data(), Form("%d",regBayes) ),"recreate");
      for (int binjet=1; binjet<= fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
        UnfProjectX[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes-1]->ProjectionX(Form("UnfProjectX_%d",binjet-1), binjet, binjet, "E");
        UnfProjectX[binjet-1] ->Write();
      }
      unfold2DoutFileRangeSys->Close();

      TFile *unfold2DoutFileRangeSys1 = new TFile(Form("%s/alljetz2D/unfold2DoutFileRegBayesSys%s.root",outDir.Data(), Form("%d",regBayes-1) ),"recreate");
      for (int binjet=1; binjet<= fUnfoldedBayes[regBayes-2]->GetNbinsY(); binjet++){
        UnfProjectX[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes-2]->ProjectionX(Form("UnfProjectX_%d",binjet-1), binjet, binjet, "E");
        UnfProjectX[binjet-1] ->Write();
      }
      unfold2DoutFileRangeSys1->Close();

      TFile *unfold2DoutFileRangeSys2 = new TFile(Form("%s/alljetz2D/unfold2DoutFileRegBayesSys%s.root",outDir.Data(), Form("%d",regBayes+1) ),"recreate");
      for (int binjet=1; binjet<= fUnfoldedBayes[regBayes+1]->GetNbinsY(); binjet++){
        UnfProjectX[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes+1]->ProjectionX(Form("UnfProjectX_%d",binjet-1), binjet, binjet, "E");
        UnfProjectX[binjet-1] ->Write();
      }
      unfold2DoutFileRangeSys2->Close();
    }


    // Unfolded/Measured/Folded spectra comparison
    // ----------------------------------
    int databinnum = 2; int unfoldbinnum = 1; //if(){}
    //
    TH2D* fUnfoldPlot = (TH2D*)fUnfoldedBayes[regBayes-1]->Clone();
    fUnfoldPlot->SetLineColor(kRed+2);
    fUnfoldPlot->SetLineWidth(4);
    TCanvas *cUnFoldedFold = new TCanvas("UnFoldedMeasuredFolded","UnFolded,Measured,Folded", 1400, 900);
    TLegend *lUnFoldedFold = new TLegend(0.25,0.25, 0.85, 0.85);
    cUnFoldedFold->Divide(3,2);
    cUnFoldedFold->SetLogz();
    for (int binjet=1; binjet< fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
        /***********************************
        ############# Unfolding ##################
        ************************************/
        cUnFoldedFold->cd(binjet);

        TH1D* hData2Dplot = (TH1D*)hData2D->ProjectionX("",binjet+2,binjet+2,"E")->Clone();
        //hData2Dplot->SetLineColor(kBlack);
        hData2Dplot->SetLineWidth(2);
        hData2Dplot->SetMarkerStyle(20);
        hData2Dplot->SetMarkerSize(1);
        hData2Dplot->Scale(1,"width");
        hData2Dplot->Draw();
        if(binjet<5){
        hData2Dplot->SetTitle(Form("Jet pt: %s GeV/c",jetpttitle[binjet].Data()));}
        else{
        hData2Dplot->SetTitle(Form("%s",jetpttitle[binjet].Data()));}
        //fUnfoldedBayes[regBayes-1]->SetLineColor(kRed+2);
        //fUnfoldedBayes[regBayes-1]->ProjectionX("",binjet+2,binjet+2,"E")->Draw("same");
    TH1D* h2 = (TH1D*)fUnfoldPlot->ProjectionX("",binjet+1,binjet+1,"E")->Clone();
    h2->SetLineColor(kRed+2);
    h2->SetLineWidth(2);
    h2->SetMarkerStyle(20);
    h2->SetMarkerColor(kRed+2);
        h2->Scale(1,"width");
        h2->Draw("same");
        //fUnfoldPlot->ProjectionX("",binjet+2,binjet+2,"E")->Draw("same");
        TH1D* hfold = (TH1D*)folded[regBayes-1]->ProjectionX("",binjet+2,binjet+2,"E")->Clone();
        hfold->SetLineColor(kGreen+2);
        hfold->Scale(1,"width");
        hfold->Draw("same");
        //kine efficiency
        TH1D* hKineUnf = (TH1D*)hKine->ProjectionX("",binjet+1,binjet+1)->Clone();
        TH1D* hUnfoldKine = (TH1D*)fUnfoldPlot->ProjectionX("",binjet+1,binjet+1)->Clone();
        hUnfoldKine->Divide(hKineUnf);
        hUnfoldKine->SetLineColor(kOrange+1);
        hUnfoldKine->SetMarkerColor(kOrange+1);
        hUnfoldKine->SetMarkerStyle(21);
        hUnfoldKine->Scale(1,"width");
        hUnfoldKine->Draw("same");

        hData2Dplot->SetMaximum(hData2Dplot->GetMaximum()*1.3);
        hData2Dplot->SetMinimum(hData2Dplot->GetMinimum()*0.5);
    }
    cUnFoldedFold->cd(fUnfoldedBayes[regBayes-1]->GetNbinsY()+1);
        TH1D* hData2Dplot = (TH1D*)hData2D->ProjectionX("",3,3,"E")->Clone();
        hData2Dplot->SetLineWidth(2);
        TH1D* hfoldD = (TH1D*)folded[regBayes-1]->ProjectionX("",3,3,"E")->Clone();
        hfoldD->SetLineColor(kGreen+2);
        hfoldD->SetLineWidth(2);
        lUnFoldedFold->AddEntry(hData2Dplot,"Measured","l");
        lUnFoldedFold->AddEntry(fUnfoldPlot,Form("Unfolded, iter=%d",regBayes),"l");
        lUnFoldedFold->AddEntry(hfoldD,"Folded","l");
        
        TH1D* hKineUnf1 = (TH1D*)hKine->ProjectionX("",2,2)->Clone();
        TH1D* hUnfoldKine1 = (TH1D*)fUnfoldPlot->ProjectionX("",2,2)->Clone();
        hUnfoldKine1->Divide(hKineUnf1);
        hUnfoldKine1->SetLineColor(kOrange+1);
        hUnfoldKine1->SetMarkerColor(kOrange+1);
        hUnfoldKine1->SetMarkerStyle(21);
        lUnFoldedFold->AddEntry(hUnfoldKine1,"Kine Eff corr");
    lUnFoldedFold->Draw("same");
    cUnFoldedFold->SaveAs(Form("%s/alljetz2D/unfoldMeasFold.pdf",outDir.Data()));
    cUnFoldedFold->SaveAs(Form("%s/alljetz2D/unfoldMeasFold.png",outDir.Data()));
    cUnFoldedFold->SaveAs(Form("%s/alljetz2D/unfoldMeasFold.svg",outDir.Data()));
    // Unfolded/Measured spectra comparison
    // ----------------------------------
    TCanvas *cUnFolded = new TCanvas("UnFoldedMeasured","UnFolded/Measured", 1400, 900);
    cUnFolded->Divide(3,2);
    cUnFolded->SetLogz();
    for (int binjet=1; binjet< fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
        /***********************************
        ############# Unfolding ##################
        ************************************/
        cUnFolded->cd(binjet);

        TH1D* hData2Dplot = (TH1D*)hData2D->ProjectionX("",binjet+databinnum,binjet+databinnum,"E")->Clone();
        //hData2Dplot->SetLineColor(kBlack);
        hData2Dplot->SetLineWidth(2);
        hData2Dplot->SetMarkerStyle(20);
        hData2Dplot->SetMarkerSize(1);
        hData2Dplot->Draw();
        hData2Dplot->SetTitle(Form("Jet pt: %s GeV/c",jetpttitle[binjet].Data()));
        for(int ivar=0; ivar<NTrials; ivar++){//changes
            fUnfoldedBayes[ivar]->SetLineColor(colortable[ivar]);
            fUnfoldedBayes[ivar]->ProjectionX("",binjet+unfoldbinnum,binjet+unfoldbinnum,"E")->Draw("same");
        }
    }
    TLegend *lUnFolded = new TLegend(0.25,0.25, 0.85, 0.85);
    cUnFolded->cd(fUnfoldedBayes[regBayes-1]->GetNbinsY()+1);
    for(int ivar=0; ivar<NTrials; ivar++){//changes
        TH1D* hRatio = (TH1D*)fUnfoldedBayes[ivar]->ProjectionX("",1,1,"E");
        hRatio->SetLineColor(colortable[ivar]);
        lUnFolded->AddEntry(hRatio,Form("Bayes iter = %d",ivar+1),"l");
    }
    lUnFolded->Draw();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
databinnum += -2; //delete
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Folded/Measured spectra comparison
    // ----------------------------------
    TCanvas *cFolded = new TCanvas("FoldedMeasured","Folded/Measured", 1400, 900);
    cFolded->Divide(3,2);
    //for (int binjet=1; binjet<= fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
    for (int binjet=1; binjet<= 6; binjet++){
        /***********************************
        ############# folding ##################
        ************************************/
        cFolded->cd(binjet);

        TH1D* hData2Dplot = (TH1D*)hData2D->ProjectionX("",binjet+databinnum,binjet+databinnum,"E")->Clone();
        //hData2Dplot->SetLineColor(kBlack);
        hData2Dplot->SetLineWidth(2);
        hData2Dplot->SetMarkerStyle(20);
        hData2Dplot->SetMarkerSize(1);
        hData2Dplot->Draw();
//        hData2Dplot->SetTitle(Form("Jet pt: %s GeV/c",jetpttitle[binjet].Data()));
        for(int ivar=0; ivar<NTrials; ivar++){//changes
            folded[ivar]->SetLineColor(colortable[ivar]);
            folded[ivar]->ProjectionX("",binjet+databinnum,binjet+databinnum,"E")->Draw("same");
        }
    }
    TLegend *lFolded2 = new TLegend(0.25,0.25, 0.85, 0.85);
    cFolded->cd(fUnfoldedBayes[regBayes-1]->GetNbinsY());
    for(int ivar=0; ivar<NTrials; ivar++){//changes
        TH1D* hRatio = (TH1D*)folded[ivar]->ProjectionX("",1,1,"E");
        hRatio->SetLineColor(colortable[ivar]);
//        lFolded2->AddEntry(hRatio,Form("Bayes iter = %d",ivar+1),"l");
    }
//    lFolded2->Draw();

    // Folded/Measured spectra comparison ratio
    // ----------------------------------------
    TLine *lm = new TLine(0.4, 1, 1.02, 1);
    lm->SetLineStyle(2);

    TCanvas *cFoldedR = new TCanvas("FoldedMeasuredRatio","Folded/Measured Ratio", 1400, 900);
    cFoldedR->Divide(3,2);
    TH1D* hRatio[fUnfoldedBayes[regBayes-1]->GetNbinsY()][NTrials];
    for (int binjet=1; binjet<= fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
        //TLegend *lFoldedR = new TLegend(0.65,0.65, 0.85, 0.85);
        /***********************************
        ############# folding ##################
        ************************************/
        cFoldedR->cd(binjet);
        TH1D* hDbase = (TH1D*)hData2D->ProjectionX("",binjet+databinnum,binjet+databinnum,"E");
        for(int ivar=0; ivar<NTrials; ivar++){//changes
            hRatio[binjet-1][ivar] = (TH1D*)folded[ivar]->ProjectionX("",binjet+databinnum,binjet+databinnum,"E");
            hRatio[binjet-1][ivar]->Divide(hDbase);
            hRatio[binjet-1][ivar]->SetLineColor(colortable[ivar]);
            for (int binN=1; binN <=hRatio[binjet-1][ivar]->GetNbinsX(); binN++){
                hRatio[binjet-1][ivar]->SetBinError(binN,0);
            }
//           hRatio[binjet-1][ivar]->SetTitle(Form("Jet pt: %s GeV/c",jetpttitle[binjet].Data()));
            if(!ivar)hRatio[binjet-1][ivar]->Draw();
            else hRatio[binjet-1][ivar]->Draw("same");
            hRatio[binjet-1][ivar]->GetYaxis()->SetRangeUser(0.92,1.08);
            //lFoldedR->AddEntry(hRatio,Form("Bayes iter = %d",ivar+1),"l");
        }
        //lFoldedR->Draw("same");
        lm->Draw("same");
    }
    TLegend *lFoldedR2 = new TLegend(0.25,0.25, 0.85, 0.85);
    cFoldedR->cd(fUnfoldedBayes[regBayes-1]->GetNbinsY());
    for(int ivar=0; ivar<NTrials; ivar++){//changes
        TH1D* hRatio = (TH1D*)folded[ivar]->ProjectionX("",1,1,"E");
        hRatio->SetLineColor(colortable[ivar]);
        lFoldedR2->AddEntry(hRatio,Form("Bayes iter = %d",ivar+1),"l");
    }
//    lFoldedR2->Draw();
    
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
