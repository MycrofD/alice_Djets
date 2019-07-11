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
const int fptbinsZGenN=10;
double fptbinsZGenA[fptbinsZGenN+1]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
bool BayesTrueRangeSys=1;
const int fJetptbinsGenN=9; TString truebinnum="";
double fJetptbinsGenA[fJetptbinsGenN+1]={0,2,3,4,5,7,10,15,50,60};
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
int dim[nDim]   = {zRec,jetRec,zGen,jetGen};//for extacting 6D info from THnSparse
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
    TH2D* hZjetGRebin = new TH2D("hGenRebin","hGenRebin", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsN, fJetptbinsA);
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
        //sparseMC[i]->GetAxis(5)->SetRangeUser(fptbinsZGenA[0],fptbinsZGenA[fptbinsZGenN]);
        sparseMC[i]->GetAxis(1)->SetRangeUser(fJetptbinsGenA[0],fJetptbinsGenA[fJetptbinsGenN]);
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

   // /**************************
   // #### was for MC closure. To show how gen level looks like
   // **************************/
   // TH2D *hZjetGen = (TH2D*)hZjetRecGen->Projection(3,2,"E");
   // //TH2D *hZjetGen = (TH2D*)hZjetRecGen->Projection(1,0,"E");
   // TH2D *hZjetGenRebin = Rebin2D("hZjetGenRebin", hZjetGen, fptbinsZGenN, fptbinsZGenA, fJetptbinsGenN, fJetptbinsGenA, 0);
   // TCanvas* cZGen = new TCanvas();
   // cZGen->SetLogz();
   // //hZjetGenRebin->Draw("colz");
   // //hZjetGenRebin->Draw("TEXT SAME");
   // //WeightMatrixY1Y2(hZjetRecGen,hZjetGen);

    /***************************
    #### unfolding settings ####
    ***************************/
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    RooUnfoldResponse response (hZjetRRebin, hZjetGRebin);
    //----------- fill 4D histo response matrix
    //int eventcount = 0;
    double binval = 0;

    TH2D *dResZ1  = new TH2D("dRes1" ,"" ,10,0,1,10,0,1 );
    dResZ1->GetYaxis()->SetTitle("Reco: jet pt 5-7 GeV");
 //   dResZ1->SetTitleSize(0.4,"Y");
    TH2D *dResZ2  = new TH2D("dRes2" ,"" ,10,0,1,10,0,1 );
    TH2D *dResZ3  = new TH2D("dRes3" ,"" ,10,0,1,10,0,1 );
    TH2D *dResZ4  = new TH2D("dRes4" ,"" ,10,0,1,10,0,1 );
    dResZ4->GetXaxis()->SetTitle("z gen");
    dResZ4->GetYaxis()->SetTitle("z reco");
    TH2D *dResZ5  = new TH2D("dRes5" ,"" ,10,0,1,10,0,1 );
    dResZ5->GetYaxis()->SetTitle("jet pt 7-10 GeV");
    TH2D *dResZ6  = new TH2D("dRes6" ,"" ,10,0,1,10,0,1 );
    TH2D *dResZ7  = new TH2D("dRes7" ,"" ,10,0,1,10,0,1 );
    TH2D *dResZ8  = new TH2D("dRes8" ,"" ,10,0,1,10,0,1 );
    TH2D *dResZ9  = new TH2D("dRes9" ,"" ,10,0,1,10,0,1 );
    dResZ9->GetYaxis()->SetTitle("jet pt 10-15 GeV");
    TH2D *dResZ10 = new TH2D("dRes10","" ,10,0,1,10,0,1 );
    TH2D *dResZ11 = new TH2D("dRes11","" ,10,0,1,10,0,1 );
    TH2D *dResZ12 = new TH2D("dRes12","" ,10,0,1,10,0,1 );
    TH2D *dResZ13 = new TH2D("dRes13","Gen: jet pt 5-7  GeV",10,0,1,10,0,1 );
    dResZ13->GetYaxis()->SetTitle("jet pt 15-50 GeV");
    TH2D *dResZ14 = new TH2D("dRes14","jet pt 7-10 GeV",10,0,1,10,0,1 );
    TH2D *dResZ15 = new TH2D("dRes15","jet pt 10-15 GeV",10,0,1,10,0,1 );
    TH2D *dResZ16 = new TH2D("dRes16","jet pt 15-50 GeV",10,0,1,10,0,1 );

    int respcount = 0;
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

            bool measurement_ok = kTRUE;
            //if (i==0 || j ==0 || k ==0 || m == 0){
            // if (i_center<0 || j_center<0){
            //     measurement_ok = kFALSE;
            //     cout<<i<<"--"<<j<<endl;//"--"<<k<<"--"<<m<<endl;
            //     //eventcount+=1;
            // }
            if(DR_center<2 || DG_center<2){
                  measurement_ok = kFALSE;
            }
            else if((DR_center<3 || DG_center<3) && jG_center>=7 ){
                  measurement_ok = kFALSE;
            }
            else if((DR_center<5 || DG_center<5) && jG_center>=10 ){
                  measurement_ok = kFALSE;
            }
            if (measurement_ok){
                response.Fill(zR_center,jR_center,zG_center,jG_center,weight);
		respcount+=1;
		//--drawing response matrices.
		//--we need 4x4 squares for each jetptxjetpt intervals. 
		//--each square should have 10x10 squares for each zxz bins.
		if(jR_center>fJetptbinsA[0] && jR_center<=fJetptbinsA[1]){
		    if(jG_center>fJetptbinsA[0] && jG_center<=fJetptbinsA[1]){
		    	binval = dResZ1->GetBinContent(dResZ1->GetXaxis()->FindBin(zG_center), dResZ1->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ1->SetBinContent(dResZ1->GetXaxis()->FindBin(zR_center),dResZ1->GetYaxis()->FindBin(zG_center),binval);
		    }
		    else if(jG_center>fJetptbinsA[1] && jG_center<=fJetptbinsA[2]){
		    	binval = dResZ2->GetBinContent(dResZ2->GetXaxis()->FindBin(zG_center), dResZ2->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ2->SetBinContent(dResZ2->GetXaxis()->FindBin(zR_center),dResZ2->GetYaxis()->FindBin(zG_center),binval);
		    }
		    else if(jG_center>fJetptbinsA[2] && jG_center<=fJetptbinsA[3]){
		    	binval = dResZ3->GetBinContent(dResZ3->GetXaxis()->FindBin(zG_center), dResZ3->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ3->SetBinContent(dResZ3->GetXaxis()->FindBin(zR_center),dResZ3->GetYaxis()->FindBin(zG_center),binval);
		    }
		    else if(jG_center>fJetptbinsA[3] && jG_center<=fJetptbinsA[4]){
		    	binval = dResZ4->GetBinContent(dResZ4->GetXaxis()->FindBin(zG_center), dResZ4->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ4->SetBinContent(dResZ4->GetXaxis()->FindBin(zR_center),dResZ4->GetYaxis()->FindBin(zG_center),binval);
		    }

		}
		else if(jR_center>fJetptbinsA[1] && jR_center<=fJetptbinsA[2]){
		    if(jG_center>fJetptbinsA[0] && jG_center<=fJetptbinsA[1]){
		    	binval = dResZ5->GetBinContent(dResZ5->GetXaxis()->FindBin(zG_center), dResZ5->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ5->SetBinContent(dResZ5->GetXaxis()->FindBin(zR_center),dResZ5->GetYaxis()->FindBin(zG_center),binval);
		    }
		    else if(jG_center>fJetptbinsA[1] && jG_center<=fJetptbinsA[2]){
		    	binval = dResZ6->GetBinContent(dResZ6->GetXaxis()->FindBin(zG_center), dResZ6->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ6->SetBinContent(dResZ6->GetXaxis()->FindBin(zR_center),dResZ6->GetYaxis()->FindBin(zG_center),binval);
		    }
		    else if(jG_center>fJetptbinsA[2] && jG_center<=fJetptbinsA[3]){
		    	binval = dResZ7->GetBinContent(dResZ7->GetXaxis()->FindBin(zG_center), dResZ7->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ7->SetBinContent(dResZ7->GetXaxis()->FindBin(zR_center),dResZ7->GetYaxis()->FindBin(zG_center),binval);
		    }
		    else if(jG_center>fJetptbinsA[3] && jG_center<=fJetptbinsA[4]){
		    	binval = dResZ8->GetBinContent(dResZ8->GetXaxis()->FindBin(zG_center), dResZ8->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ8->SetBinContent(dResZ8->GetXaxis()->FindBin(zR_center),dResZ8->GetYaxis()->FindBin(zG_center),binval);
		    }

		}
		else if(jR_center>fJetptbinsA[2] && jR_center<=fJetptbinsA[3]){
		    if(jG_center>fJetptbinsA[0] && jG_center<=fJetptbinsA[1]){
		    	binval = dResZ9->GetBinContent(dResZ9->GetXaxis()->FindBin(zG_center), dResZ9->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ9->SetBinContent(dResZ9->GetXaxis()->FindBin(zR_center),dResZ9->GetYaxis()->FindBin(zG_center),binval);
		    }
		    else if(jG_center>fJetptbinsA[1] && jG_center<=fJetptbinsA[2]){
		    	binval = dResZ10->GetBinContent(dResZ10->GetXaxis()->FindBin(zG_center), dResZ10->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ10->SetBinContent(dResZ10->GetXaxis()->FindBin(zR_center),dResZ10->GetYaxis()->FindBin(zG_center),binval);
		    }
		    else if(jG_center>fJetptbinsA[2] && jG_center<=fJetptbinsA[3]){
		    	binval = dResZ11->GetBinContent(dResZ11->GetXaxis()->FindBin(zG_center), dResZ11->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ11->SetBinContent(dResZ11->GetXaxis()->FindBin(zR_center),dResZ11->GetYaxis()->FindBin(zG_center),binval);
		    }
		    else if(jG_center>fJetptbinsA[3] && jG_center<=fJetptbinsA[4]){
		    	binval = dResZ12->GetBinContent(dResZ12->GetXaxis()->FindBin(zG_center), dResZ12->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ12->SetBinContent(dResZ12->GetXaxis()->FindBin(zR_center),dResZ12->GetYaxis()->FindBin(zG_center),binval);
		    }

		}
		else if(jR_center>fJetptbinsA[3] && jR_center<=fJetptbinsA[4]){
		    if(jG_center>fJetptbinsA[0] && jG_center<=fJetptbinsA[1]){
		    	binval = dResZ13->GetBinContent(dResZ13->GetXaxis()->FindBin(zG_center), dResZ13->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ13->SetBinContent(dResZ13->GetXaxis()->FindBin(zR_center),dResZ13->GetYaxis()->FindBin(zG_center),binval);
		    }
		    else if(jG_center>fJetptbinsA[1] && jG_center<=fJetptbinsA[2]){
		    	binval = dResZ14->GetBinContent(dResZ14->GetXaxis()->FindBin(zG_center), dResZ14->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ14->SetBinContent(dResZ14->GetXaxis()->FindBin(zR_center),dResZ14->GetYaxis()->FindBin(zG_center),binval);
		    }
		    else if(jG_center>fJetptbinsA[2] && jG_center<=fJetptbinsA[3]){
		    	binval = dResZ15->GetBinContent(dResZ15->GetXaxis()->FindBin(zG_center), dResZ15->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ15->SetBinContent(dResZ15->GetXaxis()->FindBin(zR_center),dResZ15->GetYaxis()->FindBin(zG_center),binval);
		    }
		    else if(jG_center>fJetptbinsA[3] && jG_center<=fJetptbinsA[4]){
		    	binval = dResZ16->GetBinContent(dResZ16->GetXaxis()->FindBin(zG_center), dResZ16->GetXaxis()->FindBin(zR_center));
		    	binval+=1;
			display +=1;
		    	dResZ16->SetBinContent(dResZ16->GetXaxis()->FindBin(zR_center),dResZ16->GetYaxis()->FindBin(zG_center),binval);
		    }
	        }
            }
            //else{
            //    response.Miss(k_center,m_center,weight);
            //}//cout<<i<<"-"<<j<<"-"<<k<<"-"<<m<<endl;
        }
    }
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
		
	cout<<respcount<<"------"<<display<<endl;
        
        //double total = dResZ1->Integral() + dResZ2->Integral() + dResZ3->Integral() + dResZ4->Integral()
	//       	+ dResZ5->Integral() + dResZ6->Integral() + dResZ7->Integral() + dResZ8->Integral() 
	//	+ dResZ9->Integral() + dResZ10->Integral() + dResZ11->Integral() + dResZ12->Integral() 
	//	+ dResZ13->Integral() + dResZ14->Integral() + dResZ15->Integral() + dResZ16->Integral();
	//
	//for (int ir = 1; ir<=10; ir++){
	//    for (int ic = 1; ic<= 10; ic++){
	//	dResZ1 ->SetBinContent(ir,ic,dResZ1 ->GetBinContent(ir,ic)/total);
	//	dResZ2 ->SetBinContent(ir,ic,dResZ2 ->GetBinContent(ir,ic)/total);
	//	dResZ3 ->SetBinContent(ir,ic,dResZ3 ->GetBinContent(ir,ic)/total);
	//	dResZ4 ->SetBinContent(ir,ic,dResZ4 ->GetBinContent(ir,ic)/total);
	//	dResZ5 ->SetBinContent(ir,ic,dResZ5 ->GetBinContent(ir,ic)/total);
	//	dResZ6 ->SetBinContent(ir,ic,dResZ6 ->GetBinContent(ir,ic)/total);
	//	dResZ7 ->SetBinContent(ir,ic,dResZ7 ->GetBinContent(ir,ic)/total);
	//	dResZ8 ->SetBinContent(ir,ic,dResZ8 ->GetBinContent(ir,ic)/total);
	//	dResZ9 ->SetBinContent(ir,ic,dResZ9 ->GetBinContent(ir,ic)/total);
	//	dResZ10->SetBinContent(ir,ic,dResZ10->GetBinContent(ir,ic)/total);
	//	dResZ11->SetBinContent(ir,ic,dResZ11->GetBinContent(ir,ic)/total);
	//	dResZ12->SetBinContent(ir,ic,dResZ12->GetBinContent(ir,ic)/total);
	//	dResZ13->SetBinContent(ir,ic,dResZ13->GetBinContent(ir,ic)/total);
	//	dResZ14->SetBinContent(ir,ic,dResZ14->GetBinContent(ir,ic)/total);
	//	dResZ15->SetBinContent(ir,ic,dResZ15->GetBinContent(ir,ic)/total);
	//	dResZ16->SetBinContent(ir,ic,dResZ16->GetBinContent(ir,ic)/total);
	//    }
	//	cout<<"================+"<<endl;
	//}
	//
	
	//dResZ1->Scale(1.0/dResZ1->Integral());
	//dResZ2->Scale(1.0/dResZ2->Integral());
	//dResZ3->Scale(1.0/dResZ3->Integral());
	//dResZ4->Scale(1.0/dResZ4->Integral());
	//dResZ5->Scale(1.0/dResZ5->Integral());
	//dResZ6->Scale(1.0/dResZ6->Integral());
	//dResZ7->Scale(1.0/dResZ7->Integral());
	//dResZ8->Scale(1.0/dResZ8->Integral());
	//dResZ9->Scale(1.0/dResZ9->Integral());
	//dResZ10->Scale(1.0/dResZ10->Integral());
	//dResZ11->Scale(1.0/dResZ11->Integral());
	//dResZ12->Scale(1.0/dResZ12->Integral());
	//dResZ13->Scale(1.0/dResZ13->Integral());
	//dResZ14->Scale(1.0/dResZ14->Integral());
	//dResZ15->Scale(1.0/dResZ15->Integral());
	//dResZ16->Scale(1.0/dResZ16->Integral());
	
	//h->GetZaxis()->SetRangeUser(0.9, 1.1)

	TCanvas *cResMat = new TCanvas("ResponseMatrix","ResponseMatrix",800,800);
	cResMat->Divide(4,4);
	/*
        double zmax = 1;
        dResZ1 ->GetZaxis()->SetRangeUser(0,zmax);
        dResZ2 ->GetZaxis()->SetRangeUser(0,zmax);
        dResZ3 ->GetZaxis()->SetRangeUser(0,zmax);
        dResZ4 ->GetZaxis()->SetRangeUser(0,zmax);
        dResZ5 ->GetZaxis()->SetRangeUser(0,zmax);
        dResZ6 ->GetZaxis()->SetRangeUser(0,zmax);
        dResZ7 ->GetZaxis()->SetRangeUser(0,zmax);
        dResZ8 ->GetZaxis()->SetRangeUser(0,zmax);
        dResZ9 ->GetZaxis()->SetRangeUser(0,zmax);
        dResZ10->GetZaxis()->SetRangeUser(0,zmax);
        dResZ11->GetZaxis()->SetRangeUser(0,zmax);
        dResZ12->GetZaxis()->SetRangeUser(0,zmax);
        dResZ13->GetZaxis()->SetRangeUser(0,zmax);
        dResZ14->GetZaxis()->SetRangeUser(0,zmax);
        dResZ15->GetZaxis()->SetRangeUser(0,zmax);
        dResZ16->GetZaxis()->SetRangeUser(0,zmax);
*/
	cResMat->cd(13);dResZ1 ->Draw("colz");
	cResMat->cd(14);dResZ2 ->Draw("colz");
	cResMat->cd(15);dResZ3 ->Draw("colz");
	cResMat->cd(16);dResZ4 ->Draw("colz");
	cResMat->cd(9) ;dResZ5 ->Draw("colz");
	cResMat->cd(10);dResZ6 ->Draw("colz");
	cResMat->cd(11);dResZ7 ->Draw("colz");
	cResMat->cd(12);dResZ8 ->Draw("colz");
	cResMat->cd(5) ;dResZ9 ->Draw("colz");
	cResMat->cd(6) ;dResZ10->Draw("colz");
	cResMat->cd(7) ;dResZ11->Draw("colz");
	cResMat->cd(8) ;dResZ12->Draw("colz");
	cResMat->cd(1) ;dResZ13->Draw("colz");
	cResMat->cd(2) ;dResZ14->Draw("colz");
	cResMat->cd(3) ;dResZ15->Draw("colz");
	cResMat->cd(4) ;dResZ16->Draw("colz");


	TCanvas *cResMat2 = new TCanvas("ResponseMatrix2","ResponseMatrix2",800,800);
	cResMat2->Divide(4,4);
	cResMat2->cd(13);dResZ1 ->Draw("");
	cResMat2->cd(14);dResZ2 ->Draw("");
	cResMat2->cd(15);dResZ3 ->Draw("");
	cResMat2->cd(16);dResZ4 ->Draw("");
	cResMat2->cd(9) ;dResZ5 ->Draw("");
	cResMat2->cd(10);dResZ6 ->Draw("");
	cResMat2->cd(11);dResZ7 ->Draw("");
	cResMat2->cd(12);dResZ8 ->Draw("");
	cResMat2->cd(5) ;dResZ9 ->Draw("");
	cResMat2->cd(6) ;dResZ10->Draw("");
	cResMat2->cd(7) ;dResZ11->Draw("");
	cResMat2->cd(8) ;dResZ12->Draw("");
	cResMat2->cd(1) ;dResZ13->Draw("");
	cResMat2->cd(2) ;dResZ14->Draw("");
	cResMat2->cd(3) ;dResZ15->Draw("");
	cResMat2->cd(4) ;dResZ16->Draw("");

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

    /***********************************
    // FD systematics
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

    TH2D *fUnfoldedBayesScale = (TH2D*)fUnfoldedBayes[regBayes-1]->Clone();
    TH1D *UnfProjectX[4];
    TH1D *hDataProjectX[4];
    TH1D *UnfProjectXScale[4];
    //************************************Unfolding systematics
    TH1D *UnfProjectXIterPre[4];
    TH1D *UnfProjectXIterPost[4];
    //************************************
    TH1D *UnfProjectXUp[4];
    TH1D *hDataProjectXUp[4];
    TH1D *UnfProjectXScaleUp[4];
    TH1D *UnfProjectXDo[4];
    TH1D *hDataProjectXDo[4];
    TH1D *UnfProjectXScaleDo[4];
    //************************************
    TFile *unfold2DoutFile = new TFile(Form("%s/alljetz2D/unfold2DoutFile%s.root",outDir.Data(),truebinnum.Data()),"recreate");
    for (int binjet=1; binjet<= fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
        UnfProjectX[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes-1]->ProjectionX(Form("UnfProjectX_%d",binjet), binjet, binjet, "E");
        hDataProjectX[binjet-1] = (TH1D*)hData2D->ProjectionX(Form("hDatarojectX_%d",binjet), binjet, binjet, "E");
        UnfProjectXScale[binjet-1] = (TH1D*)UnfProjectX[binjet-1]->Clone();
        UnfProjectXScale[binjet-1]->Scale(1,"width");
        UnfProjectX[binjet-1]     ->Write(Form("UnfProjectX_%d",binjet-1));
        hDataProjectX[binjet-1]   ->Write(Form("hDataProjectX_%d",binjet-1));
        //UnfProjectXScale[binjet-1]->Write(Form("UnfProjectXScale_%d",binjet-1));
        UnfProjectXIterPre[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes-2]->ProjectionX(Form("UnfProjectXIterPre_%d",binjet), binjet, binjet, "E");
        UnfProjectXIterPre[binjet-1]     ->Write(Form("UnfProjectXIterPre_%d",binjet-1));
        //************************************Unfolding systematics
        UnfProjectXIterPost[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes]->ProjectionX(Form("UnfProjectXIterPost_%d",binjet), binjet, binjet, "E");
        UnfProjectXIterPost[binjet-1]     ->Write(Form("UnfProjectXIterPost_%d",binjet-1));
        //************************************
        if(FDsys){
            UnfProjectXUp[binjet-1] = (TH1D*)fUnfoldedBayesUp[regBayes-1]->ProjectionX(Form("UnfProjectXUp_%d",binjet), binjet, binjet, "E");
            UnfProjectXDo[binjet-1] = (TH1D*)fUnfoldedBayesDo[regBayes-1]->ProjectionX(Form("UnfProjectXDo_%d",binjet), binjet, binjet, "E");
            hDataProjectXUp[binjet-1] = (TH1D*)hData2DUp->ProjectionX(Form("hDataProjectXUp_%d",binjet), binjet, binjet, "E");
            hDataProjectXDo[binjet-1] = (TH1D*)hData2DDo->ProjectionX(Form("hDataProjectXDo_%d",binjet), binjet, binjet, "E");
            UnfProjectXScaleUp[binjet-1] = (TH1D*)UnfProjectXUp[binjet-1]->Clone();
            UnfProjectXScaleDo[binjet-1] = (TH1D*)UnfProjectXDo[binjet-1]->Clone();
            UnfProjectXScaleUp[binjet-1]->Scale(1,"width");
            UnfProjectXScaleDo[binjet-1]->Scale(1,"width");
            UnfProjectXUp[binjet-1]     ->Write(Form("UnfProjectXUp_%d",binjet-1));
            hDataProjectXUp[binjet-1]   ->Write(Form("hDataProjectXUp_%d",binjet-1));
            UnfProjectXDo[binjet-1]     ->Write(Form("UnfProjectXDo_%d",binjet-1));
            hDataProjectXDo[binjet-1]   ->Write(Form("hDataProjectXDo_%d",binjet-1));
        }
        for(int binz=1;binz<=fUnfoldedBayesScale->GetNbinsX();binz++){
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

    if (BayesTrueRangeSys){
      TFile *unfold2DoutFileRangeSys = new TFile(Form("%s/alljetz2D/unfold2DoutFileRangeSys_%s.root",outDir.Data(),truebinnum.Data()),"recreate");
      for (int binjet=1; binjet<= fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
        UnfProjectX[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes-1]->ProjectionX(Form("UnfProjectX_%d",binjet), binjet, binjet, "E");
        UnfProjectX[binjet-1] ->Write(Form("UnfProjectX_%d",binjet-1));
      }
      unfold2DoutFileRangeSys->Close();
    }
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
