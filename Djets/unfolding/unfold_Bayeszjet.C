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
//Int_t randNumClos=26;
Double_t jetmaxResponse = 50;
TString jetpttitle[5] = {"2-5","5-7", "7-10","10-15","15-50"};
// to fill the response matrix
bool BayesParameterSys=1;
TString truebinnum="";//TString truebinnum="JES4";//TString truebinnum="0350";
const int DnDim   = 6;//the six dimensions
//int Ddim[DnDim]   = {0,1,5,6,2,7};//for extacting 6D info from THnSparse
int zRec = 0, jetRec = 1, DRec = 2;
int zGen = 5, jetGen = 6, DGen = 7; 
int Ddim[DnDim]   = {zRec,jetRec,zGen,jetGen,DRec,DGen};//for extacting 6D info from THnSparse
TH1D *fRawSpec2Dproj[fJetptbinsN];
TH1D *fRawSpec2DprojUp[fJetptbinsN];//FD systematics
TH1D *fRawSpec2DprojDo[fJetptbinsN];//FD systematics

TH1D *LoadRawSpec(TString fn, TString sname, TString spostfix="");
TH2D *Rebin2D(const char* name, TH2D *h, int nx, const double *binx, int ny, const double *biny, bool crop);
void SparseToTree(TString MCfile, Bool_t isPostfix, TString postfix);
////===================== for 2D unfolding
int colortable[] = {kMagenta, kViolet, kBlue, kCyan+2, kGreen+4, kGreen+1, kYellow+1, kOrange+1, kRed, kRed+2};
int linesytle[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
/*****************************************************
 ############ begining of the macro ##################
 *****************************************************/
void unfold_Bayeszjet(
  TString listName = "",
  bool ispostfix=0,
  bool isprefix=0,
  bool DjetEff=1,
  TString prompteff="",
  TString effFile = "",
  TString datafile = "file.root",
  TString outDir = "out", // output directory
  const int regBayes = 5,  // default reg. parameter for the bayes unfolding
  bool isPrior = 0,  // if to use prior different than the true spectrum from the sim
  int priorType = 1,   // if isPrior == 1, choose type of the prior prior 1 is steep, prior 2 is gentle
  //bool useDeltaPt = 1,  // if to use a separate bkg. fluctuation matrix
  bool isFDUpSpec = 0,
  bool isFDDownSpec = 0,
  const int NTrials = 10,  //number of total trials
  Int_t randNumClos=10
)
{
    cout<<"random num clos:  "<<randNumClos<<endl;
    
    bool bClosure=1;double datajets=0;TRandom3 *random1 = new TRandom3(0);
    bool FDsys = 1;
    bool DptcutTrue=1;
    bool fd1D = 0; bool fd2D = 1;
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
    //bClosure
    TH2D* hMCjetD  = new TH2D("hMCjetD","hMCjetD", fptbinsZMeasN,  fptbinsZMeasA,  fJetptbinsN, fJetptbinsA); hMCjetD->Sumw2();
    TH2D* hMCjetDTree=new TH2D("hMCjetDTree","hMCjetDTree", fptbinsZMeasN,  fptbinsZMeasA,  fJetptbinsN, fJetptbinsA); hMCjetDTree->Sumw2();

    if(fd1D){
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
                double cont  = fRawSpec2Dproj[binjet]->GetBinContent(binz); double contS = fRawSpec2DprojScale->GetBinContent(binz);
                double contErr = fRawSpec2Dproj[binjet]->GetBinError(binz); double contErrS = fRawSpec2DprojScale->GetBinError(binz);
                hData2D->SetBinContent(binz,binjet+1,cont); hData2DS->SetBinContent(binz,binjet+1,contS); 
                hData2D->SetBinError(binz,binjet+1,contErr);hData2DS->SetBinError(binz,binjet+1,contErrS);
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
        cout<<"datajets222="<<hData2DS->Integral()<<endl;
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
        hData2D->SetName("hData2DRaw");
        hData2DS->SetName("hData2DScaled");
        hData2D ->Write();//raw 2D
        hData2DS->Write();//scaled 2D
        outFile ->Close();
    }
    else if(fd2D){
        TFile *fdsub  = TFile::Open(Form("%s/outFD_%d.root",datafile.Data(),(int)DjetEff));
        TH2D* hdataTemp = (TH2D*)fdsub->Get("hData_sub");
        TH2D* hdataTeUp = (TH2D*)fdsub->Get("hData_subup");
        TH2D* hdataTeDo = (TH2D*)fdsub->Get("hData_subdo");
        if(!hdataTemp || !hdataTeUp || !hdataTeDo){ cout<<"uh, oh!! no FD sub files!"<<endl;}
        hData2D = (TH2D*)hdataTemp->Clone("hData2D");
        hData2DUp = (TH2D*)hdataTeUp->Clone("hData2DUp");
        hData2DDo = (TH2D*)hdataTeDo->Clone("hData2DDo");
    }
    else{cout<<"no FD subtracted file!"<<endl;return;}
    //===============================
    // Extracting the files for charm efficiencies
    // ==============================
    TString promptefffiles[fJetptbinsN];
    TH1D* effHists[fJetptbinsN];
    for(int i=0;i<fJetptbinsN;i++){
        promptefffiles[i]=prompteff+Form("%d_%d.root",(int)fJetptbinsA[i],(int)fJetptbinsA[i+1]);
        TFile* fileEff=new TFile(promptefffiles[i].Data(),"read");
        effHists[i]=(TH1D*)fileEff->Get("hEff_reb");
    }
    // ==============================
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
    THnSparseD *hZjetRecGenD;
    THnSparseD *hZjetRGD[NDMC];
    //---- Storing the response matrix for drawing purposes ------------
    TList *histList[NDMC];THnSparseD *sparseMC[NDMC];THnSparseD *sparsereco[NDMC];
    //--- Empty 2D histograms to define binning of response matrix at reco and gen level --------------------
    TH2D* hZjetRRebin = new TH2D("hRecRebin","hRecRebin", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    TH2D* hZjetR = new TH2D("hRec","hRec", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    hZjetRRebin->Sumw2();
    TH2D* hZjetGRebin = new TH2D("hGenRebin","hGenRebin", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsN, fJetptbinsA);
    hZjetGRebin->Sumw2();
    //--- Empty hist for reco level replacing with resolution
    TH2D* hZjetRRebin_resol = new TH2D("hRecRebinResol","hRecRebinResol", 50, -1.0,1.0, fJetptbinsN, fJetptbinsA);

    //--- Reading from the ResponseMatrix of each D meson in a loop
    for(int i=0; i<NDMC; i++){
        if(!isprefix){
            if(!ispostfix) {histList[i] =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));}
            else{cout<<"-----rewrite stuff, blah blah, check again----------------"<<endl; return;}
        }
        else{
            if(ispostfix) {histList[i] =  (TList*)dir->Get(Form("%s%sMBN%dMCrec",histName.Data(),listName.Data(),i));}
            else{cout<<"-----rewrite stuff, blah blah, check again----------------"<<endl; return;}
        }
        sparseMC[i] = (THnSparseD*)histList[i]->FindObject("ResponseMatrix");
        sparseMC[i]->GetAxis(4)->SetRangeUser(-0.9+fRpar,0.9-fRpar);//pseudo-rapidity cuts
        sparseMC[i]->GetAxis(9)->SetRangeUser(-0.9+fRpar,0.9-fRpar);//pseudo-rapidity cuts

        
        //----- getting min and max of each dimension.... //is it necessary?
        for (int idim=0; idim< sparseMC[i]->GetNdimensions(); idim++){
        	auto axis = sparseMC[i]->GetAxis(i);
        	int min = 0; int max = axis->GetNbins()+1;
        	axis->SetRange(min,max);
        }
        hZjetRGD[i] = (THnSparseD*)sparseMC[i]->Projection(DnDim,Ddim,"E");// z, jet and D components
        //-----
        if (!i){
                hZjetRecGenD = (THnSparseD*)hZjetRGD[0]->Clone("hZjetRecGenD");
        }
        else {
                hZjetRecGenD->Add(hZjetRGD[i]);
        }
    }//end of NDMC for loop //hZjetRecGenD->SaveAs("proj.root");

    cout<<hZjetRecGenD->GetAxis(0)->GetXmin()<<endl;
    cout<<hZjetRecGenD->GetAxis(1)->GetXmin()<<endl;
    cout<<hZjetRecGenD->GetAxis(2)->GetXmin()<<endl;

    //********************************************************
    // Kinematic efficiency histos
    //********************************************************/
    TH2D* kineEffCutRec  = new TH2D("hkineEffCutRec","hkineEffCutRec", fptbinsZMeasN,  fptbinsZMeasA,  fJetptbinsN, fJetptbinsA); kineEffCutRec->Sumw2();
    TH2D* kineEffFulRec  = new TH2D("hkineEffFulRec","hkineEffFulRec", fptbinsZMeasN,  fptbinsZMeasA,  fJetptbinsN, fJetptbinsA); kineEffFulRec->Sumw2();
    TH2D* kineEffCutGen  = new TH2D("hkineEffCutGen","hkineEffCutGen", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsN, fJetptbinsA); kineEffCutGen->Sumw2();
    TH2D* kineEffFulGen  = new TH2D("hkineEffFulGen","hkineEffFulGen", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsN, fJetptbinsA); kineEffFulGen->Sumw2();
    /***************************
    #### unfolding settings ####
    ***************************/
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    RooUnfoldResponse response (hZjetRRebin, hZjetGRebin);
    RooUnfoldResponse responseTree (hZjetRRebin, hZjetGRebin);
    RooUnfoldResponse resol (hZjetRRebin_resol, hZjetGRebin);
    RooUnfoldResponse resol_noeff (hZjetRRebin_resol, hZjetGRebin);
    RooUnfoldResponse responseFine(hZjetRRebin, hZjetGRebin);
    //=============================
    //CLOSURE , bClosure
    //=============================
    TH2D* hZjetRRebinClos = new TH2D("hRecRebinClos","hRecRebinClos", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    hZjetRRebinClos->Sumw2();
    TH2D* hZjetGRebinClos = new TH2D("hGenRebinClos","hGenRebinClos", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsN, fJetptbinsA);
    hZjetGRebinClos->Sumw2();
    RooUnfoldResponse responseClos (hZjetRRebinClos, hZjetGRebinClos);
    //----------- fill 4D histo response matrix
    double MCjets1 = 0;
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

            //if(zR_center>=1)zR_center=0.99999;
            //if(zG_center>=1)zG_center=0.99999;

            if(isPrior){
                if     (priorType==1) {weight *= 1+(  zG_center-  jG_center/70)/2;}
                else if(priorType==2) {weight *= 1-(  zG_center-  jG_center/70)/2;}
                else if(priorType==3) {weight *= 1+(2*zG_center-  jG_center/70)/3;}
                else if(priorType==4) {weight *= 1-(2*zG_center-  jG_center/70)/3;}
                else if(priorType==5) {weight *= 1+(3*zG_center-  jG_center/70)/4;}
                else if(priorType==6) {weight *= 1-(3*zG_center-  jG_center/70)/4;}
                else if(priorType==7) {weight *= 1+(1+zG_center-2*jG_center/70)/3;}
                else if(priorType==8) {weight *= 1-(1+zG_center-2*jG_center/70)/3;}
                else weight = weight;
            }
            //========================================================================
            double eff_c = 1;
            if    (jR_center>=fJetptbinsA[0] && jR_center<=fJetptbinsA[1]){eff_c=effHists[0]->GetBinContent(effHists[0]->GetXaxis()->FindBin(DR_center));}
            else if(jR_center>fJetptbinsA[1] && jR_center<=fJetptbinsA[2]){eff_c=effHists[1]->GetBinContent(effHists[1]->GetXaxis()->FindBin(DR_center));}
            else if(jR_center>fJetptbinsA[2] && jR_center<=fJetptbinsA[3]){eff_c=effHists[2]->GetBinContent(effHists[2]->GetXaxis()->FindBin(DR_center));}
            else if(jR_center>fJetptbinsA[3] && jR_center<=fJetptbinsA[4]){eff_c=effHists[3]->GetBinContent(effHists[3]->GetXaxis()->FindBin(DR_center));}
            else if(jR_center>fJetptbinsA[4] && jR_center<=fJetptbinsA[5]){eff_c=effHists[4]->GetBinContent(effHists[4]->GetXaxis()->FindBin(DR_center));}
            //========================================================================
            bool measurement_pre = kTRUE, measurement_pos = kTRUE;//kine eff booleans
            for(int iBin=0;iBin<5;iBin++){
                //if(jG_center>=fJetptbinsA[iBin] && jG_center<=fJetptbinsA[iBin+1]){
                //if(jG_center>=fJetptbinsA[iBin] && jG_center<=fJetptbinsA[5]){
                if(jG_center>=fJetptbinsA[0] && jG_center<=jetmaxResponse){
                    if(DG_center>=fDptRangesA[iBin] && DG_center<=fDptRangesAUp[4]){
                        //if(fRpar-0.9 <= jmatch[9] && jmatch[9] <= 0.9-fRpar){
                            kineEffFulRec->Fill(zR_center,jR_center,weight);
                            kineEffFulGen->Fill(zG_center,jG_center,weight);
                            //pre unfolding kine eff--//if(DG_center<=fDptRangesAUp[4] || jG_center<=fJetptbinsA[fJetptbinsN] || jG_center>=fJetptbinsA[0] || zG_center>=fptbinsZFinalA[0]){ kineEffCutRec->Fill(zR_center,jR_center,weight); }
                            if(DG_center>fDptRangesAUp[4] || jG_center>fJetptbinsA[fJetptbinsN] || jG_center<fJetptbinsA[0] || zG_center<fptbinsZFinalA[0]){measurement_pre = kFALSE;}
                            if (measurement_pre){ kineEffCutRec->Fill(zR_center,jR_center,weight); }
                            //post unfolding kine eff--//if(DR_center<=fDptRangesAUp[4] || jR_center<=fJetptbinsA[fJetptbinsN] || jR_center>=fJetptbinsA[0] || zR_center>=fptbinsZMeasA[0]){ kineEffCutGen->Fill(zG_center,jG_center,weight); }
                            if(DR_center>fDptRangesAUp[4] || jR_center>fJetptbinsA[fJetptbinsN] || jR_center<fJetptbinsA[0] || zR_center<fptbinsZMeasA[0]){measurement_pos = kFALSE;}
                            if (measurement_pos){ kineEffCutGen->Fill(zG_center,jG_center,weight);}
                            if(jR_center>=fJetptbinsA[iBin] && jR_center<=fJetptbinsA[iBin+1]){
                                if(DR_center>=fDptRangesA[iBin] && DR_center<=fDptRangesAUp[iBin]){
                                    //if(fRpar-0.9 <= jmatch[4] && jmatch[4] <= 0.9-fRpar){
                                        response.Fill(zR_center,jR_center,zG_center,jG_center,weight/eff_c);
                                        resol.Fill((zR_center-zG_center)/zG_center,jR_center,zG_center,jG_center,weight/eff_c);
                                        resol_noeff.Fill((zR_center-zG_center)/zG_center,jR_center,zG_center,jG_center,weight);
                                            hMCjetD->Fill(zG_center,jG_center,weight/eff_c);
                                            MCjets1 += weight/eff_c;
                                        if(bClosure){
                            //hZjetGRebinClos->Fill(zG_center,jG_center,weight);//gen level closure
                            //hZjetRRebinClos->Fill(zR_center,jR_center,weight);//gen level closure
                            //responseClos.Fill(zR_center,jR_center,zG_center,jG_center,weight);//response
                                        }
                                }
                            }
                    }
                }
            }
        }
    }
    else{cout<<"No response!!"<<endl;return;}
    //=============================
    //CLOSURE , bClosure
    //=============================
    /*
    TH2D* hZjetRRebinClos = new TH2D("hRecRebinClos","hRecRebinClos", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    hZjetRRebinClos->Sumw2();
    TH2D* hZjetGRebinClos = new TH2D("hGenRebinClos","hGenRebinClos", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsN, fJetptbinsA);
    hZjetGRebinClos->Sumw2();
    RooUnfoldResponse responseClos (hZjetRRebinClos, hZjetGRebinClos);
    */
    if(bClosure){
        //reading TTree
        TFile* fTreeSparse = nullptr;   
        Double_t jmatch[13];//Dzero     
        Double_t bincontent;
        TTree* tree_ = NULL;
        TString tmp(effFile);//original MC file
        TString treepost = "TTree";if(ispostfix) treepost+=listName;
        treepost+=".root";
        tmp.Remove(effFile.Length()-5,5).Append(treepost);
        fTreeSparse = new TFile(tmp,"read");
        if(!(fTreeSparse->IsOpen())){
            //writing TTree
            cout<<"TTree file not found, converting THnSparse to TTree right now!"<<endl;
            SparseToTree(effFile,ispostfix,listName);
            //reading TTree
            cout<<"reading TTree: "<<tmp<<endl;
            fTreeSparse = new TFile(tmp,"read");
        }
        cout<<"reading TTree done"<<tmp<<endl;
        tree_ = (TTree*)fTreeSparse->Get("ResponseMatrixSum_tree");
        tree_->SetBranchAddress("coord",&jmatch);
        tree_->SetBranchAddress("bincontent",&bincontent);
        //datajets//bClosure
        datajets = hData2D->Integral();cout<<"datajets="<<datajets<<endl;
        cout<<"datajets="<<hData2DS->Integral()<<endl;
        double MCjets = hMCjetD->Integral();
        cout<<"MC jets================================="<<hMCjetD->Integral()<<endl;
        cout<<"MC jets=================================="<<MCjets1<<endl;
        //scaling the prompt jets and response for closure
        double DJscaling = 0.3;
        double RMscaling = 1-datajets/MCjets; RMscaling = 0.8; //04: 0.878577//06:0.88274403:0.887817;02:0.872606
        cout<<"RMscaling:"<<RMscaling<<endl;
        //loop over sparse bins
        for(int i=0; i<tree_->GetEntries();i++){
            tree_->GetEntry(i);
            // loop over each jet-entry: why? to get a random value about the RMscaling
            double RMw = 0;
            double DJw = 0;
            for(int ientry = 0; ientry < bincontent; ientry++){
                if(random1->Uniform(1) <= RMscaling) RMw++;
                if(random1->Uniform(1) <= DJscaling) DJw++;
            }
            //
            Double_t zR_center=jmatch[0], jR_center=jmatch[1], DR_center=jmatch[2], zG_center=jmatch[5], jG_center=jmatch[6], DG_center=jmatch[7];
            double eff_c = 1;
            // *
            if    (jR_center>=fJetptbinsA[0] && jR_center<=fJetptbinsA[1]){eff_c=effHists[0]->GetBinContent(effHists[0]->GetXaxis()->FindBin(DR_center));}
            else if(jR_center>fJetptbinsA[1] && jR_center<=fJetptbinsA[2]){eff_c=effHists[1]->GetBinContent(effHists[1]->GetXaxis()->FindBin(DR_center));}
            else if(jR_center>fJetptbinsA[2] && jR_center<=fJetptbinsA[3]){eff_c=effHists[2]->GetBinContent(effHists[2]->GetXaxis()->FindBin(DR_center));}
            else if(jR_center>fJetptbinsA[3] && jR_center<=fJetptbinsA[4]){eff_c=effHists[3]->GetBinContent(effHists[3]->GetXaxis()->FindBin(DR_center));}
            else if(jR_center>fJetptbinsA[4] && jR_center<=fJetptbinsA[5]){eff_c=effHists[4]->GetBinContent(effHists[4]->GetXaxis()->FindBin(DR_center));}
            // * /
            //========================================================================
            for(int iBin=0;iBin<5;iBin++){
                if(jG_center>=fJetptbinsA[iBin] && jG_center<=fJetptbinsA[iBin+1]){
                //if(jR_center>=fJetptbinsA[0] && jR_center<=fJetptbinsA[5]){
                    if(DG_center>=fDptRangesA[iBin] && DG_center<=fDptRangesAUp[iBin]){
                    //if(DR_center>=fDptRangesA[0] && DR_center<=fDptRangesAUp[4]){
                        if(fRpar-0.9 <= jmatch[9] && jmatch[9] <= 0.9-fRpar){
                        //if(fRpar-0.9 <= jmatch[4] && jmatch[4] <= 0.9-fRpar){
                            hZjetGRebinClos->Fill(jmatch[5],jmatch[6],bincontent-RMw);//gen level closure
                            //hZjetGRebinClos->Fill(jmatch[5],jmatch[6],DJw);//gen level closure
                        }
                    }
                }
            }
            for(int iBin=0;iBin<5;iBin++){
                //if(jG_center>=fJetptbinsA[iBin] && jG_center<=fJetptbinsA[iBin+1]){
                //if(jG_center>=fJetptbinsA[iBin] && jG_center<=fJetptbinsA[5]){
                if(jG_center>=fJetptbinsA[0] && jG_center<=jetmaxResponse){
                    if(DG_center>=fDptRangesA[iBin] && DG_center<=fDptRangesAUp[4]){
                    //if(DG_center>=fDptRangesA[iBin] && DG_center<=fDptRangesAUp[iBin]){
                        if(fRpar-0.9 <= jmatch[9] && jmatch[9] <= 0.9-fRpar){
                            if(jR_center>=fJetptbinsA[iBin] && jR_center<=fJetptbinsA[iBin+1]){
                            //if(jR_center>=fJetptbinsA[iBin] && jR_center<=fJetptbinsA[5]){
                            //if(jR_center>=fJetptbinsA[iBin] && jR_center<=jetmaxResponse){
                                //if(DR_center>=fDptRangesA[iBin] && DR_center<=fDptRangesAUp[4]){
                                if(DR_center>=fDptRangesA[iBin] && DR_center<=fDptRangesAUp[iBin]){
                                    if(fRpar-0.9 <= jmatch[4] && jmatch[4] <= 0.9-fRpar){
                                        hZjetRRebinClos->Fill(jmatch[0],jmatch[1],(bincontent-RMw)/eff_c);//rec level 
                                        //hZjetRRebinClos->Fill(jmatch[0],jmatch[1],(DJw)/eff_c);//rec level 
                                        responseClos.Fill(jmatch[0],jmatch[1],jmatch[5],jmatch[6],RMw/eff_c);//response
                                        responseTree.Fill(jmatch[0],jmatch[1],jmatch[5],jmatch[6],bincontent/eff_c);//response
                                        hMCjetDTree->Fill(zG_center,jG_center,bincontent/eff_c);
                                        //responseClos.Fill(jmatch[0],jmatch[1],jmatch[5],jmatch[6],bincontent/eff_c);//response
                                    }
                                }
                            }
                        }
                    }
                }
            }


        }
        TCanvas* cClos = new TCanvas("cClos","cClos",800,600);
        hZjetRRebinClos->Draw();
        hZjetGRebinClos->Draw("same");

    }

        // ---------------------------------------
        cout<<"MC jets Tree================================="<<hMCjetDTree->Integral()<<endl;
        TCanvas *cMC = new TCanvas("MC2D", "M2D", 800, 600);
        cMC  ->SetLogz();
        hMCjetD->Draw("colz");
        hMCjetD->Draw("TEXT SAME");
        cMC  ->SaveAs(Form("%s/alljetz2D/hMCjetD.pdf",outDir.Data()));
        cMC  ->SaveAs(Form("%s/alljetz2D/hMCjetD.png",outDir.Data()));
        cMC  ->SaveAs(Form("%s/alljetz2D/hMCjetD.svg",outDir.Data()));
    
    //
    //
    //********************************************************
    // Kinematic efficiency begins. for pre and post unfolding
    //********************************************************/
    TCanvas* cKineRec = new TCanvas("cKineRec","cKineRec",800,600);
    TH2D* hKineRec = (TH2D*) kineEffFulRec->Clone();
    hKineRec->Divide(kineEffCutRec,kineEffFulRec,1,1,"B");
    hKineRec->SetTitle("KineEffRec");
    hKineRec->Draw("colz");
    hKineRec->Draw("TEXT same");
    cKineRec  ->SaveAs(Form("%s/alljetz2D/cKineRec.pdf",outDir.Data()));
    cKineRec  ->SaveAs(Form("%s/alljetz2D/cKineRec.png",outDir.Data()));
    cKineRec  ->SaveAs(Form("%s/alljetz2D/cKineRec.svg",outDir.Data()));
    
    TCanvas* cKineGen = new TCanvas("hKineGen","hKineGen",800,600);
    TH2D* hKineGen = (TH2D*) kineEffCutGen->Clone();
    hKineGen->Divide(kineEffCutGen,kineEffFulGen,1,1,"B");
    hKineGen->SetTitle("KineEffGen");
    hKineGen->Draw("colz");
    hKineGen->Draw("TEXT same");
    cKineGen  ->SaveAs(Form("%s/alljetz2D/cKineGen.pdf",outDir.Data()));
    cKineGen  ->SaveAs(Form("%s/alljetz2D/cKineGen.png",outDir.Data()));
    cKineGen  ->SaveAs(Form("%s/alljetz2D/cKineGen.svg",outDir.Data()));

    //setting errors to zero
    for(int i=0;i<=hKineRec->GetNbinsY();i++){for(int j=0;j<=hKineRec->GetNbinsX();j++){hKineRec->SetBinError(j,i,0);}}
    for(int i=0;i<=hKineGen->GetNbinsY();i++){for(int j=0;j<=hKineGen->GetNbinsX();j++){hKineGen->SetBinError(j,i,0);}}

    //*********************************
    // Kinematic efficiency end
    //*********************************/

    //*********************************
    // Execution of unfolding
    // *******************************/
    TH2D *fUnfoldedBayes[NTrials]; TH2D *folded[NTrials]; TString outName = "unfoldedSpectrum";
    TH2D *fUnfoldedBayesClos[NTrials];TH2D *foldedClos[NTrials]; TString outNameClos = "unfoldedSpectrumClos";
    TCanvas* cUnfolded = new TCanvas("cUnfolded","cUnfolded",800,600); cUnfolded->SetLogz();
    TLegend* leg =  new TLegend(0.15,0.5,0.30,0.85); leg->SetBorderSize(0);
    //------------ do unfolding NTrials times ------------
    hData2D->Multiply(hData2D,hKineRec,1,1,"B");
    hZjetRRebinClos->Multiply(hZjetRRebinClos,hKineRec,1,1,"B");
    for(Int_t ivar=0; ivar<NTrials; ivar++){//changes
        /*****************************************
        ############# unfolding ##################
        *****************************************/
        RooUnfoldBayes unfold (&responseTree, hData2D, ivar+1);
        fUnfoldedBayes[ivar] = (TH2D*)unfold.Hreco();
        folded[ivar] = (TH2D*)responseTree.ApplyToTruth(fUnfoldedBayes[ivar]);
        if(bClosure){
            RooUnfoldBayes unfoldClos (&responseClos, hZjetRRebinClos, ivar+1);
            fUnfoldedBayesClos[ivar] = (TH2D*)unfoldClos.Hreco();
            foldedClos[ivar] = (TH2D*)responseClos.ApplyToTruth(fUnfoldedBayesClos[ivar]);
        }

        fUnfoldedBayes[ivar]->SetTitle("unfolded z-jetpt spectrum");
        if(ivar == regBayes-1){
            fUnfoldedBayes[ivar]->Draw("colz");
            fUnfoldedBayes[ivar]->Draw("TEXT SAME");}
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
    TH2D *fUnfoldedBayesUp[NTrials]; TH2D *foldedUp[NTrials]; TH2D *fUnfoldedBayesDo[NTrials]; TH2D *foldedDo[NTrials];
    if(FDsys){
        for(int ivar=0; ivar<NTrials; ivar++){//changes
            RooUnfoldBayes unfoldUp (&response, hData2DUp, ivar+1); RooUnfoldBayes unfoldDo (&response, hData2DDo, ivar+1);
            fUnfoldedBayesUp[ivar] = (TH2D*)unfoldUp.Hreco(); fUnfoldedBayesDo[ivar] = (TH2D*)unfoldDo.Hreco();
            foldedUp[ivar] = (TH2D*)response.ApplyToTruth(fUnfoldedBayesUp[ivar]); foldedDo[ivar] = (TH2D*)response.ApplyToTruth(fUnfoldedBayesDo[ivar]);
        }
    }
    /***********************************
    // FD systematics end
    ************************************/

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
    TH2D* resolobj = (TH2D*)resol.Hresponse()->Clone();
    TH2D* resol_noeffobj = (TH2D*)resol_noeff.Hresponse()->Clone();
    TH2D* respTruth = (TH2D*)response.Htruth()->Clone();
    TH2D* respMeas = (TH2D*)response.Hmeasured()->Clone();
    respTruth->Scale(1/response.Htruth()->Integral());
    respMeas->Scale(1/response.Hmeasured()->Integral());
    TFile *responseobject = new TFile(Form("%s/alljetz2D/responseobject%s%s.root",outDir.Data(), truebinnum.Data(), isPrior?Form("%d",priorType):""  ),"recreate");
    response.Write();
    respobj->SetName("responseNow");
    resolobj->SetName("resolNow");
    resol_noeffobj->SetName("resol_noeffNow");
    respobj->Write();
    resolobj->Write();
    resol_noeffobj->Write();
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
    respobjP->SetTitle("Response");
    respobjP->GetYaxis()->SetTitle("z, jetpt, GEN");
    respobjP->GetYaxis()->SetTitleOffset(1.5);
    respobjP->GetXaxis()->SetTitle("z, jetpt, RECO");
    respobjP->GetXaxis()->SetTitleOffset(1.5);
    respobjP->Draw("colz");
    //setting axes and labels properly---------
    TGaxis *axisV1 = new TGaxis(0, 0, 0, 1, 0.4, 0.6, 1, "-L");axisV1->SetLabelFont(43);axisV1->SetLabelSize(15);axisV1->SetLabelOffset(0.03);axisV1->Draw();
    TGaxis *axisV2 = new TGaxis(0, 2, 0, 5, 0.7, 1.0, 3, "-L");axisV2->SetLabelFont(43);axisV2->SetLabelSize(15);axisV2->SetLabelOffset(0.03);axisV2->Draw();
    TGaxis *axisV3 = new TGaxis(0, 6, 0,10, 0.6, 1.0, 4, "-L");axisV3->SetLabelFont(43);axisV3->SetLabelSize(15);axisV3->SetLabelOffset(0.03);axisV3->Draw();
    TGaxis *axisV4 = new TGaxis(0,11, 0,15, 0.6, 1.0, 4, "-L");axisV4->SetLabelFont(43);axisV4->SetLabelSize(15);axisV4->SetLabelOffset(0.03);axisV4->Draw();
    TGaxis *axisV5 = new TGaxis(0,16, 0,20, 0.6, 1.0, 4, "-L");axisV5->SetLabelFont(43);axisV5->SetLabelSize(15);axisV5->SetLabelOffset(0.03);axisV5->Draw();
    TGaxis *axisV6 = new TGaxis(0,21, 0,25, 0.6, 1.0, 4, "-L");axisV6->SetLabelFont(43);axisV6->SetLabelSize(15);axisV6->SetLabelOffset(0.03);axisV6->Draw();
    TGaxis *axisV7 = new TGaxis(0, 0, 0, 5,   2,   5, 1, "-L");axisV7->SetLabelFont(43);axisV7->SetLabelSize(19);axisV7->SetLabelOffset(0.06);axisV7->Draw();
    TGaxis *axisV8 = new TGaxis(0,10, 0,15,   7,  10, 1, "-L");axisV8->SetLabelFont(43);axisV8->SetLabelSize(19);axisV8->SetLabelOffset(0.06);axisV8->Draw();
    TGaxis *axisV9 = new TGaxis(0,20, 0,25,  15,  50, 1, "-L");axisV9->SetLabelFont(43);axisV9->SetLabelSize(19);axisV9->SetLabelOffset(0.06);axisV9->Draw();

    TGaxis *axisH1 = new TGaxis( 0, 0, 1, 0, 0.4,0.6, 1, "-+L");axisH1->SetLabelFont(43);    axisH1->SetLabelSize(15);    axisH1->SetLabelOffset(0.005);    axisH1->Draw();
    TGaxis *axisH2 = new TGaxis( 2, 0, 5, 0, 0.7,1.0, 3, "-+L");axisH2->SetLabelFont(43);    axisH2->SetLabelSize(15);    axisH2->SetLabelOffset(0.005);    axisH2->Draw();
    TGaxis *axisH3 = new TGaxis( 6, 0,10, 0, 0.6,1.0, 4, "-+L");axisH3->SetLabelFont(43);    axisH3->SetLabelSize(15);    axisH3->SetLabelOffset(0.005);    axisH3->Draw();
    TGaxis *axisH4 = new TGaxis(11, 0,15, 0, 0.6,1.0, 4, "-+L");axisH4->SetLabelFont(43);    axisH4->SetLabelSize(15);    axisH4->SetLabelOffset(0.005);    axisH4->Draw();
    TGaxis *axisH5 = new TGaxis(16, 0,20, 0, 0.6,1.0, 4, "-+L");axisH5->SetLabelFont(43);    axisH5->SetLabelSize(15);    axisH5->SetLabelOffset(0.005);    axisH5->Draw();
    TGaxis *axisH6 = new TGaxis(21, 0,25, 0, 0.6,1.0, 4, "-+L");axisH6->SetLabelFont(43);    axisH6->SetLabelSize(15);    axisH6->SetLabelOffset(0.005);    axisH6->Draw();
    TGaxis *axisH7 = new TGaxis( 0, 0, 5, 0,   2,  5, 1, "-+L");axisH7->SetLabelFont(43);    axisH7->SetLabelSize(19);    axisH7->SetLabelOffset(0.04);    axisH7->Draw();
    TGaxis *axisH8 = new TGaxis(10, 0,15, 0,   7, 10, 1, "-+L");axisH8->SetLabelFont(43);    axisH8->SetLabelSize(19);    axisH8->SetLabelOffset(0.04);    axisH8->Draw();
    TGaxis *axisH9 = new TGaxis(20, 0,25, 0,  15, 50, 1, "-+L");axisH9->SetLabelFont(43);    axisH9->SetLabelSize(19);    axisH9->SetLabelOffset(0.04);    axisH9->Draw();
    //respobjP->Draw();
    cRespon->SaveAs(Form("%s/alljetz2D/responseobj%s%s.pdf",outDir.Data(), truebinnum.Data(), isPrior?Form("%d",priorType):""));
    cRespon->SaveAs(Form("%s/alljetz2D/responseobj%s%s.svg",outDir.Data(), truebinnum.Data(), isPrior?Form("%d",priorType):""));

    TH2D* respobjPClos;
    if(bClosure){
        TH2D* respobjClos = (TH2D*)responseClos.Hresponse()->Clone();
        TCanvas *cResponClos = new TCanvas("","",950,800);
        cResponClos->SetLogz();
        respobjPClos = (TH2D*)respobjClos->Clone();
        respobjPClos->GetYaxis()->SetTitle("GEN:z,jet-pt");
        respobjPClos->GetXaxis()->SetTitle("REC:z,jet-pt");
        respobjPClos->GetYaxis()->SetTitle("z, jetpt, GEN");
        respobjPClos->GetYaxis()->SetTitleOffset(1.5);
        respobjPClos->GetXaxis()->SetTitle("z, jetpt, RECO");
        respobjPClos->GetXaxis()->SetTitleOffset(1.5);
        respobjPClos->SetTitle("Response for MC closure");
        respobjPClos->Draw("colz");
        respobjPClos->GetYaxis()->SetLabelSize(0.);respobjPClos->GetXaxis()->SetLabelSize(0.);
        axisV1->Draw();axisV2->Draw();axisV3->Draw();axisV4->Draw();axisV5->Draw();axisV6->Draw();axisV7->Draw();axisV8->Draw();axisV9->Draw();
        axisH1->Draw();axisH2->Draw();axisH3->Draw();axisH4->Draw();axisH5->Draw();axisH6->Draw();axisH7->Draw();axisH8->Draw();axisH9->Draw();
        cResponClos->SaveAs(Form("%s/alljetz2D/CloseResponseobj%s.pdf",outDir.Data(),isPrior?Form("%d",priorType):""));
        cResponClos->SaveAs(Form("%s/alljetz2D/CloseResponseobj%s.png",outDir.Data(),isPrior?Form("%d",priorType):""));
    }

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
    //CLOSURE bClosure
    if(bClosure){
        TFile *MCClosure = new TFile(Form("%s/alljetz2D/MCClosure%s_%d.root",outDir.Data(), isPrior?Form("%d",priorType):"",randNumClos ),"recreate");
        respobjPClos->Write();
        MCClosure->Close();
    }
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
      for (int binjet=1; binjet<= fUnfoldedBayes[regBayes]->GetNbinsY(); binjet++){
        UnfProjectX[binjet-1] = (TH1D*)fUnfoldedBayes[regBayes]->ProjectionX(Form("UnfProjectX_%d",binjet-1), binjet, binjet, "E");
        UnfProjectX[binjet-1] ->Write();
      }
      unfold2DoutFileRangeSys2->Close();
    }

    // Unfolded/Measured/Folded spectra comparison
    // ----------------------------------
    TH2D* fUnfoldPlot = (TH2D*)fUnfoldedBayes[regBayes-1]->Clone();
    fUnfoldPlot->SetLineColor(kRed+2);
    fUnfoldPlot->SetLineWidth(4);
    TCanvas *cUnFoldedFold = new TCanvas("UnFoldedMeasuredFolded","UnFolded,Measured,Folded", 1400, 900);
    TCanvas *cMCvData = new TCanvas("cMCvData","cMCvData", 1400, 900);
    TCanvas *cMCvDataStat = new TCanvas("cMCvDataStat","cMCvDataStat", 1400, 900);
    TCanvas *cUnfoldedTruthClos= new TCanvas("DataUnfoldTruthClosRatio","Closure:Data, Unfolded, Truth", 1400, 900);
    TCanvas *cUnfoldedTruthClosRatio= new TCanvas("UnfoldTruthClosRatio","Closure:Ratio Unfolded,Truth", 1400, 900);
    TLegend *lUnFoldedFold = new TLegend(0.25,0.25, 0.85, 0.85);
    TLegend *lUnFoldedFoldClos = new TLegend(0.25,0.25, 0.85, 0.85);
    TLegend *lStat = new TLegend(0.25,0.25, 0.85, 0.85);
    cUnFoldedFold->Divide(3,2);cUnfoldedTruthClos->Divide(3,2);cUnfoldedTruthClosRatio->Divide(3,2);
    cMCvData->Divide(3,2);cMCvDataStat->Divide(3,2);
    TH1D* hData2Dplot = (TH1D*)hData2D->ProjectionX("",1,1,"E")->Clone();//Datajets
    TH1D* hData2DplotRaw = (TH1D*)hData2D->ProjectionX("",1,1,"E")->Clone();//Datajets
    TH1D* hDataStat = (TH1D*)hData2D->ProjectionX("",1,1,"E")->Clone();//Datajets
    TH1D* hMCjets2Dplot = (TH1D*)hMCjetD->ProjectionX("",1,1,"E")->Clone();//MCjets
    TH1D* hMCjets2DplotTree = (TH1D*)hMCjetDTree->ProjectionX("",1,1,"E")->Clone();//MCjets
    TH1D* hMCjetsStat = (TH1D*)hMCjetD->ProjectionX("",1,1,"E")->Clone();//MCjets
    TH1D* hMCjetsStatTree = (TH1D*)hMCjetDTree->ProjectionX("",1,1,"E")->Clone();//MCjets
    TH1D* hMCclosData = (TH1D*)hZjetRRebinClos->ProjectionX("",1,1,"E");//closure MC jets
    TH1D* hMCclosDataStat = (TH1D*)hZjetRRebinClos->ProjectionX("",1,1,"E");//closure MC jets
    TH1D* hMCclosTrue = (TH1D*)hZjetGRebinClos->ProjectionX("",1,1,"E");//closure MC jets
    TH1D* hMCclosTrueStat = (TH1D*)hZjetGRebinClos->ProjectionX("",1,1,"E");//closure MC jets
    TH1D* hUnfoldPlot = (TH1D*)fUnfoldPlot->ProjectionX("",1,1,"E")->Clone();
    TH1D* hfold = (TH1D*)folded[regBayes-1]->ProjectionX("",1,1,"E")->Clone();
    TH1D* hKineUnf = (TH1D*)hKineGen->ProjectionX("",1,1)->Clone();
    TH1D* hUnfoldKine = (TH1D*)fUnfoldPlot->ProjectionX("",1,1)->Clone();
    TH1D* hDataClos;TH1D* hUnfoClos;TH1D* hTrueClos;TH1D* unfByTruth;//bClosure
    for (int binjet=1; binjet<= fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
        /***********************************
        ############# Unfolding ##################
        ************************************/
        cUnFoldedFold->cd(binjet);

        hData2Dplot = (TH1D*)hData2D->ProjectionX(Form("hData_%d",binjet),binjet,binjet,"E")->Clone();
        //hData2Dplot->SetLineColor(kBlack);
        hData2Dplot->SetLineWidth(2);
        hData2Dplot->SetMarkerStyle(20);
        hData2Dplot->SetMarkerSize(1);
        hData2Dplot->Scale(1,"width");
        hData2Dplot->Draw();
        hData2Dplot->SetTitle(Form("Jet pt: %s GeV/c",jetpttitle[binjet-1].Data()));

        hUnfoldPlot = (TH1D*)fUnfoldPlot->ProjectionX(Form("Unfold_%d",binjet),binjet,binjet,"E")->Clone();
        hUnfoldPlot->SetLineColor(kRed+2);
        hUnfoldPlot->SetLineWidth(2);
        hUnfoldPlot->SetMarkerStyle(20);
        hUnfoldPlot->SetMarkerColor(kRed+2);
        hUnfoldPlot->Scale(1,"width");
        hUnfoldPlot->Draw("same");

        hfold = (TH1D*)folded[regBayes-1]->ProjectionX(Form("hfold_%d",binjet),binjet,binjet,"E")->Clone();
        hfold->SetLineColor(kGreen+2);
        hfold->Scale(1,"width");
        hfold->Draw("same");

        //kine efficiency
        hKineUnf = (TH1D*)hKineGen->ProjectionX("",binjet,binjet)->Clone();
        hUnfoldKine = (TH1D*)fUnfoldPlot->ProjectionX("",binjet,binjet)->Clone();
        hUnfoldKine->Divide(hKineUnf);
        hUnfoldKine->SetLineColor(kOrange+1);
        hUnfoldKine->SetMarkerColor(kOrange+1);
        hUnfoldKine->SetMarkerStyle(21);
        hUnfoldKine->Scale(1,"width");
        hUnfoldKine->Draw("same");


        hData2Dplot->SetMaximum(hData2Dplot->GetMaximum()*1.5);
        hData2Dplot->SetMinimum(hData2Dplot->GetMinimum()*0.5);
        if(bClosure){
            cUnfoldedTruthClos->cd(binjet);
            hDataClos = (TH1D*)hZjetRRebinClos->ProjectionX("",binjet,binjet,"E")->Clone();
            hUnfoClos = (TH1D*)fUnfoldedBayesClos[regBayes-1]->ProjectionX("",binjet,binjet,"E")->Clone();
            hUnfoClos->Divide(hKineUnf);
            hTrueClos = (TH1D*)hZjetGRebinClos->ProjectionX("",binjet,binjet,"E")->Clone();
            hDataClos->SetLineColor(kRed+2);
            hTrueClos->SetLineColor(kGreen+2);
            hUnfoClos->SetLineColor(kBlue+2);
            hTrueClos->SetMinimum(0);
            hTrueClos->SetMaximum(hTrueClos->GetMaximum()*1.5);
            hTrueClos->Draw();
            hDataClos->Draw("same");
            hUnfoClos->Draw("same");

            cUnfoldedTruthClosRatio->cd(binjet);
            unfByTruth = (TH1D*)hUnfoClos->Clone(Form("unfByTruth%d",binjet));
            unfByTruth->Divide(hTrueClos);
            unfByTruth->SetMinimum(0);unfByTruth->SetMaximum(2);
            unfByTruth->Draw();
                TFile *MCClosure = new TFile(Form("%s/alljetz2D/MCClosure%s_%d.root",outDir.Data(), isPrior?Form("%d",priorType):"",randNumClos ),"update");
                unfByTruth->Write();
                MCClosure->Close();
        }
        //MC jets
        cMCvData->cd(binjet);
        
        hData2DplotRaw = (TH1D*)hData2D->ProjectionX(Form("hData_%d",binjet),binjet,binjet,"E")->Clone();
        hMCjets2Dplot = (TH1D*)hMCjetD->ProjectionX("",binjet,binjet,"E")->Clone();
        hMCjets2DplotTree = (TH1D*)hMCjetDTree->ProjectionX("",binjet,binjet,"E")->Clone();
        hMCclosData = (TH1D*)hZjetRRebinClos->ProjectionX("",binjet,binjet,"E");//closure MC jets
        hMCclosTrue = (TH1D*)hZjetGRebinClos->ProjectionX("",binjet,binjet,"E");//closure MC jets
        hData2DplotRaw->SetLineColor(kBlue+2);
        if(bClosure){
            hMCclosData->SetLineColor(kGreen+2);
            hMCclosTrue->SetLineColor(kOrange+2);
        }
        hData2DplotRaw->SetMaximum(hData2DplotRaw->GetMaximum()*1.5);
        hData2DplotRaw->SetMinimum(hData2DplotRaw->GetMinimum()*0.5);

        hData2DplotRaw->Draw();
        //hMCjets2Dplot->Draw("same");
        if(bClosure){
            hMCclosData->Draw("same");
            hMCclosTrue->Draw("same");
        }

        //MC jets stat unc.
        cMCvDataStat->cd(binjet);
        hMCjetsStat = (TH1D*)hMCjetD->ProjectionX("",binjet,binjet,"E")->Clone();//MCjets
        hMCjetsStatTree = (TH1D*)hMCjetDTree->ProjectionX("",binjet,binjet,"E")->Clone();//MCjets
        hDataStat = (TH1D*)hData2D->ProjectionX("",binjet,binjet,"E")->Clone();//Datajets
        hMCclosDataStat = (TH1D*)hZjetRRebinClos->ProjectionX("",binjet,binjet,"E");//closure MC jets
        hMCclosTrueStat = (TH1D*)hZjetGRebinClos->ProjectionX("",binjet,binjet,"E");//closure MC jets
        //hMCjets2Dplot = (TH1D*)hMCjetD->ProjectionX("",binjet,binjet,"E")->Clone();
        for (int iStat = 1; iStat<=hMCjets2Dplot->GetNbinsX(); iStat++){
            hMCjetsStat->SetBinContent(iStat,hMCjets2Dplot->GetBinError(iStat)/hMCjets2Dplot->GetBinContent(iStat));
            hMCjetsStatTree->SetBinContent(iStat,hMCjets2DplotTree->GetBinError(iStat)/hMCjets2DplotTree->GetBinContent(iStat));
            hDataStat->SetBinContent(iStat,hData2Dplot->GetBinError(iStat)/hData2Dplot->GetBinContent(iStat));
            hMCjetsStat->SetBinError(iStat,0);
            hMCjetsStatTree->SetBinError(iStat,0);
            hDataStat->SetBinError(iStat,0);
            if(bClosure){
                hMCclosDataStat->SetBinContent(iStat,hMCclosData->GetBinError(iStat)/hMCclosData->GetBinContent(iStat));
                hMCclosTrueStat->SetBinContent(iStat,hMCclosTrue->GetBinError(iStat)/hMCclosTrue->GetBinContent(iStat));
                hMCclosDataStat->SetBinError(iStat,0);
                hMCclosTrueStat->SetBinError(iStat,0);
            }
        }
        hDataStat->SetMaximum(0.5);
        hDataStat->SetMinimum(0);
        hMCjetsStat->SetLineColor(kRed+2);
        hDataStat->SetLineColor(kBlue+2);
        if(bClosure){
            hMCclosDataStat->SetLineColor(kGreen+2);
            hMCclosTrueStat->SetLineColor(kOrange+2);
        }

        hDataStat->Draw();
        hMCjetsStat->Draw("same");
        hMCjetsStatTree->Draw("same");
        if(bClosure){
            hMCclosDataStat->Draw("same");
            hMCclosTrueStat->Draw("same");
        }

    }
    cMCvDataStat->cd(fUnfoldedBayes[regBayes-1]->GetNbinsY()+1);
    lStat->AddEntry(hDataStat,"Data stat","l");
    lStat->AddEntry(hMCjetsStat,"MC stat","l");
    lStat->AddEntry(hMCjetsStatTree,"MC Tree stat","l");
    if(bClosure){
        lStat->AddEntry(hMCclosDataStat,"MC clos data stat","l");
        lStat->AddEntry(hMCclosTrueStat,"MC clos true stat","l");
    }
    lStat->Draw("same");
    cMCvData->SaveAs(Form("%s/alljetz2D/MCvData.pdf",outDir.Data()));
    cMCvData->SaveAs(Form("%s/alljetz2D/MCvData.png",outDir.Data()));
    cMCvDataStat->SaveAs(Form("%s/alljetz2D/MCvDataStat.pdf",outDir.Data()));
    cMCvDataStat->SaveAs(Form("%s/alljetz2D/MCvDataStat.png",outDir.Data()));
    //legend in the last tpad
    cUnFoldedFold->cd(fUnfoldedBayes[regBayes-1]->GetNbinsY()+1);
    lUnFoldedFold->AddEntry(hData2Dplot,"Measured","l");
    lUnFoldedFold->AddEntry(hUnfoldPlot,Form("Unfolded, iter=%d",regBayes),"l");
    lUnFoldedFold->AddEntry(hfold,"Folded","l");
    lUnFoldedFold->AddEntry(hUnfoldKine,"Kine Eff corr");
    lUnFoldedFold->Draw("same");
    cUnFoldedFold->SaveAs(Form("%s/alljetz2D/unfoldMeasFold.pdf",outDir.Data()));
    cUnFoldedFold->SaveAs(Form("%s/alljetz2D/unfoldMeasFold.png",outDir.Data()));


    if(bClosure){
        cUnfoldedTruthClos->cd(fUnfoldedBayesClos[regBayes-1]->GetNbinsY()+1);
        lUnFoldedFoldClos->AddEntry(hDataClos,"MC data");
        lUnFoldedFoldClos->AddEntry(hTrueClos,"MC truth");
        lUnFoldedFoldClos->AddEntry(hUnfoClos,"MC unfold");
        lUnFoldedFoldClos->Draw("same");
        cUnfoldedTruthClos->SaveAs(Form("%s/alljetz2D/unfoldMeasTruthClos.pdf",outDir.Data()));
        cUnfoldedTruthClos->SaveAs(Form("%s/alljetz2D/unfoldMeasTruthClos.png",outDir.Data()));
        cUnfoldedTruthClosRatio->SaveAs(Form("%s/alljetz2D/unfoldMeasTruthClosRatio.pdf",outDir.Data()));
        cUnfoldedTruthClosRatio->SaveAs(Form("%s/alljetz2D/unfoldMeasTruthClosRatio.png",outDir.Data()));
    }
    // Unfolded/Measured spectra comparison
    // ------------------------------------
    TCanvas *cUnFolded = new TCanvas("UnFoldedMeasured","UnFolded/Measured", 1400, 900);
    cUnFolded->Divide(3,2);
    cUnFolded->SetLogz();
    for (int binjet=1; binjet< fUnfoldedBayes[regBayes-1]->GetNbinsY()+1; binjet++){
        /***********************************
        ############# Unfolding comparison
        ************************************/
        cUnFolded->cd(binjet);
        TH1D* hData2Dplot = (TH1D*)hData2D->ProjectionX("",binjet,binjet,"E")->Clone();
        //hData2Dplot->SetLineColor(kBlack);
        hData2Dplot->SetLineWidth(2);
        hData2Dplot->SetMarkerStyle(20);
        hData2Dplot->SetMarkerSize(1);
        hData2Dplot->Draw();
        hData2Dplot->SetTitle(Form("Jet pt: %s GeV/c",jetpttitle[binjet-1].Data()));
        for(int ivar=0; ivar<NTrials; ivar++){//changes
            fUnfoldedBayes[ivar]->SetLineColor(colortable[ivar]);
            fUnfoldedBayes[ivar]->ProjectionX("",binjet,binjet,"E")->Draw("same");
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
        TH1D* hData2Dplot = (TH1D*)hData2D->ProjectionX("",binjet,binjet,"E")->Clone();
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
    cFolded->cd(fUnfoldedBayes[regBayes-1]->GetNbinsY());
    for(int ivar=0; ivar<NTrials; ivar++){//changes
        TH1D* hRatio = (TH1D*)folded[ivar]->ProjectionX("",1,1,"E");
        hRatio->SetLineColor(colortable[ivar]);
    }
    // Folded/Measured spectra comparison ratio
    // ----------------------------------------
    TLine *lm = new TLine(0.4, 1, 1.02, 1);
    lm->SetLineStyle(2);

    TCanvas *cFoldedR = new TCanvas("FoldedMeasuredRatio","Folded/Measured Ratio", 1400, 900);
    cFoldedR->Divide(3,2);
    TH1D* hRatio[fUnfoldedBayes[regBayes-1]->GetNbinsY()][NTrials];
    for (int binjet=1; binjet<= fUnfoldedBayes[regBayes-1]->GetNbinsY(); binjet++){
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
        }
        lm->Draw("same");
    }
    TLegend *lFoldedR = new TLegend(0.25,0.25, 0.85, 0.85);
    cFoldedR->cd(fUnfoldedBayes[regBayes-1]->GetNbinsY());
    for(int ivar=0; ivar<NTrials; ivar++){//changes
        TH1D* hRatio = (TH1D*)folded[ivar]->ProjectionX("",1,1,"E");
        hRatio->SetLineColor(colortable[ivar]);
        lFoldedR->AddEntry(hRatio,Form("Bayes iter = %d",ivar+1),"l");
    }
    cout<<"Djeteff:"<<DjetEff<<endl;
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
//----
//void setRespAxis(TGaxis *gaxis){
//    gaxis->SetLabelFont(43);
//    gaxis->SetLabelSize(15);
//    gaxis->SetLabelOffset(0.03);
//    gaxis->Draw();
//}
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
void SparseToTree(TString MCfile, Bool_t isPostfix, TString postfix)
{
// Read the THnSparse
    TFile *File1 = new TFile(MCfile,"read");
    TDirectoryFile* dir=(TDirectoryFile*)File1->Get("DmesonsForJetCorrelations");
    TString histName;
    if(fDmesonSpecie) histName = "histosDStarMBN";
                    else histName = "histosD0MBN";
    TList *histList[NDMC];
    THnSparseF *sparseMC[NDMC];
    THnSparseF *h = nullptr;

    for(int i=0; i<NDMC; i++){
        histList[i] =  (TList*)dir->Get(Form("%s%d%sMCrec",histName.Data(),i,isPostfix?postfix.Data():""));
        sparseMC[i] = (THnSparseF*)histList[i]->FindObject("ResponseMatrix");
        std::cout<<sparseMC[i]->GetNbins()<<std::endl;
        if(!i)h = dynamic_cast<THnSparseF*>(sparseMC[0]->Clone("ResponseMatrixSum"));
        else h->Add(sparseMC[i]);
    }
    //copy from here: https://root.cern.ch/root/html/tutorials/tree/drawsparse.C.html
    // Creates a TTree and fills it with the coordinates of all
    // filled bins. The tree will have one branch for each dimension,
    // and one for the bin content.
    TString outName(MCfile);
    outName.Remove(outName.Length()-5,5);
    outName+="TTree";
    if(isPostfix)outName+=postfix;
    outName+=".root";
    std::cout<<"Converting THnSparse: "<<MCfile<<std::endl;
    std::cout<<"Into TTree: "<<outName<<std::endl;
   TFile *File2 = new TFile(outName,"RECREATE");
   Int_t dim = h->GetNdimensions();
   TString name(h->GetName()); name += "_tree";
   TString title(h->GetTitle()); title += " tree";

   TTree* tree = new TTree(name, title);
   Double_t* x = new Double_t[dim + 1];
   memset(x, 0, sizeof(Double_t) * (dim + 1));

   TString branchname;
   for (Int_t d = 0; d < dim; ++d) {
      if (branchname.Length())
         branchname += ":";
      TAxis* axis = h->GetAxis(d);
      branchname += axis->GetName();
      branchname += "/D";
   }
   tree->Branch("coord", x, branchname);
   tree->Branch("bincontent", &x[dim], "bincontent/D");

   const int percentPrint = 1;
   int totalWindows = h->GetNbins();
   int percent;
   int step = totalWindows / (100/percentPrint);
   int nextPrint = step;

   Int_t *bins = new Int_t[dim];
   for (Long64_t i = 0; i < h->GetNbins(); ++i) {
      if (i >= nextPrint)
       {
           percent = (100 * i) / totalWindows;
           std::cout << "\r" << std::string(percent/percentPrint , '|') << percent << "%";
           std::cout.flush();
           nextPrint += step;
       }
      x[dim] = h->GetBinContent(i, bins);
      for (Int_t d = 0; d < dim; ++d) {
         x[d] = h->GetAxis(d)->GetBinCenter(bins[d]);
      }

      tree->Fill();

   }
   std::cout<<std::endl<<"tree entries "<< tree->GetEntries()<<std::endl;
   File2->Write(nullptr,TObject::kOverwrite);
   File2->Close();
}
