//
// Macro to subtract B feed-down from inclusive z_||,ch spectrum before unfolding
//
// Author: A.Mohanty  (auro.mohanty@cern.ch)
//

#include "../DsignalExtraction/configDzero_ppz.h"
//Global
//------
const int DnDim   = 6;
int zRec = 0, jetRec = 1, DRec = 2;
int zGen = 5, jetGen = 6, DGen = 7;
int Ddim[DnDim]   = {zRec,jetRec,zGen,jetGen,DRec,DGen};//for extacting 6D info from THnSparse

TH1D* LoadRawSpec(TString fn, TString sname, TString spostfix="");
TH1* GetInputHist(TString inFile, string histName,TH1 *hh);
void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, double Msize = 1.1, Width_t Lwidth = 2, Style_t Lstyle = 1);
void setErrorsZero2D(TH2D* hh);
void setErrorsZero1D(TH1D* hh);

TH1D* fRawSpec2Dproj[fJetptbinsN];
//=====================================================================
void subtractFD_zjet(
    TString dataInFileDir="",
    TString outDir="",
    TString simDir="",
    TString dataAnalysisFile="",
    TString MCAnalysisFile="",
    bool isprefix = 0,
    bool ispostfix = 0,
    bool DjetEff=1,
    TString efffile=""
)
{
    gStyle->SetOptStat(0000);
//0.Data Luminousity
//--------------------------------------------------------------
    TFile* Filedata = new TFile(dataAnalysisFile,"read");
    TDirectoryFile* dir = (TDirectoryFile*)Filedata->Get("PWG3_D2H_DmesonsForJetCorrelationsMBN0");
    AliNormalizationCounter* c = (AliNormalizationCounter*)dir->Get("NormalizationCounter");
    double nEv = c->GetNEventsForNorm();
    double dataLum = nEv/(sigma_in*1000);//Luminosity in mbar
    double simScaling = BRDzero*dataLum;
//1.Extracting data/effC, combining into a 2D histogram
//--------------------------------------------------------------
    //need the signalextraction folder:DONE "%s/Z0to102_jetbin_%d_%d/JetPtSpectra_SB_eff.root"
    //define histos: DONE
    TH2* hData2D = new TH2D("hdata_zjet", "hdata_zjet", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    TString data2D[fJetptbinsN];
    //reading data from each 1D: DONE
    for (int binjet = 0; binjet < fJetptbinsN; binjet++){
        //files
        data2D[binjet]=Form("%s/Z0to102_jetbin_%d_%d/JetPtSpectra_SB_eff.root",dataInFileDir.Data(),(int)fJetptbinsA[binjet],(int)fJetptbinsA[binjet+1]);
        //histos
        fRawSpec2Dproj[binjet]=(TH1D*)LoadRawSpec(data2D[binjet].Data(),"hjetptspectrumReb");
    }
    //combining 1D data into 2D: DONE
    for (int binjet=0; binjet<fJetptbinsN; binjet++){
        for (int binz=0; binz<fptbinsZMeasN+1; binz++){
            double cont  = fRawSpec2Dproj[binjet]->GetBinContent(binz);
            double contE = fRawSpec2Dproj[binjet]->GetBinError(binz);
            hData2D->SetBinContent(binz,binjet+1,cont);
            hData2D->SetBinError(binz,binjet+1,contE);
        }
    }
    //TCanvas* cData = new TCanvas("EffCorr2Ddata","EffCorr2Ddata", 800, 600);TH1D* h1 = (TH1D*)hData2D->ProjectionX("h1",2,2,"E");h1->Draw();
    //now we have eff corrected 2D histogram ready for B FD subtraction
    
//2.Extracting simulation*effB/effC, combining into a 2D histogram
//--------------------------------------------------------------
    //need the simulation folder: DONE
    int cent =0;//central value, used afterwards. can be replaced with simNr
    int simNr=0;//0-central value
    int nFiles=fBsimN;
    TH1D* hFD[fJetptbinsN][fBsimN];
    TH1D* hFD_binned[fJetptbinsN][fBsimN];
    //storing the 1D histos from sim files: DONE
    for (int binjet=0; binjet<fJetptbinsN; binjet++){
        for (int nr=simNr; nr<nFiles; nr++){
            TString file = simDir;
            file += Form("/Jetpt%d_%d",(int)fJetptbinsA[binjet],(int)fJetptbinsA[binjet+1]);
            file += "/Z_"; file+= fRunB[nr]; file+="_Jetpt"; file+=(int)fJetptbinsA[binjet]; file += "_"; file += (int)fJetptbinsA[binjet+1];
            file += "_Dpt"; file += fDptRangesA[binjet]; file += "_"; file += fDptRangesAUp[binjet]; file += "_effScaled"; file += "_Dzero"; file += ".root";
            TH1D* htmp;
            htmp = (TH1D*) GetInputHist(file,"hPt",htmp);
            hFD[binjet][nr] = (TH1D*)htmp->Clone(Form("hFD_%d_%d",binjet,nr));
            hFD_binned[binjet][nr] = (TH1D*)htmp->Rebin(fptbinsZTrueN,Form("hFD_binned_%d_%d",binjet,nr),fptbinsZTrueA);
        }
    }
    //combining them to 2D histos: DONE
    TH2D* h2dFD[fBsimN];
    for (int nr=simNr; nr<nFiles; nr++){
        h2dFD[nr] = new TH2D(Form("h2dFD_%d",nr), Form("h2dFD_%d",nr), fptbinsZTrueN, fptbinsZTrueA, fJetptbinsN, fJetptbinsA);
        for (int binjet=0; binjet<fJetptbinsN; binjet++){//jetbins start with numbering 1, but here jetbin also represents the sl no, which starts with 0
            for (int binz = 0; binz < fptbinsZMeasN+1; binz++){//z bins start with numbering 0
                double cont = hFD_binned[binjet][nr]->GetBinContent(binz);
                double contE= hFD_binned[binjet][nr]->GetBinError(binz);

                h2dFD[nr]->SetBinContent(binz,binjet+1,cont);
                h2dFD[nr]->SetBinError(binz,binjet+1,contE);
            }
        }
    }

    //up and down sim values: DONE
    TH2D* hFD_up = (TH2D*)h2dFD[0]->Clone("hFD_up");
    TH2D* hFD_do = (TH2D*)h2dFD[0]->Clone("hFD_do");
    for(int j=1;j<=h2dFD[0]->GetNbinsY();j++){
        for(int i=1;i<=h2dFD[0]->GetNbinsX();i++){
            double maxfd = 0, minfd = 0;
            maxfd = minfd = h2dFD[0]->GetBinContent(i,j);
            for (int simsl=1;simsl<nFiles;simsl++){//run through other sim except 0
                if (h2dFD[simsl]->GetBinContent(i,j)>maxfd){maxfd=h2dFD[simsl]->GetBinContent(i,j);}
                else if (h2dFD[simsl]->GetBinContent(i,j)<minfd){minfd=h2dFD[simsl]->GetBinContent(i,j);}
            }
            hFD_up->SetBinContent(i,j,maxfd);
            hFD_up->SetBinError(i,j,0);
            hFD_do->SetBinContent(i,j,minfd);
            hFD_do->SetBinError(i,j,0);
        }
    }
    //central sim value: DONE
    TH2D* hFD_central = (TH2D*)h2dFD[0]->Clone("hFD_central");
    hFD_central->Scale(simScaling);
    hFD_up->Scale(simScaling);
    hFD_do->Scale(simScaling);

//3.Extracting the files for charm efficiencies: NOTDONE
//---------------------------------------------
    //access the efficiency folder
    double *efficiency = 0x0;
    TString effFiles[fJetptbinsN];
    TH1D* effHists[fJetptbinsN];
    for(int i=0;i<fJetptbinsN;i++){
        effFiles[i]=efffile+Form("%d_%d.root",(int)fJetptbinsA[i],(int)fJetptbinsA[i+1]);
        TFile* fileEff=new TFile(effFiles[i].Data(),"read");
        effHists[i]=(TH1D*)fileEff->Get("hEff_reb");
    }
//4.Getting the 2D Response with extra DptReco dimension to scale by charm efficiencies
//-------------------------------------------------------------------------------------
    //access the MC train output file to create the RM: DONE
    TFile* FileMC = new TFile(MCAnalysisFile,"read");
    TDirectoryFile* dirMC  = (TDirectoryFile*)FileMC->Get("DmesonsForJetCorrelations");
    TString histName;
    if(!isprefix){histName = "histosD0MBN";}
    else{histName = "histosD0";}
    //reading from the ResponseMatrix of each D meson in a loop: DONE
    TList* histList[NDMC];
    THnSparseD *sparseMC[NDMC];
    THnSparseD *hZjetRGD[NDMC];
    THnSparseD *hZjetRecGenD;
    for(int i=0; i<NDMC; i++){
        if(!isprefix){
            if(!ispostfix){
                histList[i]=(TList*)dirMC->Get(Form("%s%dFDMCrec",histName.Data(),i));
            }else{return;}
        }else{return;}
        sparseMC[i] = (THnSparseD*)histList[i]->FindObject("ResponseMatrix");
        sparseMC[i]->GetAxis(4)->SetRangeUser(-0.9+fRpar,0.9-fRpar);//pseudo-rapidity cuts on Reco level
        sparseMC[i]->GetAxis(9)->SetRangeUser(-0.9+fRpar,0.9-fRpar);//pseudo-rapidity cuts on Gen level

        hZjetRGD[i] = (THnSparseD*)sparseMC[i]->Projection(DnDim,Ddim,"E");// z, jet and D components
        if(!i){hZjetRecGenD = (THnSparseD*)hZjetRGD[0]->Clone("hZjetRecGenD");}
        else{hZjetRecGenD->Add(hZjetRGD[i]);}
    }

    //setting up dimensions for response matrix: DONE
    TH2D* hZjetRRebin = new TH2D("hRecRebin","hRecRebin", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA); hZjetRRebin->Sumw2();
    TH2D* hZjetGRebin = new TH2D("hGenRebin","hGenRebin", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsN, fJetptbinsA); hZjetGRebin->Sumw2();
    //setting up kinematic efficiency histos
    TH2D* kineEffCutPre  = new TH2D("hkineEffCutPre","hkineEffCutPre", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA); kineEffCutPre->Sumw2();
    TH2D* kineEffFulPre  = new TH2D("hkineEffFulPre","hkineEffFulPre", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA); kineEffFulPre->Sumw2();
    TH2D* kineEffCutPos  = new TH2D("hkineEffCutPos","hkineEffCutPos", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsN, fJetptbinsA); kineEffCutPos->Sumw2();
    TH2D* kineEffFulPos  = new TH2D("hkineEffFulPos","hkineEffFulPos", fptbinsZFinalN, fptbinsZFinalA, fJetptbinsN, fJetptbinsA); kineEffFulPos->Sumw2();

    //making the response: DONE
    RooUnfoldResponse response (hZjetRRebin, hZjetGRebin);
    for(int z=0; z<hZjetRecGenD->GetNbins();z++){
        int coord[DnDim+1] = {0,0,0,0,0,0};
        double content = hZjetRecGenD->GetBinContent(z,coord);
        int i=coord[0], j=coord[1], k=coord[2], m=coord[3], n=coord[4], p=coord[5];
        double weight = content;
        double zR_center = hZjetRecGenD->GetAxis(0)->GetBinCenter(i);
        double jR_center = hZjetRecGenD->GetAxis(1)->GetBinCenter(j);
        double zG_center = hZjetRecGenD->GetAxis(2)->GetBinCenter(k);
        double jG_center = hZjetRecGenD->GetAxis(3)->GetBinCenter(m);
        double DR_center = hZjetRecGenD->GetAxis(4)->GetBinCenter(n);
        double DG_center = hZjetRecGenD->GetAxis(5)->GetBinCenter(p);

        bool measurement_ok = kTRUE; bool measurement_pre = kTRUE; bool measurement_pos = kTRUE;
        double effresp = content;//needs to be changed by the charm efficiency NOTDONE
        if( DR_center<fDptRangesA[0] || DG_center< fJetptbinsA[0] ){measurement_ok=kFALSE;}//not particularly necessary
        else if( (DG_center<fDptRangesA[0] && jG_center>=fJetptbinsA[0]) || (DR_center<fDptRangesA[0] && jR_center>=fJetptbinsA[0]) ){measurement_ok=kFALSE;}
        else if( (DG_center<fDptRangesA[1] && jG_center>=fJetptbinsA[1]) || (DR_center<fDptRangesA[1] && jR_center>=fJetptbinsA[1]) ){measurement_ok=kFALSE;}
        else if( (DG_center<fDptRangesA[2] && jG_center>=fJetptbinsA[2]) || (DR_center<fDptRangesA[2] && jR_center>=fJetptbinsA[2]) ){measurement_ok=kFALSE;}
        else if( (DG_center<fDptRangesA[3] && jG_center>=fJetptbinsA[3]) || (DR_center<fDptRangesA[3] && jR_center>=fJetptbinsA[3]) ){measurement_ok=kFALSE;}
        else if( (DG_center<fDptRangesA[4] && jG_center>=fJetptbinsA[4]) || (DR_center<fDptRangesA[4] && jR_center>=fJetptbinsA[4]) ){measurement_ok=kFALSE;}

        if (measurement_ok){ kineEffFulPre->Fill(zR_center,jR_center,effresp);kineEffFulPos->Fill(zG_center,jG_center,weight);}
        measurement_pre=measurement_pos=measurement_ok;
        //pre unfolding kine = post folding kine
        if(DG_center>fDptRangesAUp[4] || jG_center>fJetptbinsA[fJetptbinsN] || jG_center<fJetptbinsA[0] || zG_center<fptbinsZMeasA[0]){measurement_pre = kFALSE;}
        if(measurement_pre){ kineEffCutPre->Fill(zR_center,jR_center,effresp); }
        //post unfolding kine = pre folding kine
        if(DR_center>fDptRangesAUp[4] || jR_center>fJetptbinsA[fJetptbinsN] || jR_center<fJetptbinsA[0] || zR_center<fptbinsZMeasA[0]){measurement_pos = kFALSE;}
        if(measurement_pos){ kineEffCutPos->Fill(zG_center,jG_center,weight ); }
        //original response
        if(DR_center>fDptRangesAUp[4] || jR_center>fJetptbinsA[fJetptbinsN] || zR_center<fptbinsZMeasA[0]) {measurement_ok=kFALSE;}
        //application of D-jet efficiency to the weight
        //??
        if(measurement_ok){
            if(DjetEff){
                double eff_c = 1;
                if     (jR_center>fJetptbinsA[0] && jR_center<=fJetptbinsA[1]){eff_c=effHists[0]->GetBinContent(effHists[0]->GetXaxis()->FindBin(DR_center));}
                else if(jR_center>fJetptbinsA[1] && jR_center<=fJetptbinsA[2]){eff_c=effHists[1]->GetBinContent(effHists[1]->GetXaxis()->FindBin(DR_center));}
                else if(jR_center>fJetptbinsA[2] && jR_center<=fJetptbinsA[3]){eff_c=effHists[2]->GetBinContent(effHists[2]->GetXaxis()->FindBin(DR_center));}
                else if(jR_center>fJetptbinsA[3] && jR_center<=fJetptbinsA[4]){eff_c=effHists[3]->GetBinContent(effHists[3]->GetXaxis()->FindBin(DR_center));}
                else if(jR_center>fJetptbinsA[4] && jR_center<=fJetptbinsA[5]){eff_c=effHists[4]->GetBinContent(effHists[4]->GetXaxis()->FindBin(DR_center));}
                weight = weight/eff_c;
            }
            response.Fill(zR_center,jR_center,zG_center,jG_center,weight);
        }
    }
    //getting kine efficiencies: DONE
    TCanvas* cKinePre = new TCanvas("cKinePre","cKinePre",800,600);
    TH2D* hKinePre = (TH2D*) kineEffFulPre->Clone();
    hKinePre->Divide(kineEffCutPre, kineEffFulPre,1,1,"B");
    hKinePre->SetTitle("KineEffPre");
    hKinePre->Draw("colz");
    hKinePre->Draw("TEXT same");
    cKinePre  ->SaveAs(Form("%s/cKinePre.pdf",outDir.Data()));
    cKinePre  ->SaveAs(Form("%s/cKinePre.png",outDir.Data()));
    cKinePre  ->SaveAs(Form("%s/cKinePre.svg",outDir.Data()));
    //getting kine efficiencies: DONE
    TCanvas* cKinePos = new TCanvas("cKinePos","cKinePos",800,600);
    TH2D* hKinePos = (TH2D*) kineEffFulPos->Clone();
    hKinePos->Divide(kineEffCutPos, kineEffFulPos,1,1,"B");
    hKinePos->SetTitle("KineEffPos");
    hKinePos->Draw("colz");
    hKinePos->Draw("TEXT same");
    cKinePos  ->SaveAs(Form("%s/cKinePos.pdf",outDir.Data()));
    cKinePos  ->SaveAs(Form("%s/cKinePos.png",outDir.Data()));
    cKinePos  ->SaveAs(Form("%s/cKinePos.svg",outDir.Data()));
    //setting errors to zero: DONE
    for(int j=0;j<=hKinePre->GetNbinsY();j++){for(int i=0;i<=hKinePre->GetNbinsX();i++){hKinePre->SetBinError(i,j,0);}}
    for(int j=0;j<=hKinePos->GetNbinsY();j++){for(int i=0;i<=hKinePos->GetNbinsX();i++){hKinePos->SetBinError(i,j,0);}}
    //apply kine 1 (kinepos)
    TH2D* hFD_centralK = (TH2D*)hFD_central->Clone("hFD_centralK");
    TH2D* hFD_upK      = (TH2D*)hFD_up->Clone("hFD_upK");
    TH2D* hFD_doK      = (TH2D*)hFD_do->Clone("hFD_doK");
    hFD_centralK->Multiply(hKinePos);
    hFD_upK->Multiply(hKinePos);
    hFD_doK->Multiply(hKinePos);
    //folding: DONE
    TH2D* hFolded = new TH2D("hFolded","hFolded", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    TH2D* hFoldUp = new TH2D("hFoldup","hFoldUp", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    TH2D* hFoldDo = new TH2D("hFoldup","hFoldDo", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    hFolded = (TH2D*)response.ApplyToTruth(hFD_centralK);
    hFoldUp = (TH2D*)response.ApplyToTruth(hFD_upK);
    hFoldDo = (TH2D*)response.ApplyToTruth(hFD_doK);
    //apply kine 2 (kinepre)
    TH2D* hFoldedK = (TH2D*)hFolded->Clone("hFoldedK");
    TH2D* hFoldUpK = (TH2D*)hFoldUp->Clone("hFoldUpK");
    TH2D* hFoldDoK = (TH2D*)hFoldDo->Clone("hFoldDoK");
    hFoldedK->Divide(hFoldedK,hKinePre);
    hFoldUpK->Divide(hFoldUpK,hKinePre);
    hFoldDoK->Divide(hFoldDoK,hKinePre);

//5.Subtracting folded simulation from data: DONE
//-----------------------------------------
    TH2D* hData_sub   = (TH2D*) hData2D->Clone("hData_sub");
    TH2D* hData_subup = (TH2D*) hData2D->Clone("hData_subup");
    TH2D* hData_subdo = (TH2D*) hData2D->Clone("hData_subdo");
    hData_sub->Add(hFoldedK,-1);
    hData_subup->Add(hFoldDoK,-1);
    hData_subdo->Add(hFoldUpK,-1);

    TCanvas* cFDsub = new TCanvas("FDsub","FDsub",900,600);
    cFDsub->Divide(3,2);
    TCanvas* cFDratio = new TCanvas("FDratio","FDratio",900,600);
    cFDratio->Divide(3,2);
    for (int binjet = 1; binjet<hData_sub->GetNbinsY()+1;binjet++){
        cFDsub->cd(binjet);
        TH1D* hhdata  =(TH1D*)hData2D->ProjectionX("",binjet,binjet,"E")->Clone(Form("hdata%d",binjet));

        TH1D* hhdata_sub  =(TH1D*)hData_sub  ->ProjectionX("",binjet,binjet,"E")->Clone(Form("hsub_c%d",binjet));
        TH1D* hhdata_subup=(TH1D*)hData_subup->ProjectionX("",binjet,binjet,"E")->Clone(Form("hsub_u%d",binjet));
        TH1D* hhdata_subdo=(TH1D*)hData_subdo->ProjectionX("",binjet,binjet,"E")->Clone(Form("hsub_d%d",binjet));

        TH1D* hhFD_ce  = (TH1D*)hFoldedK->ProjectionX("",binjet,binjet,"E")->Clone(Form("hFD_ce%d",binjet));
        TH1D* hhFD_up = (TH1D*)hFoldUpK->ProjectionX("",binjet,binjet,"E")->Clone(Form("hFD_up%d",binjet));
        TH1D* hhFD_do = (TH1D*)hFoldDoK->ProjectionX("",binjet,binjet,"E")->Clone(Form("hFD_do%d",binjet));

        TH1D* hhData=(TH1D*)hData2D->ProjectionX("",binjet,binjet,"E")->Clone(Form("hData%d",binjet));
        hhData->Draw();
        hhdata_sub  ->Draw("same");
        hhdata_subup->Draw("same");
        hhdata_subdo->Draw("same");
        hhFD_up->Draw("same");
        hhFD_do->Draw("same");
        
        cFDratio->cd(binjet);
        TH1D* hhratio = (TH1D*)hhFD_ce->Clone(Form("ratio_%d",binjet));setHistoDetails(hhratio,8,20);
        TH1D* hhratUp = (TH1D*)hhFD_up->Clone(Form("ratUp_%d",binjet));setHistoDetails(hhratUp,8,24,0,2,2);
        TH1D* hhratDo = (TH1D*)hhFD_do->Clone(Form("ratDo_%d",binjet));setHistoDetails(hhratDo,8,24,0,2,2);
        hhratio->Divide(hhdata);hhratUp->Divide(hhdata);hhratDo->Divide(hhdata);
        setErrorsZero1D(hhratUp);setErrorsZero1D(hhratDo);//setting errors to zero
        hhratio->GetYaxis()->SetRangeUser(0,1.15);
        hhratio->Draw();
        hhratUp->Draw("same");
        hhratDo->Draw("same");
    }

//6.Writing histograms to file
    TFile* outFile = new TFile(Form("%s/histosForFD_%d.root",outDir.Data(),(int)DjetEff),"recreate");
    hData2D->Write();
    outFile->Close();
return;
}

//========================================================================
/// load 1D raw spectrum
TH1D *LoadRawSpec(TString fn, TString sname, TString spostfix="") {
    TFile *f  = TFile::Open(fn);
    if (!f) { Error("LoadRawSpectrum","Raw spectrum file %s not found.",fn.Data()); return 0; }
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

TH1* GetInputHist(TString inFile, string histName,TH1 *hh){
    TFile *jetPtFile = new TFile(inFile,"read");
    hh = (TH1F*)jetPtFile->Get(histName.c_str());
    return hh;

}

void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, double Msize, Width_t Lwidth, Style_t Lstyle){

    hh->SetMarkerColor(color);
    hh->SetMarkerStyle(Mstyle);;
    hh->SetLineColor(color);
    hh->SetLineWidth(Lwidth);
    hh->SetMarkerSize(Msize);
    hh->SetLineStyle(Lstyle);
   // hh->SetName(name.c_str());
    hh->SetTitle("");
    //hh->GetXaxis()->SetTitle("p_{T,ch jet} (GeV/c)");
    hh->GetXaxis()->SetTitle("#it{z}_{||}");

}
void setErrorsZero2D(TH2D* hh){
    for(int j=1;j<=hh->GetNbinsY();j++){
        for(int i=1;i<=hh->GetNbinsX();i++){
            hh->SetBinError(i,j,0);
        }
    }
}
void setErrorsZero1D(TH1D* hh){
    for(int i=1;i<=hh->GetNbinsX();i++){
        hh->SetBinError(i,0);
    }
}
