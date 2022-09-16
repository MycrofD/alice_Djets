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

TString jetpttitle[5] = {"2-5","5-7", "7-10","10-15","15-50"};

TH1D* LoadRawSpec(TString fn, TString sname, TString spostfix="");
TH1* GetInputHist(TString inFile, string histName,TH1 *hh);
void setHistoDetails(TH1 *hh, Color_t color, Style_t Mstyle, double Msize = 1.1, Width_t Lwidth = 2, Style_t Lstyle = 1);
void setErrorsZero2D(TH2D* hh);
void setErrorsZero1D(TH1D* hh);

TH1D* fRawSpec2Dproj[fJetptbinsN];
//=====================================================================
void subtractFD_zjet(
    TString listName = "",
    TString dataInFileDir="",
    TString outDir="",
    TString simDir="",
    TString dataAnalysisFile="",
    TString MCAnalysisFile="",
    bool isprefix = 0,
    bool ispostfix = 0,
    bool DjetEff=1,
    TString prompteff=""
)
{
    gStyle->SetOptStat(0000);
    gSystem->Exec(Form("mkdir %s",outDir.Data()));
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
            for (int binz = 0; binz < fptbinsZTrueN+1; binz++){//z bins start with numbering 0
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

//3.Extracting the files for charm efficiencies: DONE
//---------------------------------------------
    //access the efficiency folder
    TString promptefffiles[fJetptbinsN];
    TH1D* effHists[fJetptbinsN];
    for(int i=0;i<fJetptbinsN;i++){
        promptefffiles[i]=prompteff+Form("%d_%d.root",(int)fJetptbinsA[i],(int)fJetptbinsA[i+1]);
        TFile* fileEff=new TFile(promptefffiles[i].Data(),"read");
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
            if(!ispostfix){histList[i]=(TList*)dirMC->Get(Form("%s%dFDMCrec",histName.Data(),i));}
            else{cout<<"BLAH, BLAH, BLAH....."<<endl;return;}
        }
        else{
            if(ispostfix) {histList[i]=(TList*)dirMC->Get(Form("%s%sMBN%dFDMCrec",histName.Data(),listName.Data(),i));}
            else{cout<<"BLAH, BLAH, BLAH....."<<endl;return;}
        }
        sparseMC[i] = (THnSparseD*)histList[i]->FindObject("ResponseMatrix");
        sparseMC[i]->GetAxis(4)->SetRangeUser(-0.9+fRpar,0.9-fRpar);//pseudo-rapidity cuts on Reco level
        sparseMC[i]->GetAxis(9)->SetRangeUser(-0.9+fRpar,0.9-fRpar);//pseudo-rapidity cuts on Gen level

        hZjetRGD[i] = (THnSparseD*)sparseMC[i]->Projection(DnDim,Ddim,"E");// z, jet and D components
        if(!i){hZjetRecGenD = (THnSparseD*)hZjetRGD[0]->Clone("hZjetRecGenD");}
        else{hZjetRecGenD->Add(hZjetRGD[i]);}
    }

    //setting up dimensions for response matrix: DONE
    TH2D* hZjetRRebin = new TH2D("hRecRebin","hRecRebin", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA); hZjetRRebin->Sumw2();
    TH2D* hZjetGRebin = new TH2D("hGenRebin","hGenRebin", fptbinsZTrueN, fptbinsZTrueA, fJetptbinsN, fJetptbinsA); hZjetGRebin->Sumw2();
    //setting up kinematic efficiency histos
    TH2D* kineEffCutRec  = new TH2D("hkineEffCutRec","hkineEffCutRec", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA); kineEffCutRec->Sumw2();
    TH2D* kineEffFulRec  = new TH2D("hkineEffFulRec","hkineEffFulRec", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA); kineEffFulRec->Sumw2();
    TH2D* kineEffCutGen  = new TH2D("hkineEffCutGen","hkineEffCutGen", fptbinsZTrueN, fptbinsZTrueA, fJetptbinsN, fJetptbinsA); kineEffCutGen->Sumw2();
    TH2D* kineEffFulGen  = new TH2D("hkineEffFulGen","hkineEffFulGen", fptbinsZTrueN, fptbinsZTrueA, fJetptbinsN, fJetptbinsA); kineEffFulGen->Sumw2();

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

        bool measurement_ok = kTRUE; bool measurement_rec = kTRUE; bool measurement_gen = kTRUE;
        if( DR_center<fDptRangesA[0] || DG_center< fJetptbinsA[0] ){measurement_ok=kFALSE;}//not particularly necessary
        else if( (DG_center<fDptRangesA[0] && jG_center>=fJetptbinsA[0]) || (DR_center<fDptRangesA[0] && jR_center>=fJetptbinsA[0]) ){measurement_ok=kFALSE;}
        else if( (DG_center<fDptRangesA[1] && jG_center>=fJetptbinsA[1]) || (DR_center<fDptRangesA[1] && jR_center>=fJetptbinsA[1]) ){measurement_ok=kFALSE;}
        else if( (DG_center<fDptRangesA[2] && jG_center>=fJetptbinsA[2]) || (DR_center<fDptRangesA[2] && jR_center>=fJetptbinsA[2]) ){measurement_ok=kFALSE;}
        else if( (DG_center<fDptRangesA[3] && jG_center>=fJetptbinsA[3]) || (DR_center<fDptRangesA[3] && jR_center>=fJetptbinsA[3]) ){measurement_ok=kFALSE;}
        else if( (DG_center<fDptRangesA[4] && jG_center>=fJetptbinsA[4]) || (DR_center<fDptRangesA[4] && jR_center>=fJetptbinsA[4]) ){measurement_ok=kFALSE;}

        if (measurement_ok){ kineEffFulRec->Fill(zR_center,jR_center,weight);kineEffFulGen->Fill(zG_center,jG_center,weight);}
        measurement_rec=measurement_gen=measurement_ok;//reseting kine eff booleans
        //pre folding kine
        if(DR_center>fDptRangesAUp[4] || jR_center>fJetptbinsA[fJetptbinsN] || jR_center<fJetptbinsA[0] || zR_center<fptbinsZMeasA[0]){measurement_gen = kFALSE;}
        //if(DR_center>fDptRangesAUp[4] || jR_center>fJetptbinsA[fJetptbinsN]){measurement_gen = kFALSE;}
        if(measurement_gen){ kineEffCutGen->Fill(zG_center,jG_center,weight); }
        //post folding kine
        if(DG_center>fDptRangesAUp[4] || jG_center>fJetptbinsA[fJetptbinsN] || jG_center<fJetptbinsA[0] || zG_center<fptbinsZTrueA[0]){measurement_rec = kFALSE;}
        if(measurement_rec){ kineEffCutRec->Fill(zR_center,jR_center,weight ); }
        //original response
        if(DR_center>fDptRangesAUp[4] || jR_center>fJetptbinsA[fJetptbinsN] || zR_center<fptbinsZMeasA[0]) {measurement_ok=kFALSE;}
        //application of D-jet efficiency to the weight
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
    TCanvas* cKineRec = new TCanvas("cKineRec","cKineRec",800,600);
    TH2D* hKineRec = (TH2D*) kineEffFulRec->Clone();
    hKineRec->Divide(kineEffCutRec, kineEffFulRec,1,1,"B");
    hKineRec->SetTitle("KineEffRec");
    hKineRec->Draw("colz");
    hKineRec->Draw("TEXT same");
    cKineRec  ->SaveAs(Form("%s/cKineRec.pdf",outDir.Data()));
    cKineRec  ->SaveAs(Form("%s/cKineRec.png",outDir.Data()));
    cKineRec  ->SaveAs(Form("%s/cKineRec.svg",outDir.Data()));
    //getting kine efficiencies: DONE
    TCanvas* cKineGen = new TCanvas("cKineGen","cKineGen",800,600);
    TH2D* hKineGen = (TH2D*) kineEffFulGen->Clone();
    hKineGen->Divide(kineEffCutGen, kineEffFulGen,1,1,"B");
    hKineGen->SetTitle("KineEffGen");
    hKineGen->Draw("colz");
    hKineGen->Draw("TEXT same");
    cKineGen  ->SaveAs(Form("%s/cKineGen.pdf",outDir.Data()));
    cKineGen  ->SaveAs(Form("%s/cKineGen.png",outDir.Data()));
    cKineGen  ->SaveAs(Form("%s/cKineGen.svg",outDir.Data()));
    //setting errors to zero: DONE
    for(int j=0;j<=hKineRec->GetNbinsY();j++){for(int i=0;i<=hKineRec->GetNbinsX();i++){hKineRec->SetBinError(i,j,0);}}
    for(int j=0;j<=hKineGen->GetNbinsY();j++){for(int i=0;i<=hKineGen->GetNbinsX();i++){hKineGen->SetBinError(i,j,0);}}
    //apply kine 1 (kinegen)
    TH2D* hFD_centralK = (TH2D*)hFD_central->Clone("hFD_centralK");
    TH2D* hFD_upK      = (TH2D*)hFD_up->Clone("hFD_upK");
    TH2D* hFD_doK      = (TH2D*)hFD_do->Clone("hFD_doK");
    hFD_centralK->Multiply(hKineGen);
    hFD_upK->Multiply(hKineGen);
    hFD_doK->Multiply(hKineGen);
    //folding: DONE
    TH2D* hFolded = new TH2D("hFolded","hFolded", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    TH2D* hFoldUp = new TH2D("hFoldup","hFoldUp", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    TH2D* hFoldDo = new TH2D("hFoldup","hFoldDo", fptbinsZMeasN, fptbinsZMeasA, fJetptbinsN, fJetptbinsA);
    hFolded = (TH2D*)response.ApplyToTruth(hFD_centralK);
    hFoldUp = (TH2D*)response.ApplyToTruth(hFD_upK);
    hFoldDo = (TH2D*)response.ApplyToTruth(hFD_doK);
    //apply kine 2 (kine_rec)
    TH2D* hFoldedK = (TH2D*)hFolded->Clone("hFoldedK");
    TH2D* hFoldUpK = (TH2D*)hFoldUp->Clone("hFoldUpK");
    TH2D* hFoldDoK = (TH2D*)hFoldDo->Clone("hFoldDoK");
    hFoldedK->Divide(hFoldedK,hKineRec);
    hFoldUpK->Divide(hFoldUpK,hKineRec);
    hFoldDoK->Divide(hFoldDoK,hKineRec);
    //saving the response
    TCanvas *cRespon = new TCanvas("","",950,800);
    TH2D* respobj = (TH2D*)response.Hresponse()->Clone();
    cRespon->SetLogz();
    respobj->Draw("colz");//Write
    cRespon->SaveAs(Form("%s/cRespon.pdf",outDir.Data()));
    cRespon->SaveAs(Form("%s/cRespon.svg",outDir.Data()));

//5.Subtracting folded simulation from data: DONE
//-----------------------------------------
    TH2D* hData_sub   = (TH2D*) hData2D->Clone("hData_sub");
    TH2D* hData_subup = (TH2D*) hData2D->Clone("hData_subup");
    TH2D* hData_subdo = (TH2D*) hData2D->Clone("hData_subdo");
    hData_sub->Add(hFoldedK,-1);
    hData_subup->Add(hFoldDoK,-1);
    hData_subdo->Add(hFoldUpK,-1);

    TCanvas* cFDsub = new TCanvas("FDsub","FDsub",1800,1200);
    cFDsub->Divide(3,2);
    TCanvas* cFDratio = new TCanvas("FDratio","FDratio",1800,1000);
    cFDratio->Divide(3,2);
    for (int binjet = 1; binjet<hData_sub->GetNbinsY()+1;binjet++){
        cFDsub->cd(binjet);
        cFDsub->cd(binjet)->SetLogy();
        TH1D* hhdata  =(TH1D*)hData2D->ProjectionX("",binjet,binjet,"E")->Clone(Form("hdata%d",binjet));

        TH1D* hhdata_sub  =(TH1D*)hData_sub  ->ProjectionX("",binjet,binjet,"E")->Clone(Form("hsub_c%d",binjet));
        TH1D* hhdata_subup=(TH1D*)hData_subup->ProjectionX("",binjet,binjet,"E")->Clone(Form("hsub_u%d",binjet));
        TH1D* hhdata_subdo=(TH1D*)hData_subdo->ProjectionX("",binjet,binjet,"E")->Clone(Form("hsub_d%d",binjet));

        TH1D* hhFD_ce  = (TH1D*)hFoldedK->ProjectionX("",binjet,binjet,"E")->Clone(Form("hFD_ce%d",binjet));
        TH1D* hhFD_up = (TH1D*)hFoldUpK->ProjectionX("",binjet,binjet,"E")->Clone(Form("hFD_up%d",binjet));
        TH1D* hhFD_do = (TH1D*)hFoldDoK->ProjectionX("",binjet,binjet,"E")->Clone(Form("hFD_do%d",binjet));

        TH1D* hhData=(TH1D*)hData2D->ProjectionX("",binjet,binjet,"E")->Clone(Form("hData%d",binjet));
        setHistoDetails(hhData,kGreen+2,21,0.6,1,1);
        setHistoDetails(hhdata_sub  ,2,20,0.8,1,1);
        setHistoDetails(hhdata_subup,2,20,0,1,2);setErrorsZero1D(hhdata_subup);//setting errors to zero
        setHistoDetails(hhdata_subdo,2,20,0,1,2);setErrorsZero1D(hhdata_subdo);
        setHistoDetails(hhFD_up,4,20,0,1,2);
        setHistoDetails(hhFD_do,4,20,0,1,2);
        
//        hhData->Scale(1,"width");//clone and scale
//        hhdata_sub  ->Scale(1,"width");
//        hhdata_subup->Scale(1,"width");
//        hhdata_subdo->Scale(1,"width");
//        hhFD_ce->Scale(1,"width");
//        hhFD_up->Scale(1,"width");
//        hhFD_do->Scale(1,"width");

        hhData->GetYaxis()->SetRangeUser(0.8*hhFD_do->GetMinimum(),2.5*hhData->GetMaximum());
        hhData->SetTitle(Form("jetpt %s",jetpttitle[binjet-1].Data()));
        hhData->Draw(); 
        hhdata_sub  ->Draw("same");
        hhdata_subup->Draw("same");
        hhdata_subdo->Draw("same");
        hhFD_ce->Draw("same");
        hhFD_up->Draw("same");
        hhFD_do->Draw("same");

        TLegend *leg1 = new TLegend(0.105,0.75,0.50,0.85);
        leg1->SetBorderSize(0);
        //leg1->SetStyle(0);
        leg1->AddEntry(hhData,"Uncorrected D-Jet yield","p");
        leg1->AddEntry(hhFD_ce,"B Feed-Down (POWHEG+PYTHIA)","p");
        leg1->AddEntry(hhdata_sub,"FD corrected yield","p");
        leg1->Draw("same");

        //ratio canvas cd()        
        cFDratio->cd(binjet);
        TH1D* hhratio = (TH1D*)hhFD_ce->Clone(Form("ratio_%d",binjet));setHistoDetails(hhratio,8,20);
        TH1D* hhratUp = (TH1D*)hhFD_up->Clone(Form("ratUp_%d",binjet));setHistoDetails(hhratUp,8,24,0,2,2);
        TH1D* hhratDo = (TH1D*)hhFD_do->Clone(Form("ratDo_%d",binjet));setHistoDetails(hhratDo,8,24,0,2,2);
        hhratio->Divide(hhdata);hhratUp->Divide(hhdata);hhratDo->Divide(hhdata);
        setErrorsZero1D(hhratUp);setErrorsZero1D(hhratDo);//setting errors to zero
        hhratio->GetYaxis()->SetRangeUser(0,1.15);
        hhratio->SetTitle(Form("jetpt %s",jetpttitle[binjet-1].Data()));
        hhratio->GetYaxis()->SetTitle("Bfd/data");
        hhratio->Draw();
        hhratUp->Draw("same");
        hhratDo->Draw("same");
    }
    cFDsub->SaveAs(Form("%s/cFDsub.pdf",outDir.Data()));
    cFDsub->SaveAs(Form("%s/cFDsub.png",outDir.Data()));
    cFDratio->SaveAs(Form("%s/cFDratio.pdf",outDir.Data()));
    cFDratio->SaveAs(Form("%s/cFDratio.png",outDir.Data()));

//6.Writing histograms to file
    TFile* outFile = new TFile(Form("%s/outFD_%d.root",outDir.Data(),(int)DjetEff),"recreate");
    hData2D->Write();
    hData_sub->Write();
    hData_subup->Write();
    hData_subdo->Write();

    for (int binjet = 1; binjet<hData_sub->GetNbinsY()+1;binjet++){
        TH1D* hhFD_ce  = (TH1D*)hFoldedK->ProjectionX("",binjet,binjet,"E")->Clone(Form("hFD_ce%d",binjet));//FD central
        TH1D* hhFD_up = (TH1D*)hFoldUpK->ProjectionX("",binjet,binjet,"E")->Clone(Form("hFD_up%d",binjet));//FD up 
        TH1D* hhFD_do = (TH1D*)hFoldDoK->ProjectionX("",binjet,binjet,"E")->Clone(Form("hFD_do%d",binjet));//FD down 
        TH1D* hhratio = (TH1D*)hhFD_ce->Clone(Form("ratio_%d",binjet));setHistoDetails(hhratio,8,20);//FD central
        TH1D* hhratUp = (TH1D*)hhFD_up->Clone(Form("ratUp_%d",binjet));setHistoDetails(hhratUp,8,24,0,2,2);
        TH1D* hhratDo = (TH1D*)hhFD_do->Clone(Form("ratDo_%d",binjet));setHistoDetails(hhratDo,8,24,0,2,2);

        TH1D* hhdata  =(TH1D*)hData2D->ProjectionX("",binjet,binjet,"E")->Clone(Form("hdata%d",binjet));//data
        hhratio->Divide(hhdata);hhratUp->Divide(hhdata);hhratDo->Divide(hhdata);
        setErrorsZero1D(hhratUp);setErrorsZero1D(hhratDo);//setting errors to zero
        hhratio->GetYaxis()->SetRangeUser(0,1.15);
        hhratio->SetTitle(Form("jetpt %s",jetpttitle[binjet-1].Data()));
        hhratio->GetYaxis()->SetTitle("Bfd/data");
        hhratio->Write();
        hhratUp->Write();
        hhratDo->Write();

        TH1D* hhdata_sub  =(TH1D*)hData_sub  ->ProjectionX("",binjet,binjet,"E")->Clone(Form("hsub_c%d",binjet));
        TH1D* hhdata_subup=(TH1D*)hData_subup->ProjectionX("",binjet,binjet,"E")->Clone(Form("hsub_u%d",binjet));
        TH1D* hhdata_subdo=(TH1D*)hData_subdo->ProjectionX("",binjet,binjet,"E")->Clone(Form("hsub_d%d",binjet));
        hhdata_sub  ->Write();
        hhdata_subup->Write();
        hhdata_subdo->Write();
    }
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
