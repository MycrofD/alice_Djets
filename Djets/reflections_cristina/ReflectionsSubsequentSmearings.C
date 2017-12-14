#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <TRandom.h>

#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFMassFitter.h"
#endif

enum {kBoth, kParticleOnly, kAntiParticleOnly};
enum {kExpo=0, kLinear, kPol2};
enum {kGaus=0, kDoubleGaus};


Bool_t cutsappliedondistr=kFALSE;//kTRUE;
const Int_t nPtBins=13;
Double_t ptlims[nPtBins+1]={1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,16.,24.,36.,50.};
Int_t rebin[nPtBins]=      {6, 8, 8,10, 8, 10,10,10,10,12, 16, 24,24};

Float_t fixSoverRefAt[nPtBins]={-1.,-1.,-1.,-1.,-1.,-1,-1.,-1.,-1.,-1.,-1.,-1.,-1}; //for data 0-36
Float_t minRangeFit[nPtBins]={1.72,1.72,1.72,1.72,1.72,1.72,1.72,1.72,1.72,1.72,1.68.,1.68,1.68};
Float_t maxRangeFit[nPtBins]={2.05,2.05,2.05,2.05,2.05,2.05,2.05,2.05,2.05,2.05,2.15,2.1,2.1};
Float_t sigmaRef[nPtBins]={0.009,0.010,0.011,0.011,0.013,0.013,0.13,0.014,0.014,0.017,0.021,0.022,0.24};
Float_t massRangeForCounting[nPtBins]={0.03,0.03,0.03,0.035,0.05,0.05,0.05,0.055,0.055,0.06,0.08,0.08,0.08}; // GeV 0-24



Int_t typeb=kExpo;
Int_t types=kGaus;
Int_t optPartAntiPart=kBoth;//kAntiParticleOnly;//Both;
Int_t factor4refl=0;


Int_t rebinRflAn=1;
Float_t bkgFactor=1.0;
Float_t sigFactor=1.0;

Int_t typeb=kExpo;
Int_t types=kGaus;//era gaus
Int_t optPartAntiPart=kBoth;//kAntiParticleOnly;//Both;
Int_t factor4refl=0;

const Int_t nsamples=1;

Bool_t LoadD0toKpiMCHistos(TString fileName,TString fileName1,TString fileName2,TString fileName3, TH1F** hSig, TH1F** hBkg, TH1F** hRfl, const char *CutsType);
Bool_t LoadD0toKpiHistos(TObjArray* listFiles, TH1F** hMass, const char *CutsType);
//void Reflections("AnalysisResults.root","../data15n.root","","","","_5tev_central_topo","_massD0_pass4_central_topoTRUE","_5tev_central_topo"){
void Reflections(TString fileNameMC1="",TString fileNameMC2="",TString fileNameMC3="",TString fileNameMC4="", TString fileNameb="AnalysisResults1406_b.root", TString fileNamec="AnalysisResults1407_c.root", TString fileNamed="AnalysisResults1408_d.root", TString fileNamee="AnalysisResults1409_e.root", const char *CutsTypeMC="_massD0_pass4_central3", const char *CutsType="_pass4_central_kFirst3"){

  gROOT->SetStyle("Plain");	
  gStyle->SetOptTitle(0);

  TH1F** hsig=new TH1F*[nPtBins];
  TH1F** hbkg=new TH1F*[nPtBins];
  TH1F** hrfl=new TH1F*[nPtBins];
  TH1F** hmass=new TH1F*[nPtBins];
  for(Int_t i=0;i<nPtBins;i++){
    hsig[i]=0x0;
    hbkg[i]=0x0;
    hrfl[i]=0x0;
    hmass[i]=0x0;
  }

  Bool_t retCode;
  retCode=LoadD0toKpiMCHistos(fileNameMC1.Data(),fileNameMC2.Data(),fileNameMC3.Data(),fileNameMC4.Data(),hsig,hbkg,hrfl,CutsTypeMC);

  if(!retCode){
    printf("ERROR in reading MC input files\n");
    return;
  }

  TObjArray* listFiles=new TObjArray();
  if(fileNameb!="") listFiles->AddLast(new TObjString(fileNameb.Data()));
  if(fileNamec!="") listFiles->AddLast(new TObjString(fileNamec.Data()));
  if(fileNamed!="") listFiles->AddLast(new TObjString(fileNamed.Data()));
  if(fileNamee!="") listFiles->AddLast(new TObjString(fileNamee.Data()));
  if(listFiles->GetEntries()==0){
    printf("Missing file names in input\n");
    return;
  }

  Float_t massD;
  Bool_t retCodeData;
  retCodeData=LoadD0toKpiHistos(listFiles,hmass,CutsType);
  massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  if(!retCodeData){
    printf("ERROR in reading data input files\n");
    return;
  } 
  

  TH1D* hCntSig1=new TH1D("hCntSig1","hCntSig1",nPtBins,ptlims);
  TH1D* hCntSig2=new TH1D("hCntSig2","hCntSig2",nPtBins,ptlims);
  TH1D* hNDiffCntSig1=new TH1D("hNDiffCntSig1","hNDiffCntSig1",nPtBins,ptlims);
  TH1D* hNDiffCntSig2=new TH1D("hNDiffCntSig2","hNDiffCntSig2",nPtBins,ptlims);
  TH1D* hSignal=new TH1D("hSignal","hSignal",nPtBins,ptlims);
  TH1D* hRelErrSig=new TH1D("hRelErrSig","hRelErrSig",nPtBins,ptlims);
  TH1D* hInvSignif=new TH1D("hInvSignif","hInvSignif",nPtBins,ptlims);
  TH1D* hBackground=new TH1D("hBackground","hBackground",nPtBins,ptlims);
  TH1D* hBackgroundNormSigma=new TH1D("hBackgroundNormSigma","hBackgroundNormSigma",nPtBins,ptlims);
  TH1D* hSignificance=new TH1D("hSignificance","hSignificance",nPtBins,ptlims);
  TH1D* hMass=new TH1D("hMass","hMass",nPtBins,ptlims);
  TH1D* hSigma=new TH1D("hSigma","hSigma",nPtBins,ptlims);

  Int_t nMassBins=hmass[1]->GetNbinsX();
  Double_t hmin=hmass[1]->GetBinLowEdge(3);
  Double_t hmax=hmass[1]->GetBinLowEdge(nMassBins-2)+hmass[1]->GetBinWidth(nMassBins-2);
  Float_t minBinSum1=hmass[1]->FindBin(massD-massRangeForCounting[1]);
  Float_t maxBinSum1=hmass[1]->FindBin(massD+massRangeForCounting[1]);
  Int_t iPad=1;
   
  TF1* funBckStore1=0x0;
  TF1* funBckStore2=0x0;
  TF1* funBckStore3=0x0;

  TF1** fB1Array=new TF1*[nPtBins];

  AliHFMassFitter** fitter=new AliHFMassFitter*[nPtBins];
  Double_t arrchisquare[nPtBins];
  TCanvas* c1= new TCanvas("c1","MassSpectra",1000,800);
  Int_t nx=5, ny=3;
  if(nPtBins==12){
    nx=4;
    ny=3;
  }
  
  c1->Divide(nx,ny);

  TCanvas *myCanvas[nPtBins];
  TPaveText* ptBin[nPtBins];

  Double_t sig,errsig,s,errs,b,errb;
  for(Int_t iBin=0; iBin<nPtBins; iBin++){

    minBinSum1=hmass[iBin]->FindBin(massD-massRangeForCounting[iBin]);
    maxBinSum1=hmass[iBin]->FindBin(massD+massRangeForCounting[iBin]);

    c1->cd(iPad);
    gPad->SetTicks();
    
    
    Int_t origNbins=hmass[iBin]->GetNbinsX();
    hmass[iBin]->GetXaxis()->SetTitle("Invariant Mass (K#pi) (GeV/c^{2})");
    hmass[iBin]->GetYaxis()->SetTitle(Form("Entries / %1.0f MeV/c^{2}", (1000*(hmass[iBin]->GetXaxis()->GetBinWidth(1))*rebin[iBin])));
    fitter[iBin]=new AliHFMassFitter(hmass[iBin],hmin,hmax,rebin[iBin],typeb,types);
   
    fitter[iBin]->SetRangeFit(1.72,2.05);
    
   
    if (iBin==nPtBins-1 /*|| iBin==nPtBins-2*/) {
      //default
      fitter[iBin]->SetRangeFit(1.65,2.15);
    }

    if (iBin==nPtBins-2) {
      //default
      fitter[iBin]->SetRangeFit(1.68,2.15);
    }

    rebin[iBin]=origNbins/fitter[iBin]->GetBinN();
    fitter[iBin]->SetReflectionSigmaFactor(factor4refl);
    fitter[iBin]->SetInitialGaussianMean(massD);
    Bool_t out=fitter[iBin]->MassFitter(0);
    if(!out) {
      fitter[iBin]->GetHistoClone()->Draw();
      continue;
    }
    Double_t mass=fitter[iBin]->GetMean();
    Double_t massErr=fitter[iBin]->GetMeanUncertainty();
    Double_t sigma=fitter[iBin]->GetSigma();
    Double_t sigmaErr=fitter[iBin]->GetSigmaUncertainty();
 
    Float_t halfRange=3*sigma;
    Float_t minBinSum=hmass[iBin]->FindBin(massD-halfRange);
    Float_t maxBinSum=hmass[iBin]->FindBin(massD+halfRange);
    Printf("minBinSum (5sigma)=%f\t minBinSum (100MeV)=%f", minBinSum, minBinSum1);
    Printf("maxBinSum (5sigma)=%f\t maxBinSum (100MeV)=%f", maxBinSum, maxBinSum1);

    arrchisquare[iBin]=fitter[iBin]->GetReducedChiSquare();
    TF1* fB1=fitter[iBin]->GetBackgroundFullRangeFunc();
    fB1Array[iBin]=(TF1*)fitter[iBin]->GetBackgroundFullRangeFunc();
    TF1* fB2=fitter[iBin]->GetBackgroundRecalcFunc();
    TF1* fM=fitter[iBin]->GetMassFunc();
    if(iBin==0 && fB1) funBckStore1=new TF1(*fB1);
    if(iBin==0 && fB2) funBckStore2=new TF1(*fB2);
    if(iBin==0 && fM) funBckStore3=new TF1(*fM);

    fitter[iBin]->DrawHere(gPad);
    fitter[iBin]->Signal(3,s,errs);
    fitter[iBin]->Background(3,b,errb);
    fitter[iBin]->Significance(3,sig,errsig);
    Double_t ry=fitter[iBin]->GetRawYield();
    Double_t ery=fitter[iBin]->GetRawYieldError();

    ptBin[iBin] = new TPaveText(0.6,0.55,0.8,0.6,"NDC");
    ptBin[iBin]->SetFillStyle(0);
    ptBin[iBin]->SetBorderSize(0);
    ptBin[iBin]->AddText(0.,0.,Form("%1.1f<p_{T}<%1.1f",ptlims[iBin],ptlims[iBin+1]));
    ptBin[iBin]->SetTextFont(42);
    ptBin[iBin]->SetTextSize(0.08);
    gPad->cd();
    ptBin[iBin]->Draw();
    gPad->Update();

    myCanvas[iBin] = new TCanvas(Form("myCanvas_%d",iBin),Form("Invariant mass pt bin %d",iBin));
    fitter[iBin]->DrawHere(gPad);
    gPad->SetTicks();
    gPad->cd();
    ptBin[iBin]->Draw();
    gPad->Update();

    Float_t cntSig1=0.;
    Float_t cntSig2=0.;
    Float_t cntErr=0.;
    for(Int_t iMB=minBinSum1; iMB<=maxBinSum1; iMB++){
      Float_t bkg1=fB1 ? fB1->Eval(hmass[iBin]->GetBinCenter(iMB))/rebin[iBin] : 0;
      Float_t bkg2=fB2 ? fB2->Eval(hmass[iBin]->GetBinCenter(iMB))/rebin[iBin] : 0;
      cntSig1+=(hmass[iBin]->GetBinContent(iMB)-bkg1);
      cntSig2+=(hmass[iBin]->GetBinContent(iMB)-bkg2);
      cntErr+=(hmass[iBin]->GetBinContent(iMB));
    }
    hCntSig1->SetBinContent(iBin+1,cntSig1);
    hCntSig1->SetBinError(iBin+1,TMath::Sqrt(cntErr));
    hNDiffCntSig1->SetBinContent(iBin+1,(ry-cntSig1)/ry);
    hNDiffCntSig1->SetBinError(iBin+1,TMath::Sqrt(cntErr)/ry);
    hCntSig2->SetBinContent(iBin+1,cntSig2);
    hNDiffCntSig2->SetBinContent(iBin+1,(ry-cntSig2)/ry);
    hNDiffCntSig2->SetBinError(iBin+1,TMath::Sqrt(cntErr)/ry);
    hCntSig2->SetBinError(iBin+1,TMath::Sqrt(cntErr));
    hSignal->SetBinContent(iBin+1,ry);
    hSignal->SetBinError(iBin+1,ery);
    hRelErrSig->SetBinContent(iBin+1,errs/s);
    hInvSignif->SetBinContent(iBin+1,1/sig);
    hInvSignif->SetBinError(iBin+1,errsig/(sig*sig));
    hBackground->SetBinContent(iBin+1,b); //consider sigma
    hBackground->SetBinError(iBin+1,errb);
    hBackgroundNormSigma->SetBinContent(iBin+1,b/(3*fitter[iBin]->GetSigma())*(3*0.012)); //consider sigma
    hBackgroundNormSigma->SetBinError(iBin+1,errb);
    hSignificance->SetBinContent(iBin+1,sig);
    hSignificance->SetBinError(iBin+1,errsig);
    hMass->SetBinContent(iBin+1,mass);
    hMass->SetBinError(iBin+1,massErr);
    hSigma->SetBinContent(iBin+1,sigma);
    hSigma->SetBinError(iBin+1,sigmaErr);

    iPad++;
    
  } // end loop on pt bins

  // toy invariant mass distributions for reflections study

  const Int_t nSmearing=1;
  
  Float_t intS[nPtBins];
  Float_t intSMC[nPtBins];
  Float_t intRMC[nPtBins];
  for(Int_t ipt=0;ipt<nPtBins;ipt++){ 
    intS[ipt]=0.;
    intSMC[ipt]=0.;
    intRMC[ipt]=0.;
  }
  Float_t intRMCSmeared[nPtBins][nSmearing];
  Float_t intSMCSmeared[nPtBins][nSmearing];
  for(Int_t i=0; i<nSmearing; i++){
    for(Int_t ipt=0;ipt<nPtBins;ipt++){ 
      intRMCSmeared[ipt][i]=0.;
      intSMCSmeared[ipt][i]=0.;
    }
  }

  TH1D **hSignalInjected=new TH1D*[nSmearing];
  TH1D **hRflInjected=new TH1D*[nSmearing];

  Printf("initializiang sig and rfl injected");
  for(Int_t i=0; i<nSmearing; i++){
    hSignalInjected[i]=new TH1D(Form("hSignalInjected_stepSmearing%d",i),"MC signal scaled and smeared;p_{T} (GeV/c);S injected",nPtBins, ptlims);
    hRflInjected[i]=new TH1D(Form("hReflectionsInjected_stepSmearing%d",i),"MC reflections scaled and smeared;p_{T} (GeV/c);R injected",nPtBins, ptlims);
  }
  Printf("end initialization sig and rfl injected");  

  TH1D** hSigData=new TH1D*[nPtBins];
  TH1D* hBkgDataSmeared[nPtBins][nSmearing];
  TH1D* hSigMCScaledSmeared[nPtBins][nSmearing];
  TH1D** hSigMCScaledSmearedToDraw=new TH1D*[nPtBins];
  TH1D* hRflMCScaledSmeared[nPtBins][nSmearing];
  TH1D* hMassToy[nPtBins][nSmearing];
  TH1D* hMassToyNoRfl[nPtBins][nSmearing];
  for(Int_t i=0; i<nSmearing; i++){
    for(Int_t ipt=0;ipt<nPtBins;ipt++){ 
      hMassToy[ipt][i]=new TH1D(Form("hMT_ipt%d_iSmearing%d",ipt,i),Form("hMT_ipt%d_iSmearing%d",ipt,i),1,0.,1.);
      hBkgDataSmeared[ipt][i]=new TH1D(Form("hB_ipt%d_iSmearing%d",ipt,i),Form("hB_ipt%d_iSmearing%d",ipt,i),1,0.,1.);
      hSigMCScaledSmeared[ipt][i]=new TH1D(Form("hS_ipt%d_iSmearing%d",ipt,i),Form("hS_ipt%d_iSmearing%d",ipt,i),1,0.,1.);
      hRflMCScaledSmeared[ipt][i]=new TH1D(Form("hR_ipt%d_iSmearing%d",ipt,i),Form("hR_ipt%d_iSmearing%d",ipt,i),1,0.,1.);

    }
  }

  TH1D* hRflOverSignalMC=new TH1D("hRflOverSignalMC","hRflOverSignalMC;p_{T} (GeV/c);R/S",nPtBins,ptlims);
  TH1D* hSignalMC=new TH1D("hSignalMC","hSignalMC;p_{T} (GeV/c);S",nPtBins,ptlims);
  TH1D* hReflectionsMC=new TH1D("hReflectionsMC","hReflectionsMC;p_{T} (GeV/c);R",nPtBins,ptlims);

  TCanvas* c8= new TCanvas("c8","Reflections MC",900,800);
  c8->Divide(nx,ny);
  TCanvas* c9= new TCanvas("c9","Signal MC",900,800);
  c9->Divide(nx,ny);
  TCanvas* c10= new TCanvas("c10","Signal MC from fit",900,800);
  c10->Divide(nx,ny);
  TCanvas* c11= new TCanvas("c11","Signal MC from fit normalized",900,800);
  c11->Divide(nx,ny);
  TCanvas* c12= new TCanvas("c12","Total without reflections",900,800);
  c12->Divide(nx,ny);
  TCanvas* c13= new TCanvas("c13","MC Ratio Reflections/Signal",900,800);
  TCanvas* c14= new TCanvas("c14","MC Reflections and Signal",900,800);

  Int_t zero=0;
  TFile* outf3=new TFile(Form("ToyInvariantMassDistributions%s_StepSmearing%d.root",CutsType,zero),"recreate");
  //outf3->cd();

  TCanvas* c2= new TCanvas("c2","Background",1000,800);
  TCanvas* c3= new TCanvas("c3","Signal",1000,800);
  TCanvas* c4= new TCanvas("c4","Reflections",1000,800);
  TCanvas* c5= new TCanvas("c5","Total",1000,800);
  TCanvas* c6= new TCanvas("c6","Signal Injected",900,800);
  TCanvas* c7= new TCanvas("c7","Reflections Injected",900,800);
  Printf("Nx = %d, Ny = %d", nx, ny);
  c2->Divide(nx,ny);
  c3->Divide(nx,ny);
  c4->Divide(nx,ny);
  c5->Divide(nx,ny);

  TLegend *leg = new TLegend(0.73,0.62,0.87,0.89, "", "brNDC");
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetLineColor(kWhite);

  Int_t jPad=1;

  TString formula="[0]/([2]*TMath::Sqrt(2*TMath::Pi()))*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))";

  for(Int_t ipt=0;ipt<nPtBins;ipt++){ 

    hmass[ipt]->Rebin(rebinRflAn);
    hsig[ipt]->Rebin(rebinRflAn);
    hrfl[ipt]->Rebin(rebinRflAn);

    Printf("Pt Bin %d", ipt);

    Printf("Initializiang hSigData[%d][%d]",ipt,zero);
    hSigData[ipt]=(TH1D*)hmass[ipt]->Clone(Form("hSigData_ptBin%d",ipt));
    hSigData[ipt]->SetDirectory(0);
    hSigData[ipt]->Reset("ICES");
    hSigData[ipt]->Sumw2();
    if(hSigData[ipt]) {Printf("End initialization hSigData"); hSigData[ipt]->Print();}

    Printf("Initializiang hBkgDataSmeared[%d][%d]",ipt,zero);
    hBkgDataSmeared[ipt][0]=(TH1D*)hmass[ipt]->Clone(Form("hBkgData_ptBin%d_stepSmearing%d",ipt,zero));
    hBkgDataSmeared[ipt][0]->SetDirectory(0);
    hBkgDataSmeared[ipt][0]->Reset("ICES");
    hBkgDataSmeared[ipt][0]->Sumw2();
    hBkgDataSmeared[ipt][0]->GetXaxis()->SetTitle("Invariant Mass (K#pi) (GeV/c^{2})");
    hBkgDataSmeared[ipt][0]->GetYaxis()->SetTitle(Form("Entries / %1.0f MeV/c^{2}", (1000*(hBkgDataSmeared[ipt][0]->GetXaxis()->GetBinWidth(1)))));
    if(hBkgDataSmeared[ipt][0]) {Printf("End initialization hBkgDataSmeared"); hBkgDataSmeared[ipt][0]->Print();}
    //c8->cd(ipt+1);
    //hBkgDataSmeared[ipt][0]->Draw("hist");

    Printf("Initializiang hSigMCScaledSmeared[%d][%d]",ipt,zero);
    hSigMCScaledSmeared[ipt][0]=(TH1D*)hsig[ipt]->Clone(Form("hSigMCScaled_ptBin%d_stepSmearing%d",ipt,zero));
    hSigMCScaledSmeared[ipt][0]->SetDirectory(0);
    hSigMCScaledSmeared[ipt][0]->Reset("ICES");
    hSigMCScaledSmeared[ipt][0]->Sumw2();
    hSigMCScaledSmeared[ipt][0]->GetXaxis()->SetTitle("Invariant Mass (K#pi) (GeV/c^{2})");
    hSigMCScaledSmeared[ipt][0]->GetYaxis()->SetTitle(Form("Entries / %1.0f MeV/c^{2}", (1000*(hSigMCScaledSmeared[ipt][0]->GetXaxis()->GetBinWidth(1)))));
    if(hSigMCScaledSmeared[ipt][0]) {Printf("end initialization hSigMCScaledSmeared"); hSigMCScaledSmeared[ipt][0]->Print();}
    //c8->cd(ipt+1);
    //hSigMCScaledSmeared[ipt][0]->Draw("same hist");

    Printf("Initializiang hRflMCScaledSmeared[%d][%d]",ipt,zero);
    hRflMCScaledSmeared[ipt][0]=(TH1D*)hrfl[ipt]->Clone(Form("hRflMCScaled_ptBin%d_stepSmearing%d",ipt,zero));
    hRflMCScaledSmeared[ipt][0]->SetDirectory(0);
    hRflMCScaledSmeared[ipt][0]->Sumw2();
    hRflMCScaledSmeared[ipt][0]->GetXaxis()->SetTitle("Invariant Mass (K#pi) (GeV/c^{2})");
    hRflMCScaledSmeared[ipt][0]->GetYaxis()->SetTitle(Form("Entries / %1.0f MeV/c^{2}", (1000*(hRflMCScaledSmeared[ipt][0]->GetXaxis()->GetBinWidth(1)))));
    if(hRflMCScaledSmeared[ipt][0]) {Printf("End initialization hRflMCScaledSmeared"); hRflMCScaledSmeared[ipt][0]->Print();}
    c8->cd(ipt+1);
    hrfl[ipt]->Draw();

    Int_t nmassbins=hBkgDataSmeared[ipt][0]->GetNbinsX();
    Printf("nMassBins = %d, rebin[1] = %d", nmassbins, rebin[1]);
    for (Int_t mBin=1; mBin<=nmassbins; mBin++){
      //filling BKG and S histograms - DATA
      if(!fB1Array[ipt]) Printf("NO f BKGROUND");
      Float_t bkg=fB1Array[ipt] ? (Float_t)((fB1Array[ipt]->Eval(hSigData[ipt]->GetBinCenter(mBin)))/rebin[ipt])*rebinRflAn : 0.;
      Float_t bkgErr=0.;
      Float_t signal=0.;
      Float_t signalErr=0.;
      if(bkg!=0) {
	signal=(hmass[ipt]->GetBinContent(mBin))-bkg;
	signalErr=TMath::Sqrt(signal);
	bkgErr=TMath::Sqrt(bkg);
      }
      hSigData[ipt]->SetBinContent(mBin,signal);
      hSigData[ipt]->SetBinError(mBin,signalErr);
      hBkgDataSmeared[ipt][0]->SetBinContent(mBin,bkg);
      hBkgDataSmeared[ipt][0]->SetBinError(mBin,bkgErr);
      }
    if(hBkgDataSmeared[ipt][0]){Printf("End filling hSigData and hBkGDataSmeared[%d][0]",ipt); hBkgDataSmeared[ipt][0]->Print();}

    //fitting MC signal histo antake the shape as template
    TF1 *gaussMCSignal=new TF1("gaussMCSig",formula.Data(),1.65,2.15);
    gaussMCSignal->SetParName(0,"IntegralSgn");
    gaussMCSignal->SetParName(1,"Mean");
    gaussMCSignal->SetParName(2,"Sigma");
    gaussMCSignal->SetParameter(0,1);
    gaussMCSignal->SetParameter(1,1.864);
    gaussMCSignal->SetParameter(2,0.010);
    gaussMCSignal->SetLineColor(kOrange+2);
    c9->cd(ipt+1);
    gStyle->SetOptFit(11111);
    hsig[ipt]->Fit("gaussMCSig","RI","",1.65,2.15);
    hsig[ipt]->Draw();

    Int_t nmassbinsSig=hSigMCScaledSmeared[ipt][0]->GetNbinsX();
    Printf("nMassBins = %d, rebin = %d", nmassbinsSig, rebinRflAn);
    for (Int_t mBin=1; mBin<=nmassbinsSig; mBin++){
      if(!gaussMCSignal) Printf("NO f Signal");
      Float_t signal=gaussMCSignal ? (Float_t)(gaussMCSignal->Eval(hSigMCScaledSmeared[ipt][0]->GetBinCenter(mBin)))/rebinRflAn : 0.;
      Float_t signalErr=TMath::Sqrt(signal);
      hSigMCScaledSmeared[ipt][0]->SetBinContent(mBin,signal);
      hSigMCScaledSmeared[ipt][0]->SetBinError(mBin,signalErr);
      }
    hSigMCScaledSmearedToDraw[ipt]=(TH1D*)hSigMCScaledSmeared[ipt][0]->Clone(Form("hSigMCFromFitNotScaled_ptBin%d_stepSmearing%d",ipt,zero));
    hSigMCScaledSmearedToDraw[ipt]->SetDirectory(0);
    hSigMCScaledSmearedToDraw[ipt]->SetLineColor(kRed);
    c10->cd(ipt+1);
    hSigMCScaledSmearedToDraw[ipt]->Draw("hist");

   
    /*if(ipt==nPtBins-3){
      hSigData[ipt]->GetXaxis()->SetRangeUser(1.68,2.15);
      hBkgDataSmeared[ipt][0]->GetXaxis()->SetRangeUser(1.68,2.15);
      hSigMCScaledSmeared[ipt][0]->GetXaxis()->SetRangeUser(1.68,2.15);
      hRflMCScaledSmeared[ipt][0]->GetXaxis()->SetRangeUser(1.68,2.15);
      }
      else */if(ipt==nPtBins-1 || ipt==nPtBins-2){
      hSigData[ipt]->GetXaxis()->SetRangeUser(1.65,2.15);
      hBkgDataSmeared[ipt][0]->GetXaxis()->SetRangeUser(1.65,2.15);
      hSigMCScaledSmeared[ipt][0]->GetXaxis()->SetRangeUser(1.65,2.15);
      hRflMCScaledSmeared[ipt][0]->GetXaxis()->SetRangeUser(1.65,2.15);
    }
    else{
      hSigData[ipt]->GetXaxis()->SetRangeUser(1.72,2.05);
      hBkgDataSmeared[ipt][0]->GetXaxis()->SetRangeUser(1.72,2.05);
      hSigMCScaledSmeared[ipt][0]->GetXaxis()->SetRangeUser(1.72,2.05);
      hRflMCScaledSmeared[ipt][0]->GetXaxis()->SetRangeUser(1.72,2.05);
    }
    Printf("End setting ranges");

    hBkgDataSmeared[ipt][0]->Scale(bkgFactor);
    hSigData[ipt]->Scale(sigFactor);

    //signal - DATA
    intS[ipt]=hSigData[ipt]->Integral();
    Printf("Integral signal data pt bin %d = %f", ipt, intS[ipt]);

    //signal - MC
    if(ipt==nPtBins-1 || ipt==nPtBins-2) {
      intSMC[ipt]=hsig[ipt]->Integral(hsig[ipt]->FindBin(1.65),hsig[ipt]->FindBin(2.15));
      intRMC[ipt]=hrfl[ipt]->Integral(hrfl[ipt]->FindBin(1.65),hrfl[ipt]->FindBin(2.15));
    }
    /*else if (ipt==nPtBins-3) {
      intSMC[ipt]=hsig[ipt]->Integral(hsig[ipt]->FindBin(1.68),hsig[ipt]->FindBin(2.15));
      intRMC[ipt]=hrfl[ipt]->Integral(hrfl[ipt]->FindBin(1.68),hrfl[ipt]->FindBin(2.15));
      }*/
    else{
      intSMC[ipt]=hsig[ipt]->Integral(hsig[ipt]->FindBin(1.72),hsig[ipt]->FindBin(2.05));
      intRMC[ipt]=hrfl[ipt]->Integral(hrfl[ipt]->FindBin(1.72),hrfl[ipt]->FindBin(2.05));
    }
    Printf("Integral signal MC pt bin %d = %f", ipt, intSMC[ipt]);

    //reflections over signal MC ratio vs. pT
    hRflOverSignalMC->SetBinContent(ipt+1,intRMC[ipt]/intSMC[ipt]);
    hSignalMC->SetBinContent(ipt+1,intSMC[ipt]);
    hReflectionsMC->SetBinContent(ipt+1,intRMC[ipt]);
    
    //scaling signal histogram - MC
    intSMCSmeared[ipt][0]=hSigMCScaledSmeared[ipt][0]->Integral();
    hSigMCScaledSmeared[ipt][0]->Scale(intS[ipt]/intSMCSmeared[ipt][0]);

    c11->cd(ipt+1);
    hSigMCScaledSmeared[ipt][0]->Draw("hist");
    
    hRflMCScaledSmeared[ipt][0]->Scale(intS[ipt]/intSMCSmeared[ipt][0]);

    intRMCSmeared[ipt][0]=hRflMCScaledSmeared[ipt][0]->Integral();
    Printf("Integral rfl MC pt bin %d = %f", ipt, intRMCSmeared[ipt][0]);
    
    hSignalInjected[0]->SetBinContent(ipt+1,intS[ipt]);
    hRflInjected[0]->SetBinContent(ipt+1,intRMCSmeared[ipt][0]);
    hMassToy[ipt][0]=(TH1D*)hBkgDataSmeared[ipt][0]->Clone(Form("hMass_ptBin%d_stepSmearing%d",ipt,zero));
    hMassToy[ipt][0]->SetDirectory(0);
    hMassToy[ipt][0]->Add(hSigMCScaledSmeared[ipt][0],1);
    hMassToy[ipt][0]->Add(hRflMCScaledSmeared[ipt][0],1);
    hMassToyNoRfl[ipt][0]=(TH1D*)hBkgDataSmeared[ipt][0]->Clone(Form("hMassNoRfl_ptBin%d_stepSmearing%d",ipt,zero));
    hMassToyNoRfl[ipt][0]->SetDirectory(0);
    hMassToyNoRfl[ipt][0]->Add(hSigMCScaledSmeared[ipt][0],1);

    Printf("START DRAWING");
    Printf("jPAD = %d", jPad);
   
    hBkgDataSmeared[ipt][0]->SetMarkerStyle(20);
    hBkgDataSmeared[ipt][0]->SetMarkerColor(kBlack);
    hBkgDataSmeared[ipt][0]->SetLineColor(kBlack);
    hBkgDataSmeared[ipt][0]->SetLineWidth(1);
    if(ipt==0) hBkgDataSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*1400.);
    if(ipt==1) hBkgDataSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*350.);
    if(ipt==2) hBkgDataSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*120.);
    if(ipt==3) hBkgDataSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*220.);
    if(ipt==4) hBkgDataSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*140.);
    if(ipt==5) hBkgDataSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*120.);
    if(ipt==6) hBkgDataSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*40.);
    if(ipt==7) hBkgDataSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*70.);
    if(ipt==8) hBkgDataSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*50.);
    if(ipt==9) hBkgDataSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*20.);
    c2->cd(jPad);
    gPad->SetTicks();
    Printf("c2");
    hBkgDataSmeared[ipt][0]->Draw("hist");

    hSigMCScaledSmeared[ipt][0]->SetMarkerStyle(20);
    hSigMCScaledSmeared[ipt][0]->SetMarkerColor(kBlack);
    hSigMCScaledSmeared[ipt][0]->SetLineColor(kBlack);
    hSigMCScaledSmeared[ipt][0]->SetLineWidth(1);
    if(ipt==0) hSigMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*230.);
    if(ipt==1) hSigMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*280.);
    if(ipt==2) hSigMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*230.);
    if(ipt==3) hSigMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*230.);
    if(ipt==4) hSigMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*140.);
    if(ipt==5) hSigMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*100.);
    if(ipt==6) hSigMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*50.);
    if(ipt==7) hSigMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*90.);
    if(ipt==8) hSigMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*40.);
    if(ipt==9) hSigMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*15.);
    c3->cd(jPad);
    gPad->SetTicks();
    Printf("c3");
    hSigMCScaledSmeared[ipt][0]->Draw("hist");

    hRflMCScaledSmeared[ipt][0]->SetMarkerStyle(20);
    hRflMCScaledSmeared[ipt][0]->SetMarkerColor(kBlack);
    hRflMCScaledSmeared[ipt][0]->SetLineColor(kBlack);
    hRflMCScaledSmeared[ipt][0]->SetLineWidth(1);
    if(ipt==0) hRflMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*10.);
    if(ipt==1) hRflMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*15.);
    if(ipt==2) hRflMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*15.);
    if(ipt==3) hRflMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*20.);
    if(ipt==4) hRflMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*15.);
    if(ipt==5) hRflMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*10.);
    if(ipt==6) hRflMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*10.);
    if(ipt==7) hRflMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*20.);
    if(ipt==8) hRflMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*10.);
    if(ipt==9) hRflMCScaledSmeared[ipt][0]->GetYaxis()->SetRangeUser(0.,sigFactor*rebinRflAn*8.);
    c4->cd(jPad);
    gPad->SetTicks();
    Printf("c4");
    hRflMCScaledSmeared[ipt][0]->Draw("hist");

    hMassToy[ipt][0]->SetMarkerStyle(20);
    hMassToy[ipt][0]->SetMarkerColor(kBlack);
    hMassToy[ipt][0]->SetLineColor(kBlack);
    hMassToy[ipt][0]->SetLineWidth(1);
    if(ipt==0) hMassToy[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*1400.);
    if(ipt==1) hMassToy[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*500.);
    if(ipt==2) hMassToy[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*300.);
    if(ipt==3) hMassToy[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*400.);
    if(ipt==4) hMassToy[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*250.);
    if(ipt==5) hMassToy[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*160.);
    if(ipt==6) hMassToy[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*70.);
    if(ipt==7) hMassToy[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*140.);
    if(ipt==8) hMassToy[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*60.);
    if(ipt==9) hMassToy[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*20.);
    c5->cd(jPad);
    gPad->SetTicks();
    Printf("c5");
    hMassToy[ipt][0]->Draw("hist");

    hMassToyNoRfl[ipt][0]->SetMarkerStyle(20);
    hMassToyNoRfl[ipt][0]->SetMarkerColor(kBlack);
    hMassToyNoRfl[ipt][0]->SetLineColor(kBlack);
    hMassToyNoRfl[ipt][0]->SetLineWidth(1);
    if(ipt==0) hMassToyNoRfl[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*1400.);
    if(ipt==1) hMassToyNoRfl[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*500.);
    if(ipt==2) hMassToyNoRfl[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*300.);
    if(ipt==3) hMassToyNoRfl[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*400.);
    if(ipt==4) hMassToyNoRfl[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*250.);
    if(ipt==5) hMassToyNoRfl[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*160.);
    if(ipt==6) hMassToyNoRfl[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*70.);
    if(ipt==7) hMassToyNoRfl[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*140.);
    if(ipt==8) hMassToyNoRfl[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*60.);
    if(ipt==9) hMassToyNoRfl[ipt][0]->GetYaxis()->SetRangeUser(0.,bkgFactor*rebinRflAn*20.);
    c12->cd(jPad);
    gPad->SetTicks();
    Printf("c12");
    hMassToyNoRfl[ipt][0]->Draw("hist");
    
    jPad++;

    Printf("START WRITING");
    outf3->cd();    
    hBkgDataSmeared[ipt][0]->Write();
    hSigMCScaledSmeared[ipt][0]->Write();
    hRflMCScaledSmeared[ipt][0]->Write();
    hMassToy[ipt][0]->Write();
    hMassToyNoRfl[ipt][0]->Write();

  } //end loop on pt bins
  c6->cd();
  c6->SetTicks();
  Printf("c6");
  hSignalInjected[0]->SetMarkerStyle(20);
  hSignalInjected[0]->SetMarkerColor(kBlack);
  hSignalInjected[0]->SetLineColor(kBlack);
  hSignalInjected[0]->SetLineWidth(2);
  hSignalInjected[0]->Draw("hist");


  c7->cd();
  c7->SetTicks();
  Printf("c7");
  hRflInjected[0]->SetMarkerStyle(20);
  hRflInjected[0]->SetMarkerColor(kBlack);
  hRflInjected[0]->SetLineColor(kBlack);
  hRflInjected[0]->SetLineWidth(2);
  hRflInjected[0]->Draw("hist");

  c13->cd();
  c13->SetTicks();
  hRflOverSignalMC->SetLineColor(kRed);
  hRflOverSignalMC->SetLineWidth(2);
  hRflOverSignalMC->Draw("hist");
  c14->cd();
  c14->SetTicks();
  hSignalMC->SetLineColor(kBlue);
  hSignalMC->SetLineWidth(2);
  hSignalMC->Draw("hist");
  hReflectionsMC->SetLineColor(kGreen+1);
  hReflectionsMC->SetLineWidth(2);
  hReflectionsMC->Draw("same hist");

  Printf("END LOOP TO FILL STEP SMEARING = 0 HISTOS");
  
  outf3->cd();
  hSignalInjected[0]->Write();
  hRflInjected[0]->Write();
  hRflOverSignalMC->Write();
  hSignalMC->Write();
  hReflectionsMC->Write();
  outf3->Close();
  
  Printf("STEP SMEARING = 0 HISTOS WRITTEN");  
  
  Printf("Loop on smearing steps");
  for(Int_t iSmearing=1; iSmearing<nSmearing; iSmearing++){ //loop on smearing steps
    
    TFile* outf2=new TFile(Form("ToyInvariantMassDistributions%s_StepSmearing%d.root",CutsType,iSmearing),"recreate");
    if(!outf2) Printf("FAILED FILE CREATION");

    Printf("Loop on pt bins");
    for(Int_t iptBin=0; iptBin<nPtBins; iptBin++){ //loop on pt bins

      TRandom smearingGenerator;

      //smearing BKG - DATA 
      Int_t nMassBinsBkg=hBkgDataSmeared[iptBin][0]->GetNbinsX();
      Printf("First loop on mass bins, %d", nMassBinsBkg);
      Printf("Smearing hBkgDataSmeared[%d][%d]",iptBin,iSmearing);
      hBkgDataSmeared[iptBin][iSmearing]=(TH1D*)hBkgDataSmeared[iptBin][0]->Clone(Form("hBkgData_ptBin%d_stepSmearing%d",iptBin,iSmearing));
      hBkgDataSmeared[iptBin][iSmearing]->SetDirectory(0);
      if(hBkgDataSmeared[iptBin][iSmearing]) hBkgDataSmeared[iptBin][iSmearing]->Print();
      for (Int_t mBin=1; mBin<=nMassBinsBkg; mBin++){
	Float_t bkg=hBkgDataSmeared[iptBin][iSmearing]->GetBinContent(mBin);
	//if(mBin%20==0) Printf("bkg %f",bkg);
	// TRandom smearingGeneratorB;
	// Float_t bkgSmeared=smearingGeneratorB.Poisson(bkg);
	Float_t bkgSmeared=gRandom->Poisson(bkg);
	//if(mBin%20==0) Printf("bkgSmeared %f",bkgSmeared);
	Float_t bkgSmearedErr=TMath::Sqrt(bkgSmeared);
	hBkgDataSmeared[iptBin][iSmearing]->SetBinContent(mBin,bkgSmeared);
	hBkgDataSmeared[iptBin][iSmearing]->SetBinError(mBin,bkgSmearedErr);
      }
      Printf("End first loop on mass bins");
      
      //smearing signal histogram - MC 
      Int_t nMassBinsMC=hSigMCScaledSmeared[iptBin][0]->GetNbinsX();
      Printf("nMassBinsMC %d", nMassBinsMC);
      Printf("Second loop on mass bins, %d", nMassBinsMC);
      Printf("Smearing hSigMCScaled[%d][%d]",iptBin,iSmearing);
      hSigMCScaledSmeared[iptBin][iSmearing]=(TH1D*)hSigMCScaledSmeared[iptBin][0]->Clone(Form("hSigMCScaled_ptBin%d_stepSmearing%d",iptBin,iSmearing));
      hSigMCScaledSmeared[iptBin][iSmearing]->SetDirectory(0);
      if(hSigMCScaledSmeared[iptBin][iSmearing]) hSigMCScaledSmeared[iptBin][iSmearing]->Print();
      for(Int_t im=0; im<nMassBinsMC; im++){
	Float_t sigMC=hSigMCScaledSmeared[iptBin][iSmearing]->GetBinContent(im+1);
	//if(im%20==0)Printf("sigMC %f",sigMC);
	// TRandom smearingGeneratorS;
	// Float_t sigMCSmeared=smearingGeneratorS.Poisson(sigMC);
	Float_t sigMCSmeared=gRandom->Poisson(sigMC);
	//if(im%20==0)Printf("sigMCSmeared %f",sigMCSmeared);
	Float_t sigMCSmearedErr=TMath::Sqrt(sigMCSmeared);
	hSigMCScaledSmeared[iptBin][iSmearing]->SetBinContent(im+1,sigMCSmeared);
	hSigMCScaledSmeared[iptBin][iSmearing]->SetBinError(im+1,sigMCSmearedErr);
      }
      Printf("End second loop on mass bins");
      intSMCSmeared[iptBin][iSmearing]=hSigMCScaledSmeared[iptBin][iSmearing]->Integral();
     
     
      //smearing R histogram - MC
      Int_t nMassBinsRfl=hRflMCScaledSmeared[iptBin][0]->GetNbinsX();
      Printf("Third loop on mass bins, %d", nMassBinsRfl);
      Printf("Smearing hRflMCScaledSmeared[%d][%d]",iptBin,iSmearing);
      hRflMCScaledSmeared[iptBin][iSmearing]=(TH1D*)hRflMCScaledSmeared[iptBin][0]->Clone(Form("hRflMCScaled_ptBin%d_stepSmearing%d",iptBin,iSmearing));
      hRflMCScaledSmeared[iptBin][iSmearing]->SetDirectory(0);
      if(hRflMCScaledSmeared[iptBin][iSmearing]) hRflMCScaledSmeared[iptBin][iSmearing]->Print();
      //scaling R histogram - MC
      hRflMCScaledSmeared[iptBin][iSmearing]->Scale(intSMCSmeared[iptBin][iSmearing]/intS[iptBin]);
      for(Int_t im=0; im<nMassBinsRfl; im++){
	Float_t rflMC=hRflMCScaledSmeared[iptBin][iSmearing]->GetBinContent(im+1);
	//if(im%20==0)Printf("rflMC %f",rflMC);
	// TRandom smearingGeneratorR;
	// Float_t rflMCSmeared=smearingGeneratorR.Poisson(rflMC);
	Float_t rflMCSmeared=gRandom->Poisson(rflMC);
	//if(im%20==0)Printf("rflMCSmeared %f",rflMCSmeared);
	Float_t rflMCSmearedErr=TMath::Sqrt(rflMCSmeared);      
	hRflMCScaledSmeared[iptBin][iSmearing]->SetBinContent(im+1,rflMCSmeared);
	hRflMCScaledSmeared[iptBin][iSmearing]->SetBinError(im+1,rflMCSmearedErr);
      }
      Printf("End third loop on mass bins");
      intRMCSmeared[iptBin][iSmearing]=hRflMCScaledSmeared[iptBin][iSmearing]->Integral();
    
      hSignalInjected[iSmearing]->SetBinContent(iptBin+1,intSMCSmeared[iptBin][iSmearing]);
      hRflInjected[iSmearing]->SetBinContent(iptBin+1,intRMCSmeared[iptBin][iSmearing]);
      hMassToy[iptBin][iSmearing]=(TH1D*)hBkgDataSmeared[iptBin][iSmearing]->Clone(Form("hMass_ptBin%d_stepSmearing%d",iptBin,iSmearing));
      hMassToy[iptBin][iSmearing]->SetDirectory(0);
      hMassToy[iptBin][iSmearing]->Add(hSigMCScaledSmeared[iptBin][iSmearing],1);
      hMassToy[iptBin][iSmearing]->Add(hRflMCScaledSmeared[iptBin][iSmearing],1);
      hMassToyNoRfl[iptBin][iSmearing]=(TH1D*)hBkgDataSmeared[iptBin][iSmearing]->Clone(Form("hMassNoRfl_ptBin%d_stepSmearing%d",iptBin,iSmearing));
      hMassToyNoRfl[iptBin][iSmearing]->SetDirectory(0);
      hMassToyNoRfl[iptBin][iSmearing]->Add(hSigMCScaledSmeared[iptBin][iSmearing],1);

      outf2->cd();
      hBkgDataSmeared[iptBin][iSmearing]->Write();
      hSigMCScaledSmeared[iptBin][iSmearing]->Write();
      hRflMCScaledSmeared[iptBin][iSmearing]->Write();
      hMassToy[iptBin][iSmearing]->Write();
      hMassToyNoRfl[iptBin][iSmearing]->Write();

    } //end loop on pt bins

    outf2->cd();
    hSignalInjected[iSmearing]->Write();
    hRflInjected[iSmearing]->Write();
    outf2->Close();

  }// end loop on smearing steps

  
  leg->AddEntry(hSignalInjected[0],"step 0","lep");
 
  for(Int_t i=1; i<nSmearing; i++){
    Int_t kPad=1;
    Int_t color=i+1;
    if(i==9) color=12; 
    for(Int_t j=0; j<nPtBins; j++){

      hBkgDataSmeared[j][i]->SetMarkerStyle(20+i);
      hBkgDataSmeared[j][i]->SetMarkerColor(color);
      hBkgDataSmeared[j][i]->SetLineColor(color);
      hBkgDataSmeared[j][i]->SetLineWidth(1);
      c2->cd(kPad);
      hBkgDataSmeared[j][i]->Draw("same hist");

      hSigMCScaledSmeared[j][i]->SetMarkerStyle(20+i);
      hSigMCScaledSmeared[j][i]->SetMarkerColor(color);
      hSigMCScaledSmeared[j][i]->SetLineColor(color);
      hSigMCScaledSmeared[j][i]->SetLineWidth(1);
      c3->cd(kPad);
      hSigMCScaledSmeared[j][i]->Draw("same hist");

      hRflMCScaledSmeared[j][i]->SetMarkerStyle(20+i);
      hRflMCScaledSmeared[j][i]->SetMarkerColor(color);
      hRflMCScaledSmeared[j][i]->SetLineColor(color);
      hRflMCScaledSmeared[j][i]->SetLineWidth(1);
      c4->cd(kPad);
      hRflMCScaledSmeared[j][i]->Draw("same hist");

      hMassToy[j][i]->SetMarkerStyle(20+i);
      hMassToy[j][i]->SetMarkerColor(color);
      hMassToy[j][i]->SetLineColor(color);
      hMassToy[j][i]->SetLineWidth(1);
      c5->cd(kPad);
      hMassToy[j][i]->Draw("same hist");

      hMassToyNoRfl[j][i]->SetMarkerStyle(20+i);
      hMassToyNoRfl[j][i]->SetMarkerColor(color);
      hMassToyNoRfl[j][i]->SetLineColor(color);
      hMassToyNoRfl[j][i]->SetLineWidth(1);
      c12->cd(kPad);
      hMassToyNoRfl[j][i]->Draw("same hist");

      kPad++;
    }
    c6->cd();
    hSignalInjected[i]->SetMarkerStyle(20+i);
    hSignalInjected[i]->SetMarkerColor(color);
    hSignalInjected[i]->SetLineColor(color);
    hSignalInjected[i]->SetLineWidth(2);
    hSignalInjected[i]->Draw("same hist");

    leg->AddEntry(hSignalInjected[i],Form("step %d",i),"lep");

    c7->cd();
    hRflInjected[i]->SetMarkerStyle(20+i);
    hRflInjected[i]->SetMarkerColor(color);
    hRflInjected[i]->SetLineColor(color);
    hRflInjected[i]->SetLineWidth(2);
    hRflInjected[i]->Draw("same hist");
  }

  /*c2->cd(1);
  leg->Draw();
  c3->cd(1);
  leg->Draw();
  c4->cd(1);
  leg->Draw();
  c5->cd(1);
  leg->Draw();
  c6->cd();
  leg->Draw();
  c7->cd();
  leg->Draw();*/


  // end reflections


  

  TCanvas *cpar=new TCanvas("cpar","Fit params",1200,600);
  cpar->Divide(2,1);
  cpar->cd(1);
  hMass->SetMarkerStyle(20);
  hMass->Draw("PE");
  hMass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hMass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
  cpar->cd(2);
  hSigma->SetMarkerStyle(20);
  hSigma->Draw("PE");
  hSigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hSigma->GetYaxis()->SetTitle("Sigma (GeV/c^{2})");

  TCanvas* csig=new TCanvas("csig","Results",1200,600);
  csig->Divide(3,1);
  csig->cd(1);
  hSignal->SetMarkerStyle(20);
  hSignal->SetMarkerColor(4);
  hSignal->SetLineColor(4);
  hSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hSignal->GetYaxis()->SetTitle("Signal");
  hSignal->Draw("P");
  hCntSig1->SetMarkerStyle(26);
  hCntSig1->SetMarkerColor(2);
  hCntSig1->SetLineColor(2);
  hCntSig1->Draw("PSAME");
  hCntSig2->SetMarkerStyle(29);
  hCntSig2->SetMarkerColor(kGray+1);
  hCntSig2->SetLineColor(kGray+1);
  hCntSig2->Draw("PSAME");
  TLegend* legg=new TLegend(0.4,0.7,0.89,0.89);
  legg->SetFillColor(0);
  TLegendEntry* ent=legg->AddEntry(hSignal,"From Fit","PL");
  ent->SetTextColor(hSignal->GetMarkerColor());
  ent=legg->AddEntry(hCntSig1,"Bin Count. only bkg","PL");
  ent->SetTextColor(hCntSig1->GetMarkerColor());
  ent=legg->AddEntry(hCntSig2,"Bon Count. sign+bkg","PL");
  ent->SetTextColor(hCntSig2->GetMarkerColor());
  legg->Draw();
  csig->cd(2);
  hBackground->SetMarkerStyle(20);
  hBackground->Draw("P");
  hBackground->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hBackground->GetYaxis()->SetTitle("Background");
  csig->cd(3);
  hSignificance->SetMarkerStyle(20);
  hSignificance->Draw("P");
  hSignificance->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hSignificance->GetYaxis()->SetTitle("Significance");

  TCanvas* cDiffS=new TCanvas("cDiffS","Difference",1200,600);
  cDiffS->Divide(2,1);
  cDiffS->cd(1);
  hRelErrSig->SetMarkerStyle(20); //fullcircle
  hRelErrSig->SetTitleOffset(1.2);  
  hRelErrSig->SetTitle("Rel Error from Fit;p_{T} (GeV/c);Signal Relative Error");
  hRelErrSig->Draw("P");
  hInvSignif->SetMarkerStyle(21); //fullsquare
  hInvSignif->SetMarkerColor(kMagenta+1);
  hInvSignif->SetLineColor(kMagenta+1);
  hInvSignif->Draw("PSAMES");
  TLegend* leg2=new TLegend(0.4,0.7,0.89,0.89);
  leg2->SetFillColor(0);
  TLegendEntry* ent2=leg2->AddEntry(hRelErrSig,"From Fit","P");
  ent2->SetTextColor(hRelErrSig->GetMarkerColor());
  ent2=leg2->AddEntry(hInvSignif,"1/Significance","PL");
  ent2->SetTextColor(hInvSignif->GetMarkerColor());
  leg2->Draw();

  cDiffS->cd(2);
  hNDiffCntSig1->SetMarkerStyle(26);
  hNDiffCntSig1->SetMarkerColor(2);
  hNDiffCntSig1->SetLineColor(2);
  hNDiffCntSig1->SetTitle("Cmp Fit-Count;p_{T} (GeV/c);(S_{fit}-S_{count})/S_{fit}");
  hNDiffCntSig1->Draw("P");
  hNDiffCntSig2->SetMarkerStyle(29);
  hNDiffCntSig2->SetMarkerColor(kGray+1);
  hNDiffCntSig2->SetLineColor(kGray+1);
  hNDiffCntSig2->Draw("PSAME");
  TLegend* leg1=new TLegend(0.4,0.7,0.89,0.89);
  leg1->SetFillColor(0);
  TLegendEntry* ent1=leg1->AddEntry(hNDiffCntSig1,"Bin Count. only bkg","PL");
  ent1->SetTextColor(hNDiffCntSig1->GetMarkerColor());
  ent1=leg1->AddEntry(hNDiffCntSig2,"Bin Count. sign+bkg","PL");
  ent1->SetTextColor(hNDiffCntSig2->GetMarkerColor());
  leg1->Draw();

  TCanvas *cDiffonly = new TCanvas("SystErrBinCount","SystErrBinCount");
  cDiffonly->cd();
  hNDiffCntSig1->SetMarkerStyle(22);
  hNDiffCntSig1->SetMarkerColor(2);
  hNDiffCntSig1->SetMarkerSize(1.3);
  hNDiffCntSig1->SetLineColor(2);
  hNDiffCntSig1->SetLineWidth(2);
  hNDiffCntSig1->SetTitle("Cmp Fit-Count;p_{T} (GeV/c);(S_{fit}-S_{count})/S_{fit}");
  hNDiffCntSig1->Draw("P");
  hNDiffCntSig2->SetMarkerStyle(20);
  hNDiffCntSig2->SetMarkerColor(kBlue);
  hNDiffCntSig2->SetLineColor(kBlue);
  hNDiffCntSig2->SetLineWidth(2);
  hNDiffCntSig2->SetMarkerSize(1.3);
  hNDiffCntSig2->Draw("PSAME");
  TLegend* legOnly=new TLegend(0.4,0.7,0.89,0.89);
  legOnly->SetFillColor(0);
  TLegendEntry* entOnly=legOnly->AddEntry(hNDiffCntSig1,"Bin Count. only bkg","PL");
  entOnly->SetTextColor(hNDiffCntSig1->GetMarkerColor());
  entOnly=legOnly->AddEntry(hNDiffCntSig2,"Bin Count. sign+bkg","PL");
  entOnly->SetTextColor(hNDiffCntSig2->GetMarkerColor());
  TPaveText* tsyst=new TPaveText(0.2,0.90,0.8,0.95,"NDC");
  tsyst->SetFillStyle(0);
  tsyst->SetBorderSize(0);
  tsyst->AddText(0.,0.,"D^{0} #rightarrow K#pi Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV, 20-40% centrality");
  tsyst->SetTextFont(42);
  tsyst->Draw();

  TPaveText* tsystlabel=new TPaveText(0.8,0.90,1.,0.95,"NDC");
  tsystlabel->SetFillStyle(0);
  tsystlabel->SetBorderSize(0);
  tsystlabel->AddText(0.,0.,"ALICE Performance");
  tsystlabel->SetTextFont(42);
  tsystlabel->Draw();
  
  legOnly->Draw();
    
  TGraph* grReducedChiSquare=new TGraph(nPtBins,ptlims,arrchisquare);
  grReducedChiSquare->SetName("grReducedChiSquare");
  grReducedChiSquare->SetTitle("Reduced Chi2;p_{T} (GeV/c);#tilde{#chi}^{2}");
  TCanvas *cChi2=new TCanvas("cChi2","reduced chi square",600,600);
  cChi2->cd();
  grReducedChiSquare->SetMarkerStyle(21);
  grReducedChiSquare->Draw("AP");

  TCanvas* cbkgNormSigma=new TCanvas("cbkgNormSigma","Background normalized to sigma",400,600);
  cbkgNormSigma->cd();
  hBackgroundNormSigma->SetMarkerStyle(20);
  hBackgroundNormSigma->Draw("P");
  hBackgroundNormSigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hBackgroundNormSigma->GetYaxis()->SetTitle("Background #times 3 #times 0.012/ (3 #times #sigma)");
  hBackgroundNormSigma->Draw();
 


  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) {
    /*if(analysisType==kD0toKpi)*/ partname="D0";
    //if(analysisType==kDplusKpipi) partname="Dplus";
  }
  if(optPartAntiPart==kAntiParticleOnly) {
    /*if(analysisType==kD0toKpi)*/ partname="D0bar";
    //if(analysisType==kDplusKpipi) partname="Dminus";
  }

  TFile* outf=new TFile(Form("Reflections%s%s.root",partname.Data(),CutsType),"recreate");
  outf->cd();
  hMass->Write();
  hSigma->Write();
  hCntSig1->Write();
  hCntSig2->Write();
  hNDiffCntSig1->Write();
  hNDiffCntSig2->Write();
  hSignal->Write();
  hRelErrSig->Write();
  hInvSignif->Write();
  hBackground->Write();
  hBackgroundNormSigma->Write();
  hSignificance->Write();
  grReducedChiSquare->Write();
  funBckStore1->Write();
  funBckStore2->Write();
  funBckStore3->Write();
  outf->Close();

//cri
// TFile *miofile =new TFile("massfile_central.root","recreate");
 // miofile->cd();
 // for(Int_t i=0;i<nPtBins;i++){
  //prendi histo e write
  //cout<< " i " << i << endl;
  //cout<< " hMass " << hMass[i] << endl;
 // hmass[i]->Write();
 // }
 //miofile->Close();
//cri

  TFile* outfMC=new TFile(Form("BackgroundSignalReflectionsMC%s%s.root",partname.Data(),CutsType),"recreate");
  outfMC->cd();
  for(Int_t i=0;i<nPtBins;i++){
    hsig[i]->Write(Form("histSgn_%d",i));
    hbkg[i]->Write(Form("histBkg_%d",i));
    hrfl[i]->Write(Form("histRfl_%d",i));
  }
  outfMC->Close();

  return;

}
  
Bool_t LoadD0toKpiMCHistos(TString fileName1,TString fileName2,TString fileName3,TString fileName4, TH1F** hSig, TH1F** hBkg, TH1F** hRfl, const char *CutsType){
  TFile *f[4];
  f[0]=TFile::Open(fileName1.Data());
  f[1]=TFile::Open(fileName2.Data());
  f[2]=TFile::Open(fileName3.Data());
  f[3]=TFile::Open(fileName4.Data());
 for(Int_t i=0;i<4;i++){ 
 if(!f[i]){
    printf("ERROR: file %s does not exist\n",fileName1.Data());
    return kFALSE;
  }
}
cout<< " f " << f[0] << " 1 " << f[1] << " 2 " << f[2] << " f3 " << f[3] << endl;
  cout<< " prima di open file " << endl;
//  printf("Open File %s\n",f[i]->GetName());
cout<< " dopo di open file " << endl;
  
  
  TString dirname="PWG3_D2H_D0InvMass";
  if(optPartAntiPart==kParticleOnly) dirname+="D0";
  else if(optPartAntiPart==kAntiParticleOnly) dirname+="D0bar";
  if(cutsappliedondistr) dirname+="C";
  dirname += CutsType;
  //dirname += "MC";
  
  TDirectory *dir[4];
  for(Int_t i=0;i<4;i++){
 dir[i] = (TDirectory*)f[i]->Get(dirname);
  if(!dir[i]){
    printf("ERROR: directory %s not found in %s\n",dirname.Data(),fileName1.Data());
    return kFALSE;
  }
}
cout<< " dir " << dir[0] << " dir " << dir[1] << " dir " << dir[2] << " " << dir[3] << endl;

  TString listmassname="coutputmassD0Mass";
  if(optPartAntiPart==kParticleOnly) listmassname+="D0";
  else if(optPartAntiPart==kAntiParticleOnly) listmassname+="D0bar";
  if(cutsappliedondistr) listmassname+="C";
  listmassname += CutsType;
  listmassname += "0100";
  printf("List name %s\n",listmassname.Data());
cout<< " prima di prendere hlist " << endl;
  TList *hlist[4];
  for(Int_t i=0;i<4;i++)
  hlist[i]=(TList*)dir[i]->Get(listmassname);

  cout<<  " hlist " << hlist[0] << " hlist " << hlist[1] << " hlist " << hlist[2]<< " hlist " << hlist[3] << endl;
  
  TString cutsobjname="cutsD0";
  if(optPartAntiPart==kParticleOnly) cutsobjname+="D0";
  else if(optPartAntiPart==kAntiParticleOnly) cutsobjname+="D0bar";
  if(cutsappliedondistr) cutsobjname+="C";
  cutsobjname += CutsType;
  //cutsobjname += "MC";
  cutsobjname += "0100";
  printf("Cuts name %s\n",cutsobjname.Data());
cout<< " prima del check sul cut " << endl;
  AliRDHFCutsD0toKpi *dcuts=(AliRDHFCutsD0toKpi*)dir[0]->Get(cutsobjname);
  if(!dcuts) {
    printf("ERROR: Cut objects do not match\n");
    return kFALSE;
  }
 
cout<< " dopo il check sul cuts " << endl;
  Int_t nPtBinsCuts=dcuts->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;

  printf("Loop over the pt bins\n");

  Int_t iFinBin=0;
  for(Int_t i=0;i<nPtBinsCuts;i++){
cout<< " pt bin " << endl;
    if(ptlimsCuts[i]>=ptlims[iFinBin+1]) iFinBin+=1; 
    if(iFinBin>nPtBins) break;
    if(ptlimsCuts[i]>=ptlims[iFinBin] && 
       ptlimsCuts[i+1]<=ptlims[iFinBin+1]){
     //cri
for(Int_t iFile=0; iFile<4; iFile++){
cout<< " loop files " <<iFile << endl;
      TH1F *htempS;
      Printf(" Looking for histo histSig_%d in pt bin %d",i,i);
      cout<< " hlist " << hlist[iFile] <<endl;
      htempS=(TH1F*)hlist[iFile]->FindObject(Form("histSgn_%d",i));
      cout<< " htempS " << htempS << endl;
      if(!hSig[iFinBin]){
         cout<<" sono nel caso in cui non c'e' hsig e lo creo " << endl;
	hSig[iFinBin]=new TH1F(*htempS);
      }else{
        cout<< " hSign e' riempito con Add " << endl;
	hSig[iFinBin]->Add(htempS);
      }
      cout<< " ora passo al histo backgroud " << endl;
      TH1F *htempB;
      Printf(" Looking for histo histBkg_%d in pt bin %d",i,i);
      htempB=(TH1F*)hlist[iFile]->FindObject(Form("histBkg_%d",i));
	
      if(!hBkg[iFinBin]){
	hBkg[iFinBin]=new TH1F(*htempB);
      }else{
	hBkg[iFinBin]->Add(htempB);
      }

      TH1F *htempR;
      Printf(" Looking for histo histRfl_%d in pt bin %d",i,i);
      htempR=(TH1F*)hlist[iFile]->FindObject(Form("histRfl_%d",i));
	
      if(!hRfl[iFinBin]){
	hRfl[iFinBin]=new TH1F(*htempR);
      }else{
	hRfl[iFinBin]->Add(htempR);
      }
    }
  }
 }
  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="D0";
  if(optPartAntiPart==kAntiParticleOnly) partname="D0bar";
  
  TFile* outf=new TFile(Form("Reflections%s%s.root",partname.Data(),CutsType),"recreate");
  outf->cd();
  dcuts->Write();
  outf->Close();
  return kTRUE;
}

Bool_t LoadD0toKpiHistos(TObjArray* listFiles, TH1F** hMass, const char *CutsType){
  //
  Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
  AliRDHFCutsD0toKpi** dcuts=new AliRDHFCutsD0toKpi*[nFiles];

  Int_t nReadFiles=0;
  for(Int_t iFile=0; iFile<nFiles; iFile++){
    TString fName=((TObjString*)listFiles->At(iFile))->GetString();    
    TFile *f=TFile::Open(fName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fName.Data());
      continue;
    }
    printf("Open File %s\n",f->GetName());

    TString dirname="PWG3_D2H_D0InvMass";
    if(optPartAntiPart==kParticleOnly) dirname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) dirname+="D0bar";
    if(cutsappliedondistr) dirname+="C";
    dirname += CutsType;
    TDirectory *dir = (TDirectory*)f->Get(dirname);
    if(!dir){
      printf("ERROR: directory %s not found in %s\n",dirname.Data(),fName.Data());
      continue;
    }
    TString listmassname="coutputmassD0Mass";
    if(optPartAntiPart==kParticleOnly) listmassname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) listmassname+="D0bar";
    if(cutsappliedondistr) listmassname+="C";
    listmassname += CutsType;
    listmassname += "0100";
    printf("List name %s\n",listmassname.Data());
    hlist[nReadFiles]=(TList*)dir->Get(listmassname);

    TString cutsobjname="cutsD0";
    if(optPartAntiPart==kParticleOnly) cutsobjname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) cutsobjname+="D0bar";
    if(cutsappliedondistr) cutsobjname+="C";
    cutsobjname += CutsType;
    cutsobjname += "0100";
    printf("Cuts name %s\n",cutsobjname.Data());

    dcuts[nReadFiles]=(AliRDHFCutsD0toKpi*)dir->Get(cutsobjname);
    if(!dcuts[nReadFiles]) {
      printf("ERROR: Cut objects do not match\n");
      return kFALSE;
    }
    if(nReadFiles>0){
      Bool_t sameCuts=dcuts[nReadFiles]->CompareCuts(dcuts[0]);
      if(!sameCuts){
	printf("ERROR: Cut objects do not match\n");
	return kFALSE;
      }
    }
    nReadFiles++;
  }
  if(nReadFiles<nFiles){
    printf("WARNING: not all requested files have been found\n");
    if (nReadFiles==0) {
      printf("ERROR: Any file/dir found\n");
      return kFALSE;
    }
  }
  //  printf("Cuts type %s, particle/antipart %d\n",CutsType,optPartAntiPart);

  Int_t nPtBinsCuts=dcuts[0]->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;

  printf("Loop over the pt bins\n");

  Int_t iFinBin=0;
  for(Int_t i=0;i<nPtBinsCuts;i++){
    if(ptlimsCuts[i]>=ptlims[iFinBin+1]) iFinBin+=1; 
    if(iFinBin>nPtBins) break;
    if(ptlimsCuts[i]>=ptlims[iFinBin] && 
       ptlimsCuts[i+1]<=ptlims[iFinBin+1]){
      for(Int_t iFile=0; iFile<nReadFiles; iFile++){
	TH1F * htemp; 
	printf(" Looking for histo histMass_%d in pt bin %d",i,i);
	htemp=(TH1F*)hlist[iFile]->FindObject(Form("histMass_%d",i));
	
	if(!hMass[iFinBin]){
	  hMass[iFinBin]=new TH1F(*htemp);
	}else{
	  hMass[iFinBin]->Add(htemp);
	}
      }
    }
  }
  return kTRUE;
}
