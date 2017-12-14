/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//
//  Class to extract D-jet pT or z spectrum
//
//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------

#include "AliDJetsPtSpectrum.h"

//___________________________________________________________________________________________
AliDJetsPtSpectrum::AliDJetsPtSpectrum():
fFileOutput(0x0),
fFileNameOutput("out"),
fOutputDir("out"),
fDmesonSpecie(kDStarD0pi),
fDmesonLabel("Dstar"),
fYieldApproach(kEffScale),
fpTmin(0),
fpTmax(99),
fzmin(-2),
fzmax(2),
fjetetamin(-0.5),
fjetetamax(0.5),
fnDbins(0),			
fDbinpTedges(0x0), 
fnJetbins(0),
fJetbinpTedges(0x0),
fDEffValues(0x0), 
fnSigmaSign(3),

fInvMass3D(0x0),
fRawPtSpectrum(0x0),
fEffCorrPtSpectrum(0x0)

{

  if(fDmesonSpecie==kD0toKpi){
    fnSigmaBkg = [-8,-5,5,8];
    fDmass = ;
    fDsigma = ;
  }
  else if(fDmesonSpecie==kDStarD0pi){
    fnSigmaBkg = [-8,-5,5,13];
    fDmass = 0.1455;
    fDsigma = 0.00065;
    
  }
  
  

}

//___________________________________________________________________________________________
AliDJetsPtSpectrum::AliDJetsPtSpectrum(const AliDJetsPtSpectrum &source):
// copy constructor
fFileInput(source.fFileInput),
fFileRawOutput(source.fFileRawOutput),
fDmesonSpecie(source.fDmesonSpecie),
fDmesonLabel(source.fDmesonLabel),
fYieldApproach(source.fYieldApproach),
fpTmin(source.fpTmin),
fpTmax(source.fpTmax),
fnDbins(source.fnDbins),			
fDbinpTedges(source.fDbinpTedges),    		
fDEffValues(source.fDEffValues),
fMassPlot(source.fMassPlot)
{


}

//___________________________________________________________________________________________
AliDJetsPtSpectrum::~AliDJetsPtSpectrum() {
//destructor

  
  delete[] fnSigmaBkg;
  if(fDbinpTedges)    delete[] fDbinpTedges;   
  if(fJetbinpTedges)  delete[] fJetbinpTedges;   
  if(fDEffValues)     delete[] fDEffValues; 
  if(fMeanSigmaVar)   delete[] fMeanSigmaVar; 
  
  

}

//___________________________________________________________________________________________
Bool_t AliDJetsPtSpectrum::SetDmesonSpecie(DMesonSpecies k){

  if(k<0 || k>1) {
    printf("Error! D meson specie not correctly set!\n");
    return kFALSE;
  } else if(k==0) fDmesonLabel="Dzero";
  else fDmesonLabel="Dstar";

  fDmesonSpecie=k;
  return kTRUE;

}

//___________________________________________________________________________________________
Bool_t AliDJetsPtSpectrum::ExtractInputMassPlot(const Int_t nFiles, TString *datafiles, Bool_t postfix, TString listName){

  std::cout << "===== Configuration:\nD meson: " << fDmesonLabel << "\nMethod (eff.scale/sideband): " << fYieldApproach << std::endl;

  Bool_t success = kTRUE;

  TFile *File;
  TDirectoryFile* dir;
  TList *histList;
  THnSparseF *sparse;
  
  TString branchHistName;
  if(fDmesonSpecie==kD0toKpi)         branchHistName = "histosDMBN";      // !!!! #TODO
  else if(fDmesonSpecie==kDStarD0pi)  branchHistName = "histosDStarMBN";

  for (int j=0;j<nFiles;j++){
        File = new TFile(datafiles[j],"read");
        dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
        if(!File){
          std::cout << "File " << File << " cannot be opened! check your file path!" << std::endl; 
          return kFALSE;
        }
    
        for(int i=0;i<fNDdata; i++){
            if(postfix) histList =  (TList*)dir->Get(Form("%s%d%s",branchHistName.Data(),i,listName.Data()));
            else histList =  (TList*)dir->Get(Form("%s%d",branchHistName.Data(),i));
            sparse = (THnSparseF*)histList->FindObject("hsDphiz");
            sparse->GetAxis(0)->SetRangeUser(fzmin,fzmax); 
            //sparse->GetAxis(1)->SetRangeUser(jetmin,jetmax);
            sparse->GetAxis(5)->SetRangeUser(fjetetamin,fjetetamax);
            if(j==0 && i==0) fInvMass3D=(TH3D*)sparse->Projection(3,1,2);
            else fInvMass3D->Add((TH3D*)sparse->Projection(3,1,2));
        }
    
  }

  if(!fInvMass3D) {
      std::cout << "Error in extracting 3D raw mass plot! Exiting..." << std::endl;
      return kFALSE;
  }

  if(fYieldApproach==kEffScale) success = ExtractInputMassPlotEffScale();
  if(fYieldApproach==kSideband) success = ExtractInputMassPlotSideband();  
 
 
  if(success) std::cout << "Extracted mass spectrum" << std::endl;


  return success;

}



//___________________________________________________________________________________________
Bool_t AliDJetsPtSpectrum::ExtractInputMassPlotDstarEffScale() {


   std::cout << "===== Extracting raw/eff. corrected jet pT spectrum for Eff. Scaled method =====" << std::endl;  

    if(!fMassPlot) {
      std::cout << "Error in extracting the mass plot! Exiting..." << std::endl;
      return kFALSE;
    }     

    return kTRUE;
       
}

//___________________________________________________________________________________________
Bool_t AliDJetsPtSpectrum::ExtractInputMassPlotDstarSideband() {

    std::cout << "===== Extracting raw/eff. corrected jet pT spectrum for SB method =====" << std::endl;    

    TH1F *hmean = new TH1F("hmean","hmean",fnDbins,fDbinpTedges);
    TH1F *hsigma = new TH1F("hsigma","hsigma",fnDbins,fDbinpTedges);
    TH1F *hrelErr = new TH1F("hrelErr","hrelErr",fnDbins,fDbinpTedges);
    TH1F *hsign = new TH1F("hsign","hsign",fnDbins,fDbinpTedges);
    TH1F *hsb = new TH1F("hsb","hsb",fnDbins,fDbinpTedges);
    TH1F *hSignal = new TH1F("hSignal","hSignal",fnDbins,fDbinpTedges);
    hSignal->Sumw2();

    TH1F *hmass[fnDbins];
    TH1F *hjetpt[fnDbins];
    TH1F *hjetpt_s[fnDbins];
    TH1F *hjetptsub[fnDbins];
    TH1F *hjetptcorr[fnDbins];
    TF1 *fullfit[fnDbins];

     for(int i=0; i<fnDbins; i++){
                 
        TH1F *hh=(TH1F*)hInvMassptD->ProjectionX(Form("hh_%d",i),hInvMassptD->GetYaxis()->FindBin(jetmin), hInvMassptD->GetYaxis()->FindBin(jetmax)-1,hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[i]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[i+1])-1);
        hh->Rebin(2);

        TH1F *hmassfit = (TH1F*)hh->Clone("hmassfit");
        float hmin=TMath::Max(minf,hmassfit->GetBinLowEdge(2));
        float hmax=TMath::Min(maxf,hmassfit->GetBinLowEdge(hmassfit->GetNbinsX()));
       // AliHFMassFitter* fitterp=new AliHFMassFitter((TH1F*)hmassfit,hmin,hmax,1,bkgtype,0);
        AliHFInvMassFitter* fitterp=new AliHFInvMassFitter((TH1F*)hmassfit,hmin,hmax,fbkgtype,0);
        fitterp->SetInitialGaussianMean(fDmass);
        fitterp->SetInitialGaussianSigma(fDsigma);
        fitterp->MassFitter(kFALSE);
        
        TH1F* h=fitterp->GetHistoClone();
        fullfit[i]=h->GetFunction("funcmass");
        fullfit[i]->SetName(Form("fullfit_%d",i));
        //((TList*)h->GetListOfFunctions())->RemoveAt(2);
        hmass[i] = (TH1F*)h->Clone(Form("hmass_%",i));
        hmass[i]->SetName(Form("hmass_%d",i));
         
        float Dsigma = fitterp->GetSigma();
        float Dmean = fitterp->GetMean();
        float signal_l_min = Dmean+fnSigmaBkg[0]*Dsigma;
        float signal_l_max = Dmean+fnSigmaBkg[1]*Dsigma;
        float signal_u_min = Dmean+fnSigmaBkg[2]*Dsigma;
        float signal_u_max = Dmean+fnSigmaBkg[3]*Dsigma;
        float signal_c_min = Dmean-fnSigmaSign*Dsigma;
        float signal_c_max = Dmean+fnSigmaSign*Dsigma;
        
        int binmin = hmass[i]->GetXaxis()->FindBin(signal_c_min);
        int binmax = hmass[i]->GetXaxis()->FindBin(signal_c_max);
        double binwidth = hmass[i]->GetXaxis()->GetBinWidth(1)*0.5;
        
        Double_t s,serr,b,berr,signf,signferr;
        fitterp->Signal(fnSigmaSign,s,serr);
        fitterp->Background(hmass[i]->GetXaxis()->GetBinCenter(binmin)-binwidth, hmass[i]->GetXaxis()->GetBinCenter(binmax-1)+binwidth ,b,berr);
        fitterp->Significance(fnSigmaSign,signf,signferr);
        Double_t sob=s/b, soberr;
        soberr=TMath::Sqrt((serr/b)*(serr/b) + (s/b/b*berr)*(s/b/b*berr));
       
      
        //-------- jet pt spectrum - signal
        hjetpt[i]=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_%d",i),hInvMassptD->GetXaxis()->FindBin(signal_c_min), hInvMassptD->GetXaxis()->FindBin(signal_c_max)-1,hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[i]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[i+1])-1);
        
        //------ jet pt spectrum - side bands
        TH1F* hjetpt_s1=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_s1%d",i),hInvMassptD->GetXaxis()->FindBin(signal_l_min), hInvMassptD->GetXaxis()->FindBin(signal_l_max)-1,hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[i]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[i+1])-1);
        hjetpt_s1->Sumw2();
        TH1F* hjetpt_s2=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_s2%d",i),hInvMassptD->GetXaxis()->FindBin(signal_u_min), hInvMassptD->GetXaxis()->FindBin(signal_u_max)-1,hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[i]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[i+1])-1);
        hjetpt_s2->Sumw2();
        hjetpt_s[i] = (TH1F*)hjetpt_s1->Clone(Form("hjetpt_s_%d",i));
        hjetpt_s->Sumw2();
        hjetpt_s[i]->Add(hjetpt_s2);
        
        // scale background from side bands to the background under the peak
        double scaling = b / hjetpt_s[i]->Integral();
        hjetpt_s[i]->Scale(scaling);
        //------- subtract background from signal jet
        hjetptsub[i] = (TH1F*)hjetpt[i]->Clone(Form("hjetptsub_%d",i));
        hjetptsub[i]->Add(hjetpt_s[i],-1);
        
        if(!i) fRawPtSpectrum = (TH1F*)hjetptsub[i]->Clone("fRawPtSpectrum");
        else fRawPtSpectrum->Add(hjetptsub[i]);
        
        if(fEffCorr){
        //------- correct for D* efficiency
          hjetptcorr[i] = (TH1F*)hjetptsub[i]->Clone(Form("hjetptcorr_%d",i));
          hjetptcorr[i]->Scale(1./fDEffValues[i]); // D efficiency
          //----- add corrected jet pt spectrum distributions in each D pt bin into a one distribution
          if(!i) fEffCorrPtSpectrum = (TH1F*)hjetptcorr[i]->Clone("fEffCorrPtSpectrum");
          else fEffCorrPtSpectrum->Add(hjetptcorr[i]);
        }
        
        // ---------------- fitting results
        hmean->SetBinContent(i+1,fitterp->GetMean()*1000);
        hmean->SetBinError(i+1, fitterp->GetMeanUncertainty()*1000);
     
        hsigma->SetBinContent(i+1,fitterp->GetSigma()*1000);
        hsigma->SetBinError(i+1, fitterp->GetSigmaUncertainty()*1000);
        
        hrelErr->SetBinContent(i+1,serr/s);
        hsign->SetBinContent(i+1,signf);
        hsign->SetBinError(i+1,signferr);
        hsb->SetBinContent(i+1,sob);
        hsb->SetBinError(i+1,soberr);
        
        hSignal->SetBinContent(i+1,s);
        hSignal->SetBinError(i+1,serr);
        
    }
    
    if(!fRawPtSpectrum) {
      std::cout << "Error in extracting raw jet pT spectrum plot! Exiting..." << std::endl;
      return kFALSE;
    }   
    if(fEffCorr && !fEffCorrPtSpectrum) {
      std::cout << "Error in extracting efficiency corrected jet pT spectrum plot! Exiting..." << std::endl;
      return kFALSE;
    }   
 
    TFile *ofile = new TFile(Form("%s/JetPtSpectra_SB_%s_%s_ptD%d.root",outdir.Data(),prod.Data(), bEff ? "eff" : "noEff", (int)ptDbins[0]),"RECREATE");
    hmean->Write();
    hsigma->Write();
    hsign->Write();
    hsb->Write();
    hrelErr->Write();
    
    hjetptspectrum->Write();
    hSignal->Write();
    hjetptspectrumReb->Write();
    
    for(int i=0; i<ptbinsDN; i++){
        hjetpt[i]->Write();
        hjetpt_s[i]->Write();
        hjetptsub[i]->Write();
        hmass[i]->Write();
        hmass_l[i]->Write();
        hmass_u[i]->Write();
        hmass_c[i]->Write();
        fullfit[i]->Write();
        hjetptcorr[i]->Write();
    }
    
    ofile->Close();
 
 
    return kTRUE;
       
}

Bool_t AliDJetsPtSpectrum::SetSigmaBkg(Int_t nbins, Double_t *sigma){
  
    //if( (sizeof(sigma)/sizeof(*sigma)) != 4){
    if( nbins != 4){
        std::cout << "!!!! sigma bkg requires 4 values !!!!!" << std::endl;
        return;
    }
    
    for(int i=0;i<4;i++) fnSigmaBkg[i] = sigma[i];
    
}

//___________________________________________________________________________________________
void AliDJetsPtSpectrum::SetDmesonPtBins(Int_t nbins, Double_t* ptedges) {

  if(!nbins) return;
  fnDbins=nbins;
  fDbinpTedges = new Double_t[fnDbins+1];
  for(int i=0;i<fnDbins+1;i++) {
    fDbinpTedges[i]=ptedges[i];
  }

  return;
}

//___________________________________________________________________________________________
void AliDJetsPtSpectrum::SetJetPtBins(Int_t nbins, Double_t* ptedges) {

  if(!nbins) return;
  fnJetbins=nbins;
  fJetbinpTedges = new Double_t[fnJetbins+1];
  for(int i=0;i<fnJetbins+1;i++) {
    fJetbinpTedges[i]=ptedges[i];
  }

  return;
}


void AliDJetsPtSpectrum::setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Size_t size, Width_t width, int scale){
    
    if(scale)h->Scale(1,"width");
    h->SetMarkerStyle(Mstyle);
    h->SetMarkerColor(color);
    h->SetMarkerSize(size);
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetTitle(0);
 
    return;

}

void AliDJetsPtSpectrum::SaveCanvas(TCanvas *c, TString name){
    
    c->SaveAs(Form("%s.png",name.Data());
    c->SaveAs(Form("%s.pdf",name.Data());
}

void AliDJetsPtSpectrum::setStyle(){
    
    gStyle->SetOptStat(000);
    gStyle->SetLegendFont(42);
    //gStyle->SetLegendTextSize(0.05);
    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetPadRightMargin(0.02);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetTitleOffset(1.,"x");
    gStyle->SetTitleOffset(0.9.,"y");
    gStyle->SetTitleSize(0.04,"xyz");
    gStyle->SetLabelSize(0.03,"xyz");
    
}

