 //-----------------------------------------------------------------------
 //  Author B.Trzeciak
 //  Utrecht University
 //  barbara.antonina.trzeciak@cern.ch
 //-----------------------------------------------------------------------

 const int     fptbinsDN = 12;
 double        fptbinsDA[fptbinsDN+1] = { 1,2,3,4,5,6,7,8,10,12,16,24,36 };

void signalExtraction_refl(
  TString data = "$HOME/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_fast_R03_D0MC_def.root",
  TString out = "$HOME/Work/alice/analysis/pp5TeV/D0jet/outMC/reflections",
  bool postfix = 1, TString listName = "cut2",
  bool isMoreFiles = 0,
  TString prod = "kl"   // for more than 1 file, for one file leave it empty)
)
{

  int fRebinMass = 2;
  double jetmin = -10, jetmax = 50;
  double zmin = -2, zmax = 2;
  double minf = 1.65, maxf = 2.1;

    TString plotsDir = "/plots";
    TString outdir = out;
    gSystem->Exec(Form("mkdir %s",outdir.Data()));
    gSystem->Exec(Form("mkdir %s%s",outdir.Data(),plotsDir.Data()));

    if(!isMoreFiles) prod="";
    int nFiles = (int)prod.Length();

    TString histName = "histosD0MBN";
    // get analysis output file
    TString datafile;
    TFile *File;
    TDirectoryFile* dir;
    TList *histList;
    THnSparseF *sparse;

    TH2D *hInvMassptDSig, *hInvMassptDRefl;
    if(!isMoreFiles) {
      datafile = data;
      File = new TFile(datafile,"read");
      if(!File) { cout << "==== WRONG FILE WITH DATA =====\n\n"; return ;}
      dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");

      for(int i=0;i<3; i++){
          if(postfix) histList =  (TList*)dir->Get(Form("%s%d%sMCrec",histName.Data(),i,listName.Data()));
          else histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
          sparse = (THnSparseF*)histList->FindObject("hsDphiz");
          sparse->GetAxis(0)->SetRangeUser(zmin,zmax);
          sparse->GetAxis(1)->SetRangeUser(jetmin,jetmax);
          if(i==0) {
            hInvMassptDSig = (TH2D*)sparse->Projection(2,6);
            hInvMassptDRefl = (TH2D*)sparse->Projection(2,7);
          }
          else {
            hInvMassptDSig->Add((TH2D*)sparse->Projection(2,6));
            hInvMassptDRefl->Add((TH2D*)sparse->Projection(2,7));
          }
      }
    }
    else {
      for (int j=0;j<nFiles;j++){
          datafile = data;
          datafile += prod.Data()[j];
          datafile += ".root";
          File = new TFile(datafile,"read");
          if(!File) { cout << "==== WRONG FILE WITH DATA =====\n\n"; return ;}
          dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");

          for(int i=0;i<2; i++){
              if(postfix) histList =  (TList*)dir->Get(Form("%s%d%sMCrec",histName.Data(),i,listName.Data()));
              else histList =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
              sparse = (THnSparseF*)histList->FindObject("hsDphiz");
              sparse->GetAxis(0)->SetRangeUser(zmin,zmax);
              sparse->GetAxis(1)->SetRangeUser(jetmin,jetmax);
              if(i==0 && j==0) {
                hInvMassptDSig = (TH2D*)sparse->Projection(2,6);
                hInvMassptDRefl = (TH2D*)sparse->Projection(2,7);
              }
              else {
                hInvMassptDSig->Add((TH2D*)sparse->Projection(2,6));
                hInvMassptDRefl->Add((TH2D*)sparse->Projection(2,7));
              }
          }
      }
    }

    TH1D *hsig[fptbinsDN], *hrefl[fptbinsDN];
    TH1F *hFitReflNewTemp[fptbinsDN], *ratio[fptbinsDN];
    TString formulaSig = "[0]/([2]*TMath::Sqrt(2*TMath::Pi()))*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))";
    TString formulaRef = "[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]/( TMath::Sqrt(2.*TMath::Pi())*[5])*TMath::Exp(-(x-[4])*(x-[4])/(2.*[5]*[5]))";

    TCanvas *cSig = new TCanvas("cSig","cSig",1600,1800);
    cSig->Divide(4,3);
    TCanvas *cRefl = new TCanvas("cRefl","cRefl",1600,1800);
    cRefl->Divide(4,3);
    TCanvas *cRefl2 = new TCanvas("cRefl2","cRefl2",1600,1800);
    cRefl2->Divide(4,3);
    TCanvas *cRatio = new TCanvas("cRatio","cRatio",1600,1800);
    cRatio->Divide(4,3);

    for(int i=0; i<fptbinsDN; i++){
      hsig[i] = (TH1D*)hInvMassptDSig->ProjectionX(Form("histSgn_%d_%d",(int)fptbinsDA[i],(int)fptbinsDA[i+1]),hInvMassptDSig->GetYaxis()->FindBin(fptbinsDA[i]), hInvMassptDSig->GetYaxis()->FindBin(fptbinsDA[i+1])-1);
      hrefl[i] = (TH1D*)hInvMassptDRefl->ProjectionX(Form("hrefl_%d_%d",(int)fptbinsDA[i],(int)fptbinsDA[i+1]),hInvMassptDRefl->GetYaxis()->FindBin(fptbinsDA[i]), hInvMassptDRefl->GetYaxis()->FindBin(fptbinsDA[i+1])-1);
      hsig[i]->Rebin(fRebinMass);
      hrefl[i]->Rebin(fRebinMass);

      hsig[i]->GetXaxis()->SetRangeUser(minf,maxf);
      hrefl[i]->GetXaxis()->SetRangeUser(minf,maxf);
      hsig[i]->SetTitle(Form("%.1lf < pt^{D} < %.1lf",fptbinsDA[i],fptbinsDA[i+1]));
      hrefl[i]->SetTitle(Form("%.1lf < pt^{D} < %.1lf",fptbinsDA[i],fptbinsDA[i+1]));

      TF1 *gaussMCSignal=new TF1("gaussMCSig",formulaSig.Data(),minf,maxf);
      gaussMCSignal->SetParName(0,"IntegralSgn");
      gaussMCSignal->SetParName(1,"Mean");
      gaussMCSignal->SetParName(2,"Sigma");
      gaussMCSignal->SetParameter(0,1);
      gaussMCSignal->SetParameter(1,1.864);
      gaussMCSignal->SetParameter(2,0.010);
      gaussMCSignal->SetLineColor(kOrange+2);
      cSig->cd(i+1);
      gStyle->SetOptFit(11111);
      hsig[i]->Fit("gaussMCSig","RI","",1.65,2.15);
      hsig[i]->Draw();

      cRefl->cd(i+1);
      TF1 *doublegaussMCRefl=new TF1("doublegaussMCRefl",formulaRef.Data(),minf,maxf);
      doublegaussMCRefl->SetParName(0,"IntegralRefl");
    //  doublegaussMCRefl->SetParName(1,"Mean");
    //  doublegaussMCRefl->SetParName(2,"Sigma");
    //  doublegaussMCRefl->SetParameter(0,1);
    //  doublegaussMCRefl->SetParameter(1,1.864);
    //  doublegaussMCRefl->SetParameter(2,0.010);
      doublegaussMCRefl->SetLineColor(kRed+2);
      doublegaussMCRefl->SetParameter(0, 1);
      doublegaussMCRefl->SetParameter(1, 1);
      doublegaussMCRefl->SetParameter(2, 1);
      doublegaussMCRefl->SetParameter(3, 1);
      doublegaussMCRefl->SetParameter(4, 1);
      doublegaussMCRefl->SetParameter(5, 1);
      hrefl[i]->Fit("doublegaussMCRefl", "MLFR");
      hrefl[i]->Draw();

      cRefl2->cd(i+1);
      TF1 *fFitRefl = hrefl[i]->GetFunction("doublegaussMCRefl");
      //fFitReflection->cd();
      hFitReflNewTemp[i] = (TH1F*)hrefl[i]->Clone(Form("histRflFittedDoubleGaus_pt%d_%d",(int)fptbinsDA[i],(int)fptbinsDA[i+1]));
      ratio[i] = (TH1F*)hrefl[i]->Clone(Form("ratioRelDistr_pt%d_%d", (int)fptbinsDA[i],(int)fptbinsDA[i+1]));

      for(Int_t iBin2=1; iBin2<=hrefl[i]->GetNbinsX(); iBin2++){
        hFitReflNewTemp[i]->SetBinContent(iBin2, 0.);
        ratio[i]->SetBinContent(iBin2, 0.);

        hFitReflNewTemp[i]->SetBinContent(iBin2, fFitRefl->Eval(hrefl[i]->GetBinCenter(iBin2)));
        ratio[i]->SetBinContent(iBin2, (hrefl[i]->GetBinContent(iBin2) / fFitRefl->Eval(hrefl[i]->GetBinCenter(iBin2))));

      }
      cRefl2->cd(i+1);
      hFitReflNewTemp[i]->Draw();

      cRatio->cd(i+1);
      ratio[i]->GetYaxis()->SetRangeUser(-1.5, 3.);
      ratio[i]->Draw();
      ratio[i]->Fit("pol0", "FM");

    }


    // --------------------------------------------------------
    // ----------- write to output file
    TFile *ofile = new TFile(Form("%s/reflectionTemplates_%s.root",outdir.Data(), postfix ? listName.Data() : "pPb" ),"RECREATE");
    for(int i=0; i<fptbinsDN; i++){
      hsig[i]->Write();
      hrefl[i]->Write();
      hFitReflNewTemp[i]->Write();
      ratio[i]->Write();
    }
    ofile->Close();
    // --------------------------------------------------------

}
