//
//  Execute with:
//  gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ANALYSIS/macros -I$ROOTSYS/include");
//  gROOT->LoadMacro("AliDJetRawYieldUncertaintyLocal.cxx++")
//  .L ExtractDJetRawYieldUncertainty.C
//  EvaluateBinPerBinUncertainty(...) //to be done for each pT bin in which you have a mass spectrum
//  ExtractDJetRawYieldUncertainty(...) //to build the uncertainty for the various bins of the jet pT spectrum
// 


//double sigmaD[] = {0.000557321, 0.00056017, 0.000548255, 0.00055629, 0.000582418, 0.000598629, 0.000598187, 0.000635046, 0.000704751, 0.000772872 };
double sigmaD[] = { 0.000555598, 0.000559867, 0.000551223, 0.000557854, 0.000580121, 0.000584797, 0.000598181, 0.000627724, 0.000693677, 0.000739915 };
double sigmajet[] = {0.000523, 0.000501, 0.000503, 0.000525, 0.000524, 0.000564, 0.000564, 0.000564};

const int ptbinsDN = 10;
double ptDbins[ptbinsDN+1] = { 3,4,5,6,7,8,10,12,16,24,36 };
  
  
void EvaluateBinPerBinUncertainty( 
   Int_t bin = 0,
   Int_t specie=AliDJetRawYieldUncertaintyLocal::kDStarD0pi,  //D-meson decay channel
   Int_t method=AliDJetRawYieldUncertaintyLocal::kSideband, // kSideband,  //kEffScale,  //yield extraction methodliDJetRawYieldUncertainty::kDStarD0pi
   Double_t zmin=-1.,   //lower z edge
   Double_t zmax=2.,   //upper z edge
   Bool_t refl=kFALSE
   )
{

  double ptmin = ptDbins[bin];
  double ptmax = ptDbins[bin+1];


  AliDJetSignalExtraction *interface = new AliDJetSignalExtraction();
  Bool_t flagSpecie = interface->SetDmesonSpecie(specie);
  if(!flagSpecie) return;
  interface->SetYieldMethod(method);
  interface->SetPtBinEdgesForMassPlot(ptmin,ptmax);
  interface->SetZedges(zmin,zmax);
  interface->SetFitReflections(refl);
  

  if(specie==0) SetInputParametersDzero(interface);  // check the names and the values in the method!!
  else if(specie==1) SetInputParametersDstar(interface);  // check the names and the values in the method!!
  else if {printf("Error in setting the D-meson specie! Exiting...\n"); return kFALSE;}

  interface->SetDebugLevel(2); //0 = just do the job; 1 = additional printout; 2 = print individual fits

  Bool_t extract = interface->ExtractInputMassPlot();
  if(!extract) {
    printf("Error in extracting the mass plot! Exiting...\n");
    return;
  }
  
  

  interface->ClearObjects();

  return;
}


//________________________________________
void SetInputParametersDstar(AliDJetRawYieldUncertaintyLocal *interface){


  Int_t nDbins = 10;
  Double_t ptDbins[11] = {3,4,5,6,7,8,10,12,16,24,36};
  
  Int_t nJetbins = 9;
 // Double_t ptJetbins[9] = {2,4,6,8,10,12,16,24,40};
  Double_t ptJetbins[10] = {3,4,5,6,8,10,14,20,30,50};

  Double_t DMesonEff[10]; 
  TFile *FileEff = new TFile("/home/basia/Work/alice/analysis/pPb_run2/Efficiencies/out_806Preliminary/DjetEff_prompt_jetpt2_50.root");
	TH1F *hEff = (TH1F*)FileEff->Get("hEff_reb");
	for(int i=0;i<ptbinsDN;i++){
		double pt = (ptDbins[i]+ptDbins[i+1]) / 2.;
		double eff = hEff->GetBinContent(hEff->GetXaxis()->FindBin(pt));
		DMesonEff[i] = eff;
	}

    
  interface->SetInputFilename("/home/basia/Work/alice/analysis/pPb_run2/outData/cuts806_preliminary/AnalysisResults_FASTwoSDD.root");
  interface->SetInputDirname("DmesonsForJetCorrelations");
  interface->SetInputListname("histosDStarMBN");
  interface->SetInputObjectname("hsDphiz");
   
  interface->SetDmesonPtBins(nDbins,ptDbins);
  interface->SetJetPtBins(nJetbins,ptJetbins);
  interface->SetDmesonEfficiency(DMesonEff);
  //interface->SetSigmaToFixDPtBins(sigmaD);
  //interface->SetSigmaToFixJetPtBins(sigmajet);


  return;
}


void setStyle(){
    
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
