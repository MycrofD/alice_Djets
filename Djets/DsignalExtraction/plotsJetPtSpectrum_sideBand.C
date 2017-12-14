#include "style.C"
#include <string>
#include <sstream>
#include <iostream>
    
setHistoDetails(TH1 *h, Color_t color, Style_t Mstyle, Width_t width);
    
    double zmin = -2, zmax = 2.;
    double jetmin = -5, jetmax = 50;
    double plotmin = 0, plotmax = 30;
    double jetplotmin = 0, jetplotmax = 30;
    
    //------- D pT bins
  
    const int ptbinsDN = 10;
    float ptDbins[ptbinsDN+1] = { 3,4,5,6,7,8,10,12,16,24,36 };
    
  //  const int ptbinsDN = 11;
   // float ptDbins[ptbinsDN+1] = { 2,3,4,5,6,7,8,10,12,16,24,36 };
   
    
    //------- jet pT bins
   /* const int ptbinsN = 9;
    double ptbins[ptbinsN+1] = { 0,2,4,6,8,10,12,16,24,40 };
    //float ptbins[ptbinsN+1] = { 0,2,4,6,8,10,12,15,30,50 };
    
    const int ptbinsN = 8;
    double ptbins[ptbinsN+1] = { 2,4,6,8,10,12,16,24,40 };
    */
    
    const int ptbinsN = 9;
    const double ptbins[ptbinsN+1] = { 3,4,5,6,8,10,14,20,30,50 };
    
   // const int ptbinsN = 13;
   // const double ptbins[ptbinsN+1] = {-5,0,1,2,3,4,5,6,8,10,14,20,30,50 };
    
    
    int ND = 4; // number of D mesons to analyse
    
    bool savePlots = 1;
    bool bEff = 0;
    bool isptcut = 1;
    bool isdetails = 0;
    
    
    //------- efficiency file
    double efficiency[ptbinsDN];
    
    
    // mass fit params
    Float_t minf = 0.140;
    Float_t maxf = 0.160;
    Int_t   bkgtype = 5;
    Float_t mass = 0.145;
    Float_t sigma = 0.00065;
   
 
    TH1F* hmass[ptbinsDN];
    TH1F *hmass_l[ptbinsDN];
    TH1F *hmass_u[ptbinsDN];
    TH1F *hmass_c[ptbinsDN];
    TF1* fullfit[ptbinsDN];
    
    TH1F* hjetpt[ptbinsDN];
    TH1F *hjetpt_s[ptbinsDN];
    TH1F *hjetptsub[ptbinsDN];
    TH1F *hjetptcorr[ptbinsDN];
    
    TH1F *hjetptspectrum;
    TH1F *hrawjetptspectrum;
    
    TH1F *hmean;
    TH1F *hsigma;
    TH1F *hrelErr;
    TH1F *hsign;
    TH1F *hsb;
    TH1F *hSignal;
    
    TH3D* hInvMassptD;

    
    bool isEta = 0;
    double jetEta = 0.5;
    

void plotsJetPtSpectrum_sideBand(TString efffile = "/home/basia/Work/alice/analysis/pPb_run2/Efficiencies/out_806BaseNewEffRM_oldImpr/DjetEff_prompt_jetpt2_50.root", TString dataFile = "../outData/",  TString outdir="", bool isEff = 1, bool save = 1, bool postfix = 0, TString listName = "Cut",  int Rpar = 4 )
{


    setStyle();
    savePlots = save;
    bEff = isEff;
    
    // get analysis output file
    TString datafile = dataFile;
  
    
	TFile *File = new TFile(datafile,"read");
	TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
    TList *histList;
    THnSparseF *sparse;
    for(int i=0;i<ND; i++){
        if(postfix) histList =  (TList*)dir->Get(Form("histosDStarMBN%d%s",i,listName.Data()));
        else histList =  (TList*)dir->Get(Form("histosDStarMBN%d",i));
        sparse = (THnSparseF*)histList->FindObject("hsDphiz");
        sparse->GetAxis(0)->SetRangeUser(zmin,zmax); 
        sparse->GetAxis(1)->SetRangeUser(jetmin,jetmax);
        if(isEta) sparse->GetAxis(5)->SetRangeUser(-jetEta,jetEta);
        if(!i) hInvMassptD=(TH3D*)sparse->Projection(3,1,2);
        else {
                TH3D *h3 = (TH3D*)sparse->Projection(3,1,2);
                hInvMassptD->Add(h3);
        }
    }
    
      
    if(bEff){
	TFile *FileEff = new TFile(efffile.Data(),"read");
	TH1F *hEff = (TH1F*)FileEff->Get("hEff_reb");
	for(int i=0;i<ptbinsDN;i++){
		double pt = (ptDbins[i]+ptDbins[i+1]) / 2.;
		double eff = hEff->GetBinContent(hEff->GetXaxis()->FindBin(pt));
		efficiency[i] = eff;
	}
    }

    // --------------------------------------------------------
    // fit inv. mass in D pT bins and get raw jet pT spectra in the signal and side-bands regions
    rawJetSpectra(outdir,prod);
    
   
   //------------------- draw output histos -------------------------------
   
   // draw and save fit params
   saveFitParams(outdir,prod);
   
   
    TPaveText *pvProd = new TPaveText(0.75,0.65,0.9,0.7,"brNDC");
    pvProd->SetFillStyle(0);
    pvProd->SetBorderSize(0);
    pvProd->AddText(Form("%s",prod.Data()));
    
    TPaveText *pvCuts = new TPaveText(0.75,0.6,0.9,0.65,"brNDC");
    pvCuts->SetFillStyle(0);
    pvCuts->SetBorderSize(0);
    pvCuts->AddText(Form("%s",outdir.Data()));
    
    
     // before eff. correction
    TCanvas *cRawSpectrum = new TCanvas("cRawSpectrum","cRawSpectrum",800,600);
    setHistoDetails(hrawjetptspectrum,0,kBlue+1,20,1.2);
    hrawjetptspectrum->GetXaxis()->SetTitle("p_{T}^{ch,jet} (GeV/c)");
    hrawjetptspectrum->GetYaxis()->SetTitle("Raw dN/dp_{T}");
    hrawjetptspectrum->GetXaxis()->SetRangeUser(plotmin,plotmax);
    hrawjetptspectrum->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");

    
    TPaveText *pv3 = new TPaveText(0.65,0.65,0.92,0.85,"brNDC");
    pv3->SetFillStyle(0);
    pv3->SetBorderSize(0);
    if(!isptcut) pv3->AddText(Form("charged jets, p-Pb, R=0.%d",Rpar));
    else pv3->AddText(Form("#splitline{charged jets, p-Pb, R=0.%d}{p_{T,D*}>%dGeV/c}",Rpar,(Int_t)ptDbins[0]));
     

    hjetptspectrum->GetYaxis()->SetTitle("dN/dp_{T}");
    hjetptspectrum->GetXaxis()->SetRangeUser(plotmin,plotmax);
    
    TH1F *hjetptspectrumReb_tmp = (TH1F*)hjetptspectrum->Clone("hjetptspectrumReb_tmp");
    TH1F *hjetptspectrumReb = (TH1F*)hjetptspectrumReb_tmp->Rebin(ptbinsN,"hjetptspectrumReb",ptbins);
     
    TCanvas *cSpectrum = new TCanvas("cSpectrum","cSpectrum",800,600);
    cSpectrum->SetLogy();
    setHistoDetails(hjetptspectrum,0,kRed,20,1.2);
    hjetptspectrum->GetXaxis()->SetTitle("p_{T,ch jet} (GeV/c)");
    hjetptspectrum->Draw();
    pv3->Draw("same");
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    if(savePlots) SaveCanvas(cSpectrum,outdir+"/plots/jetPtSpectrum_SB_"+prod);
     
    setHistoDetails(hjetptspectrumReb,1,kBlue,20,1.2);
    hjetptspectrumReb->GetXaxis()->SetTitle("p_{T,ch jet} (GeV/c)");
    TCanvas *cSpectrumRebin = new TCanvas("cSpectrumRebin","cSpectrumRebin",800,600);
    cSpectrumRebin->SetLogy();
    hjetptspectrumReb->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    pv3->Draw("same");
    
     if(savePlots) SaveCanvas(cSpectrumRebin,outdir+"/plots/jetPtSpectrum_SB_Rebin_"+prod);
    
    TH1D *hjetptspectrumRebUnc = (TH1D*)hjetptspectrumReb->Clone("hjetptspectrumRebUnc");
    hjetptspectrumRebUnc->GetYaxis()->SetTitle("Rel. unc.");
    
     for(int j=0; j<=hjetptspectrumReb->GetNbinsX();j++){
                double err;
                if(hjetptspectrumReb->GetBinContent(j)) err = hjetptspectrumReb->GetBinError(j)/hjetptspectrumReb->GetBinContent(j);
                else err = 0;
                hjetptspectrumRebUnc->SetBinContent(j,err);
                hjetptspectrumRebUnc->SetBinError(j,0);
    }
    
    TCanvas *cSpectrumRebinUnc = new TCanvas("cSpectrumRebinUnc","cSpectrumRebinUnc",800,600);

    hjetptspectrumRebUnc->GetXaxis()->SetRangeUser(jetplotmin,jetplotmax);
    hjetptspectrumRebUnc->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    pv3->Draw("same");
    
    
     if(savePlots) SaveCanvas(cSpectrumRebinUnc,outdir+"/plots/jetPtSpectrumUnc_SB_Rebin_"+prod);
    
   
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
   
}


void rawJetSpectra(TString outdir, TString prod){
    
    hmean = new TH1F("hmean","hmean",ptbinsDN,ptDbins);
    hsigma = new TH1F("hsigma","hsigma",ptbinsDN,ptDbins);
    hrelErr = new TH1F("hrelErr","hrelErr",ptbinsDN,ptDbins);
    hsign = new TH1F("hsign","hsign",ptbinsDN,ptDbins);
    hsb = new TH1F("hsb","hsb",ptbinsDN,ptDbins);
    hSignal = new TH1F("hSignal","hSignal",ptbinsDN,ptDbins);
    hSignal->Sumw2();
    
    hInvMassptD->GetXaxis()->SetTitle("m^{D*}");
    hInvMassptD->GetYaxis()->SetTitle("p_{T}^{jet}");
    hInvMassptD->SetTitle();

    setHistoDetails(hmean,0,2,20);
    hmean->GetYaxis()->SetTitle("signal #mu (MeV/c^{2})");
    setHistoDetails(hsigma,0,4,20);
    hsigma->GetYaxis()->SetTitle("signal #sigma (MeV/c^{2})");
    setHistoDetails(hsign,0,4,20);
    hsign->GetYaxis()->SetTitle("significance");
    setHistoDetails(hsb,0,8,20);
    hsb->GetYaxis()->SetTitle("S/B");
    setHistoDetails(hrelErr,0,6,20);
    hrelErr->GetYaxis()->SetTitle("rel. unc.");
    
    TPaveText *pvProd = new TPaveText(0.8,0.25,0.98,0.3,"brNDC");
    pvProd->SetFillStyle(0);
    pvProd->SetBorderSize(0);
    pvProd->AddText(Form("%s",prod.Data()));
    
    TPaveText *pvCuts = new TPaveText(0.8,0.3,0.98,0.35,"brNDC");
    pvCuts->SetFillStyle(0);
    pvCuts->SetBorderSize(0);
    pvCuts->AddText(Form("%s",outdir.Data()));
    
    
    TCanvas *c2 = new TCanvas("c2","c2",1200,1200);
    c2->Divide(3,4);
    TCanvas *c2jet = new TCanvas("c2jet","c2jet",1200,1200);
    c2jet->Divide(3,4);
    TCanvas *c2jetcorr = new TCanvas("c2jetcorr","c2jetcorr",1200,1200);
    c2jetcorr->Divide(3,4);
    
    
    for(int i=0; i<ptbinsDN; i++){
                 
        TH1F *hh=(TH1F*)hInvMassptD->ProjectionX(Form("hh_%d",i),hInvMassptD->GetYaxis()->FindBin(jetmin), hInvMassptD->GetYaxis()->FindBin(jetmax)-1,hInvMassptD->GetZaxis()->FindBin(ptDbins[i]), hInvMassptD->GetZaxis()->FindBin(ptDbins[i+1])-1);
          //hh->Rebin(2);
        hh->GetXaxis()->SetRangeUser(0.14,0.155);
        hh->SetTitle(Form("%.1lf < pt^{D*} < %.1lf",ptDbins[i],ptDbins[i+1]));
      

        c2->cd(i+1);
        //hh->Draw();

        TH1F *hmassfit = (TH1F*)hh->Clone("hmassfit");

        float hmin=TMath::Max(minf,hmassfit->GetBinLowEdge(2));
        float hmax=TMath::Min(maxf,hmassfit->GetBinLowEdge(hmassfit->GetNbinsX()));
        AliHFMassFitter* fitterp=new AliHFMassFitter((TH1F*)hmassfit,hmin,hmax,1,bkgtype,0);
        fitterp->SetInitialGaussianMean(mass);
        fitterp->SetInitialGaussianSigma(sigma);
        fitterp->MassFitter(kFALSE);
        //fitterp->Draw("same");
        TH1F* h=fitterp->GetHistoClone();
        fullfit[i]=h->GetFunction("funcmass");
        fullfit[i]->SetName(Form("fullfit_%d",i));
        ((TList*)h->GetListOfFunctions())->RemoveAt(2);
        hmass[i] = (TH1F*)h->Clone(Form("hmass_%",i));
        hmass[i]->SetName(Form("hmass_%d",i));
        hmass[i]->Draw("PE");
        hmass[i]->GetYaxis()->SetTitle("Entries");
        hmass[i]->GetYaxis()->SetTitleOffset(1.4);
        hmass[i]->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
        fullfit[i]->Draw("same");
        
         
        float Dsigma = fitterp->GetSigma();
        float Dmean = fitterp->GetMean();
        float signal_l_min = Dmean-8*Dsigma;
        float signal_l_max = Dmean-5*Dsigma;
        float signal_u_min = Dmean+5*Dsigma;
        float signal_u_max = Dmean+13*Dsigma;
        float signal_c_min = Dmean-3*Dsigma;
        float signal_c_max = Dmean+3*Dsigma;
        
        int binmin = hmass[i]->GetXaxis()->FindBin(signal_c_min);
        int binmax = hmass[i]->GetXaxis()->FindBin(signal_c_max);
        double binwidth = hmass[i]->GetXaxis()->GetBinWidth(1)*0.5;
        
        TPaveText *pv=new TPaveText(0.48,0.7,0.9,0.9,"brNDC");
        pv->SetFillStyle(0);
        pv->SetBorderSize(0);
        Double_t s,serr,b,berr,signf,signferr;
        fitterp->Signal(3,s,serr);
        cout << "Signal = " << s << "Error = " << serr << endl;
        //fitterp->Background(3,b,berr);
        fitterp->Background(hmass[i]->GetXaxis()->GetBinCenter(binmin)-binwidth, hmass[i]->GetXaxis()->GetBinCenter(binmax-1)+binwidth ,b,berr);
      
        fitterp->Significance(3,signf,signferr);
        Double_t sob=s/b, soberr;
        soberr=TMath::Sqrt((serr/b)*(serr/b) + (s/b/b*berr)*(s/b/b*berr));
        pv->AddText(Form("Signif.(3#sigma) = (%.1f #pm %.1f)", signf,signferr));
        pv->AddText(Form("#mu = (%.2f #pm %.2f) MeV/#it{c}^{2}", fitterp->GetMean()*1000,fitterp->GetMeanUncertainty()*1000));
        pv->AddText(Form("#sigma = (%.2f #pm %.2f) MeV/#it{c}^{2}", fitterp->GetSigma()*1000,fitterp->GetSigmaUncertainty()*1000));

        Bool_t twodigits=kTRUE;
        if(soberr*100. > 35.) twodigits=kFALSE;
        //if(twodigits) pv->AddText(Form("S/B (3#sigma) = (%.2f #pm %.2f)", sob,soberr));
        //else pv->AddText(Form("S/B (3#sigma) = (%.1f #pm %.1f)", sob,soberr));
        if(twodigits) pv->AddText(Form("S (3#sigma) = (%.2f #pm %.2f)", s,serr));
        else pv->AddText(Form("S (3#sigma) = (%.1f #pm %.1f)", s,serr));
        pv->Draw("same");
    
        if(isdetails) pvProd->Draw("same");
        if(isdetails) pvCuts->Draw("same");
       
        hmass_l[i] = (TH1F*)hmass[i]->Clone("hmass_l");
        hmass_l[i]->GetXaxis()->SetRangeUser(signal_l_min,signal_l_max);
        hmass_l[i]->SetName(Form("hmass_l_%d",i));
        hmass_u[i] = (TH1F*)hmass[i]->Clone("hmass_u");
        hmass_u[i]->GetXaxis()->SetRangeUser(signal_u_min,signal_u_max);
        hmass_u[i]->SetName(Form("hmass_u_%d",i));
        hmass_c[i] = (TH1F*)hmass[i]->Clone("hmass_c");
        hmass_c[i]->GetXaxis()->SetRangeUser(signal_c_min,signal_c_max);
        hmass_c[i]->SetName(Form("hmass_c_%d",i));
        
        hmass_l[i]->SetFillColor(kBlue+2);
        hmass_u[i]->SetFillColor(kBlue+2);
        hmass_c[i]->SetFillColor(kRed+2);
        hmass_l[i]->SetFillStyle(3004);
        hmass_u[i]->SetFillStyle(3004);
        hmass_c[i]->SetFillStyle(3005);
        
        hmass_l[i]->Draw("hsame");
        hmass_u[i]->Draw("hsame");
        hmass_c[i]->Draw("hsame");
        
        //-------- jet pt spectrum - signal
        hjetpt[i]=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_%d",i),hInvMassptD->GetXaxis()->FindBin(signal_c_min), hInvMassptD->GetXaxis()->FindBin(signal_c_max)-1,hInvMassptD->GetZaxis()->FindBin(ptDbins[i]), hInvMassptD->GetZaxis()->FindBin(ptDbins[i+1])-1);
        hjetpt[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
        hjetpt[i]->SetTitle(Form("%.1lf < pt^{D*} < %.1lf",ptDbins[i],ptDbins[i+1]));
        hjetpt[i]->SetMarkerColor(kRed+2);
        hjetpt[i]->SetLineColor(kRed+2);
        
        //------ jet pt spectrum - side bands
        TH1F* hjetpt_s1=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_s1%d",i),hInvMassptD->GetXaxis()->FindBin(signal_l_min), hInvMassptD->GetXaxis()->FindBin(signal_l_max)-1,hInvMassptD->GetZaxis()->FindBin(ptDbins[i]), hInvMassptD->GetZaxis()->FindBin(ptDbins[i+1])-1);
        hjetpt_s1->SetTitle(Form("%.1lf < pt^{D*} < %.1lf",ptDbins[i],ptDbins[i+1]));
        hjetpt_s1->SetMarkerColor(kBlue+2);
        hjetpt_s1->SetLineColor(kBlue+2);
        TH1F* hjetpt_s2=(TH1F*)hInvMassptD->ProjectionY(Form("hjetpt_s2%d",i),hInvMassptD->GetXaxis()->FindBin(signal_u_min), hInvMassptD->GetXaxis()->FindBin(signal_u_max)-1,hInvMassptD->GetZaxis()->FindBin(ptDbins[i]), hInvMassptD->GetZaxis()->FindBin(ptDbins[i+1])-1);
        hjetpt_s2->SetTitle(Form("%.1lf < pt^{D*} < %.1lf",ptDbins[i],ptDbins[i+1]));
        hjetpt_s2->SetMarkerColor(kBlue+2);
        hjetpt_s2->SetLineColor(kBlue+2);
        hjetpt_s[i] = (TH1F*)hjetpt_s1->Clone(Form("hjetpt_s_%d",i));
        hjetpt_s[i]->Add(hjetpt_s2);
        
        
        // scale background from side bands to the background under the peak
        double scaling = b / hjetpt_s[i]->Integral();
        hjetpt_s[i]->Scale(scaling);
        
        
        c2jet->cd(i+1);
        gPad->SetLogy();
        
        hjetpt[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
        hjetpt_s[i]->GetXaxis()->SetRangeUser(plotmin,plotmax);
        hjetpt[i]->GetXaxis()->SetTitle("p_{T,ch jet} (GeV/c)");
        
        hjetpt[i]->Draw("ep");
        hjetpt_s[i]->Draw("epsame");
        
        TLegend *l1 = new TLegend(0.75,0.6,0.9,0.8);
        l1->AddEntry(hjetpt[i],"signal","l");
        l1->AddEntry(hjetpt_s[i],"SB","l");
        l1->Draw("same");
        
        
        //------- subtract background from signal jet
        hjetptsub[i] = (TH1F*)hjetpt[i]->Clone(Form("hjetptsub_%d",i));
        hjetptsub[i]->Add(hjetpt_s[i],-1);
        hjetptsub[i]->SetMarkerColor(kGreen+3);
        hjetptsub[i]->SetLineColor(kGreen+3);
        
        //c2jetsub->cd(i+1);
        hjetptsub[i]->Draw("epsame");
        if(isdetails) pvProd->Draw("same");
        if(isdetails) pvCuts->Draw("same");
        
        TLegend *l2 = new TLegend(0.6,0.75,0.9,0.9);
        l2->AddEntry(hjetptsub[i],"bkg subtracted","l");
        l1->AddEntry(hjetptsub[i],"sig-SB","l");
        //l2->Draw("same");
        
        if(!i) hrawjetptspectrum = (TH1F*)hjetptsub[i]->Clone("hrawjetptspectrum");
        else hrawjetptspectrum->Add(hjetptsub[i]);
        
        //------- correct for D* efficiency
        hjetptcorr[i] = (TH1F*)hjetptsub[i]->Clone(Form("hjetptcorr_%d",i));
        if(bEff) hjetptcorr[i]->Scale(1./efficiency[i]); // D efficiency
        hjetptcorr[i]->SetMarkerColor(kBlue+3);
        hjetptcorr[i]->SetLineColor(kBlue+3);
        //hjetptcorr->GetXaxis()->SetRangeUser(plotmin,plotmax);
        
        c2jetcorr->cd(i+1);
        hjetptcorr[i]->Draw("ep");
        if(isdetails) pvProd->Draw("same");
        if(isdetails) pvCuts->Draw("same");
        
        TLegend *l3 = new TLegend(0.6,0.7,0.9,0.9);
        l3->AddEntry(hjetptcorr[i],"#splitline{corrected jet p_{T}}{spectrum}","l");
        l3->Draw("same");
        
        //----- add corrected jet pt spectrum distributions in each D pt bin into a one distribution
        
        if(!i) hjetptspectrum = (TH1F*)hjetptcorr[i]->Clone("hjetptspectrum");
        else hjetptspectrum->Add(hjetptcorr[i]);
        
        
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
  
    if(savePlots) SaveCanvas(c2,outdir+"/plots/invMass_"+prod);
    if(savePlots) SaveCanvas(c2jet,outdir+"/plots/jetRawSpectrum"+prod);
    
    
}

void  saveFitParams(TString outdir,TString prod){
    
    TPaveText *pvProd = new TPaveText(0.78,0.75,0.9,0.8,"brNDC");
    pvProd->SetFillStyle(0);
    pvProd->SetBorderSize(0);
    pvProd->AddText(Form("%s",prod.Data()));
    
    TPaveText *pvCuts = new TPaveText(0.78,0.8,0.9,0.85,"brNDC");
    pvCuts->SetFillStyle(0);
    pvCuts->SetBorderSize(0);
    pvCuts->AddText(Form("%s",outdir.Data()));
    
    
     hmean->GetYaxis()->SetTitleOffset(1.4);
    hmean->GetYaxis()->SetTitleSize(0.04);
    hmean->GetYaxis()->SetLabelSize(0.03);
    
    hsigma->GetYaxis()->SetTitleOffset(1.4);
    hsigma->GetYaxis()->SetTitleSize(0.04);
    hsigma->GetYaxis()->SetLabelSize(0.03);
    
    hmean->GetYaxis()->SetRangeUser(144.8,146);
    hsigma->GetYaxis()->SetRangeUser(0.35,0.8);
    
    hSignal->SetName("hSignal");
    setHistoDetails(hSignal,0,2,20);
    hSignal->GetYaxis()->SetTitle("yield");
    hSignal->GetYaxis()->SetTitleOffset(1.2);
    hSignal->GetYaxis()->SetTitleSize(0.04);
    hSignal->GetYaxis()->SetLabelSize(0.03);
    hsign->GetYaxis()->SetTitleOffset(1.2);
    hsign->GetYaxis()->SetTitleSize(0.04);
    hsign->GetYaxis()->SetLabelSize(0.03);
    hrelErr->GetYaxis()->SetTitleOffset(1.2);
    hrelErr->GetYaxis()->SetTitleSize(0.04);
    hrelErr->GetYaxis()->SetLabelSize(0.03);
    
    
    double mean = 145.421;
    TLine *lm = new TLine(ptDbins[0], mean, ptDbins[ptbinsDN], mean);
    lm->SetLineStyle(2);
    
    TCanvas *cMassFit = new TCanvas("cMassFit","cMassFit",1200,600);
    cMassFit->Divide(2,1);
    cMassFit->cd(1);
    gPad->SetLeftMargin(0.15);
    hmean->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    lm->Draw("same");
    cMassFit->cd(2);
    hsigma->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    

     if(savePlots)  SaveCanvas(cMassFit,outdir+"/plots/gaussianParams_"+prod);
    
    TCanvas *cSignal = new TCanvas("cSignal","cSignal",1200,1200);
    cSignal->Divide(2,2);
    
    cSignal->cd(1);
    hSignal->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    cSignal->cd(2);
    hsign->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    cSignal->cd(3);
    hrelErr->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    cSignal->cd(4);
    hsb->Draw();
    if(isdetails) pvProd->Draw("same");
    if(isdetails) pvCuts->Draw("same");
    
     if(savePlots) SaveCanvas(cSignal,outdir+"/plots/signalParams_"+prod);
    
}
void setHistoDetails(TH1 *h, int scale, Color_t color, Style_t Mstyle, Size_t size = 0.9, Width_t width=2){
    
    if(scale)h->Scale(1,"width");
    h->SetMarkerStyle(Mstyle);
    h->SetMarkerColor(color);
    h->SetMarkerSize(size);
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetTitle(0);
    h->GetXaxis()->SetTitle("p_{T,D*}(GeV/c)");
    
    return;

}

void SaveCanvas(TCanvas *c, TString name = "tmp"){
    
    c->SaveAs(Form("%s_pTD%d.png",name.Data(),(int)ptDbins[0]));
    c->SaveAs(Form("%s_pTD%d.pdf",name.Data(),(int)ptDbins[0]));
   
   
   
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
    //style();
    
}
