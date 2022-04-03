import ROOT
ROOT.gROOT.SetBatch()
plotmin, plotmax = -1.0, 1.0

ptbinsJetN = 10;
ptJetbins] = [2,3,4,5,6,8,10,14,20,30,50]

"""
void drawResZ()
{

    style();
    gStyle->SetOptStat(000);

    gStyle->SetLegendFont(42);


   //TFile *inFile = new TFile("/home/kvapil/work/analysis/pp_run2/D0jet/BaseCuts/Default_AnalysisResults_Run2.root/FDsubtraction/JetPtSpectrum_FDsub.root","read");

   //TFile *File = new TFile("/mnt/hgfs/vmware/data_R04_050219/MC/AnalysisResults_Run2w18b.root","read");
   TFile *File = new TFile("/eos/user/a/amohanty/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_1061_R04_ppMC_5cuts.root","read");
   if(!File) { std::cout << "==== WRONG FILE WITH DATA =====\n\n"; return ;}
   TDirectoryFile *dir=dynamic_cast<TDirectoryFile*>(File->Get("DmesonsForJetCorrelations"));
   TString histName = "histosD0MBN";
   TList *histList;
   THnSparseF *sparse;

   TH1D* hbin1;
   TH1D* hbin2;
   TH1D* hbin3;


   for(int i=0;i<2; i++){
       histList =  dynamic_cast<TList*>(dir->Get(Form("%s%dMCrec",histName.Data(),i)));
       sparse = dynamic_cast<THnSparseF*>(histList->FindObject("ResponseMatrix"));
       sparse->GetAxis(5)->SetRangeUser(0.4,0.6);
       sparse->GetAxis(6)->SetRangeUser(7,10); //jet
       sparse->GetAxis(2)->SetRangeUser(3,10);
       sparse->GetAxis(7)->SetRangeUser(3,10);
       if(i==0) {
         hbin1 = dynamic_cast<TH1D*>(sparse->Projection(11));
       }
       else {
         hbin1->Add(dynamic_cast<TH1D*>(sparse->Projection(11)));
       }
     //  sparse->GetAxis(2)->SetRangeUser(5,15);
   //    sparse->GetAxis(7)->SetRangeUser(5,15);
       sparse->GetAxis(5)->SetRangeUser(0.7,0.8);
       if(i==0) {
         hbin2 = dynamic_cast<TH1D*>(sparse->Projection(11));
       }
       else {
         hbin2->Add(dynamic_cast<TH1D*>(sparse->Projection(11)));
       }
   //    sparse->GetAxis(2)->SetRangeUser(5,36);
    //   sparse->GetAxis(7)->SetRangeUser(5,36);
       sparse->GetAxis(5)->SetRangeUser(0.9,1.05);
       if(i==0) {
         hbin3 = dynamic_cast<TH1D*>(sparse->Projection(11));
       }
       else {
         hbin3->Add(dynamic_cast<TH1D*>(sparse->Projection(11)));
       }
   }

 //   TH1D *hFD_ratio = (TH1D*)inFile->Get("hFD_ratio");





            hbin1->SetTitle("");
            hbin1->SetMarkerColor(kRed+1);
            hbin1->SetLineColor(kRed+1);
            hbin1->SetMarkerStyle(20);
            hbin1->SetMarkerSize(1.2);
            hbin1->GetXaxis()->SetTitle("#Delta#it{z}^{ch}_{||}");
            hbin1->GetYaxis()->SetTitle("Probability Density");
            hbin1->GetXaxis()->SetLabelSize(0.04);
            hbin1->GetXaxis()->SetTitleSize(0.05);
            hbin1->GetXaxis()->SetTitleOffset(1.);
            hbin1->GetYaxis()->SetTitleOffset(1.1);
            hbin1->GetYaxis()->SetLabelSize(0.045);
            hbin1->GetYaxis()->SetTitleSize(0.05);
            hbin1->GetXaxis()->SetRangeUser(plotmin,plotmax);

            //hbin1->SetMaximum(hbin1->GetMaximum()*3);
            hbin1->SetMaximum(hbin1->GetMaximum()*10);
            hbin1->Scale(1./hbin1->GetEntries());
            hbin1->Scale(1,"width");


            hbin2->SetTitle("");
            hbin2->SetMarkerColor(kBlue+1);
            hbin2->SetLineColor(kBlue+1);
            hbin2->SetMarkerStyle(21);
            hbin2->SetMarkerSize(1.2);
            hbin2->GetXaxis()->SetTitle("#Delta#it{z}^{ch}_{||}");
            hbin2->GetYaxis()->SetTitle("Probability Density");
            hbin2->GetXaxis()->SetLabelSize(0.04);
            hbin2->GetXaxis()->SetTitleSize(0.05);
            hbin2->GetXaxis()->SetTitleOffset(1.);
            hbin2->GetYaxis()->SetTitleOffset(1.1);
            hbin2->GetYaxis()->SetLabelSize(0.045);
            hbin2->GetYaxis()->SetTitleSize(0.05);
            hbin2->GetXaxis()->SetRangeUser(plotmin,plotmax);
            hbin2->SetMaximum(hbin1->GetMaximum()*10);
            hbin2->Scale(1./hbin2->GetEntries());
            hbin2->Scale(1,"width");


            hbin3->SetTitle("");
            hbin3->SetMarkerColor(kGreen+1);
            hbin3->SetLineColor(kGreen+1);
            hbin3->SetMarkerStyle(22);
            hbin3->SetMarkerSize(1.2);
            hbin3->GetXaxis()->SetTitle("#Delta#it{z}^{ch}_{||}");
            hbin3->GetYaxis()->SetTitle("Probability Density");
            hbin3->GetXaxis()->SetLabelSize(0.04);
            hbin3->GetXaxis()->SetTitleSize(0.05);
            hbin3->GetXaxis()->SetTitleOffset(1.);
            hbin3->GetYaxis()->SetTitleOffset(1.1);
            hbin3->GetYaxis()->SetLabelSize(0.045);
            hbin3->GetYaxis()->SetTitleSize(0.05);
            hbin3->GetXaxis()->SetRangeUser(plotmin,plotmax);
            hbin3->SetMaximum(hbin1->GetMaximum()*10);
            hbin3->Scale(1./hbin3->GetEntries());
            hbin3->Scale(1,"width");;

        //    hbin1->GetYaxis()->SetRangeUser(1E-3,70);
        //    hbin2->GetYaxis()->SetRangeUser(1E-3,70);
        //    hbin3->GetYaxis()->SetRangeUser(1E-3,70);


   // TLegend *leg = new TLegend(0.15,0.55,0.5,0.70);
    TLegend *leg = new TLegend(0.58,0.2,0.92,0.35);
    //leg->SetTextSize(0.045);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->AddEntry(hbin1,"0.4 < #it{Z}^{ch,gen}_{||} < 0.6","p");
    leg->AddEntry(hbin2,"0.7 < #it{Z}^{ch,gen}_{||} < 0.8","p");
    leg->AddEntry(hbin3,"0.9 < #it{Z}^{ch,gen}_{||} < 1.0","p");



    TPaveText *pt = new TPaveText(0.15,0.68,0.5,0.95,"NB NDC");
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(13);
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
    //TText *text = pt->AddText("ALICE Preliminary");
    TText *text = new TText;
    text = pt->AddText("This Thesis"); //uncomment
    text = pt->AddText("PYTHIA6+GEANT3, pp, #sqrt{#it{s}} = 13 TeV");
    text = pt->AddText(Form("charged jets, anti-#it{k}_{T}, #it{R} = 0.%d, |#it{#eta}_{lab}^{jet}| < 0.%d",4,5));
    text = pt->AddText(Form ("with D^{0}, %d < #it{p}_{T,D^{0}} < %d GeV/#it{c}",3,10));
    text = pt->AddText(Form ("with D^{0}, %d < #it{p}^{gen}_{T,ch jet}  < %d GeV/#it{c}",7,10));

    hbin1->GetYaxis()->SetRangeUser(1E-4,1300);

    TCanvas *cEff = new TCanvas("cEff","cEff",1000,800);
    cEff->SetLogy();
    hbin1->Draw();
    hbin2->Draw("same");
    hbin3->Draw("same");
    pt->Draw("same");



    leg->Draw("same");


   // cEff->SaveAs("DjetpTres.png");
    cEff->SaveAs("DjetpTres_Z_thesis.pdf");
   // cEff->SaveAs("DjetpTres.eps");




}
"""
