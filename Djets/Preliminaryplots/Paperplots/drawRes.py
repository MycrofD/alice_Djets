import os, os.path, sys
import ROOT
import style_settings
import array
ROOT.TH1.AddDirectory(False)


plotmin = -0.7
plotmax =0.4;

# ptbinsJetN = 11;
# ptJetbins = [2,3,4,5,6,8,10,12,14,20,30,50]

ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetPadBottomMargin(0.1)
ROOT.gStyle.SetPadTopMargin(1)
ROOT.gStyle.SetPadRightMargin(0.13)

def drawRes():
    gStyle->SetOptStat(000);

    gStyle->SetLegendFont(42);



   TFile *File = new TFile(
           //"/home/kvapil/work/analysis/pp_run2/D0jet/data_200519/MC/AnalysisResults_Run2.root"
           //"/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outData/AnalysisResults_745_R04_pp_5cuts.root"
           "/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/AnalysisResults_1061_R04_ppMC_5cuts.root"
           ,"read");
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
       sparse->GetAxis(0)->SetRangeUser(-2,2);
       sparse->GetAxis(1)->SetRangeUser(5,6);
       sparse->GetAxis(2)->SetRangeUser(2,36);
       if(i==0) {
         hbin1 = dynamic_cast<TH1D*>(sparse->Projection(10));
       }
       else {
         hbin1->Add(dynamic_cast<TH1D*>(sparse->Projection(10)));
       }

       sparse->GetAxis(1)->SetRangeUser(10,12);
       if(i==0) {
         hbin2 = dynamic_cast<TH1D*>(sparse->Projection(10));
       }
       else {
         hbin2->Add(dynamic_cast<TH1D*>(sparse->Projection(10)));
       }
       sparse->GetAxis(1)->SetRangeUser(20,30);
       if(i==0) {
         hbin3 = dynamic_cast<TH1D*>(sparse->Projection(10));
       }
       else {
         hbin3->Add(dynamic_cast<TH1D*>(sparse->Projection(10)));
       }
   }



            hbin1->SetTitle("");
            hbin1->SetMarkerColor(kRed+2);
            hbin1->SetLineColor(kRed+2);
            hbin1->SetMarkerStyle(24);//20
            hbin1->SetMarkerSize(1.2);
            hbin1->GetXaxis()->SetTitle("#Delta_{#it{p}_{T}}");
            hbin1->GetYaxis()->SetTitle("Probability Density");
            hbin1->GetXaxis()->SetLabelSize(0.04);
            hbin1->GetXaxis()->SetTitleSize(0.04);
            hbin1->GetXaxis()->SetTitleOffset(1.);
            hbin1->GetYaxis()->SetTitleOffset(1.3);
            hbin1->GetYaxis()->SetLabelSize(0.04);
            hbin1->GetYaxis()->SetTitleSize(0.04);
            hbin1->GetXaxis()->SetRangeUser(plotmin,plotmax);

            hbin1->SetMaximum(hbin1->GetMaximum()*3);
            hbin1->Scale(1./hbin1->GetEntries());
            hbin1->Scale(1,"width");


            hbin2->SetTitle("");
            hbin2->SetMarkerColor(kBlue+2);
            hbin2->SetLineColor(kBlue+2);
            hbin2->SetMarkerStyle(25);//21
            hbin2->SetMarkerSize(1.2);
            hbin2->GetXaxis()->SetTitle("#Delta#it{p}_{T}");
            hbin2->GetYaxis()->SetTitle("Probability Density");
            //hbin2->GetXaxis()->SetLabelSize(0.04);
            //hbin2->GetXaxis()->SetTitleSize(0.05);
            //hbin2->GetXaxis()->SetTitleOffset(1.);
            //hbin2->GetYaxis()->SetTitleOffset(1.1);
            //hbin2->GetYaxis()->SetLabelSize(0.045);
            //hbin2->GetYaxis()->SetTitleSize(0.05);
            //hbin2->GetXaxis()->SetRangeUser(plotmin,plotmax);
            //hbin2->SetMaximum(hbin1->GetMaximum()*3);
            hbin2->Scale(1./hbin2->GetEntries());
            hbin2->Scale(1,"width");


            hbin3->SetTitle("");
            hbin3->SetMarkerColor(kGreen+2);
            hbin3->SetLineColor(kGreen+2);
            hbin3->SetMarkerStyle(27);//22
            hbin3->SetMarkerSize(1.4);
            hbin3->GetXaxis()->SetTitle("#Delta#it{p}_{T}");
            hbin3->GetYaxis()->SetTitle("Probability Density");
            //hbin3->GetXaxis()->SetLabelSize(0.04);
            //hbin3->GetXaxis()->SetTitleSize(0.05);
            //hbin3->GetXaxis()->SetTitleOffset(1.);
            //hbin3->GetYaxis()->SetTitleOffset(1.1);
            //hbin3->GetYaxis()->SetLabelSize(0.045);
            //hbin3->GetYaxis()->SetTitleSize(0.05);
            //hbin3->GetXaxis()->SetRangeUser(plotmin,plotmax);
            //hbin3->SetMaximum(hbin1->GetMaximum()*3);
            hbin3->Scale(1./hbin3->GetEntries());
            hbin3->Scale(1,"width");;

        //    hbin1->GetYaxis()->SetRangeUser(1E-3,70);
        //    hbin2->GetYaxis()->SetRangeUser(1E-3,70);
        //    hbin3->GetYaxis()->SetRangeUser(1E-3,70);


    TLegend *leg = new TLegend(0.55,0.65,0.82,0.85);
   // TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
    //leg->SetTextSize(0.045);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->AddEntry(hbin1,"5 < #it{p}^{gen}_{T,ch jet} < 6 GeV/#it{c}","p");
    leg->AddEntry(hbin2,"10 < #it{p}^{gen}_{T,ch jet} < 12 GeV/#it{c}","p");
    leg->AddEntry(hbin3,"20 < #it{p}^{gen}_{T,ch jet} < 30 GeV/#it{c}","p");



    TPaveText *pt = new TPaveText(0.15,0.7,0.5,0.9,"NB NDC");
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(13);
    pt->SetTextFont(42);
    pt->SetTextSize(0.032);
    //TText *text = pt->AddText("ALICE Preliminary");
    TText *text = new TText;
    text = pt->AddText("ALICE PYTHIA6, pp, #sqrt{#it{s}} = 5.02 TeV"); //uncomment
    text = pt->AddText("Prompt D^{0} #rightarrow K^{-}#pi^{+} and charge conj.");
    text = pt->AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4");
    text = pt->AddText("|#it{#eta}_{lab}^{jet}| < 0.5");
    text = pt->AddText("#it{p}_{T, D^{0}} > 2 GeV/#it{c}");
    //text = pt->AddText("ALICE Simulation"); //uncomment
    //text = pt->AddText("PYTHIA6, pp, #sqrt{#it{s}} = 13 TeV");
    //text = pt->AddText(Form("charged jets, anti-#it{k}_{T}, #it{R} = 0.%d, |#it{#eta}_{lab}^{jet}| < 0.%d",4,5));
    //text = pt->AddText(Form ("with D^{0}, %d < #it{p}_{T,D^{0}} < %d GeV/#it{c}",2,36));

    TCanvas *cEff = new TCanvas("cEff","cEff",1150,800);
    cEff->SetLogy();
    hbin1->Draw();
    hbin2->Draw("same");
    hbin3->Draw("same");
    pt->Draw("same");



    leg->Draw("same");


    cEff->SaveAs("DjetpTres.png");
    cEff->SaveAs("DjetpTres.pdf");
    cEff->SaveAs("DjetpTres.eps");




