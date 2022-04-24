//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------

#include "config.h"

TTree* toTree(THnSparse* h);
void DetRM(
        bool isPrompt = 1, 
        TString datafile = "../outMC/AnalysisResults_fast_D0MCPythia_SMQcorr2.root", 
        TString outDir = "plots",
        bool postfix = 0, 
        TString listName = "FD", 
        bool isprefix=0,
        TString effFilePnP=""
        )
//const int DimResN = 3;int jetR=1;int jetG=6;int DptG=7;
//int DimRes = {jetR,jetG,DptG}; 
{
    gStyle->SetOptStat(0000); //Mean and RMS shown
    gStyle->SetPadRightMargin(0.1);
    gSystem->Exec(Form("mkdir %s",outDir.Data()));
    gSystem->Exec(Form("mkdir %s/plots",outDir.Data()));

    TFile *File = new TFile(datafile,"read");
    TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
    TString histName;
    if(!isprefix){
            if(fDmesonSpecie) histName = "histosDStarMBN";
            else histName = "histosD0MBN";}
    else{
            if(fDmesonSpecie) histName = "histosDStarMBN";
            else histName = "histosD0";}
    //closure
    TString datafile_resp=Form("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/closure_%d/AnalysisResults_input_100.root",Rpar);
    //TString datafile_input=Form("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/closure_%d/secondhalf/AnalysisResults_test2.root",Rpar);
    TString datafile_input=Form("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/outMC/closure_%d/test/AnalysisResults_test.root",Rpar);
    TFile *File_resp = new TFile(datafile_resp,"read");
    TFile *File_input = new TFile(datafile_input,"read");
    TDirectoryFile* dir_resp=(TDirectoryFile*)File_resp->Get("DmesonsForJetCorrelations");
    TDirectoryFile* dir_input=(TDirectoryFile*)File_input->Get("DmesonsForJetCorrelations");
    TList* histList_resp[NDMC];
    TList* histList_input[NDMC];
    //------
    float jetmin = 0, jetmax = 60;
    float Dptmin = fptbinsDA[0], Dptmax = fptbinsDA[fptbinsDN];

    TPaveText *pv2 = new TPaveText(0.15,0.8,0.25,0.9,"brNDC");
    pv2->SetFillStyle(0);
    pv2->SetBorderSize(0);
    pv2->AddText(Form("R=0.%d",Rpar));

    TH1D *hMCpt;
	TH1D *hMCpt_reco;
    TH2D *hPtJet[NDMC];
    THnSparseD *hPtJetD[NDMC];//to store Djet eff corrected response
    THnSparseD *hPtJetD_resp[NDMC];//to store Djet eff corrected response for closure
    TH1D *hPtG_clos[NDMC];//for closure, true level input
    THnSparseD *hPtR_cl2d[NDMC];//for closure, reco level input
    TH1D *hPtG[NDMC];//matched, not required
    TH1D *hPtR[NDMC];//matched, not required

	TList *histList[NDMC];
	THnSparseD *sparseMC[NDMC];
	THnSparseD *RespMC[NDMC];
	THnSparseD *GMC[NDMC];
	THnSparseD *RMC[NDMC];

    TH2D *hPtJet2d_noeff; 
    THnSparseD *hPtJetD3d;
    THnSparseD *hPtJetD3d_resp;//for closure response
    TH1D *hPtJetGen;
    TH1D *hPtJetGen_clos;
    TH1D *hPtJetRec;
    TH1D *hPtJetRec_clos;
    THnSparseD *hPtJetRec_cl2d;

    TTree* sparsetree[NDMC];

    for(int i=0; i<NDMC; i++){
        if(!isprefix){
            if(postfix){histList[i] =  (TList*)dir->Get(Form("%s%d%sMCrec",histName.Data(),i,listName.Data())); }
            else {
                if(isPrompt) {
                    histList[i] =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
                    histList_resp[i] =  (TList*)dir_resp->Get(Form("%s%dMCrec",histName.Data(),i));
                    histList_input[i] =  (TList*)dir_input->Get(Form("%s%dMCrec",histName.Data(),i));
                }
                else {
                    histList[i] =  (TList*)dir->Get(Form("%s%dFDMCrec",histName.Data(),i));
                    histList_resp[i] =  (TList*)dir_resp->Get(Form("%s%dFDMCrec",histName.Data(),i));
                    histList_input[i] =  (TList*)dir_input->Get(Form("%s%dFDMCrec",histName.Data(),i));
                }
            }
        }
        else{
             if(postfix) {
                     if(isPrompt){ histList[i] =  (TList*)dir->Get(Form("%s%sMBN%dMCrec",histName.Data(),listName.Data(),i)); }
                     else{    histList[i] =  (TList*)dir->Get(Form("%s%sMBN%dFDMCrec",histName.Data(),listName.Data(),i)); }
             }
             else { cout<<"-----postfix has to be true if prefix is true!! check again----------------"<<endl; return;       }
        }


        sparseMC[i] = (THnSparseD*)histList[i]->FindObject("ResponseMatrix");
        RespMC[i] = (THnSparseD*)histList_resp[i]->FindObject("ResponseMatrix");

        THnSparseD* gMC = (THnSparseD*)histList_input[i]->FindObject("ResponseMatrix");
        //sparsetree[i] = toTree(gMC);
        GMC[i] = (THnSparseD*)gMC->Clone(Form("GMC_%d",i));
        THnSparseD* rMC = (THnSparseD*)histList_input[i]->FindObject("ResponseMatrix");
        RMC[i] = (THnSparseD*)rMC->Clone(Form("RMC_%d",i));

        //Dpt and jetpt cuts on Reco level
        sparseMC[i]->GetAxis(2)->SetRangeUser(Dptmin,Dptmax);
        sparseMC[i]->GetAxis(1)->SetRangeUser(jetmin,jetmax);
        RespMC[i]->GetAxis(2)->SetRangeUser(Dptmin,Dptmax);//for closure, response
        RespMC[i]->GetAxis(1)->SetRangeUser(jetmin,jetmax);//for closure, response
        RMC[i]->GetAxis(2)->SetRangeUser(Dptmin,Dptmax);//for closure, data input
        RMC[i]->GetAxis(1)->SetRangeUser(jetmin,jetmax);//for closure, data input

        //Dpt and jetpt cuts on Gen level
        if(fDmesonSpecie){// Dstar 
            sparseMC[i]->GetAxis(6)->SetRangeUser(Dptmin,Dptmax);
            sparseMC[i]->GetAxis(5)->SetRangeUser(jetmin,jetmax); 
        }
        else {// Dzero 
            sparseMC[i]->GetAxis(7)->SetRangeUser(Dptmin,Dptmax);
            sparseMC[i]->GetAxis(6)->SetRangeUser(jetmin,jetmax);
            RespMC[i]->GetAxis(7)->SetRangeUser(Dptmin,Dptmax);//closure, response
            RespMC[i]->GetAxis(6)->SetRangeUser(jetmin,jetmax);//closure, response
            GMC[i]->GetAxis(7)->SetRangeUser(Dptmin,Dptmax);//closure, true input
            GMC[i]->GetAxis(6)->SetRangeUser(jetmin,jetmax);//closure, true input
        }

        //eta cuts on Gen and Reco levels
    	if(!fDmesonSpecie) {//Dzero
    		sparseMC[i]->GetAxis(4)->SetRangeUser(-(0.9-fRpar),0.9-fRpar);
    		sparseMC[i]->GetAxis(9)->SetRangeUser(-(0.9-fRpar),0.9-fRpar);
    		RespMC[i]->GetAxis(4)->SetRangeUser(-(0.9-fRpar),0.9-fRpar);//closure, response
    		RespMC[i]->GetAxis(9)->SetRangeUser(-(0.9-fRpar),0.9-fRpar);//closure, response
    		RMC[i]->GetAxis(4)->SetRangeUser(-(0.9-fRpar),0.9-fRpar);//closure, data input
    		GMC[i]->GetAxis(9)->SetRangeUser(-(0.9-fRpar),0.9-fRpar);//closure, true input
    	}

        const int DNresDim = 3;
        if(fDmesonSpecie){//Projection(DnDim,Ddim,"E");// z, j
            hPtG[i] = (TH1D*)sparseMC[i]->Projection(5); //Dstar tmp
            hPtR[i] = (TH1D*)sparseMC[i]->Projection(1);
            int DresDim[DNresDim] = {1,5,2};
            hPtJet[i] = (TH2D*)sparseMC[i]->Projection(5,1,"E"); //Dstar tmp
            hPtJetD[i] = (THnSparseD*)sparseMC[i]->Projection(DNresDim,DresDim,"E");//additional Dpt gen axis
            hPtJetD_resp[i] = (THnSparseD*)RespMC[i]->Projection(DNresDim,DresDim,"E");//additional Dpt gen axis
        }
        else{
            hPtG[i] = (TH1D*)sparseMC[i]->Projection(6);
            hPtR[i] = (TH1D*)sparseMC[i]->Projection(1);
            int DresDim[DNresDim] = {1,6,2};
            hPtJet[i] = (TH2D*)sparseMC[i]->Projection(6,1,"E");
            hPtJetD[i] = (THnSparseD*)sparseMC[i]->Projection(DNresDim,DresDim,"E");//additional Dpt gen axis
            hPtJetD_resp[i] = (THnSparseD*)RespMC[i]->Projection(DNresDim,DresDim,"E");//additional Dpt gen axis
            hPtG_clos[i] = (TH1D*)GMC[i]->Projection(6,"E");//
            hPtR_cl2d[i] = (THnSparseD*)RMC[i]->Projection(DNresDim,DresDim,"E");//
        }



        hPtJet[i]->Sumw2(); 
        hPtJetD[i]->Sumw2();//
        hPtJetD_resp[i]->Sumw2();//closure response
        hPtG_clos[i]->Sumw2();//closure true
        hPtR_cl2d[i]->Sumw2();//closure data input
        hPtJet[i]->SetName(Form("hPtJet_%d",i)); 
        hPtJetD[i]->SetName(Form("hPtJetD_%d",i));
        hPtG_clos[i]->SetName(Form("hPtG_clos_%d",i));
        hPtR_cl2d[i]->SetName(Form("hPtR_cl2d_%d",i));


        hPtG[i]->SetName(Form("hPtG_%d",i));
        hPtR[i]->SetName(Form("hPtR_%d",i));
        hPtG[i]->Sumw2();
        hPtR[i]->Sumw2();

		if (!i){
  		    hPtJet2d_noeff = (TH2D*)hPtJet[0]->Clone("hPtJet2d_noeff");
  		    hPtJetD3d = (THnSparseD*)hPtJetD[0]->Clone("hPtJetD3d");
  		    hPtJetD3d_resp = (THnSparseD*)hPtJetD_resp[0]->Clone("hPtJetD3d_resp");
  		    hPtJetGen = (TH1D*)hPtG[0]->Clone("hPtJetGen");
  		    hPtJetRec = (TH1D*)hPtR[0]->Clone("hPtJetRec");
            hPtJetRec_cl2d=(THnSparseD*)hPtR_cl2d[0]->Clone("hPtR_cl2d");
            hPtJetGen_clos=(TH1D*)hPtG_clos[0]->Clone("hPtG_clos");
        }
        else{
            hPtJet2d_noeff->Add(hPtJet[i]);
            hPtJetD3d->Add(hPtJetD[i]);
            hPtJetD3d_resp->Add(hPtJetD_resp[i]);
      		hPtJetGen->Add(hPtG[i]);
      		hPtJetRec->Add(hPtR[i]);
            hPtJetRec_cl2d->Add(hPtR_cl2d[i]);
            hPtJetGen_clos->Add(hPtG_clos[i]);
        }
	}

    //--- Djet eff scaling and filling a new response
    TFile *fileEff = new TFile(effFilePnP.Data(),"read");
    TH1D* effHist; effHist = (TH1D*)fileEff->Get("hEff_reb");
    TH2D *hPtJet2dDjet = new TH2D("hjetRecGen","hjetRecGen", 60, 0, 60, 60, 0, 60); 
    TH2D *hPtJet2dDjet_clos = new TH2D("hjetRecGen_resp","hjetRecGen_resp", 60, 0, 60, 60, 0, 60); //for closure response
    TH1D *hGenNum = new TH1D("hGenNum","hGenNum",fptbinsJetMeasN,fptbinsJetMeasA);
    TH1D *hGenDen = new TH1D("hGenDen","hGenDen",fptbinsJetMeasN,fptbinsJetMeasA);//gen level dist, for kine eff
    TH1D *hRecNum = new TH1D("hRecNum","hRecNum",fptbinsJetMeasN,fptbinsJetMeasA);
    TH1D *hRecDen = new TH1D("hRecDen","hRecDen",fptbinsJetMeasN,fptbinsJetMeasA);
    cout<<"Number of bins in the response: "<<hPtJetD3d->GetNbins()<<endl;
    for (int z = 0; z < hPtJetD3d->GetNbins()+1; z++){
        int coord[3]={0,0,0};
        double content = hPtJetD3d->GetBinContent(z,coord);
        int i = coord[0], j = coord[1], k = coord[2];
        double weight = content;
        double jR_center = hPtJetD3d->GetAxis(0)->GetBinCenter(i);
        double jG_center = hPtJetD3d->GetAxis(1)->GetBinCenter(j);
        double DR_center = hPtJetD3d->GetAxis(2)->GetBinCenter(k);
        double eff_c = 1;
        eff_c = effHist->GetBinContent(effHist->GetXaxis()->FindBin(DR_center));
        weight = weight/eff_c;
        hPtJet2dDjet->Fill(jR_center, jG_center, weight);

        //kine eff stuff
        hGenDen->Fill(jG_center, content); hRecDen->Fill(jR_center, weight);
        if(jG_center>=fptbinsJetMeasA[0] || jG_center <= fptbinsJetMeasA[fptbinsJetMeasN]){hGenNum->Fill(jG_center, content);}
        if(jR_center>=fptbinsJetMeasA[0] || jR_center <= fptbinsJetMeasA[fptbinsJetMeasN]){hRecNum->Fill(jR_center, weight);}
    }
    //kine eff, setting errors zero
    TH1D* hGenEff = (TH1D*)hGenNum->Clone("hGenEff");
    TH1D* hRecEff = (TH1D*)hRecNum->Clone("hRecEff");
    hGenEff->Divide(hGenDen);
    hRecEff->Divide(hRecDen);
    for(int i=0; i<hGenEff->GetXaxis()->GetNbins();i++){
        hGenEff->SetBinError(i+1,0);
        hRecEff->SetBinError(i+1,0);
    }
    TCanvas* cc_kineEff = new TCanvas("cc_kf","cc_kf", 600,500);
    //hGenEff->GetYaxis()->SetRangeUser(0.99,1.01);
    hRecEff->Draw(); //hRecEff->Draw("same");
    cc_kineEff->SaveAs(Form("%s/kineEff_%d.png",outDir.Data(),isPrompt));
    // closure------- response
    for (int z = 0; z < hPtJetD3d_resp->GetNbins()+1; z++){
        int coord[3]={0,0,0};
        double content = hPtJetD3d_resp->GetBinContent(z,coord);
        int i = coord[0], j = coord[1], k = coord[2];
        double weight = content;
        double jR_center = hPtJetD3d_resp->GetAxis(0)->GetBinCenter(i);
        double jG_center = hPtJetD3d_resp->GetAxis(1)->GetBinCenter(j);
        double DR_center = hPtJetD3d_resp->GetAxis(2)->GetBinCenter(k);
        double eff_c = 1;
        eff_c = effHist->GetBinContent(effHist->GetXaxis()->FindBin(DR_center));
        weight = weight/eff_c;
        hPtJet2dDjet_clos->Fill(jR_center, jG_center, weight);
    }
    // closure--------- reco input D-jet eff input
    hPtJetRec_clos = new TH1D("hPtJetRec_clos","hPtJetRec_clos", 60, 0, 60); 
    for (int z = 0; z < hPtJetRec_cl2d->GetNbins()+1; z++){
        int coord[3]={0,0,0};
        double content = hPtJetRec_cl2d->GetBinContent(z,coord);
        double DR_center = hPtJetRec_cl2d->GetAxis(2)->GetBinCenter(coord[2]);
        double jR_center = hPtJetRec_cl2d->GetAxis(0)->GetBinCenter(coord[0]);
        double eff_c = 1;
        eff_c = effHist->GetBinContent(effHist->GetXaxis()->FindBin(DR_center));
        double weight = content/eff_c;
        hPtJetRec_clos->Fill(jR_center, weight);
    }
    // --------------
    hPtJet2d_noeff->SetTitle(""); 
    hPtJet2dDjet->SetTitle("");
    
    hPtJet2d_noeff->SetName("hPtJet2d_noeff"); 
    hPtJet2dDjet->SetName("hPtJet2d");
    hPtJet2d_noeff->GetXaxis()->SetTitle("p_{T,ch jet}^{rec.} (GeV/#it{c})");
    hPtJet2d_noeff->GetYaxis()->SetTitle("p_{T,ch jet}^{gen.} (GeV/#it{c})");
    hPtJet2dDjet->GetXaxis()->SetTitle("p_{T,ch jet}^{rec.} (GeV/#it{c})");
    hPtJet2dDjet->GetYaxis()->SetTitle("p_{T,ch jet}^{gen.} (GeV/#it{c})");

    hPtJetGen->SetName("hPtJetGen");
    hPtJetRec->SetName("hPtJetRec");
    hPtJetGen->SetLineColor(kBlue+2);
    hPtJetRec->SetLineColor(kRed+2);

    hPtJetGen_clos->SetName("hPtJetGen_clos");
    hPtJetRec_clos->SetName("hPtJetRec_clos");
    hPtJet2dDjet_clos->SetName("hPtJet2d_clos");
//    hPtJetGen->SetLineColor(kBlue+2);
//    hPtJetRec->SetLineColor(kRed+2);

    TCanvas *cjetPt = new TCanvas("cjetPt","cjetPt",800,600);
    hPtJetRec_clos->SetLineColor(kRed+2);
    hPtJetGen_clos->SetLineWidth(2);
    cjetPt->SetLogy();
    hPtJetRec_clos->Draw();
    hPtJetGen_clos->Draw("same");
    cjetPt->SaveAs(Form("%s/true_effRec_%d.png",outDir.Data(), isPrompt));

    TCanvas *cjetPt2d = new TCanvas("cjetPt2d","cjetPt2d",800,600);
    cjetPt2d->SetLogz();
    hPtJet2d_noeff->Draw("colz");
    pv2->Draw("same");

    TCanvas *cjetPt2dDjet = new TCanvas("cjetPt2dDjet","cjetPt2dDjet",800,600);
    cjetPt2dDjet->SetLogz();
    hPtJet2dDjet->Draw("colz");
    pv2->Draw("same");

    cjetPt2d->SaveAs(Form("%s/plots/DetMatrix_%s_Dpt%d_%d.png",outDir.Data(), isPrompt ? "prompt" : "nonPrompt", (int)Dptmin, (int)Dptmax));
    cjetPt2dDjet->SaveAs(Form("%s/plots/DetMatrix_%s_Dpt%d_%d_Djet.png",outDir.Data(), isPrompt ? "prompt" : "nonPrompt", (int)Dptmin, (int)Dptmax));

    //TFile *ofile = new TFile(Form("%s/DetMatrix_Dpt%d_%d.root",outDir, (int)Dptmin, (int)Dptmax),"RECREATE");
    //if(Djet)...//introduce flag for efficiency scaling of the response matrix as an option
    TFile *ofile = new TFile(Form("%s/DetMatrix_%s.root",outDir.Data(), isPrompt ? "prompt" : "nonPrompt" ),"RECREATE");

    hPtJetGen->Write();
    hPtJetRec->Write();
    hPtJet2d_noeff->Write();
    hPtJet2dDjet->Write();//D-jet eff scaled response
    hPtJetGen_clos->Write();//3. 25% MC no eff, true
    hPtJetRec_clos->Write();//2. 25% MC Djet eff scaled
    hPtJet2dDjet_clos->Write();//1. 75% response for closure, Djet eff scaled
    ofile->Close();

   return;
}

//change THnSparse to TTree
TTree* toTree(THnSparse* h)
 {
    // Creates a TTree and fills it with the coordinates of all
    // filled bins. The tree will have one branch for each dimension,
    // and one for the bin content.

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

    Int_t *bins = new Int_t[dim];
    for (Long64_t i = 0; i < h->GetNbins(); ++i) {
       x[dim] = h->GetBinContent(i, bins);
       for (Int_t d = 0; d < dim; ++d) {
          x[d] = h->GetAxis(d)->GetBinCenter(bins[d]);
       }

       tree->Fill();
    }

    delete [] bins;
    //delete [] x;
    return tree;
 }
