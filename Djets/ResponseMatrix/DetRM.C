//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------

#include "config.h"

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

    float jetmin = 0, jetmax = 60;
    float Dptmin = fptbinsDA[0], Dptmax = fptbinsDA[fptbinsDN];

    TPaveText *pv2 = new TPaveText(0.15,0.8,0.25,0.9,"brNDC");
    pv2->SetFillStyle(0);
    pv2->SetBorderSize(0);
    pv2->AddText(Form("R=0.%d",Rpar));

    TH1D *hMCpt;
	TH1D *hMCpt_reco;
    TH2D *hPtJet[NDMC];
    THnSparseD *hPtJetD[NDMC];//to store Djet eff corrected spectra
    TH1D *hPtG[NDMC];//matched, not required
    TH1D *hPtR[NDMC];//matched, not required

	TList *histList[NDMC];TList *histList2[NDMC];
	THnSparseD *sparseMC[NDMC];
	THnSparseD *GMC[NDMC];
	THnSparseD *RMC[NDMC];

    TH2D *hPtJet2d_noeff; 
    THnSparseD *hPtJetD3d;
    TH1D *hPtJetGen;
    TH1D *hPtJetRec;

    TH1D* hGMC[NDMC];//closure
    TH1D* hRMC[NDMC];
    TH1D* hGMCall;
    TH1D* hRMCall;

    for(int i=0; i<NDMC; i++){
        if(!isprefix){
            if(!postfix) { 
                if(isPrompt){
                    histList[i] =  (TList*)dir->Get(Form("%s%dMCrec",histName.Data(),i));
                }
                else{
                    histList[i] =  (TList*)dir->Get(Form("%s%dFDMCrec",histName.Data(),i));
                }
            }
            else {
                histList[i] =  (TList*)dir->Get(Form("%s%d%sMCrec",histName.Data(),i,listName.Data())); 
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
        GMC[i] = (THnSparseD*)sparseMC[i]->Clone(Form("GMC_%d",i));
        RMC[i] = (THnSparseD*)sparseMC[i]->Clone(Form("RMC_%d",i));
        //Dpt and jetpt cuts on Reco level
        sparseMC[i]->GetAxis(2)->SetRangeUser(Dptmin,Dptmax);
        sparseMC[i]->GetAxis(1)->SetRangeUser(jetmin,jetmax);
        RMC[i]->GetAxis(2)->SetRangeUser(Dptmin,Dptmax);//for closure
        RMC[i]->GetAxis(1)->SetRangeUser(jetmin,jetmax);//for closure
        //RMC[i]->GetAxis(6)->SetRangeUser(jetmin,jetmax);
        //Dpt and jetpt cuts on Gen level
        if(fDmesonSpecie){// Dstar 
            sparseMC[i]->GetAxis(6)->SetRangeUser(Dptmin,Dptmax);
            sparseMC[i]->GetAxis(5)->SetRangeUser(jetmin,jetmax); 
        }
        else {// Dzero 
            sparseMC[i]->GetAxis(7)->SetRangeUser(Dptmin,Dptmax);
            sparseMC[i]->GetAxis(6)->SetRangeUser(jetmin,jetmax);
            GMC[i]->GetAxis(7)->SetRangeUser(Dptmin,Dptmax);
            GMC[i]->GetAxis(6)->SetRangeUser(jetmin,jetmax);
        }

        //eta cuts on Gen and Reco levels
    	if(!fDmesonSpecie){//Dzero
    		sparseMC[i]->GetAxis(4)->SetRangeUser(-(0.9-fRpar),0.9-fRpar);
    		sparseMC[i]->GetAxis(9)->SetRangeUser(-(0.9-fRpar),0.9-fRpar);
    	    RMC[i]->GetAxis(4)->SetRangeUser(-(0.9-fRpar),0.9-fRpar);
    	    GMC[i]->GetAxis(9)->SetRangeUser(-(0.9-fRpar),0.9-fRpar);
    	}

        const int DNresDim = 3;
        if(fDmesonSpecie){//Projection(DnDim,Ddim,"E");// z, j
            int DresDim[DNresDim] = {1,5,2};
            hPtJet[i] = (TH2D*)sparseMC[i]->Projection(5,1,"E"); //Dstar tmp
            hPtJetD[i] = (THnSparseD*)sparseMC[i]->Projection(DNresDim,DresDim,"E");//additional Dpt gen axis
        }
        else{
            int DresDim[DNresDim] = {1,6,2};
            hPtJet[i] = (TH2D*)sparseMC[i]->Projection(6,1,"E");
            hPtJetD[i] = (THnSparseD*)sparseMC[i]->Projection(DNresDim,DresDim,"E");//additional Dpt gen axis
            hGMC[i] = (TH1D*)GMC[i]->Projection(7,"E");//closure
            hRMC[i] = (TH1D*)RMC[i]->Projection(2,"E");//closure
        }

        hPtJet[i]->Sumw2(); 
        hPtJetD[i]->Sumw2();//not happening for THnSparse
        hPtJet[i]->SetName(Form("hPtJet_%d",i)); 
        hPtJetD[i]->SetName(Form("hPtJetD_%d",i));

        hGMC[i]->Sumw2();//closure
        hRMC[i]->Sumw2();//closure
        hGMC[i]->SetName(Form("hGMC_%d",i));
        hRMC[i]->SetName(Form("hRMC_%d",i));

        //gen level jetpt projection
        if(fDmesonSpecie) hPtG[i] = (TH1D*)sparseMC[i]->Projection(5); //Dstar tmp
        else hPtG[i] = (TH1D*)sparseMC[i]->Projection(6);

        //reco level jetpt projection
        hPtR[i] = (TH1D*)sparseMC[i]->Projection(1);

        hPtG[i]->SetName(Form("hPtG_%d",i));
        hPtR[i]->SetName(Form("hPtR_%d",i));
        hPtG[i]->Sumw2();
        hPtR[i]->Sumw2();

		if (!i){
  		    hPtJet2d_noeff = (TH2D*)hPtJet[0]->Clone("hPtJet2d_noeff");
  		    hPtJetD3d = (THnSparseD*)hPtJetD[0]->Clone("hPtJetD3d");
  		    hPtJetGen = (TH1D*)hPtG[0]->Clone("hPtJetGen");
  		    hPtJetRec = (TH1D*)hPtR[0]->Clone("hPtJetRec");
            hGMCall = (TH1D*)hGMC[0]->Clone("hGMCall");
            hRMCall = (TH1D*)hRMC[0]->Clone("hRMCall");
        }
        else{
            hPtJet2d_noeff->Add(hPtJet[i]);
            hPtJetD3d->Add(hPtJetD[i]);
      		hPtJetGen->Add(hPtG[i]);
      		hPtJetRec->Add(hPtR[i]);
            hGMCall->Add(hGMC[i]);
            hRMCall->Add(hRMC[i]);
        }
	}
    //--- Djet eff scaling and filling a new response
    TFile *fileEff = new TFile(effFilePnP.Data(),"read");
    TH1D* effHist; effHist = (TH1D*)fileEff->Get("hEff_reb");

    TCanvas* ceff = new TCanvas("ceff","ceff",600,400);
    TH1D* heffRec = (TH1D*)hRMCall->Clone("heffRec");
    TH1D* heffRecReb = (TH1D*)heffRec->Rebin(fptbinsDN,"heffRecReb",fptbinsDA);
    TH1D* heffGen = (TH1D*)hGMCall->Clone("heffGen");
    TH1D* heffGenReb = (TH1D*)heffGen->Rebin(fptbinsDN,"heffGenReb",fptbinsDA);
    heffRecReb->Divide(heffGenReb);
    heffRecReb->Draw();
    ceff->SaveAs(Form("%s/ceff_%d.png",outDir.Data(),isPrompt));

    TH2D *hPtJet2dDjet = new TH2D("hjetRecGen","hjetRecGen", 60, 0, 60, 60, 0, 60); 
    TH1D *hGenNum = new TH1D("hGenNum","hGenNum",fptbinsJetMeasN,fptbinsJetMeasA);
    TH1D *hGenDen = new TH1D("hGenDen","hGenDen",fptbinsJetMeasN,fptbinsJetMeasA);
    TH1D *hRecNum = new TH1D("hRecNum","hRecNum",fptbinsJetMeasN,fptbinsJetMeasA);
    TH1D *hRecDen = new TH1D("hRecDen","hRecDen",fptbinsJetMeasN,fptbinsJetMeasA);
    cout<<"Number of bins in the response: "<<hPtJetD3d->GetNbins()<<endl;
    for (int z = 0; z < hPtJetD3d->GetNbins(); z++){
        int coord[3]={0,0,0};
        double content = hPtJetD3d->GetBinContent(z,coord);
        int i = coord[0], j = coord[1], k = coord[2];
        double weight = content;
        double jR_center = hPtJetD3d->GetAxis(0)->GetBinCenter(i);
        double jG_center = hPtJetD3d->GetAxis(1)->GetBinCenter(j);
        double DR_center = hPtJetD3d->GetAxis(2)->GetBinCenter(k);
        //if(jR_center<0){ cout<<jR_center<<"=="<<jG_center<<endl;}

        double eff_c = 1;
        eff_c = effHist->GetBinContent(effHist->GetXaxis()->FindBin(DR_center));
        weight = weight/eff_c;
        hPtJet2dDjet->Fill(jR_center, jG_center, weight);

        hGenDen->Fill(jG_center, content); hRecDen->Fill(jR_center, weight);
        if(jR_center>=fptbinsJetMeasA[0] || jR_center <= fptbinsJetMeasA[fptbinsJetMeasN]){hGenNum->Fill(jG_center, content);}
        if(jG_center>=fptbinsJetMeasA[0] || jG_center <= fptbinsJetMeasA[fptbinsJetMeasN]){hRecNum->Fill(jR_center, weight);}
    }
    //kine eff, setting errors zero
    TH1D* hGenEff=(TH1D*)hGenNum->Clone("hGenEff");
    TH1D* hRecEff=(TH1D*)hRecNum->Clone("hRecEff");
    hGenEff->Divide(hGenDen);hRecEff->Divide(hRecDen);
    for(int i=0; i<hGenEff->GetXaxis()->GetNbins();i++){
        hGenEff->SetBinError(i+1,0);
        hRecEff->SetBinError(i+1,0);
    }
    TCanvas* cc_kineEff = new TCanvas("cc_kf","cc_kf", 600,500);
    hGenNum->GetYaxis()->SetRangeUser(0.9,1.1);
    hGenNum->Draw(); hRecNum->Draw("same");
    cc_kineEff->SaveAs(Form("%s/kineEff_%d.png",outDir.Data(),isPrompt));

if(isPrompt){
    //closure===========================from sparse
    const int DnDim=6;
    int zRec = 0, jetRec = 1, DRec = 2;
    int zGen = 5, jetGen = 6, DGen = 7;
    int Ddim[DnDim]   = {zRec,jetRec,zGen,jetGen,DRec,DGen};
    THnSparseD* GMC_0 = (THnSparseD*)GMC[0]->Projection(DnDim,Ddim,"E");//Clone("GMC_all");
    THnSparseD* GMC_1 = (THnSparseD*)GMC[1]->Projection(DnDim,Ddim,"E");//Clone("GMC_all");
    THnSparseD* RMC_0 = (THnSparseD*)RMC[0]->Projection(DnDim,Ddim,"E");//Clone("RMC_all");
    THnSparseD* RMC_1 = (THnSparseD*)RMC[1]->Projection(DnDim,Ddim,"E");//Clone("RMC_all");
    THnSparseD* GMC_all = (THnSparseD*)GMC_0->Clone("GMC_all");
    THnSparseD* RMC_all = (THnSparseD*)RMC_0->Clone("RMC_all");
    GMC_all->Add(GMC_1);
    RMC_all->Add(RMC_1);
    TH2D* h_gen;// = new TH2D("hGMC_gen","hGMC_gen",60,0,60,60,0,60);
    TH2D* h_rec = new TH2D("hGMC_rec","hGMC_rec",60,0,60,60,0,60);
    TH2D* h_eff = new TH2D("hGMC_eff","hGMC_eff",60,0,60,60,0,60);
    cout<<"====="<<GMC_all->GetNbins()<<endl;
    cout<<"====="<<RMC_all->GetNbins()<<endl;

    THnSparseD* sparse_gen=(THnSparseD*)GMC_all->Clone("sparse_gen");//Projection(DnDim,Ddim,"E");
    THnSparseD* sparse_rec=(THnSparseD*)RMC_all->Clone("sparse_rec");//Projection(DnDim,Ddim,"E");
    
    
    for (int i=1; i<sparse_rec->GetNbins()+1;i++){
        int coord[DnDim] = {0,0,0,0,0,0};
        double content = sparse_rec->GetBinContent(i,coord); 
        double DR_center = sparse_rec->GetAxis(2)->GetBinCenter(coord[2]);
        double jR_center = sparse_rec->GetAxis(1)->GetBinCenter(coord[1]);
        double eff_c = 1;
        eff_c = effHist->GetBinContent(effHist->GetXaxis()->FindBin(DR_center));
        double weight = content/eff_c;
        h_eff->Fill(DR_center,jR_center,weight);
        h_rec->Fill(DR_center,jR_center,content);
    }
    TCanvas* ccheck = new TCanvas("ccheck","ccheck",600,400);
//    //hGMC_eff1->SetLineColor(kRed+2);
    TH1D* hrec = (TH1D*)h_rec->ProjectionX("hrec",0,-1);
    TH1D* heff = (TH1D*)h_eff->ProjectionX("heff",0,-1);
    hrec->SetLineColor(kGreen+2);
//    ccheck->SetLogy();
    heff->Draw();
//    hrec->Draw("same");
    ccheck->SaveAs(Form("%s/ccheck_%d.png",outDir.Data(), isPrompt));
}
    // closure-------
    TH2D *hPtJetClosure = new TH2D("hjetRGClosure","hjetRGClosure", 60, 0, 60, 60, 0, 60); 
    int full_resp_size = hPtJetD3d->GetNbins();
    int clos_resp_size = (int)(full_resp_size*0.8);

    //random stuff
    random_device dev;
    mt19937 rng(dev());
    uniform_int_distribution<mt19937::result_type> dist6(1,full_resp_size);

    //storing the random number distribution: a check
    TH1D* h_test = new TH1D("test","test",clos_resp_size/20,0,clos_resp_size);
    for (int i = 0; i < clos_resp_size; i++){h_test->Fill(dist6(rng));}
    TCanvas* cc_test = new TCanvas("rr","rr", 600,500);h_test->Draw();cc_test->SaveAs(Form("%s/TEST_%d.png",outDir.Data(),isPrompt));
    
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


    TCanvas *cjetPt = new TCanvas("cjetPt","cjetPt",800,600);
    cjetPt->SetLogy();
    hPtJetGen->Draw();
    hPtJetRec->Draw("same");
    //cjetPt->SaveAs(Form("%s/pTdist_Dpt%d_%d.png",outDir, (int)Dptmin, (int)Dptmax));

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

    ofile->Close();

   return;
}
