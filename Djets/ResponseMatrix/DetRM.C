//-----------------------------------------------------------------------
//  Author B.Trzeciak
//  Utrecht University
//  barbara.antonina.trzeciak@cern.ch
//-----------------------------------------------------------------------

#include "config.h"

void SparseToTree(TString MCfile, Bool_t isPostfix, TString postfix);

void DetRM(
        bool isPrompt = 1
        ,TString datafile = "../outMC/AnalysisResults_fast_D0MCPythia_SMQcorr2.root"
        ,TString outDir = "plots"
        ,bool postfix = 0
        ,TString listName = "FD"
        ,bool isprefix=0
        ,TString effFilePnP=""
        ,TString FDsubfile="FDsub.root"
        ,bool bClosure = 0
        )
{
//gStyle settings
    gStyle->SetOptStat(0000); //Mean and RMS shown
    gStyle->SetPadRightMargin(0.1);
    gSystem->Exec(Form("mkdir %s",outDir.Data()));
    gSystem->Exec(Form("mkdir %s/plots",outDir.Data()));
//reading MC train output file
    TFile *File = new TFile(datafile,"read");
    TDirectoryFile* dir=(TDirectoryFile*)File->Get("DmesonsForJetCorrelations");
    TString histName;
    if(!isprefix){
            if(fDmesonSpecie) histName = "histosDStarMBN";
            else histName = "histosD0MBN";}
    else{
            if(fDmesonSpecie) histName = "histosDStarMBN";
            else histName = "histosD0";}
//minor settings
    float jetmin = 0, jetmax = 60;
    float Dptmin = fptbinsDA[0], Dptmax = fptbinsDA[fptbinsDN];
//minor settings
    TPaveText *pv2 = new TPaveText(0.15,0.8,0.25,0.9,"brNDC");
    pv2->SetFillStyle(0);
    pv2->SetBorderSize(0);
    pv2->AddText(Form("R=0.%d",Rpar));
//Histograms initiation
    TH1D *hMCpt;
	TH1D *hMCpt_reco;
    TH2D *hPtJet[NDMC];
    THnSparseD *hPtJetD[NDMC];//to store Djet eff corrected spectra
    TH1D *hPtG[NDMC];//matched, not required; gives the hPtJetGen in root file
    TH1D *hPtR[NDMC];//matched, not required; gives the hPtJetRec in root file

	TList *histList[NDMC];
	THnSparseD *sparseMC[NDMC];
	THnSparseD *GMC[NDMC];
	THnSparseD *RMC[NDMC];
    TRandom3 *random1 = new TRandom3(0);

    TH2D *hPtJet2d_noeff; 
    THnSparseD *hPtJetD3d;
    TH1D *hPtJetGen;
    TH1D *hPtJetRec;

    TH1D* hGMC[NDMC];//efficiency Dpt dependent
    TH1D* hRMC[NDMC];
    TH1D* hGMCall;
    TH1D* hRMCall;

    for(int i=0; i<NDMC; i++){
        if(!isprefix){
            if(!postfix) { //default
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
             if(postfix) {//variation
                     if(isPrompt){ histList[i] =  (TList*)dir->Get(Form("%s%sMBN%dMCrec",histName.Data(),listName.Data(),i)); }
                     else{    histList[i] =  (TList*)dir->Get(Form("%s%sMBN%dFDMCrec",histName.Data(),listName.Data(),i)); }
             }
             else { cout<<"-----postfix has to be true if prefix is true!! check again----------------"<<endl; return;       }
        }

        sparseMC[i] = (THnSparseD*)histList[i]->FindObject("ResponseMatrix");
        GMC[i] = (THnSparseD*)histList[i]->FindObject("ResponseMatrix");//GMC[i] = (THnSparseD*)histList[i]->Clone(Form("GMC_%d",i));
        RMC[i] = (THnSparseD*)histList[i]->FindObject("ResponseMatrix");//RMC[i] = (THnSparseD*)histList[i]->Clone(Form("RMC_%d",i));
        //Dpt and jetpt cuts on Reco level
        sparseMC[i]->GetAxis(2)->SetRangeUser(Dptmin,Dptmax);
        sparseMC[i]->GetAxis(1)->SetRangeUser(jetmin,jetmax);
        RMC[i]->GetAxis(2)->SetRangeUser(Dptmin,Dptmax);//for closure
        RMC[i]->GetAxis(1)->SetRangeUser(jetmin,jetmax);//for closure
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
            hGMC[i] = (TH1D*)GMC[i]->Projection(7,"E");//Dpt gen proj
            hRMC[i] = (TH1D*)RMC[i]->Projection(2,"E");//Dpt rec proj
        }

        hPtJet[i]->Sumw2(); 
        hPtJetD[i]->Sumw2();//not happening for THnSparse
        hPtJet[i]->SetName(Form("hPtJet_%d",i)); 
        hPtJetD[i]->SetName(Form("hPtJetD_%d",i));

        hGMC[i]->Sumw2();//
        hRMC[i]->Sumw2();//
        hGMC[i]->SetName(Form("hGMC_%d",i));
        hRMC[i]->SetName(Form("hRMC_%d",i));

        //gen level jetpt projection, matched
        if(fDmesonSpecie) hPtG[i] = (TH1D*)sparseMC[i]->Projection(5); //Dstar tmp
        else hPtG[i] = (TH1D*)sparseMC[i]->Projection(6);
        hPtG[i]->SetName(Form("hPtG_%d",i));
        hPtG[i]->Sumw2();
        //reco level jetpt projection, matched
        hPtR[i] = (TH1D*)sparseMC[i]->Projection(1);
        hPtR[i]->SetName(Form("hPtR_%d",i));
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

    // Dpt dependent eff checks.
    TCanvas* ceff = new TCanvas("ceff","ceff",600,400);
    TH1D* heffRec = (TH1D*)hRMCall->Clone("heffRec");
    TH1D* heffRecReb = (TH1D*)heffRec->Rebin(fptbinsDN,"heffRecReb",fptbinsDA);
    TH1D* heffGen = (TH1D*)hGMCall->Clone("heffGen");
    TH1D* heffGenReb = (TH1D*)heffGen->Rebin(fptbinsDN,"heffGenReb",fptbinsDA);
    heffRecReb->Divide(heffGenReb);
    heffRecReb->Draw();
    ceff->SaveAs(Form("%s/ceff_%d.png",outDir.Data(),isPrompt));

    //Creating the detector responses
    TH2D *hPtJet2dDjet = new TH2D("hjetRecGen","hjetRecGen", 60, 0, 60, 60, 0, 60); 
    TH1D *hGenNum = new TH1D("hGenNum","hGenNum",fptbinsJetMeasN,fptbinsJetMeasA);
    TH1D *hGenDen = new TH1D("hGenDen","hGenDen",fptbinsJetMeasN,fptbinsJetMeasA);
    TH1D *hRecNum = new TH1D("hRecNum","hRecNum",fptbinsJetMeasN,fptbinsJetMeasA);
    TH1D *hRecDen = new TH1D("hRecDen","hRecDen",fptbinsJetMeasN,fptbinsJetMeasA);
    //resolution
    TH2D *hPtJet2d_resol = new TH2D("hjetRGResol","hjetRGResol", 50, -1.0, 1.0, 60, 0, 60); 
    TH2D *hPtJet2dDjet_resol = new TH2D("hjetRecGenResol","hjetRecGenResol", 50, -1.0, 1.0, 60, 0, 60); 
    TH1D *hR_resol = new TH1D("hR_resol","hR_resol",50,-1.0,1.0);
    TH1D *hG_resol = new TH1D("hG_resol","hG_resol",fptbinsJetMeasN,fptbinsJetMeasA);
    RooUnfoldResponse resol (hR_resol, hG_resol);
    
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
//        hPtJet2dDjet->Fill(jR_center, jG_center, weight);
//        hPtJet2dDjet_resol->Fill((jR_center-jG_center)/jG_center,jG_center,weight);
//        resol.Fill((jR_center-jG_center)/jG_center,jG_center,weight);
//        hPtJet2d_resol->Fill((jR_center-jG_center)/jG_center,jG_center,content);

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


    //--Preparing ttree method for full response and closure test info
    TH1D* hMCjetClosData = new TH1D("hMCjetClosData","hMCjetClosData",60,jetmin,jetmax);
    TH1D* hMCjetClosTrue = new TH1D("hMCjetClosTrue","hMCjetClosTrue",60,jetmin,jetmax);
    TH2D* hPtJet2dDjet_clos = new TH2D("hPtJet2dDjet_clos","hPtJet2dDjet_clos", 60, 0, 60, 60, 0, 60);
    //==
    TH1D* hMCjetD = (TH1D*)hPtJet2dDjet->ProjectionX("hMCjet");//reco level proj of eff scaled response: eff scaled MC jets
    double MCjets = hMCjetD->Integral();cout<<"eff MC jets="<<MCjets<<endl;
    hMCjetD->GetXaxis()->SetRangeUser(fptbinsJetMeasA[0],fptbinsJetMeasA[fptbinsJetMeasN]);
    MCjets = hMCjetD->Integral();cout<<"eff MC jets="<<MCjets<<endl;//eff-scaled MC jets in range
    //--
    //get TTree. Create if it doesn't exist
    TFile* fTreeSparse = nullptr;
    Double_t jmatch[13];//Dzero
    Double_t bincontent;
    TTree* tree_ = NULL;
    TString tmp(datafile);
    TString treepost = "TTree";if(postfix) treepost+=listName;
    treepost+=".root";
    tmp.Remove(datafile.Length()-5,5).Append(treepost);
    fTreeSparse = new TFile(tmp,"read");
    if(!(fTreeSparse->IsOpen())){
        //writing TTree
        cout<<"TTree file not found, converting THnSparse to TTree right now!"<<endl;
        SparseToTree(datafile, postfix,listName);
        //reading TTree
        cout<<"reading TTree: "<<tmp<<endl;
        fTreeSparse = new TFile(tmp,"read");
    }
    cout<<"reading TTree done"<<tmp<<endl;
    tree_ = (TTree*)fTreeSparse->Get("ResponseMatrixSum_tree");
    tree_->SetBranchAddress("coord",&jmatch);
    tree_->SetBranchAddress("bincontent",&bincontent);

    //datajets
    double datajets = 0;
//    TFile *fileFD = new TFile(FDsubfile,"read");
//    TH1D *bFDspectrum = (TH1D*)fileFD->Get("hData_binned_sub");
//    datajets = bFDspectrum->Integral();
//    cout<<"datajets="<<datajets<<endl;
    //scaling the prompt jets and response for closure
//    double RMscaling = 1-datajets/MCjets; //this is about RM 0.9, MCFD_reco 0.1
    double RMscaling = 0.8;
    cout<<RMscaling<<" :RMscaling"<<endl;
    //loop over sparse bins
    for(int i=0; i<tree_->GetEntries();i++){
        tree_->GetEntry(i);
        // loop over each jet-entry: why? to get a random value about the RMscaling
        double RMw = 0;
        for(int ientry = 0; ientry < bincontent; ientry++){
            if(random1->Uniform(1) <= RMscaling) RMw++;
        }
        //efficiency
        double eff_c = 1;
        eff_c = effHist->GetBinContent(effHist->GetXaxis()->FindBin(jmatch[2]));
        if(jetmin <= jmatch[6] && jmatch[6] <= jetmax){
            if(Dptmin <= jmatch[7] && jmatch[7] <= Dptmax){
                if(fRpar-0.9 <= jmatch[9] && jmatch[9] <= 0.9-fRpar){
                    hMCjetClosTrue->Fill(jmatch[6],bincontent-RMw);//gen level closure hist
                    if(jetmin <= jmatch[1] && jmatch[1] <= jetmax){
                        if(Dptmin <= jmatch[2] && jmatch[2] <= Dptmax){
                            if(fRpar-0.9 <= jmatch[4] && jmatch[4] <= 0.9-fRpar){
                                hMCjetClosData->Fill(jmatch[1],(bincontent-RMw)/eff_c);
                                hPtJet2dDjet_clos->Fill(jmatch[1],jmatch[6],RMw/eff_c);
                                hPtJet2dDjet->Fill(jmatch[1], jmatch[6],bincontent/eff_c);
                                hPtJet2dDjet_resol->Fill((jmatch[1]-jmatch[6])/jmatch[6],jmatch[6],bincontent/eff_c);
                                hPtJet2d_resol->Fill((jmatch[1]-jmatch[6])/jmatch[6],jmatch[6],bincontent);
                                resol.Fill((jmatch[1]-jmatch[6])/jmatch[6],jmatch[6],bincontent/eff_c);
                                //KINE EFF to be implemented later
//                                hGenDen->Fill(jmatch[6], bincontent); hRecDen->Fill(jmatch[1], bincontent/eff_c); //kine eff
//        if(jR_center>=fptbinsJetMeasA[0] || jR_center <= fptbinsJetMeasA[fptbinsJetMeasN]){hGenNum->Fill(jG_center, content);}
//        if(jG_center>=fptbinsJetMeasA[0] || jG_center <= fptbinsJetMeasA[fptbinsJetMeasN]){hRecNum->Fill(jR_center, weight);}
                            }
                        }
                    }
                }
            }
        }
    }
    TCanvas *cMCjetTD = new TCanvas("cMCjetTD","cMCjetTD",800,600);
    cMCjetTD->SetLogy();
    hMCjetClosTrue->Draw();
    hMCjetClosData->Draw("same");
    cMCjetTD->SaveAs(Form("%s/hMCjetTD_%s.png",outDir.Data(), isPrompt ? "prompt" : "nonPrompt"));
    
    // Main response stuff
    hPtJet2d_noeff->SetTitle(""); 
    hPtJet2dDjet->SetTitle("");
    
    hPtJet2d_noeff->SetName("hPtJet2d_noeff"); 
    hPtJet2dDjet->SetName("hPtJet2d");
    TH2D* resolobj = (TH2D*)resol.Hresponse()->Clone();
    resolobj->SetName("resolNow");
    hPtJet2d_noeff->GetXaxis()->SetTitle("p_{T,ch jet}^{rec.} (GeV/#it{c})");
    hPtJet2d_noeff->GetYaxis()->SetTitle("p_{T,ch jet}^{gen.} (GeV/#it{c})");
    hPtJet2dDjet->GetXaxis()->SetTitle("p_{T,ch jet}^{rec.} (GeV/#it{c})");
    hPtJet2dDjet->GetYaxis()->SetTitle("p_{T,ch jet}^{gen.} (GeV/#it{c})");

    //drawing the responses:
    //1. normal response
    TCanvas *cjetPt2d = new TCanvas("cjetPt2d","cjetPt2d",800,600);
    cjetPt2d->SetLogz();
    hPtJet2d_noeff->Draw("colz");
    pv2->Draw("same");
    //2. eff scaled response
    TCanvas *cjetPt2dDjet = new TCanvas("cjetPt2dDjet","cjetPt2dDjet",800,600);
    cjetPt2dDjet->SetLogz();
    hPtJet2dDjet->Draw("colz");
    pv2->Draw("same");
    //saving the two responses
    cjetPt2d->SaveAs(Form("%s/plots/DetMatrix_%s_Dpt%d_%d.png",outDir.Data(), isPrompt ? "prompt" : "nonPrompt", (int)Dptmin, (int)Dptmax));
    cjetPt2dDjet->SaveAs(Form("%s/plots/DetMatrix_%s_Dpt%d_%d_Djet.png",outDir.Data(), isPrompt ? "prompt" : "nonPrompt", (int)Dptmin, (int)Dptmax));
    //matched gen are reco level jetpt distributions
    hPtJetGen->SetName("hPtJetGen");
    hPtJetRec->SetName("hPtJetRec");
    hPtJetGen->SetLineColor(kBlue+2);
    hPtJetRec->SetLineColor(kRed+2);
    TCanvas *cjetPt = new TCanvas("cjetPt","cjetPt",800,600);
    cjetPt->SetLogy();
    hPtJetGen->Draw();
    hPtJetRec->Draw("same");
    //cjetPt->SaveAs(Form("%s/pTdist_Dpt%d_%d.png",outDir, (int)Dptmin, (int)Dptmax));
    //---------------------------------
    //TFile *ofile = new TFile(Form("%s/DetMatrix_Dpt%d_%d.root",outDir, (int)Dptmin, (int)Dptmax),"RECREATE");
    //if(Djet)...//introduce flag for efficiency scaling of the response matrix as an option
    TFile *ofile = new TFile(Form("%s/DetMatrix_%s.root",outDir.Data(), isPrompt ? "prompt" : "nonPrompt" ),"RECREATE");

    hPtJetGen->Write();
    hPtJetRec->Write();
    hPtJet2d_noeff->Write();
    hPtJet2dDjet->Write();//D-jet eff scaled response
    hPtJet2dDjet_resol->Write();//resol
    hPtJet2d_resol->Write();//resol
    resolobj->Write();//resol
    if(bClosure){
        hPtJet2dDjet_clos->Write();//closure response
        //hMCjetD->Write();
        //hMCjetT->Write();
        hMCjetClosData->Write();// reco level data spec for closure
        hMCjetClosTrue->Write();// gen level true spec unmatched for closure
    }

    ofile->Close();

   return;
}

void SparseToTree(TString MCfile, Bool_t isPostfix, TString postfix)
{
// Read the THnSparse
    TFile *File1 = new TFile(MCfile,"read");
    TDirectoryFile* dir=(TDirectoryFile*)File1->Get("DmesonsForJetCorrelations");
    TString histName;
    if(fDmesonSpecie) histName = "histosDStarMBN";
                    else histName = "histosD0MBN";
    TList *histList[NDMC];
    THnSparseF *sparseMC[NDMC];
    THnSparseF *h = nullptr;

    for(int i=0; i<NDMC; i++){
        histList[i] =  (TList*)dir->Get(Form("%s%d%sMCrec",histName.Data(),i,isPostfix?postfix.Data():""));
        sparseMC[i] = (THnSparseF*)histList[i]->FindObject("ResponseMatrix");
        std::cout<<sparseMC[i]->GetNbins()<<std::endl;
        if(!i)h = dynamic_cast<THnSparseF*>(sparseMC[0]->Clone("ResponseMatrixSum"));
        else h->Add(sparseMC[i]);
    }
    //copy from here: https://root.cern.ch/root/html/tutorials/tree/drawsparse.C.html
    // Creates a TTree and fills it with the coordinates of all
    // filled bins. The tree will have one branch for each dimension,
    // and one for the bin content.
    TString outName(MCfile);
    outName.Remove(outName.Length()-5,5);
    outName+="TTree";
    if(isPostfix)outName+=postfix;
    outName+=".root";
    std::cout<<"Converting THnSparse: "<<MCfile<<std::endl;
    std::cout<<"Into TTree: "<<outName<<std::endl;
   TFile *File2 = new TFile(outName,"RECREATE");
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

   const int percentPrint = 1;
   int totalWindows = h->GetNbins();
   int percent;
   int step = totalWindows / (100/percentPrint);
   int nextPrint = step;

   Int_t *bins = new Int_t[dim];
   for (Long64_t i = 0; i < h->GetNbins(); ++i) {
      if (i >= nextPrint)
       {
           percent = (100 * i) / totalWindows;
           std::cout << "\r" << std::string(percent/percentPrint , '|') << percent << "%";
           std::cout.flush();
           nextPrint += step;
       }
      x[dim] = h->GetBinContent(i, bins);
      for (Int_t d = 0; d < dim; ++d) {
         x[d] = h->GetAxis(d)->GetBinCenter(bins[d]);
      }

      tree->Fill();

   }
   std::cout<<std::endl<<"tree entries "<< tree->GetEntries()<<std::endl;
   File2->Write(nullptr,TObject::kOverwrite);
   File2->Close();
}
