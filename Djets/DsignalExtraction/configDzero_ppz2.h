    const int     fptbinsDN = sizeof(fptbinsDA)/sizeof(fptbinsDA[0])-1;

    //====== jet pT bins --- set up your jet pT bins ======
    //const int 	  fptbinsJetN = 3;
    //double	  fptbinsJetA[fptbinsJetN+1] = {5.0, 7.0, 10.0, 16.0, 36.0, 50.0};
    double	  fptbinsJetA[] = {5.0, 7.0, 10.0, 15.0, 36.0,
         5.0, 15.0, 30.0, 15.0 , 50.0,
        10.0, 16.0, 36.0,  5.0 , 50.0, 
         3.0,  5.0,  3.0,  4.25,  5.0,
         2  ,  4  ,  5  ,  2   ,  5  ,
         1  , 50  };
    
    const int     fptbinsJetN = sizeof(fptbinsJetA)/sizeof(fptbinsJetA[0])-1;

    //---------//---------//---------//---------//---------//---------
    //====== z bins ---- set up your z (momentum fraction) bins ======
    const int     fptbinsZTrueN = 5;
    double        fptbinsZTrueA[fptbinsZTrueN+1] = {0.4, 0.6, 0.7, 0.8, 0.9, 1.02};
    const int     fptbinsZMeasN = 5;
    double        fptbinsZMeasA[fptbinsZMeasN+1] = {0.4, 0.6, 0.7, 0.8, 0.9, 1.02};
    const int     fptbinsZFinalN = 5;
    double        fptbinsZFinalA[fptbinsZFinalN+1] = {0.4,0.6, 0.7, 0.8, 0.9, 1.02};
    const int     fptbinsZSimN = 5;
    double        fptbinsZSimA[fptbinsZSimN+1] = {0.4, 0.6, 0.7, 0.8, 0.9, 1.02};
    //====== z range ---- set up your min, max z ======
    double        zmin = -2., zmax = 2.; // for D-jet pT spectrum

    //====== signal extraction config ======
    Bool_t        fUseRefl = 1;                      //-----! if to use reflections (for D0, you must have reflections files ready)
    Int_t         fbkgtype = 0;                      //-----! kExpo=0, kLin=1, kPol2=2, kNoBk=3, kPow=4, kPowEx=5
    Int_t         fsigmaSignal = 2;                  //-----! sigma for the signal region
    Float_t       fsigmaBkg[] = {-9,-4,4,9};         //-----! sigma for the SB region (both left and right side from the fit)
    Float_t       fDmass = 1.86484, fDsigma = 0.010;   //-----! initial values for D mass and sigma
    //Float_t       fDsigmafix[fptbinsDN] = {0.01,0.011,0.01175,0.0125,0.013,0.0135,0.0145,0.016,0.0175,0.0185};   //-----! initial values for D mass and sigma
    double        minf = 1.71, maxf = 2.1;           //-----! min/mass of the inv. mass distributions
    Int_t         fRebinMass = 2;                    //-----! rebining of the inv. mass distributions

    Int_t fColors[] = {1,2,8,4,kOrange-1,6,kGray+1,kCyan+1,kMagenta+2,kGreen+3,kViolet+5,kYellow+2};
    Int_t fMarkers[] = {20,21,22,23,24,25,26,27,28,29,30,32,33,34};


    //============== POWHEG simulations ============================
    //======= set up here names of your simulation files =======
    TString fRunB[] = {
//Event Gen
       "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1579613250"//evtgen central R02346
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1579857760"//mb4.5
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1579795653"//mb5
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1579792457"//F2R2
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1579698756"//F1R05
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1579622199"//uF=0.5 , uR=0.5
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1580118856"//uF=0.5 , uR=1
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1580141878"//uF= 2, uR=1
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1580122386"//uF=1 , uR=2
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1580201087"//PDF=21200
    };

    TString fDescB[] = {
      "central"
,      "muR=1,muF=0.5" 
,      "muR=0.5,muF=1"
,      "muR=2,muF=1"
,      "muR=1,muF=2"
,      "muR=0.5,muF=0.5"
,      "muR=2,muF=2"
,      "m_{b}=5"
,      "m_{b}=4.5"
,      "PDF 21200"
,      "PDF 10800"
    };

    TString fRunC[] = {
// POWHEG+Pythia6 hvq
//      "AnalysisResults_FastSim_powheg+pythia6_charm_1580461941"
////=pow+pythia8
 "AnalysisResults_FastSim_powheg+pythia8_charm_central_1593017599"//"AnalysisResults_FastSim_powheg+pythia6_charm_1580461941" //central for R02
,"AnalysisResults_FastSim_powheg+pythia8_charm_m13_1593091148"
,"AnalysisResults_FastSim_powheg+pythia8_charm_m17_1593090378"
,"AnalysisResults_FastSim_powheg+pythia8_charm_F2R2_1593354950"
,"AnalysisResults_FastSim_powheg+pythia8_charm_F1R2_1593339352"
,"AnalysisResults_FastSim_powheg+pythia8_charm_F2R1_1593341637"
,"AnalysisResults_FastSim_powheg+pythia8_charm_F05R05_1593334134"
,"AnalysisResults_FastSim_powheg+pythia8_charm_F1R05_1593338244"
,"AnalysisResults_FastSim_powheg+pythia8_charm_F05R1_1593337904"
,"AnalysisResults_FastSim_pythia6_charm_1594235514"
,"AnalysisResults_FastSim_pythia8_charm_1594137862"
,"AnalysisResults_FastSim_pythia8_charm_soft2color_1598639645"
    };
    TString fDescC[] = {
//      "central"
//,      "m_{c}=1.3"
//,      "m_{c}=1.7"
//,      "muR=2,muF=2"
//,      "muR=1,muF=2"
//,      "muR=2,muF=1"
//,      "muR=0.5,muF=0.5"
//,      "muR=1,muF=0.5"
//,      "muR=0.5,muF=1"
       "PP8:central"
,      "PP8:m_{c}=1.3"
,      "PP8:m_{c}=1.7"
,      "PP8:muR=2,muF=2"
,      "PP8:muR=1,muF=2"
,      "PP8:muR=2,muF=1"
,      "PP8:muR=0.5,muF=0.5"
,      "PP8:muR=1,muF=0.5"
,      "PP8:muR=0.5,muF=1"
,      "Pythia6"
,      "Pythia8"
,      "Pythia8:Soft-mode-2"
    };
