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
//       "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1536649348"//evtgen central R0346
       "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1579613250"//evtgen central R02346
////POWHEG+PYTHIA6
////,      "AnalysisResults_FastSim_powheg+pythia6_beauty_150593961473"	//central R=0.3
// ,      "AnalysisResults_FastSim_powheg+pythia6_beauty_1554981915"	//central R=0.3,0.4,0.6
//,      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504284947"	//R1F05
//,      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504259653"	//R05F1
//,      "AnalysisResults_FastSim_powheg+pythia6_beauty_1506803374"	//R2F1
//,      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504296768"	//R1F2
//,      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504212024"	//R05F05
//,      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504318569"	//R2F2
//,      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504202511"	//mb5
//,      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504197966"	//mb4.5
//,      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504197460"	//pdf 21200
//,      "AnalysisResults_FastSim_powheg+pythia6_beauty_1504199953"	//pdf 10800
//Event Gen
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1579857760"//mb4.5
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1579795653"//mb5
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1579792457"//F2R2
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1579698756"//F1R05
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1579622199"//uF=0.5 , uR=0.5
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1580118856"//uF=0.5 , uR=1
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1580141878"//uF= 2, uR=1
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1580122386"//uF=1 , uR=2
,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1580201087"//PDF=21200
// ,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1549115009"//mb4.5
// ,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1549110628"//mb5
// ,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1549275841"//F2R2
// ,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1549119403"//F1R05
// ,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1549019046"//uF=0.5 , uR=0.5
// ,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1549035201"//uF=0.5 , uR=1
// ,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1549219132"//uF= 2, uR=1
// ,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1549203018"//uF=1 , uR=2
// ,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1550591180"//PDF=21200
//,      "AnalysisResults_FastSim_powheg+pythia6+evtgen_beauty_1550780125"//PDF=10042
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
//,      "Evt gen"
    };

    TString fRunC[] = {
// POWHEG+Pythia6 hvq
      "AnalysisResults_FastSim_powheg+pythia6_charm_1580461941"
//-----      "AnalysisResults_FastSim_powheg+pythia6_charm_central"
//-----,"AnalysisResults_FastSim_powheg+pythia6_charm_m13_1536595965"
//-----,"AnalysisResults_FastSim_powheg+pythia6_charm_m17_1536655729"
//-----,"AnalysisResults_FastSim_powheg+pythia6_charm_F2R2_1535895146"
//-----,"AnalysisResults_FastSim_powheg+pythia6_charm_F1R2_1536594271"
//-----,"AnalysisResults_FastSim_powheg+pythia6_charm_F2R1_1535916012"
//-----,"AnalysisResults_FastSim_powheg+pythia6_charm_F05R05_1535894261"
//-----,"AnalysisResults_FastSim_powheg+pythia6_charm_F1R05_1536598175"
//-----,"AnalysisResults_FastSim_powheg+pythia6_charm_F05R1_1536604800"
//// POWHEG+Pythia6 dijets
//,"AnalysisResults_FastSim_powheg+pythia6_dijet_bb160_1553275404"
//,"AnalysisResults_FastSim_powheg+pythia6_dijet_central_1553968405"
//,"AnalysisResults_FastSim_powheg+pythia6_dijet_F1R2_1557481772"
//,"AnalysisResults_FastSim_powheg+pythia6_dijet_F1R05_1556385744"
//,"AnalysisResults_FastSim_powheg+pythia6_dijet_F2R1_1557918376"
//,"AnalysisResults_FastSim_powheg+pythia6_dijet_F2R2_1560328228"
//,"AnalysisResults_FastSim_powheg+pythia6_dijet_F05R1_1558465083"
//,"AnalysisResults_FastSim_powheg+pythia6_dijet_F05R05_1556195307"
//,"AnalysisResults_FastSim_powheg+pythia6_dijet_mc13_1554197229"
//,"AnalysisResults_FastSim_powheg+pythia6_dijet_mc17_1558464708"
//,"AnalysisResults_FastSim_powheg+pythia6_dijet_PDF212_1561024691"
//Pythia6
//-----,"AnalysisResults_FastSim_pythia6_charm_1552224706"
//,"AnalysisResults_FastSim_pythia6_dijet_1552223980"
//,"AnalysisResults_FastSim_pythia6_mb_1563284959"
//Pythia8
//----- ,"AnalysisResults_FastSim_pythia8_charm_1552144258"
//,"AnalysisResults_FastSim_pythia8_dijet_1552144719"
//,"AnalysisResults_FastSim_pythia8_mb_1563297890"
////Herwig
//,"AnalysisResults_FastSim_herwig_charm_lo_1548973692"
//,"AnalysisResults_FastSim_herwig_dijet_lo_1563298593"
//,"AnalysisResults_FastSim_herwig_mb_1551799518"
    };
    TString fDescC[] = {
      "central"
,      "m_{c}=1.3"
,      "m_{c}=1.7"
,      "muR=2,muF=2"
,      "muR=1,muF=2"
,      "muR=2,muF=1"
,      "muR=0.5,muF=0.5"
,      "muR=1,muF=0.5"
,      "muR=0.5,muF=1"
    };
