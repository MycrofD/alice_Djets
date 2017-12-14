/// \file AddTaskJetUEStudies.C
/// \brief AddTask macro for the AliAnalysisTaskJetUEStudies class.
///
/// AddTask macro for the AliAnalysisTaskJetUEStudies class.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date June 28, 2017

AliAnalysisTaskJetUEStudies* AddTaskJetUEStudies(
//AliAnalysisTaskJetUEStudiesMy* AddTaskJetUEStudies(
    TString        nTracks                        = "usedefault",
    TString        nClusters                      = "usedefault",
    Double_t       trackPtCut                     = 0.15,
    Double_t       clusECut                       = 0.30,
    TString        suffix                         = ""
)
{  
// return AliAnalysisTaskJetUEStudiesMy::AddTaskJetUEStudies(nTracks, nClusters, trackPtCut, clusECut, suffix);
 return AliAnalysisTaskJetUEStudies::AddTaskJetUEStudies(nTracks, nClusters, trackPtCut, clusECut, suffix);
}
