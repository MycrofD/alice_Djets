import ROOT

defP=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR02_paperCuts/Default/efficiency/DjetEff_prompt_jetpt5_50.root","read")
ct1P=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/CutSys2Final_DzeroR02_paperCuts/Default/efficiency/DjetEff_prompt_jetpt5_50.root","read")
ct2P=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/CutSys3Final_DzeroR02_paperCuts/Default/efficiency/DjetEff_prompt_jetpt5_50.root","read")
ct3P=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/CutSys5Final_DzeroR02_paperCuts/Default/efficiency/DjetEff_prompt_jetpt5_50.root","read")
ct4P=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/CutSys6Final_DzeroR02_paperCuts/Default/efficiency/DjetEff_prompt_jetpt5_50.root","read")


defN=ROOT.TFile("/home/jackbauer/Work/alice/analysis/pp5TeV/D0jet/results_APW/Final_DzeroR02_paperCuts/Default/efficiency/DjetEff_nonPrompt_jetpt5_50.root","read")

cP = ROOT.TCanvas('cP','cP',500,400)
h0 = defP.Get('hEff_reb')
h1 = ct1P.Get('hEff_reb')
h2 = ct2P.Get('hEff_reb')
h3 = ct3P.Get('hEff_reb')
h4 = ct4P.Get('hEff_reb')

h1.Divide(h0)
h2.Divide(h0)
h3.Divide(h0)
h4.Divide(h0)

h1.SetLineColor(ROOT.kBlue+2)
h2.SetLineColor(ROOT.kGreen+2)
h3.SetLineColor(ROOT.kOrange+1)
h4.SetLineColor(ROOT.kViolet+1)


h1.Draw()
h2.Draw('same')
h3.Draw('same')
h4.Draw('same')


wait=input()
