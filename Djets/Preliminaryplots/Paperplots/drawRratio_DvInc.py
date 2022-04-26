import ROOT

# this is used after drawDvIncJet.py

radii = [2,3,4,6]
files, hs = [], []
for i in range(len(radii)):
    files.append(ROOT.TFile("DvInc_R0"+str(radii[i])+".root","read"))
    hs.append(files[i].Get("hR_Data_R0"+str(radii[i])))

hR2by4 = hs[0].Clone("hR2by4")
hR2by4.Divide(hs[2])
hR2by4.SetLineColor(ROOT.kRed+2)

hR2by6 = hs[0].Clone("hR2by6")
hR2by6.Divide(hs[3])
hR2by6.SetLineColor(ROOT.kBlue+2)
hR2by6.SetMarkerColor(ROOT.kBlue+2)
#################
c1 = ROOT.TCanvas("R02by04-DbyInc-jets","c1",600,400)
hR2by4.GetYaxis().SetRangeUser(0,4)
hR2by4.GetYaxis().SetTitle("#frac{0.2 D-jets}{0.4 D-jets}/#frac{0.2 inc.jets}{0.4 inc.jets}")
hR2by4.Draw()
c1.Update()
wait=input()
#################
c2 = ROOT.TCanvas("R02by06-DbyInc-jets","c2",600,400)
hR2by6.GetYaxis().SetRangeUser(0,4)
hR2by6.GetYaxis().SetTitle("#frac{0.2 D-jets}{0.6 D-jets}/#frac{0.2 inc.jets}{0.6 inc.jets}")
hR2by6.Draw()
c2.Update()
wait=input()

