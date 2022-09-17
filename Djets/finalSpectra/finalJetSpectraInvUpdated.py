import os, os.path, sys
import yaml
import ROOT 
import array
import math

##----------------------------------------------------------------
from ROOT import TCanvas, TLegend, TLine, TPad, TLatex, TF1, TH1D
from funcsettings import *
from style import *
##----------------------------------------------------------------
##--- SANITY CHECK
if len(sys.argv)!=8:
    print("""   === Usage example: python file.py type_ R z sysGlobal 5/13 5/13Yaml 5/13sysUnc
                    type_ = 0: jet-pt x-sec, 1: z x-sec, 2: z PDF
                    R = jet radius. 02, 03, 04, 06. (03 only for type:0 x-section)
                    z = which jetpt interval for z: 1,2,3,4,5
                    sysGlobal = 0: (default). add LmiUnc, BRUnc, DmesonTrackingUnc, CUTvar, JES; for x-section
                                1: (R comparison). no global Unc; cut var separate. no JES. only for type:0 x-section
                                2: (energy comparison) add LumiUnc and DmesonTrackingUnc, CUTvar and JES
                    5/13TeV: 5, 13
                    5/13TeV yaml: YAML config file for corresponding energy
                    sysUnc yaml: YAML config with systematics
    python file.py 0/1/2 02/03/04/06 2/3/4/5 0/1/2 5/13 5yaml 5sysuncYaml
    """
    )
    exit()

## ----------------------------------------------------------------
type_ = int(sys.argv[1])
R = sys.argv[2]
zBin = int(sys.argv[3]);
if(type_ == 0): zBin = 0;
sysGlobal = int(sys.argv[4])
energy_int = int(sys.argv[5])
dy = 2*(9.0-int(R))/10.
pdf = {0:False,1:False,2:True}[type_]
#plotRanges=[0,2,10e-9,1]
plotRanges=[0,2.1,7e-7,1]
## ----------------------------------------------------------------
##--- LOAD YAML CONFIG, Read YAML file
yaml_config_file = sys.argv[6]
yaml_unc_config_file = sys.argv[7]
with open(yaml_config_file, 'r') as stream:
    data_loaded = yaml.safe_load(stream)
with open(yaml_unc_config_file, 'r') as unc_stream:
    unc_loaded = yaml.safe_load(unc_stream)
energy=data_loaded['energy']
outSpectraDir=data_loaded['outSpectraDir'][int(bool(type_))][R]
dataFile=data_loaded['dataFile'][int(bool(type_))][R]
dataAnalysisFile=data_loaded['dataAnalysisFile'][R]
histBase=data_loaded['histBase'][int(bool(type_))]
if(type_):
    histBase+=str(zBin-1)
## ----------------------------------------------------------------
## --------------- SET PARAMETERS + BINNING HERE ------------------
## ----------------------------------------------------------------
sigma_in=data_loaded['sigma_in'] 
BRDzero=data_loaded['BRDzero'] 
BRDzeroUnc=data_loaded['BRDzeroUnc'] 
DtrackingUnc=data_loaded['DtrackingUnc'] 
LumiUnc=data_loaded['LumiUnc'] 
##--- INITIAL VALUES
xAxisBins = data_loaded['xAxisBins'][int(bool(type_))]
xAxis = array.array('d',data_loaded['xAxis'][int(bool(type_))])
xAxisC = array.array('d',data_loaded['xAxisC'][int(bool(type_))])
jetpTbins = array.array('d',data_loaded['jetpTbins'][int(bool(type_))])
DpTbins=data_loaded['DpTbins'][int(bool(type_))][R]
##--- create folders
ROOT.gSystem.Exec("mkdir "+outSpectraDir)
if(energy_int==5 and type_ != 0):
    outSpectraDir+=data_loaded['simDirZbin'][zBin]
    ROOT.gSystem.Exec("mkdir "+outSpectraDir)
outPlotDir=outSpectraDir+"/plots/"
ROOT.gSystem.Exec("mkdir "+outPlotDir)
##--- Uncertainties
# unc part 1
#systUncD_up,systUncD_down,systUncD_JES,systUncD_CUTS = np.zeros(xAxisBins),np.zeros(xAxisBins),np.zeros(xAxisBins),np.zeros(xAxisBins)
systUncD_up,systUncD_down,systUncD_JES,systUncD_CUTS =array.array('d',[0]*xAxisBins), array.array('d',[0]*xAxisBins), array.array('d',[0]*xAxisBins), array.array('d',[0]*xAxisBins)
if(type_==0):
    systUncD_up,systUncD_down,systUncD_JES,systUncD_CUTS = unc_loaded[type_][R]['systUncD_up'],unc_loaded[type_][R]['systUncD_down'],unc_loaded[type_][R]['systUncD_JES'],unc_loaded[type_][R]['systUncD_CUTS']
else:
    #systUncD_up,systUncD_down,systUncD_JES,systUncD_CUTS = unc_loaded[1][R][zBin]['systUncD_up'],unc_loaded[1][R][zBin]['systUncD_down'],unc_loaded[1][R][zBin]['systUncD_JES'],unc_loaded[1][R][zBin]['systUncD_CUTS']
    systUncD_up,systUncD_down = unc_loaded[1][R][zBin]['systUncD_up'],unc_loaded[1][R][zBin]['systUncD_down']

# unc part 2
addBRDzeroUnc,addDtrackingUnc,addLumiUnc,separateCUTUnc,addJESUnc=0,0,0,0,0
if(type_==0 or type_==1):
    if(sysGlobal==0):
        addBRDzeroUnc,addDtrackingUnc,addLumiUnc,separateCUTUnc,addJESUnc=1,1,1,0,1
    if(sysGlobal==1):
        addBRDzeroUnc,addDtrackingUnc,addLumiUnc,separateCUTUnc,addJESUnc=0,0,0,1,0
    if(sysGlobal==2):
        addBRDzeroUnc,addDtrackingUnc,addLumiUnc,separateCUTUnc,addJESUnc=0,1,1,0,1
elif(type_==2):
    addBRDzeroUnc,addDtrackingUnc,addLumiUnc,separateCUTUnc,addJESUnc=0,0,0,0,1
    plotRanges[2],plotRanges[3]=-0.3,10
    #plotRanges[0],plotRanges[1]=0,2.1

# unc part 3
bit = -1
print("sys bitmap: ",addBRDzeroUnc,addDtrackingUnc,addLumiUnc,addJESUnc,not separateCUTUnc)
for bin_ in range(xAxisBins):
    totalUnc = 0
    if(addBRDzeroUnc): totalUnc+= BRDzeroUnc**2
    if(addDtrackingUnc): totalUnc+= DtrackingUnc**2
    if(addLumiUnc): totalUnc+= LumiUnc**2
    if(addJESUnc): totalUnc += systUncD_JES[bin_]**2
    if(separateCUTUnc==0): totalUnc += systUncD_CUTS[bin_]**2
    #systUncD_up[bin_] = np.sqrt(systUncD_up[bin_]**2 + totalUnc)
    systUncD_up[bin_] = math.sqrt(systUncD_up[bin_]**2 + totalUnc)
    #systUncD_down[bin_] = np.sqrt(systUncD_down[bin_]**2 + totalUnc)
    systUncD_down[bin_] = math.sqrt(systUncD_down[bin_]**2 + totalUnc)
    #print(systUncD_up[bin_]," ",systUncD_down[bin_])
# unc part 4: 
#           add JES manually for ratio plots.
#           create separate plot for separate CUT SYS
# unc part 5
gDataCutSys = ROOT.TGraph()
if(separateCUTUnc): gDataCutSys = ROOT.TGraph(xAxisBins,xAxisC, array.array('d',systUncD_CUTS))

##--- GET DATA and EVENTS
print("get Data and Events")
File = ROOT.TFile(dataAnalysisFile,"read");
dire = File.Get("PWG3_D2H_DmesonsForJetCorrelationsMBN0")
c = dire.Get("NormalizationCounter")
nEv = c.GetNEventsForNorm()
dataLum = nEv/(sigma_in*1000) ;
simScaling = 0.5
dataScaling = 1. /(BRDzero * dataLum)/2.;
##--- INITIATE/PREPARE CANVAS
canvas, upPad, downPad,placeholder_up,placeholder_down = PrepareCanvas(xAxisBins, xAxis, zBin, pdf, plotRanges)
##--- DATA
hDataSys, hData_binned, hDataSysRatio, hData_binned_ratio = GetData(dataFile,histBase,dataScaling,xAxisBins,xAxis,systUncD_down,systUncD_up,pdf,dy)
PlaceOnPadData(upPad,hDataSys,hData_binned,1.2);
PlaceOnPadData(downPad,hDataSysRatio,hData_binned_ratio,0);
##--- PROMPT SIMULATION
if(type_):
    simDir=data_loaded['simDir'][int(bool(type_))][R]+data_loaded['simDirZbin'][zBin]
else:
    simDir=data_loaded['simDir'][int(bool(type_))][R]
simPrefix=data_loaded['simPrefix'][int(bool(type_))]
ePowhegPythia6=data_loaded['ePowhegPythia6']
ePowhegPythia8=data_loaded['ePowhegPythia8']
ePythia6=data_loaded['ePythia6']
ePythia8=data_loaded['ePythia8']
ePythia8SoftMode2=data_loaded['ePythia8SoftMode2']
ePowhegPythia6dijet=data_loaded['ePowhegPythia6dijet']
ePowhegPythia8dijet=data_loaded['ePowhegPythia8dijet']
#--------------------
fPowhegPythia6=data_loaded['fPowhegPythia6']
fPowhegPythia8=data_loaded['fPowhegPythia8']
fPythia6=data_loaded['fPythia6']
fPythia8=data_loaded['fPythia8']
fPythia8SoftMode2=data_loaded['fPythia8SoftMode2']
fPowhegPythia6dijet=data_loaded['fPowhegPythia6dijet']
fPowhegPythia8dijet=data_loaded['fPowhegPythia8dijet']
##--- powheg pythia 6
#--------------------
if(ePowhegPythia6):
    print("get POWHEG+PYTHIA6")
    simPowhegPythia6, simPowhegPythia6_cent, simPowhegPythia6_up, simPowhegPythia6_down,simPowhegPythia6var = GetSim("simPowhegPythia6",type_, 9, fPowhegPythia6, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf,dy,simPrefix,energy_int)
    simPowhegPythia6_R, simPowhegPythia6_cent_R, simPowhegPythia6_up_R, simPowhegPythia6_down_R = GetDataSimRatio("simPowhegPythia6",hData_binned,simPowhegPythia6_cent, simPowhegPythia6_up, simPowhegPythia6_down, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf,dy)
    PlaceOnPadSim(upPad,simPowhegPythia6,ROOT.TColor.GetColor("#000099"),24,1)
    PlaceOnPadSim(dowmPad,simPowhegPythia6_R,ROOT.TColor.GetColor("#000099"),24,1)
##--powheg pythia 8
#--------------------
if(ePowhegPythia8):
    print("get POWHEG+PYTHIA8")
    simPowhegPythia8, simPowhegPythia8_cent, simPowhegPythia8_up, simPowhegPythia8_down,simPowhegPythia8var = GetSim("simPowhegPythia8",type_, 9, fPowhegPythia8, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf,dy,simPrefix,energy_int)
    simPowhegPythia8_R, simPowhegPythia8_cent_R, simPowhegPythia8_up_R, simPowhegPythia8_down_R = GetDataSimRatio("simPowhegPythia8",hData_binned,simPowhegPythia8_cent, simPowhegPythia8_up, simPowhegPythia8_down, xAxisBins, xAxis)
    PlaceOnPadSim(upPad,simPowhegPythia8,ROOT.TColor.GetColor("#000099"),24,1)
    PlaceOnPadSim(downPad,simPowhegPythia8_R,ROOT.TColor.GetColor("#000099"),24,1)
##-- pythia 6
#--------------------
if(ePythia6):
    print("get PYTHIA6")
    simPythia6, simPythia6_cent, simPythia6_up, simPythia6_down,hBlackHole = GetSim("simPythia6",type_, 1,fPythia6, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf, dy,simPrefix,energy_int)
    simPythia6_R, simPythia6_cent_R, simPythia6_up_R, simPythia6_down_R = GetDataSimRatio("simPythia6",hData_binned,simPythia6_cent, simPythia6_up, simPythia6_down, xAxisBins, xAxis)
    PlaceOnPadSim(upPad,simPythia6,ROOT.TColor.GetColor("#009933"),25,2);
    PlaceOnPadSim(downPad,simPythia6_R,ROOT.TColor.GetColor("#009933"),25,2);
##-- pythia 8
#--------------------
if(ePythia8):
    print("get PYTHIA8")
    simPythia8, simPythia8_cent, simPythia8_up, simPythia8_down,hBlackHole = GetSim("simPythia8",type_, 1,fPythia8, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf, dy, simPrefix,energy_int)
    simPythia8_R, simPythia8_cent_R, simPythia8_up_R, simPythia8_down_R = GetDataSimRatio("simPythia8",hData_binned,simPythia8_cent, simPythia8_up, simPythia8_down, xAxisBins, xAxis)
    #PlaceOnPadSim(upPad,simPythia8,ROOT.TColor.GetColor("#009933"),27,2);
    #PlaceOnPadSim(downPad,simPythia8_R,ROOT.TColor.GetColor("#009933"),27,2);
    #PlaceOnPadSim(upPad,simPythia8,ROOT.kViolet+2,27,6);
    #PlaceOnPadSim(downPad,simPythia8_R,ROOT.kViolet+2,27,6);
    PlaceOnPadSim(upPad,simPythia8,ROOT.kOrange+7,27,6);
    PlaceOnPadSim(downPad,simPythia8_R,ROOT.kOrange+7,27,6);
##-- pythia 8 soft mode 2
#------------------------
if(ePythia8SoftMode2):
    print("get PYTHIA8 soft mode2")
    simPythia8Soft2, simPythia8Soft2_cent, simPythia8Soft2_up, simPythia8Soft2_down,hBlackHole = GetSim("simPythia8Soft2",type_, 1,fPythia8SoftMode2, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf, dy, simPrefix, energy_int)
    simPythia8Soft2_R, simPythia8Soft2_cent_R, simPythia8Soft2_up_R, simPythia8Soft2_down_R = GetDataSimRatio("simPythia8Soft2",hData_binned,simPythia8Soft2_cent, simPythia8Soft2_up, simPythia8Soft2_down, xAxisBins, xAxis)
    #PlaceOnPadSim(upPad,simPythia8Soft2,ROOT.kOrange+2,28,4)
    #PlaceOnPadSim(downPad,simPythia8Soft2_R,ROOT.kOrange+2,28,4);
    PlaceOnPadSim(upPad,simPythia8Soft2,ROOT.kGreen+3,28,2)
    PlaceOnPadSim(downPad,simPythia8Soft2_R,ROOT.kGreen+3,28,2);
##-- powheg pythia 6 dijet
#-------------------------
if(ePowhegPythia6dijet):
    print("get POWHEG+PYTHIA6 dijet")
    simPowhegPythia6dijet, simPowhegPythia6dijet_cent, simPowhegPythia6dijet_up, simPowhegPythia6dijet_down,simPowhegPythia6dijetvar = GetSim("simPowhegPythia6dijet",type_, 1, fPowhegPythia6dijet, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf, dy, simPrefix, energy_int)
    simPowhegPythia6dijet_R, simPowhegPythia6dijet_cent_R, simPowhegPythia6dijet_up_R, simPowhegPythia6dijet_down_R = GetDataSimRatio("simPowhegPythia6dijet",hData_binned,simPowhegPythia6dijet_cent, simPowhegPythia6dijet_up, simPowhegPythia6dijet_down, xAxisBins, xAxis)
    PlaceOnPadSim(upPad,simPowhegPythia6dijet,ROOT.TColor.GetColor("#000099"),28,1)
    PlaceOnPadSim(downPad,simPowhegPythia6dijet_R,ROOT.TColor.GetColor("#000099"),28,1)
##-- powheg pythia 8 dijet
#-------------------------
if(ePowhegPythia8dijet):
    print("get POWHEG+PYTHIA8 dijet")
    simPowhegPythia8dijet, simPowhegPythia8dijet_cent, simPowhegPythia8dijet_up, simPowhegPythia8dijet_down,simPowhegPythia8dijetvar = GetSim("simPowhegPythia8dijet",type_, 1, fPowhegPythia8dijet, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf, dy, simPrefix, energy_int)
    simPowhegPythia8dijet_R, simPowhegPythia8dijet_cent_R, simPowhegPythia8dijet_up_R, simPowhegPythia8dijet_down_R = GetDataSimRatio("simPowhegPythia8dijet",hData_binned,simPowhegPythia8dijet_cent, simPowhegPythia8dijet_up, simPowhegPythia8dijet_down, xAxisBins, xAxis)
    PlaceOnPadSim(upPad,simPowhegPythia8dijet,ROOT.TColor.GetColor("#009999"),28,1)
    PlaceOnPadSim(downPad,simPowhegPythia8dijet_R,ROOT.TColor.GetColor("#009999"),28,1)
"""
"""
###############################################
#  -------- Legend and text PAD------
###############################################
leg = ROOT.TLegend()
shift = 0.06*(ePowhegPythia6+ePowhegPythia8+ePythia6+ePythia8+ePythia8SoftMode2+ePowhegPythia6dijet)
if(type_==0):leg = ROOT.TLegend(0.35,0.45,0.65,0.7,"","NB NDC");
#if(type_==1 or type_ ==2):leg = ROOT.TLegend(0.22,0.7-shift,0.5,0.7,"","NB NDC");
if(type_==1 or type_ ==2):leg = ROOT.TLegend(0.22,0.65-shift,0.5,0.7,"","NB NDC");
leg.SetBorderSize(0)
leg.SetTextFont(43)
leg.SetTextSize(21)
leg.SetLineColor(1)
leg.SetLineStyle(1)
leg.SetLineWidth(1)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.AddEntry(hDataSys,"Data","fp")
if(ePowhegPythia6):leg.AddEntry(simPowhegPythia6,"POWHEG hvq + PYTHIA 6","pf")
if(ePowhegPythia8):leg.AddEntry(simPowhegPythia8,"POWHEG hvq + PYTHIA 8","pf")
if(ePythia6):leg.AddEntry(simPythia6,"PYTHIA 6 Perugia 2011","l")
#if(ePythia8):leg.AddEntry(simPythia8,"PYTHIA 8 Monash 2013","l")
if(ePythia8):leg.AddEntry(simPythia8,"PYTHIA 8 HardQCD Monash 2013","l")
if(ePythia8SoftMode2):leg.AddEntry(simPythia8Soft2,"PYTHIA 8 SoftQCD Mode 2","l")
if(ePowhegPythia6dijet):leg.AddEntry(simPowhegPythia6dijet,"POWHEG dijet + PYTHIA 8 hadi","pf")
if(ePowhegPythia8dijet):leg.AddEntry(simPowhegPythia8dijet,"POWHEG dijet + PYTHIA 8","l")
pt = [ROOT.TPaveText(),ROOT.TPaveText()]
pt[0] = ROOT.TPaveText(0.2,0.9,0.85,0.95,"NB NDC")
if(type_==0):pt[1] = ROOT.TPaveText(0.2,0.75,0.85,0.9,"NB NDC")
if(type_==1 or type_==2):pt[1] = ROOT.TPaveText(0.2,0.75,0.85,0.9,"NB NDC")
for s in range(len(pt)):
    pt[s].SetBorderSize(0)
    pt[s].SetFillStyle(0)
    pt[s].SetTextAlign(13)
    pt[s].SetTextFont(43)
    pt[s].SetTextSize(22)

#pt[0]->AddText("ALICE Preliminary"); //uncomment
pt[0].AddText("ALICE, pp, #sqrt{#it{s}} = "+energy+" TeV")
pt[1].AddText("charged jets, anti-#it{k}_{T}, #it{R} = 0.%d, |#it{#eta}_{ch jet}| < 0.%d"%(int(R),9-int(R)))
if(type_==0):pt[1].AddText("with D^{0}, %d < #it{p}_{T,D^{0}} < %d GeV/#it{c}"%(DpTbins[0],DpTbins[1]))
#if(type_==1 or type_==2):pt[1].AddText("%d < #it{p}_{T,ch. jet} < %d GeV/#it{c}"%(jetpTbins[zBin-1],jetpTbins[zBin]));
#if(type_==1 or type_==2):pt[1].AddText("with D^{0}, %d < #it{p}_{T,D^{0}} < %d GeV/#it{c}"%(DpTbins[0][zBin-1],DpTbins[1][zBin-1]))
if(type_==1 or type_==2):pt[1].AddText("%d < #it{p}_{T,ch jet}< %d GeV/#it{c} with D^{0}, %d < #it{p}_{T,D^{0}} < %d GeV/#it{c}"%(jetpTbins[zBin-1],jetpTbins[zBin],DpTbins[0][zBin-1],DpTbins[1][zBin-1]))
#if(type_==1 or type_==2):pt[1].AddText(""%(DpTbins[0][zBin-1],DpTbins[1][zBin-1]))
upPad.cd()
leg.Draw()
pt[0].Draw()
pt[1].Draw()
#################################################
#  -----------  TERMINATE CANVAS  ---------------
#################################################
TerminateCanvas(upPad, downPad, placeholder_up, placeholder_down)
sysmode = ""
if(type_==0 or type_==1):
    if(sysGlobal==0): sysmode = "_fullGlobal_addedCUTandJES"
    elif(sysGlobal==1): sysmode = "_noneGlobal_separateCUTnoneJES"
    elif(sysGlobal==2): sysmode = "_noBRUnc_addedCUTandJES"
elif(type_==2):
    sysmode = "_PDF_noneGlobal_addedCUTandJES"

canvas.SaveAs(outPlotDir+"/test0_finalSpectra_"+sysmode+".png")
canvas.SaveAs(outPlotDir+"/test0_finalSpectra_"+sysmode+".pdf")
canvas.SaveAs(outPlotDir+"/test0_finalSpectra_"+sysmode+".root")
canvas.SaveAs(outPlotDir+"/test0_finalSpectra_"+sysmode+".C")
##########
for i in range(hData_binned.GetNbinsX()+1):
    print(hData_binned.GetBinContent(i+1))
    print(hData_binned.GetBinCenter(i+1))
    print("okay==")
##########
ofile = ROOT.TFile(outSpectraDir+"JetPtSpectrum_final"+sysmode+".root","recreate")
hData_binned.Write()
hDataSys.Write()
if(separateCUTUnc): gDataCutSys.Write("gDataCutSys")
hData_binned_ratio.Write()
hDataSysRatio.Write()
if(ePowhegPythia6):
    simPowhegPythia6_cent.Write()
    simPowhegPythia6_up.Write()
    simPowhegPythia6_down.Write()
    simPowhegPythia6.Write()
    simPowhegPythia6_cent_R.Write()
    simPowhegPythia6_up_R.Write()
    simPowhegPythia6_down_R.Write()
    simPowhegPythia6_R.Write()
    for ivar in range(1,9):
        simPowhegPythia6var[ivar].Write()

if(ePowhegPythia8):
    simPowhegPythia8_cent.Write()
    simPowhegPythia8_up.Write()
    simPowhegPythia8_down.Write()
    simPowhegPythia8.Write()
    simPowhegPythia8_cent_R.Write()
    simPowhegPythia8_up_R.Write()
    simPowhegPythia8_down_R.Write()
    simPowhegPythia8_R.Write()
    for ivar in range(1,9):
        simPowhegPythia8var[ivar].Write()

if(ePythia6):
    simPythia6_cent.Write()
    simPythia6_up.Write()
    simPythia6_down.Write()
    simPythia6.Write()
    simPythia6_cent_R.Write()
    simPythia6_up_R.Write()
    simPythia6_down_R.Write()
    simPythia6_R.Write()

if(ePythia8):
    simPythia8_cent.Write()
    simPythia8_up.Write()
    simPythia8_down.Write()
    simPythia8.Write()
    simPythia8_cent_R.Write()
    simPythia8_up_R.Write()
    simPythia8_down_R.Write()
    simPythia8_R.Write()

if(ePythia8SoftMode2):
    simPythia8Soft2_cent.Write()
    simPythia8Soft2_up.Write()
    simPythia8Soft2_down.Write()
    simPythia8Soft2.Write()
    simPythia8Soft2_cent_R.Write()
    simPythia8Soft2_up_R.Write()
    simPythia8Soft2_down_R.Write()
    simPythia8Soft2_R.Write()

if(ePowhegPythia6dijet):
    simPowhegPythia6dijet_cent.Write()
    simPowhegPythia6dijet_up.Write()
    simPowhegPythia6dijet_down.Write()
    simPowhegPythia6dijet.Write()
    simPowhegPythia6dijet_cent_R.Write()
    simPowhegPythia6dijet_up_R.Write()
    simPowhegPythia6dijet_down_R.Write()
    simPowhegPythia6dijet_R.Write()
    for ivar in range(1,9):
        simPowhegPythia6dijetvar[ivar].Write()

if(ePowhegPythia8dijet):
    simPowhegPythia8dijet_cent.Write()
    #simPowhegPythia8dijet_up.Write()
    #simPowhegPythia8dijet_down.Write()
    simPowhegPythia8dijet.Write()
    simPowhegPythia8dijet_cent_R.Write()
    #simPowhegPythia8dijet_up_R.Write()
    #simPowhegPythia8dijet_down_R.Write()
    simPowhegPythia8dijet_R.Write()
    #for ivar in range(1,9):
    #    simPowhegPythia8dijetvar[ivar].Write()

ofile.Close()

""" 
"""
##--- ---- ---------- ----------------- -------------------- ENDGAME
print("press ENTER to exit");wait=input()
