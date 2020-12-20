import os, os.path, sys
import yaml
import matplotlib.pyplot as plt
import ROOT as RT
import ROOT 
import rootpy as rp
import numpy as np
import scipy as sp
from rootpy.io import root_open
import array
#import root_numpy as rtnp
from root_numpy import hist2array
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

from matplotlib import colors as mcolors
import seaborn as sns
##----------------------------------------------------------------
from ROOT import TCanvas, TLegend, TLine, TPad, TLatex, TF1, TH1D
from funcsettings import *
from style import *
##----------------------------------------------------------------
##--- SANITY CHECK
if len(sys.argv)!=7:
    print("""   === Usage example: python file.py type_ R z sysGlobal 5/13TeV sysUnc
                    type_ = 0: jet-pt x-sec, 1: z x-sec, 2: z PDF
                    R = jet radius. 02, 03, 04, 06. (03 only for type:0 x-section)
                    z = which jetpt interval for z: 1,2,3,4,5
                    sysGlobal = 0: (default). add LmiUnc, BRUnc, DmesonTrackingUnc, CUTvar, JES; for x-section
                                1: (R comparison). no global Unc; cut var separate. no JES. only for type:0 x-section
                                2: (energy comparison) add LumiUnc and DmesonTrackingUnc, CUTvar and JES
                    5/13TeV: YAML config file for corresponding energy
                    sysUnc: YAML config with systematics
    """
    )
    exit()

## ----------------------------------------------------------------
type_=int(sys.argv[1])
R=sys.argv[2]
zBin=int(sys.argv[3]);
if(type_ == 0): zBin = 0;
sysGlobal=int(sys.argv[4])
dy=2*(9.0-int(R))/10.
pdf = {0:False,1:False,2:True}[type_]
plotRanges=[0,2,10e-9,1]
## ----------------------------------------------------------------
##--- LOAD YAML CONFIG, Read YAML file
yaml_config_file = sys.argv[5]
yaml_unc_config_file = sys.argv[6]
with open(yaml_config_file, 'r') as stream:
    data_loaded = yaml.safe_load(stream)
with open(yaml_unc_config_file, 'r') as unc_stream:
    unc_loaded = yaml.safe_load(unc_stream)
energy=data_loaded['energy']
outSpectraDir=data_loaded['outSpectraDir'][R]
dataFile=data_loaded['dataFile'][R]
dataAnalysisFile=data_loaded['dataAnalysisFile'][R]
histBase="unfoldedSpectrum";
## ----------------------------------------------------------------
## --------------- SET PARAMETERS + BINNING HERE ------------------
## ----------------------------------------------------------------
sigma_in=data_loaded['sigma_in'] #sigma_in = 0.0578;
BRDzero=data_loaded['BRDzero'] #0.0389;
BRDzeroUnc=data_loaded['BRDzeroUnc'] #0.0004;
DtrackingUnc=data_loaded['DtrackingUnc'] #DtrackingUnc = 0.021;
LumiUnc=data_loaded['LumiUnc'] 
##--- INITIAL VALUES
xAxisBins = data_loaded['xAxisBins'][type_]
xAxis = array.array('d',data_loaded['xAxis'][type_])
xAxisC = array.array('d',data_loaded['xAxisC'][type_])
jetpTbins = array.array('d',data_loaded['jetpTbins'][type_])
DpTbins=array.array('d',data_loaded['DpTbins'][type_][R])
##--- Uncertainties
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
    plotRanges[2],plotRanges[3]=0,10
systUncD_up,systUncD_down = np.zeros(xAxisBins),np.zeros(xAxisBins)
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
PlaceOnPadData(upPad,hDataSys,hData_binned,0.9);
PlaceOnPadData(downPad,hDataSysRatio,hData_binned_ratio,0);
##--- PROMPT SIMULATION
simDir=data_loaded['simDir']
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
    print("get POWHWG+PYTHIA6")
    simPowhegPythia6, simPowhegPythia6_cent, simPowhegPythia6_up, simPowhegPythia6_down,simPowhegPythia6var = GetSim("simPowhegPythia6",type_, 9, fPowhegPythia6, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf,dy)
    simPowhegPythia6_R, simPowhegPythia6_cent_R, simPowhegPythia6_up_R, simPowhegPythia6_down_R = GetDataSimRatio("simPowhegPythia6",hData_binned,simPowhegPythia6_cent, simPowhegPythia6_up, simPowhegPythia6_down, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf,dy)
    PlaceOnPadSim(upPad,simPowhegPythia6,ROOT.TColor.GetColor("#000099"),24,1)
    PlaceOnPadSim(dowmPad,simPowhegPythia6_R,ROOT.TColor.GetColor("#000099"),24,1)
##--powheg pythia 8
#--------------------
if(ePowhegPythia8):
    print("get POWHWG+PYTHIA8")
    simPowhegPythia8, simPowhegPythia8_cent, simPowhegPythia8_up, simPowhegPythia8_down,simPowhegPythia8var = GetSim("simPowhegPythia8",type_, 9, fPowhegPythia8, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf,dy)
    simPowhegPythia8_R, simPowhegPythia8_cent_R, simPowhegPythia8_up_R, simPowhegPythia8_down_R = GetDataSimRatio("simPowhegPythia8",hData_binned,simPowhegPythia8_cent, simPowhegPythia8_up, simPowhegPythia8_down, xAxisBins, xAxis)#, DpTbins, jetpTbins, zBin, pdf,dy)
    PlaceOnPadSim(upPad,simPowhegPythia8,ROOT.TColor.GetColor("#000099"),24,1)
    PlaceOnPadSim(downPad,simPowhegPythia8_R,ROOT.TColor.GetColor("#000099"),24,1)
##-- pythia 6
#--------------------
if(ePythia6):
    print("get PYTHIA6")
    simPythia6, simPythia6_cent, simPythia6_up, simPythia6_down,hBlackHole = GetSim("simPythia6",type_, 1,fPythia6, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf, dy)
    simPythia6_R, simPythia6_cent_R, simPythia6_up_R, simPythia6_down_R = GetDataSimRatio("simPythia6",hData_binned,simPythia6_cent, simPythia6_up, simPythia6_down, xAxisBins, xAxis)#, DpTbins, jetpTbins, zBin, pdf, dy)
    PlaceOnPadSim(upPad,simPythia6,ROOT.TColor.GetColor("#009933"),25,2);
    PlaceOnPadSim(downPad,simPythia6_R,ROOT.TColor.GetColor("#009933"),25,2);
##-- pythia 8
#--------------------
if(ePythia8):
    print("get PYTHIA8")
    simPythia8, simPythia8_cent, simPythia8_up, simPythia8_down,hBlackHole = GetSim("simPythia8",type_, 1,fPythia8, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf, dy)
    simPythia8_R, simPythia8_cent_R, simPythia8_up_R, simPythia8_down_R = GetDataSimRatio("simPythia8",hData_binned,simPythia8_cent, simPythia8_up, simPythia8_down, xAxisBins, xAxis)#, DpTbins, jetpTbins, zBin, pdf, dy)
    PlaceOnPadSim(upPad,simPythia8,ROOT.TColor.GetColor("#009933"),27,2);
    PlaceOnPadSim(downPad,simPythia8_R,ROOT.TColor.GetColor("#009933"),27,2);
##-- pythia 8 soft mode 2
#------------------------
if(ePythia8SoftMode2):
    print("get PYTHIA8 soft mode2")
    simPythia8Soft2, simPythia8Soft2_cent, simPythia8Soft2_up, simPythia8Soft2_down,hBlackHole = GetSim("simPythia8Soft2",type_, 1,fPythia8SoftMode2, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf, dy)
    simPythia8Soft2_R, simPythia8Soft2_cent_R, simPythia8Soft2_up_R, simPythia8Soft2_down_R = GetDataSimRatio("simPythia8Soft2",hData_binned,simPythia8Soft2_cent, simPythia8Soft2_up, simPythia8Soft2_down, xAxisBins, xAxis)
    PlaceOnPadSim(upPad,simPythia8Soft2,ROOT.kOrange+2,28,4)
    PlaceOnPadSim(downPad,simPythia8Soft2_R,ROOT.kOrange+2,28,4);
##-- powheg pythia 6 dijet
#-------------------------
if(ePowhegPythia6dijet):
    print("get POWHWG+PYTHIA6 dijet")
    simPowhegPythia6dijet, simPowhegPythia6dijet_cent, simPowhegPythia6dijet_up, simPowhegPythia6dijet_down,simPowhegPythia6dijetvar = GetSim("simPowhegPythia6dijet",type_, 1, fPowhegPythia6dijet, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf, dy)
    simPowhegPythia6dijet_R, simPowhegPythia6dijet_cent_R, simPowhegPythia6dijet_up_R, simPowhegPythia6dijet_down_R = GetDataSimRatio("simPowhegPythia6dijet",hData_binned,simPowhegPythia6dijet_cent, simPowhegPythia6dijet_up, simPowhegPythia6dijet_down, xAxisBins, xAxis)
    PlaceOnPadSim(upPad,simPowhegPythia6dijet,ROOT.TColor.GetColor("#000099"),28,1)
    PlaceOnPadSim(downPad,simPowhegPythia6dijet_R,ROOT.TColor.GetColor("#000099"),28,1)
##-- powheg pythia 8 dijet
#-------------------------
if(ePowhegPythia8dijet):
    print("get POWHWG+PYTHIA8 dijet")
    simPowhegPythia8dijet, simPowhegPythia8dijet_cent, simPowhegPythia8dijet_up, simPowhegPythia8dijet_down,simPowhegPythia8dijetvar = GetSim("simPowhegPythia8dijet",type_, 1, fPowhegPythia8dijet, simDir, simScaling, xAxisBins, xAxis, DpTbins, jetpTbins, zBin, pdf, dy)
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
if(type_==1 or type_ ==2):leg = ROOT.TLegend(0.22,0.7-shift,0.5,0.7,"","NB NDC");
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
if(ePythia8):leg.AddEntry(simPythia8,"PYTHIA 8 Monash 2013","l")
if(ePythia8SoftMode2):leg.AddEntry(simPythia8Soft2,"PYTHIA 8 Monash 2013 mode 2","l")
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
pt[1].AddText("Charged Jets, anti-#it{k}_{T}, #it{R} = 0.%d, |#it{#eta}_{lab}^{jet}| < 0.%d"%(int(R),9-int(R)))
if(type_==0):pt[1].AddText("with D^{0}, %d < #it{p}_{T,D^{0}} < %d GeV/#it{c}"%(DpTbins[0],DpTbins[1]))
#if(type_==1 or type_==2):pt[1].AddText("%d < #it{p}_{T,jet} < %d GeV/#it{c}"%(jetpTbins[zBin-1],jetpTbins[zBin]));
#if(type_==1 or type_==2):pt[1].AddText("with D^{0}, %d < #it{p}_{T,D^{0}} < %d GeV/#it{c}"%(DpTbins[0][zBin-1],DpTbins[1][zBin-1]))
if(type_==1 or type_==2):pt[1].AddText("%d < #it{p}_{T,jet} < %d GeV/#it{c} with D^{0}, %d < #it{p}_{T,D^{0}} < %d GeV/#it{c}"%(jetpTbins[zBin-1],jetpTbins[zBin],DpTbins[0][zBin-1],DpTbins[1][zBin-1]))
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

#canvas.SaveAs()
""" 
"""
##--- ---- ---------- ----------------- -------------------- ENDGAME
print("press ENTER to exit");wait=input()
