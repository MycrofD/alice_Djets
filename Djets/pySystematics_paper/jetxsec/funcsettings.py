import numpy as np
from scipy.optimize import curve_fit
import ROOT

def HistoStyle(histo,linestyle,linewidth,markerstyle,markersize,color):
    histo.SetLineStyle(linestyle)
    histo.SetLineWidth(linewidth)
    histo.SetMarkerSize(markersize)
    histo.SetLineColor(color)
    histo.SetMarkerColor(color)

def GetDigitTextFromRatio(histo,bincenters):
    listx=""
    for i in bincenters:
        x = (histo.GetBinContent(histo.FindBin(i)) - 1)*100
        listx+="%.1f"%(abs(x))
        listx+=", "
    print(listx)
    return listx[:-2]

def GetDigitText(histo,bincenters):
    listx=""
    for i in bincenters:
        x = (histo.GetBinContent(histo.FindBin(i)))*100
        listx+="%.1f"%(abs(x))
        listx+=", "
    listx = listx[:-2]+" in %"
    print(listx)
    return listx

def rootRMS(histpylist): #takes histos in list
    if(len(histpylist))==1: print('ERROR: not enough hists'); return histpylist[0]
    hrmslist = histpylist[0].Clone('hrmslist')
    for binx in range(1,hrmslist.GetNbinsX()+1):
        bincontentx2 = 0
        hist_rejects = 0
        for i in range(1,len(histpylist)):
            diff_bin = np.sqrt((histpylist[i].GetBinContent(binx)-histpylist[0].GetBinContent(binx))**2)
            if (diff_bin>0.4): hist_rejects += 1; bincontentx2 += 0
            else: bincontentx2 += (histpylist[i].GetBinContent(binx)-histpylist[0].GetBinContent(binx))**2
        bincontentx = np.sqrt((bincontentx2)/(len(histpylist)-1-hist_rejects))
        hrmslist.SetBinContent(binx,bincontentx)
        hrmslist.SetBinError(binx,0)
            
    return hrmslist

def rootMEAN(histpylist): #takes histos in list
    if(len(histpylist))==1: print('ERROR: not enough hists'); return histpylist[0]
    hmeanlist = histpylist[0].Clone('hmeanlist')
    for binx in range(hmeanlist.GetNbinsX()+1):
        bincontentxsum = 0
        for i in range(1,len(histpylist)):
            bincontentxsum += (histpylist[i].GetBinContent(binx)-histpylist[0].GetBinContent(binx))
        bincontentx = ((bincontentxsum)/(len(histpylist)-1))
        hmeanlist.SetBinContent(binx,abs(bincontentx))
        hmeanlist.SetBinError(binx,0)
            
    return hmeanlist

def yieldPtrials(histpylist,binpt):
    hYPT = ROOT.TH1D('hYPT'+str(int(binpt*10)),'hYPT'+str(int(binpt*10)),len(histpylist),0,len(histpylist))
    for i in range(len(histpylist)):
        yieldbin = histpylist[i].GetBinContent(histpylist[i].GetXaxis().FindBin(binpt))
        errorbin = histpylist[i].GetBinError(histpylist[i].GetXaxis().FindBin(binpt))
        hYPT.SetBinContent(i+1,yieldbin);hYPT.SetBinError(i+1,errorbin)
        hYPT.GetYaxis().SetMaxDigits(3);hYPT.GetYaxis().SetRangeUser(0,yieldbin+errorbin*10)
        hYPT.GetYaxis().SetTitleOffset(1.4);
        hYPT.GetYaxis().SetTitle('Raw Yield')
        hYPT.GetXaxis().SetTitle('Trial #')

    return hYPT

def yieldDist(histpylist,binpt):
    yieldbin0 = int(histpylist[0].GetBinContent(histpylist[0].GetXaxis().FindBin(binpt)))
    errorbin0 = int(histpylist[0].GetBinError(histpylist[0].GetXaxis().FindBin(binpt)))
    rangemax=yieldbin0*2
    hYDT = ROOT.TH1D('hYDT'+str(int(binpt*10)),'hYDT'+str(int(binpt*10)),100,0,yieldbin0+5*errorbin0)
    for i in range(len(histpylist)):
        yieldbin = histpylist[i].GetBinContent(histpylist[i].GetXaxis().FindBin(binpt))
        hYDT.Fill(yieldbin)
        hYDT.GetXaxis().SetRangeUser(0,yieldbin0+5*errorbin0)

    return hYDT


def rootRMS_mt(histpylist,sigUp, sigDo): #takes histos in list
    if(len(histpylist))==1: print('ERROR: not enough hists'); return histpylist[0]
    hrmslist = histpylist[0].Clone('hrmslist')
    for binx in range(hrmslist.GetNbinsX()+1):
        bincontentx2 = 0
        outliers = 0
        total = 0
        for i in range(1,len(histpylist)):
            binval = histpylist[i].GetBinContent(binx)
            errval = histpylist[i].GetBinError(binx)
            binval0= histpylist[0].GetBinContent(binx)
            errval0= histpylist[0].GetBinError(binx)

            #if (binval>0 and binval>binval0-0.5*errval0 and binval<binval0+0.5*errval0):
            if (binval>0 and binval>binval0-sigDo*errval0 and binval<binval0+sigUp*errval0):
                bincontentx2 += (histpylist[i].GetBinContent(binx)-histpylist[0].GetBinContent(binx))**2
                total += 1
            else:
                outliers+=1
        try:
            #bincontentx = np.sqrt((bincontentx2)/(len(histpylist)-1-outliers))
            bincontentx = np.sqrt((bincontentx2)/(total))
            hrmslist.SetBinContent(binx,bincontentx)
            hrmslist.SetBinError(binx,0)
        except:
            print('outliers=',outliers)
            
    return hrmslist


def funcLine(x, a, b):
    return a*x + b

def funcBola(x, a, b, c):
    return a + np.sqrt(b*(x-c))
def funcBola2(x, a, b, c):
    return a - np.sqrt(b*(x-c))
