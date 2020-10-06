import numpy as np
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
    for binx in range(hrmslist.GetNbinsX()+1):
        bincontentx2 = 0
        for i in range(1,len(histpylist)):
            bincontentx2 += (histpylist[i].GetBinContent(binx)-histpylist[0].GetBinContent(binx))**2
        bincontentx = np.sqrt((bincontentx2)/(len(histpylist)-1))
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
        hYPT.SetBinContent(i+1,yieldbin)
        hYPT.SetBinError(i+1,errorbin)
        hYPT.GetYaxis().SetMaxDigits(3)
        hYPT.GetYaxis().SetRangeUser(yieldbin-errorbin*10,yieldbin+errorbin*5)
        hYPT.GetYaxis().SetTitleOffset(1.4)
        hYPT.GetYaxis().SetTitle('dN/dp_{T}   ')
        hYPT.GetXaxis().SetTitle('Trial number')

    #hYPT.GetYaxis().SetRangeUser(0.75*yieldbin,1.25*yieldbin)

    return hYPT
