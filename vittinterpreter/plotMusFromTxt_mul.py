#!/usr/bin/env python2.7
###!/usr/bin/env python
# 
# --------------------------------------------
# Standard python import
from optparse import OptionParser, make_option
import fnmatch, glob, os, sys, json, itertools, array
import re
#sys.argv.append( '-b' )
from array import array
## ------------------------------------------------------------------------------------------------------------------------------------------------------


#from templates_maker import buildRooDataSet
import ROOT
from ROOT import TH2D, TH1D, TFile, TProfile, TCanvas, TGraphAsymmErrors
from ROOT import RooWorkspace
from ROOT import RooAbsData
from ROOT import RooDataSet
from ROOT import *

import os.path


from ROOT import gROOT
gROOT.ForceStyle()
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetHatchesSpacing(1)
gStyle.SetHatchesLineWidth(3)

def getListFromBuffer(buf, size):
    buf.SetSize(size)
    return list(buf)

def setYerrorsToZero(g):
    newg = TGraphAsymmErrors(g.GetN(), g.GetX(), g.GetY(), g.GetEXlow(), g.GetEXhigh(), array('d', [0.]*g.GetN() ), array('d', [0.]*g.GetN() ) )
    return newg
def setXerrorsToZero(g):
    newg = TGraphAsymmErrors(g.GetN(), g.GetX(), g.GetY(), array('d', [0.]*g.GetN() ), array('d', [0.]*g.GetN() ), g.GetEYlow(), g.GetEYhigh()  )
    return newg

def convertToHist(g):
    nbins = g.GetN()
    x = getListFromBuffer(g.GetX(),g.GetN())
    edx = getListFromBuffer(g.GetEXlow(),g.GetN())
    eux = getListFromBuffer(g.GetEXhigh(),g.GetN())    
    bins = map(lambda ix,iex: ix - iex, x,edx)
    bins.append(x[-1]+eux[-1])
    hist = TH1D(str(g.GetName())+"_hist", g.GetTitle(), nbins, array('d',bins))
    y = getListFromBuffer(g.GetY(),g.GetN())
    for ix, iy in zip(x,y):
        hist.Fill(ix,iy)
    hist.SetLineColor( g.GetLineColor() )    
    hist.SetLineWidth( g.GetLineWidth() )    
    hist.SetFillColor( g.GetFillColor() )    
    hist.SetFillStyle( g.GetFillStyle() )
    return hist

def AddGraphs(g1,g2):
    x1 = getListFromBuffer(g1.GetX(),g1.GetN())
    x2 = getListFromBuffer(g2.GetX(),g2.GetN())
    if not x1 == x2:
        print 'ERRROR: cannot add TGraphs with different x-coordinates!'
        return 0
    edx1 = getListFromBuffer(g1.GetEXlow(),g1.GetN())
    edx2 = getListFromBuffer(g2.GetEXlow(),g2.GetN())
    eux1 = getListFromBuffer(g1.GetEXhigh(),g1.GetN())
    eux2 = getListFromBuffer(g2.GetEXhigh(),g2.GetN())

    y1 = getListFromBuffer(g1.GetY(),g1.GetN())
    y2 = getListFromBuffer(g2.GetY(),g2.GetN())
    edy1 = getListFromBuffer(g1.GetEYlow(),g1.GetN())
    edy2 = getListFromBuffer(g2.GetEYlow(),g2.GetN())
    euy1 = getListFromBuffer(g1.GetEYhigh(),g1.GetN())
    euy2 = getListFromBuffer(g2.GetEYhigh(),g2.GetN())

    yS   = map(lambda a,b: a+b, y1,y2)
    edyS = map(lambda a,b: np.sqrt(a**2+b**2), edy1,edy2)
    euyS = map(lambda a,b: np.sqrt(a**2+b**2), euy1,euy2)
    xS   = x1
    edxS = edx1
    euxS = edx2
    gsum = TGraphAsymmErrors(g1.GetN(), array('d',xS), array('d',yS), array('d',edxS), array('d',euxS), array('d',edyS), array('d',euyS))
    return gsum


ROOT.myColorA0   = ROOT.TColor.GetColor("#ff8000")
ROOT.myColorA1   = ROOT.TColor.GetColor("#ffbf80")

ROOT.myColorB1   = ROOT.TColor.GetColor("#538cc6")

import numpy as np



observables={}
observables["Pt"]                           = dict( xlabel="p_{T}^{#gamma#gamma} (GeV)"                      , ylabel=    "d#sigma_{fid}/dp_{T}^{#gamma#gamma} (fb/GeV)"                      , rescaleLastToSecondLast=True  , name="p_{T}^{#gamma#gamma}",                            unit=" GeV",                                       legend="r", yspeed=1.0, extratext1='', extratext2='')
observables["AbsRapidity"]                  = dict( xlabel="|y_{#gamma#gamma}|"                              , ylabel=    "d#sigma_{fid}/d|y_{#gamma#gamma}| (fb)"                            , rescaleLastToSecondLast=False  , name="|y_{#gamma#gamma}|",                              unit="",                                          legend="r", yspeed=1.0, extratext1='', extratext2='')
observables["CosThetaStar"]                 = dict( xlabel="|cos(#theta*)|"                                  , ylabel=    "d#sigma_{fid}/d|cos(#theta*)| (fb)"                                , rescaleLastToSecondLast=False  , name="|cos(#theta*)|",                                  unit="",                                          legend="r", yspeed=1.0, extratext1='', extratext2='')
                                                                                                                                                                                              
observables["Jet2p5Pt0"]                    = dict( xlabel="p_{T}^{j_{1}} (GeV)"                             , ylabel=    "d#sigma_{fid}/dp_{T}^{j_{1}} (fb/GeV)"                             , rescaleLastToSecondLast=True  , name="p_{T}^{j_{1}}",                                   unit=" GeV",                                       legend="r", yspeed=1.0, extratext1='p_{T}^{j_{1}} > 30 GeV, |#eta^{j_{1}}| < 2.5', extratext2='')
observables["Jet2p5AbsRapidity0"]           = dict( xlabel="|y^{j_{1}}|"                                     , ylabel=    "d#sigma_{fid}/d|y^{j_{1}}| (fb)"                                   , rescaleLastToSecondLast=False  , name="|y^{j_{1}}|",                                     unit="",                                          legend="r", yspeed=1.0, extratext1='p_{T}^{j_{1}} > 30 GeV, |#eta^{j_{1}}| < 2.5', extratext2='')
observables["AbsDeltaPhiGgJet0Eta2p5"]      = dict( xlabel="|#Delta #phi^{#gamma#gamma,j_{1}}|"              , ylabel=    "d#sigma_{fid}/d|#Delta #phi^{#gamma#gamma,j_{1}}| (fb)"            , rescaleLastToSecondLast=False  , name="|#Delta #phi^{#gamma#gamma,j_{1}}|",              unit="",                                          legend="r", yspeed=1.5, extratext1='p_{T}^{j_{1}} > 30 GeV, |#eta^{j_{1}}| < 2.5', extratext2='')
observables["AbsDeltaRapidityGgJet0Eta2p5"] = dict( xlabel="|#Delta y^{#gamma#gamma,j_{1}}|"                 , ylabel=    "d#sigma_{fid}/d|#Delta y^{#gamma#gamma,j_{1}}| (fb)"               , rescaleLastToSecondLast=True  , name="|#Delta y^{#gamma#gamma,j_{1}}|",                 unit="",                                          legend="r", yspeed=1.0, extratext1='p_{T}^{j_{1}} > 30 GeV, |#eta^{j_{1}}| < 2.5', extratext2='')
                                                                                                                                                                                              
observables["PtNjets2p5"]                   = dict( xlabel="p_{T}^{#gamma#gamma} x n_{jet}"                  , ylabel=    "d#sigma_{fid}/d^{2}p_{T}^{#gamma#gamma} n_{jet} (fb/GeV)"          , rescaleLastToSecondLast=True  , name="p_{T}^{#gamma#gamma}",                            unit=" GeV",                                       legend="r", yspeed=1.0, extratext1='p_{T}^{jet} > 30 GeV, |#eta^{jet}| < 2.5', extratext2='')
observables["PtNjets2p5_0"]                 = dict( xlabel="p_{T}^{#gamma#gamma}, n_{jet}=0 (GeV)"           , ylabel=    "d#sigma_{fid}/dp_{T}^{#gamma#gamma}, n_{jet}=0 (fb/GeV)"           , rescaleLastToSecondLast=True  , name="p_{T}^{#gamma#gamma}",                            unit=" GeV",                                       legend="r", yspeed=1.0, extratext1='p_{T}^{jet} > 30 GeV, |#eta^{jet}| < 2.5', extratext2='')
observables["PtNjets2p5_1"]                 = dict( xlabel="p_{T}^{#gamma#gamma}, n_{jet}=1 (GeV)"           , ylabel=    "d#sigma_{fid}/dp_{T}^{#gamma#gamma}, n_{jet}=1 (fb/Gev)"           , rescaleLastToSecondLast=True  , name="p_{T}^{#gamma#gamma}",                            unit=" GeV",                                       legend="r", yspeed=1.0, extratext1='p_{T}^{jet} > 30 GeV, |#eta^{jet}| < 2.5', extratext2='')
observables["PtNjets2p5_1plus"]             = dict( xlabel="p_{T}^{#gamma#gamma}, n_{jet}>1 (GeV)"           , ylabel=    "d#sigma_{fid}/dp_{T}^{#gamma#gamma}, n_{jet}>1 (fb/GeV)"           , rescaleLastToSecondLast=True  , name="p_{T}^{#gamma#gamma}",                            unit=" GeV",                                       legend="r", yspeed=1.0, extratext1='p_{T}^{jet} > 30 GeV, |#eta^{jet}| < 2.5', extratext2='')
                                                                                                                                                                                              
observables["AbsDeltaEtaJJEta4p7"]          = dict( xlabel="|#Delta #eta^{j_{1},j_{2}}|"                     , ylabel=    "d#sigma_{fid}/d|#Delta #eta^{j_{1},j_{2}}| (fb)"                   , rescaleLastToSecondLast=True  , name="|#Delta #eta^{j_{1},j_{2}}|",                     unit="",                                          legend="r", yspeed=1.5, extratext1='p_{T}^{j_{1,2}} > 30 GeV, |#eta^{j_{1,2}}| < 4.7', extratext2='')
observables["AbsDeltaPhiGgJjEta4p7"]        = dict( xlabel="|#Delta #phi^{#gamma#gamma,j_{1}j_{2}}|"         , ylabel=    "d#sigma_{fid}/d|#Delta #phi^{#gamma#gamma,j_{1}j_{2}}| (fb)"       , rescaleLastToSecondLast=False  , name="|#Delta #phi^{#gamma#gamma,j_{1}j_{2}}|",         unit="",                                          legend="r", yspeed=1.0, extratext1='p_{T}^{j_{1,2}} > 30 GeV, |#eta^{j_{1,2}}| < 4.7', extratext2='')
observables["AbsDeltaPhiJjEta4p7"]          = dict( xlabel="|#Delta #phi^{j_{1},j_{2}}|"                     , ylabel=    "d#sigma_{fid}/d|#Delta #phi^{j_{1},j_{2}}| (fb)"                   , rescaleLastToSecondLast=False  , name="|#Delta #phi^{j_{1},j_{2}}|",                     unit="",                                          legend="r", yspeed=2.5, extratext1='p_{T}^{j_{1,2}} > 30 GeV, |#eta^{j_{1,2}}| < 4.7', extratext2='')
observables["Jet4p7Pt1"]                    = dict( xlabel="p_{T}^{j_{2}} (GeV)"                             , ylabel=    "d#sigma_{fid}/dp_{T}^{j_{2}} (fb/GeV)"                             , rescaleLastToSecondLast=True  , name="p_{T}^{j_{2}}",                                   unit=" GeV",                                       legend="r", yspeed=2.5, extratext1='p_{T}^{j_{1,2}} > 30 GeV, |#eta^{j_{1,2}}| < 4.7', extratext2='')
observables["Jet4p7AbsRapidity1"]           = dict( xlabel="|y^{j_{2}}|"                                     , ylabel=    "d#sigma_{fid}/d|y^{j_{2}}| (fb)"                                   , rescaleLastToSecondLast=False  , name="|y^{j_{2}}|",                                     unit="",                                          legend="r", yspeed=1.0, extratext1='p_{T}^{j_{1,2}} > 30 GeV, |#eta^{j_{1,2}}| < 4.7', extratext2='')
observables["ZeppenfeldEta4p7"]             = dict( xlabel="|#bar{#eta}^{j_{1}j_{2}} - #eta^{#gamma#gamma}|" , ylabel=    "d#sigma_{fid}/d|#bar{#eta}^{j_{1}j_{2}} - #eta^{#gamma#gamma}| (fb)", rescaleLastToSecondLast=True , name="|#bar{#eta}^{j_{1}j_{2}} - #eta^{#gamma#gamma}|", unit="",                                          legend="r", yspeed=1.0, extratext1='p_{T}^{j_{1,2}} > 30 GeV, |#eta^{j_{1,2}}| < 4.7', extratext2='')
observables["MjjEta4p7"]                    = dict( xlabel="m^{j_{1}j_{2}} (GeV)"                            , ylabel=    "d#sigma_{fid}/dm^{j_{1}j_{2}} (fb/GeV)"                            , rescaleLastToSecondLast=True  , name="m^{j_{1}j_{2}}",                                  unit=" GeV",                                       legend="r", yspeed=2.0, extratext1='p_{T}^{j_{1,2}} > 30 GeV, |#eta^{j_{1,2}}| < 4.7', extratext2='')
                                                                                                                                                                                              
                                                                                                                                                                                              
observables["Njets2p5"]                     = dict( xlabel="n_{jet}"                                         , ylabel=    "d#sigma_{fid}/dn_{jet} (fb)"                                       , rescaleLastToSecondLast=False  , name="n_{jet}",                                         unit="",                                          legend="r", yspeed=1.0, extratext1='p_{T}^{jet} > 30 GeV, |#eta^{jet}| < 2.5', extratext2='', customLabels=['0','1','2','3','>3'])
observables["MET"]                          = dict( xlabel="p_{T}^{miss} (GeV)"                              , ylabel=    "d#sigma_{fid}/dp_{T}^{miss} (fb/GeV)"                              , rescaleLastToSecondLast=True  , name="p_{T}^{miss}",                                    unit=" GeV",                                       legend="r", yspeed=1.0, extratext1='', extratext2='')
observables["NjetsBflavorTight2p5"]         = dict( xlabel="n_{jet}^{b}"                                     , ylabel=    "d#sigma_{fid}/dn_{jet}^{b} (fb)"                                   , rescaleLastToSecondLast=False  , name="n_{jet}^{b}",                                     unit="",                                          legend="r", yspeed=1.0, extratext1='p_{T}^{jet} > 30 GeV, |#eta^{jet}| < 2.5', extratext2='', customLabels=['0','1','>1'])
observables["Nleptons"]                     = dict( xlabel="n_{lepton}"                                      , ylabel=    "d#sigma_{fid}/dn_{lepton} (fb)"                                    , rescaleLastToSecondLast=False  , name="n_{lepton}",                                      unit="",                                          legend="r", yspeed=1.0, extratext1='', extratext2='', customLabels=['0','1','>1'])
observables["1Lepton1Bjet"]                 = dict( xlabel="n_{lepton} = 1, N_{b} = 1"                       , ylabel=    "#sigma_{fid}(n_{lepton} = 1, N_{b} = 1 (fb)"                       , rescaleLastToSecondLast=True  , name="n_{lepton} = 1, N_{b} = 1",                       unit="",                                          legend="r", yspeed=1.0, extratext1='', extratext2='')
observables["1LeptonHighMET"]               = dict( xlabel="n_{lepton} = 1, high p_{T}^{miss}"               , ylabel=    "#sigma_{fid}(n_{lepton} = 1, high p_{T}^{miss} (fb)"               , rescaleLastToSecondLast=True  , name="n_{lepton} = 1, high p_{T}^{miss}",               unit="",                                          legend="r", yspeed=1.0, extratext1='', extratext2='')
observables["1LeptonLowMET"]                = dict( xlabel="n_{lepton} = 1, low p_{T}^{miss}"                , ylabel=    "#sigma_{fid}(n_{lepton} = 1, low p_{T}^{miss} (fb)"                , rescaleLastToSecondLast=True  , name="n_{lepton} = 1, low p_{T}^{miss}",                unit="",                                          legend="r", yspeed=1.0, extratext1='', extratext2='')
                                                                                                                                                                                              
                                                                                                                                                                                              
observables["Jet4p7Pt1VBFlike"]             = dict( xlabel="p_{T}^{j_{2}} (GeV)"                             , ylabel=    "d#sigma_{fid}/dp_{T}^{j_{2}} (fb/GeV)"                             , rescaleLastToSecondLast=True  , name="p_{T}^{j_{2}}",                                   unit=" GeV",                                       legend="r", yspeed=2.5, extratext1='p_{T}^{j_{1,2}} > 30 GeV, |#eta^{j_{1,2}}| < 4.7', extratext2='m^{j_{1}j_{2}} > 200 GeV, |#Delta#eta^{j_{1}j_{2}}| > 3.5')
observables["AbsDeltaPhiGgJjEta4p7VBFlike"] = dict( xlabel="|#Delta #phi^{#gamma#gamma,j_{1}j_{2}}|"         , ylabel=    "d#sigma_{fid}/d|#Delta #phi^{#gamma#gamma,j_{1}j_{2}}| (fb)"       , rescaleLastToSecondLast=False  , name="|#Delta #phi^{#gamma#gamma,j_{1}j_{2}}|",         unit="",                                          legend="r", yspeed=1.0, extratext1='p_{T}^{j_{1,2}} > 30 GeV, |#eta^{j_{1,2}}| < 4.7', extratext2='m^{j_{1}j_{2}} > 200 GeV, |#Delta#eta^{j_{1}j_{2}}| > 3.5')
observables["AbsDeltaPhiJjEta4p7VBFlike"]   = dict( xlabel="|#Delta #phi^{j_{1},j_{2}}|"                     , ylabel=    "d#sigma_{fid}/d|#Delta #phi^{j_{1},j_{2}}| (fb)"                   , rescaleLastToSecondLast=False  , name="|#Delta #phi^{j_{1},j_{2}}|",                     unit="",                                          legend="r", yspeed=2.5, extratext1='p_{T}^{j_{1,2}} > 30 GeV, |#eta^{j_{1,2}}| < 4.7', extratext2='m^{j_{1}j_{2}} > 200 GeV, |#Delta#eta^{j_{1}j_{2}}| > 3.5')









def getBinBoundariesFromDataset(dname):
    nameSplit = dname.split("_")
    genbinl=0.
    genbinh=0.
    recobinl=0.
    recobinh=0.
    for ip in range(len(nameSplit)):
        if "gen" in nameSplit[ip] or "Gen" in nameSplit[ip]:
            if(len(nameSplit)>ip+1):
                genbinl = float(nameSplit[ip+1].replace("m","-").replace("p","."))
                genbinh = float(nameSplit[ip+2].replace("m","-").replace("p","."))
        if "reco" in nameSplit[ip] or "Reco" in nameSplit[ip]:
            if(len(nameSplit)>ip+1):
                recobinl = float(nameSplit[ip+1].replace("m","-").replace("p","."))
                recobinh = float(nameSplit[ip+2].replace("m","-").replace("p","."))
    return [genbinl,genbinh],[recobinl,recobinh]

def getBinBoundariesFromProcess(dname):
    nameSplit = dname.split("_")
    genbinl=0.
    genbinh=0.
    for ip in range(len(nameSplit)):
        if "gen" in nameSplit[ip] or "Gen" in nameSplit[ip]:
            if(len(nameSplit)>ip+1):
                genbinl = float(nameSplit[ip+1].replace("m","-").replace("p","."))
                genbinh = float(nameSplit[ip+2].replace("m","-").replace("p","."))
    return [genbinl,genbinh]

def getVarsName(dname):
    nameSplit = dname.split("_")
    genVar=""
    recoVar=""
    for ip in range(len(nameSplit)):
        if "gen" in nameSplit[ip] or "Gen" in nameSplit[ip]:
            genVar=nameSplit[ip]
        if "reco" in nameSplit[ip] or "Reco" in nameSplit[ip]:
            recoVar=nameSplit[ip]
    return genVar,recoVar
    
def mapPOItoProcess(line):
    POItoProc={}
    for s in filter(None,line.strip(' ').split("--PO")):
        print s.strip()
        print s.strip().split(":")
        print s.strip().split(":")[0].split("/")
        print s.strip().split(":")[1].split("[")
        POItoProc[s.strip().split(":")[1].split("[")[0]] = getBinBoundariesFromProcess( s.strip().split(":")[0].split("/")[1] )
    return POItoProc

def getBestFit(POIs, lines):
    BF={}
    for POI in POIs:
        for line in lines:
            if POI in line:
                print "".join(line.split())
#                BF[POI]= re.split('[\+-]', (("".join(line.split())).split(POI+":")[1].strip()))
                BF[POI]= line.split()[-3:]
                print 'getBF debug'
                print BF[POI]
                BF[POI][2] = abs(float(BF[POI][2]))
    for POI in POIs:
        if POI not in BF.keys():
            BF[POI]=[0.0,1.0,1.0]
    return BF
    
    
            
savefmts=['.png','.root','.pdf','.jpg']
# Main routine
def main(o,args):

    print options.files
    with open(options.files) as f:
        content = f.readlines()
    POItoProc = mapPOItoProcess(content[-1])
    print "POItoProc"
    print POItoProc

    BF = getBestFit(POItoProc.keys(),content[:-1])
    print "BF"
    print BF

    print options.filesFreezeNuis
    with open(options.filesFreezeNuis) as f:
        content = f.readlines()
    POItoProcStatOnly = mapPOItoProcess(content[-1])
    print "POItoProcStatOnly"
    print POItoProcStatOnly
    BF_StatOnly = getBestFit(POItoProcStatOnly.keys(),content[:-1])
    print "BF_StatOnly"
    print BF_StatOnly


    if not set(POItoProc.keys()) == set(POItoProcStatOnly.keys()):
        raise ValueError("input files for full unc and stat. only scans do not have the same POIs definitions and can not be matched!")

    central={}
    predictions={}
    predictions["SM"]=dict(vals={}, unc={},label='NNLOPS', shift=-2, color=myColorA0, fill=3235, leg='SM', marker=20, uncXS=0.049) #not drawn
    predictions["ggHamcatnloNNLOPS"]=dict(vals={}, unc={},label='NNLOPS_ggH_amcatnlo', shift=-2, color=629, fill=3265, leg='HX + ggH aMC@NLO, #scale[0.9]{NNLOPS}', marker=26, uncXS=0.049)  #kRed-3
    predictions["ggHamcatnlo"]=dict(vals={}, unc={},label='_ggH_amcatnlo', shift=4, color=410, fill=3256, leg='HX + ggH aMC@NLO', marker=27, uncXS=0.049) #kGreen-6
    predictions["ggHpowheg"]=dict(vals={}, unc={},label='_ggH_powheg', shift=2, color=4, fill=3253, leg='HX + ggH POWHEG', marker=32, uncXS=0.049)
#    predictions["ggHpowhegNNLOPS"]=dict(vals={}, unc={},label='NNLOPS_ggH_powheg', shift=0, color=6)
    predictions["HX"]=dict(vals={}, unc={},label='NNLOPS_HX_amcatnlo', shift=0, color=801, fill=3235, leg='HX = VBF + VH + ttH', marker=4, uncXS=0.0260) #kOrange+1##myColorA0
    centralSM={}
    ggHamcatnloNNLOPS={}
    ggHamcatnlo={}
    ggHpowheg={}
    HX={}
    up={}
    down={}
    up_StatOnly={}
    down_StatOnly={}
    xerr={}
    data=[]
    for POI in POItoProc.keys():
        data.append( {'x': 0.5*(POItoProc[POI][0]+POItoProc[POI][1]), 'y': float(BF[POI][0]), 'ySM' : float(1.), 'errx': float( 0.5*(POItoProc[POI][0]+POItoProc[POI][1]) - POItoProc[POI][0]) , 'erryup' : float(BF[POI][1]), 'errydown': float(BF[POI][2]),  'erryup_statOnly' : float(BF_StatOnly[POI][1]), 'errydown_statOnly': float(BF_StatOnly[POI][2])    })
        central[0.5*(POItoProc[POI][0]+POItoProc[POI][1])] = float(BF[POI][0])
        for pr in predictions.keys():
            predictions[pr]['vals'][0.5*(POItoProc[POI][0]+POItoProc[POI][1])] = 1.
#        centralSM[0.5*(POItoProc[POI][0]+POItoProc[POI][1])] = 1.
#        ggHamcatnloNNLOPS[0.5*(POItoProc[POI][0]+POItoProc[POI][1])] = 1.
#        ggHamcatnlo[0.5*(POItoProc[POI][0]+POItoProc[POI][1])] = 1.
#        ggHpowheg[0.5*(POItoProc[POI][0]+POItoProc[POI][1])] = 1.
#        HX[0.5*(POItoProc[POI][0]+POItoProc[POI][1])] = 1.
        xerr[ 0.5*(POItoProc[POI][0]+POItoProc[POI][1]) ] = float( 0.5*(POItoProc[POI][0]+POItoProc[POI][1]) - POItoProc[POI][0]) 
        up[0.5*(POItoProc[POI][0]+POItoProc[POI][1])] = float(BF[POI][1])
        down[0.5*(POItoProc[POI][0]+POItoProc[POI][1])] = float(BF[POI][2])

        up_StatOnly[0.5*(POItoProcStatOnly[POI][0]+POItoProcStatOnly[POI][1])] = float(BF_StatOnly[POI][1])
        down_StatOnly[0.5*(POItoProcStatOnly[POI][0]+POItoProcStatOnly[POI][1])] = float(BF_StatOnly[POI][2])

        #binBound.append((POItoProc[POI][0]+POItoProc[POI][1]))
        #binBound.append(POItoProc[POI][1])
        
    #    binBound = list(sorted(set(binBound)))
    sortedData = sorted(data, key=lambda k: k['x'])
    if options.resizeFirst != -1000:
        sortedData[0]['x'] = 0.5*(options.resizeFirst + sortedData[1]['x'] - sortedData[1]['errx'])
        sortedData[0]['errx'] = sortedData[1]['x'] - sortedData[1]['errx'] - sortedData[0]['x'] 

    if options.resizeLast != -1:
        sortedData[-1]['x'] = 0.5*(options.resizeLast + sortedData[-2]['x'] + sortedData[-2]['errx'])
        sortedData[-1]['errx'] =   sortedData[-1]['x'] - (sortedData[-2]['x'] + sortedData[-2]['errx'])

    print sortedData
    print central
    for pr in predictions.keys():
        filename = '/mnt/t3nfs01/data01/shome/vtavolar/jupyter/CMSSW_8_0_28/src/higgs_model_dep/spectrum%s_%s.npz'% (predictions[pr]['label'] ,options.variable)
        if pr in ['SM', 'ggHamcatnloNNLOPS', 'ggHamcatnlo', 'ggHpowheg']:
            filenameunc = '/mnt/t3nfs01/data01/shome/vtavolar/jupyter/CMSSW_8_0_28/src/higgs_model_dep/uncertainty_ggH_amcatnlo_%s.npz'% (options.variable)
            uncXs=0.049
        elif pr in ['HX']:
            filenameunc = '/mnt/t3nfs01/data01/shome/vtavolar/jupyter/CMSSW_8_0_28/src/higgs_model_dep/uncertainty_HX_amcatnlo_%s.npz'% (options.variable)
            uncXs=0.0260
        print filename
        if os.path.isfile(filename) and options.spectrum:
            spectrum = np.load(filename)
            print spectrum.keys()
            print spectrum
            xsecs = spectrum['spectrum'] 
            print xsecs
            binwidth = spectrum['binwidth'] 
            binBoundaries = spectrum['obs_bins'] 
            if options.resizeLast != -1:
                binBoundaries[-1]=options.resizeLast
                binwidth[-1]=binBoundaries[-1]-binBoundaries[-2]
            #divide last bin for same binwidth as second-to-last
##            if options.variable != 'CosThetaStar':
            if options.variable in observables.keys():
                if observables[options.variable]['rescaleLastToSecondLast']:
                    binwidth[-1]=binwidth[-2]
            print binwidth
            xsecOverBw = map(lambda  s, b: s/b,  xsecs, binwidth)
            unc=[0.2]*len(xsecOverBw)
            if os.path.isfile(filenameunc):
                uncSpectrum = np.load(filenameunc)
                unc = uncSpectrum['uncetainty']
            ic=0
            for c in sorted(central.keys()):
                print c
                print central[c]
                predictions[pr]['vals'][c] = predictions[pr]['vals'][c]*xsecOverBw[ic]
                predictions[pr]['unc'][c] = np.sqrt(unc[ic]**2 + uncXs**2)
                ic+=1
            print pr
            print predictions[pr]['label']
            print predictions[pr]['vals']
            print predictions[pr]['unc']

                

    for c in sorted(central.keys()):
            central[c] = central[c]*predictions['SM']['vals'][c] 

    for pr in predictions.keys():
        print predictions[pr]['label']
        print predictions[pr]['vals']

    
    print sortedData
    sortedCentral = sorted(central.keys())
    print 'sortedCentral'
    print sortedCentral
    if options.resizeFirst != -1000:
        theKey = (0.5*(options.resizeFirst + sortedData[1]['x'] - sortedData[1]['errx']))
        central[theKey] = central[sortedCentral[0]]
        del central[sortedCentral[0]]
        for pr in predictions.keys():
            predictions[pr]['vals'][theKey] = predictions[pr]['vals'][sortedCentral[0]]
            del predictions[pr]['vals'][sortedCentral[0]]
            predictions[pr]['unc'][theKey] = predictions[pr]['unc'][sortedCentral[0]]
            del predictions[pr]['unc'][sortedCentral[0]]

    if options.resizeLast != -1:
        theKey = (0.5*(options.resizeLast + sortedData[-2]['x'] + sortedData[-2]['errx']))
        central[theKey] = central[sortedCentral[-1]]
        del central[sortedCentral[-1]]
        for pr in predictions.keys():
            predictions[pr]['vals'][theKey] = predictions[pr]['vals'][sortedCentral[-1]]
            del predictions[pr]['vals'][sortedCentral[-1]]
            predictions[pr]['unc'][theKey] = predictions[pr]['unc'][sortedCentral[-1]]
            del predictions[pr]['unc'][sortedCentral[-1]]

    print "central"
    print central

    for pr in predictions.keys():
        print pr
        print predictions[pr]['label']
        print predictions[pr]['vals']
        print predictions[pr]['unc']
    
    xsecs=[]
    for dt in sortedData:
        dt['y']=central[dt['x']]
        for pr in predictions.keys():
            dt['y%s'%pr] = predictions[pr]['vals'][dt['x']]
            dt['y%sUnc'%pr] = predictions[pr]['unc'][dt['x']]
        dt['erryup'] = dt['erryup']*predictions['SM']['vals'][dt['x']]
        dt['errydown'] = dt['errydown']*predictions['SM']['vals'][dt['x']]
        dt['erryup_statOnly'] = dt['erryup_statOnly']*predictions['SM']['vals'][dt['x']]
        dt['errydown_statOnly'] = dt['errydown_statOnly']*predictions['SM']['vals'][dt['x']]
        xsecs.append(dt['y'])

    print "dt['y']"
    print dt['y']
    print 'xsecs'
    print xsecs

    print "scale up last bin to avg xsec if needed"
    lastBinXsec = sortedData[-1]['y']
    avgXsec = sum(xsecs[:-1])/float(len(xsecs))
    scaleLastBin=1.
#    if lastBinXsec < avgXsec:
#        scaleLastBin = avgXsec/lastBinXsec
    sortedData[-1]['y'] = sortedData[-1]['y']*scaleLastBin
    for pr in predictions.keys():
        sortedData[-1]['y%s'%pr] = sortedData[-1]['y%s'%pr]*scaleLastBin
##    sortedData[-1]['ySM'] = sortedData[-1]['ySM']*scaleLastBin
    sortedData[-1]['erryup'] = sortedData[-1]['erryup']*scaleLastBin
    sortedData[-1]['errydown'] = sortedData[-1]['errydown']*scaleLastBin
    sortedData[-1]['erryup_statOnly'] = sortedData[-1]['erryup_statOnly']*scaleLastBin
    sortedData[-1]['errydown_statOnly'] = sortedData[-1]['errydown_statOnly']*scaleLastBin

    xsecs[-1] = xsecs[-1]*scaleLastBin
    print 'sortedData[-1][y]'
    print sortedData[-1]['y']
    print 'xsecs[-1]'
    print xsecs[-1]
    print "search min and max for plotting"
    xsecMin = min(xsecs)
    xsecMax = max(xsecs)
    if options.hideFirstBin:
        xsecMin = min(xsecs[1:])
        xsecMax = max(xsecs[1:])
    print 'xsecMin'
    print xsecMin
    print 'xsecMax'
    print xsecMax
    print up
    print down
    print xerr
    print central.keys()
    print central.values()
##    graph = TGraphAsymmErrors( len(central.keys()), array('d', central.keys()), array('d', central.values()), array('d', xerr.values() ),array('d', xerr.values()), array('d',down.values()), array('d',up.values()) )
    newSortedData={}
    if options.variable=="PtNjets2p5":
        newSortedData["PtNjets2p5_0"]=sortedData[0:3]
        newSortedData["PtNjets2p5_0"][-1]['x'] = 65.0
        newSortedData["PtNjets2p5_0"][-1]['errx'] = 20.0
        adjustLastBin = (13000.0 - 45.0)/40.0
        newSortedData["PtNjets2p5_0"][-1]['y'] = newSortedData["PtNjets2p5_0"][-1]['y']*adjustLastBin
        for pr in predictions.keys():
            newSortedData["PtNjets2p5_0"][-1]['y%s'%pr] = newSortedData["PtNjets2p5_0"][-1]['y%s'%pr]*adjustLastBin
#        newSortedData["PtNjets2p5_0"][-1]['ySM'] = newSortedData["PtNjets2p5_0"][-1]['ySM']*adjustLastBin
        newSortedData["PtNjets2p5_0"][-1]['erryup'] = newSortedData["PtNjets2p5_0"][-1]['erryup']*adjustLastBin
        newSortedData["PtNjets2p5_0"][-1]['errydown'] = newSortedData["PtNjets2p5_0"][-1]['errydown']*adjustLastBin
        newSortedData["PtNjets2p5_0"][-1]['erryup_statOnly'] = newSortedData["PtNjets2p5_0"][-1]['erryup_statOnly']*adjustLastBin
        newSortedData["PtNjets2p5_0"][-1]['errydown_statOnly'] = newSortedData["PtNjets2p5_0"][-1]['errydown_statOnly']*adjustLastBin
        

        newSortedData["PtNjets2p5_1"]=sortedData[3:6]
        for t in  newSortedData["PtNjets2p5_1"]:
            t['x'] = t['x']-13000.0
        newSortedData["PtNjets2p5_1"][-1]['x'] = 150.0
        newSortedData["PtNjets2p5_1"][-1]['errx'] = 30.0
        adjustLastBin = (13000.0 - 120.0)/60.0
        newSortedData["PtNjets2p5_1"][-1]['y'] = newSortedData["PtNjets2p5_1"][-1]['y']*adjustLastBin
        for pr in predictions.keys():
            newSortedData["PtNjets2p5_1"][-1]['y%s'%pr] = newSortedData["PtNjets2p5_1"][-1]['y%s'%pr]*adjustLastBin
#        newSortedData["PtNjets2p5_1"][-1]['ySM'] = newSortedData["PtNjets2p5_1"][-1]['ySM']*adjustLastBin
        newSortedData["PtNjets2p5_1"][-1]['erryup'] = newSortedData["PtNjets2p5_1"][-1]['erryup']*adjustLastBin
        newSortedData["PtNjets2p5_1"][-1]['errydown'] = newSortedData["PtNjets2p5_1"][-1]['errydown']*adjustLastBin
        newSortedData["PtNjets2p5_1"][-1]['erryup_statOnly'] = newSortedData["PtNjets2p5_1"][-1]['erryup_statOnly']*adjustLastBin
        newSortedData["PtNjets2p5_1"][-1]['errydown_statOnly'] = newSortedData["PtNjets2p5_1"][-1]['errydown_statOnly']*adjustLastBin


        newSortedData["PtNjets2p5_1plus"]=sortedData[6:9]
        for t in  newSortedData["PtNjets2p5_1plus"]:
            t['x'] = t['x']-26000.0
        newSortedData["PtNjets2p5_1plus"][-1]['x'] = 400.0
        newSortedData["PtNjets2p5_1plus"][-1]['errx'] = 50.0
        print 'newSortedData'
        print newSortedData
    else:
        newSortedData[options.variable]=sortedData
    for obs in newSortedData.keys():
        sortedData=newSortedData[obs]
        if options.hideFirstBin:
            sortedData=sortedData[1:]
        if obs == "AbsDeltaPhiGgJjEta4p7VBFlike" or obs == "AbsDeltaPhiGgJjEta4p7":
            for pr in predictions.keys():
                sortedData[-1][pr]['y%s'%pr] = 0.
                sortedData[-1][pr]['y%sUc'%pr] = 0.
        print "newSortedData"
        print sortedData
        c1= TCanvas("c1","",2)
#        c1.SetBottomMargin(0.15)
#        c1.SetLeftMargin(0.15)
        pads=splitCanvas(c1)
        pads[0].cd()

        graph = TGraphAsymmErrors( len(sortedData), array('d', [a['x'] for a in sortedData]), array('d', [a['y'] for a in sortedData]), array('d', [a['errx'] for a in sortedData] ),array('d', [a['errx'] for a in sortedData]), array('d',[a['errydown'] for a in sortedData]), array('d', [a['erryup'] for a in sortedData]) )
        ErrDown = [a['errydown'] for a in sortedData]
        ErrUp = [a['erryup'] for a in sortedData]



  
        lastBinBoundary = sortedData[-1]['x'] - sortedData[-1]['errx']
        lastX = sortedData[-1]['x']
        lastMaxY = max(sortedData[-1]['y'] + 1.5* sortedData[-1]['erryup'],sortedData[-1]['ySM'] + 3*sortedData[-1]['ySM']*sortedData[-1]['ySMUnc'])
        print 'lastMaxY'
        print lastMaxY
        print sortedData[-1]['y']
        print 1.5* sortedData[-1]['erryup']
        print sortedData[-1]['ySM'] 
        print 3*sortedData[-1]['ySMUnc']
        if options.logy:
            lastMaxY=float(lastMaxY)*3
        else:
            lastMaxY=float(lastMaxY)*1.4
        print 'lastMaxY'
        print lastMaxY
        secondLastBinWidth = sortedData[-2]['errx']*2
        
        

        graph.SetName("graph_"+obs)
        graph.SetTitle("")
        graph.Print("all")

        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(kBlack)
        graph.SetMarkerSize(0.6)
        graph.SetLineColor(1)
        graph.SetFillColor(1)

        graph.GetXaxis().SetTitle("")
        graph.GetYaxis().SetTitle("d #sigma_{fid} / d x")
        if (obs in observables.keys()):
            if 'ylabel' in observables[obs].keys():
                print "variable "+str(obs)+" has a label for yaxis"
                print "setting ylabel to "+str(observables[obs]['ylabel'])
                graph.GetYaxis().SetTitle(observables[obs]['ylabel'])
        graph.GetXaxis().SetTitleFont(43)
        graph.GetYaxis().SetTitleFont(43)
        graph.GetXaxis().SetTitleSize(20)
        graph.GetYaxis().SetTitleSize(20)
        graph.GetYaxis().SetTitleOffset(1.5)
        graph.GetXaxis().SetTitleOffset(1.14)
        graph.GetYaxis().SetLabelFont(43)
        graph.GetXaxis().SetLabelFont(43)
        graph.GetYaxis().SetLabelSize(16)
        graph.GetXaxis().SetLabelSize(16)
        graph.GetYaxis().SetLabelOffset(0.007)
        graph.GetXaxis().SetLabelOffset(0.007)

        ##range user if not logy: find min HX, find max data and associated uncertainty --> range for actual plot
        ##scale up that range to create room for legends and texts
        minHX = min(a['yHX'] for a in sortedData)/2.
        maxData = max(a['y'] for a in sortedData)
        erryMaxData = filter(lambda d: d['y']==maxData, sortedData)[0]['erryup']
        maxPreds=[]
        for pr in predictions.keys():
            maxPreds.append(  max(a['y%s' %pr] for a in sortedData)  )
        maxXsecRange = max(maxData + erryMaxData, max(maxPreds) )
        minData=0.0
        erryMinData=1.0
        if any (a['y'] > 0.0 for a in sortedData):
            minData = min(a['y'] for a in sortedData if a['y']>0.)
            erryMinData = filter(lambda d: d['y']==minData, sortedData)[0]['errydown']
        print 'minHX'
        print minHX
        minXsecRange = min(minHX, minData - erryMinData)
        if options.logy:
            if minData - erryMinData < minHX/10000.:
                    minXsecRange = minData/8.
            if (minData - erryMinData) <=0.:
                minXsecRange = min(minHX, minData/8.)
        if minXsecRange <= 0. and not options.logy:
            minXsecRange = 0.
        print 'setting yaxis limits, using min and max:'
        print minXsecRange
        print maxXsecRange

        yspeed = 1.0
        if (obs in observables.keys()):
            if 'yspeed' in observables[obs].keys():
                yspeed = observables[obs]['yspeed']
        print 'yspeed'
        print yspeed

        graph.GetYaxis().SetRangeUser(minXsecRange, (maxXsecRange  - minXsecRange)*2.0*yspeed )
        if options.logy:
            pads[0].SetLogy()
            ##similarly, if logy: find min HX, find max data -> count decades -> scale up by yspeed in the exponent
            print 'setting y limits to'
            print 'min %d' %minXsecRange
            print 'max %d' %(xsecMax*(10**(  yspeed*( abs(np.log10(maxXsecRange) - np.log10( minXsecRange) ) )  )   ))
            graph.GetYaxis().SetRangeUser(minXsecRange, xsecMax*(10**(  1.1*yspeed*( abs(np.log10(maxXsecRange) - np.log10( minXsecRange) ) )  )   )    )  #xsecMax*1023. )
###            graph.GetYaxis().SetRangeUser(minHX, xsecMax*(10**(  0.8*( np.log10(xsecMax) + abs(np.log10( minXsecRange)) )  )   )    )  #xsecMax*1023. )
            ###graph.GetYaxis().SetRangeUser(xsecMin/8., xsecMax*(10**(  0.8*( np.log10(xsecMax) + abs(np.log10( xsecMin/8.)) )  )   )    )  #xsecMax*1023. )

    #    graph.GetYaxis().SetRangeUser(0.0,2.0)
    #    graph.GetYaxis().SetTitleSize(1.2)
    #    graph.GetXaxis().SetTitleSize(1.2)
        x_rg = getListFromBuffer(graph.GetX(), graph.GetN())
        xEu_rg = getListFromBuffer(graph.GetEXhigh(), graph.GetN())
        xEd_rg = getListFromBuffer(graph.GetEXlow(), graph.GetN())
        print 'xaxis range'
        print float(x_rg[0]) - float(xEd_rg[0])
        print float(x_rg[-1]) + float(xEu_rg[-1])
        graph.GetXaxis().SetLimits(float(x_rg[0]) - float(xEd_rg[0]), float(x_rg[-1]) + float(xEu_rg[-1])  )
        graph.Draw("ap0")
        pads[0].Update()
        if 'customLabels' in observables[obs].keys():
            hlabels = convertToHist(graph)
            hlabels.Reset()
            hlabels.GetXaxis().SetTitle( graph.GetXaxis().GetTitle() )
            hlabels.GetYaxis().SetTitle( graph.GetYaxis().GetTitle() )

            hlabels.GetXaxis().SetTitleFont( graph.GetXaxis().GetTitleFont() )
            hlabels.GetYaxis().SetTitleFont( graph.GetYaxis().GetTitleFont() )

            hlabels.GetXaxis().SetTitleSize( graph.GetXaxis().GetTitleSize() )
            hlabels.GetYaxis().SetTitleSize( graph.GetYaxis().GetTitleSize() )

            hlabels.GetXaxis().SetTitleOffset( graph.GetXaxis().GetTitleOffset() )
            hlabels.GetYaxis().SetTitleOffset( graph.GetYaxis().GetTitleOffset() )

            hlabels.GetXaxis().SetLabelFont( graph.GetXaxis().GetLabelFont() )
            hlabels.GetYaxis().SetLabelFont( graph.GetYaxis().GetLabelFont() )

            hlabels.GetXaxis().SetLabelSize( graph.GetXaxis().GetLabelSize() )
            hlabels.GetYaxis().SetLabelSize( graph.GetYaxis().GetLabelSize() )

            hlabels.GetXaxis().SetLabelOffset( graph.GetXaxis().GetLabelOffset() )
            hlabels.GetYaxis().SetLabelOffset( graph.GetYaxis().GetLabelOffset() )

            print 'min max hist'
            print pads[0].GetUymin(), pads[0].GetUymax()
            print graph.GetYaxis().GetXmin(),graph.GetYaxis().GetXmax()
            
            hlabels.GetYaxis().SetRangeUser( graph.GetYaxis().GetXmin() , graph.GetYaxis().GetXmax()  )
            
            for (idx,label) in enumerate(observables[obs]['customLabels']):
                hlabels.GetXaxis().SetBinLabel(idx+1, str(label))
            hlabels.Draw('axis')
            graph.Draw('p0same')

        def safeSquareSubtraction(x1,x2):
            if x1**2 - x2**2 >= 0:
                return np.sqrt(x1**2 -x2**2)
            else:
                return 1e-06
        graphSysts = TGraphAsymmErrors( len(sortedData), array('d', [a['x'] for a in sortedData]), array('d', [a['y'] for a in sortedData]), array('d', [a['errx'] for a in sortedData] ),array('d', [ a['errx'] for a in sortedData]), array('d', [safeSquareSubtraction(a['errydown'], a['errydown_statOnly']) for a in sortedData]), array('d', [safeSquareSubtraction(a['erryup'], a['erryup_statOnly']) for a in sortedData]) )
        print "graphSysts"
        graphSysts.Print()
        shadedBlue= TColor.GetColorTransparent(myColorB1, 0.6)
        graphSysts.SetLineColor(shadedBlue)
        graphSysts.SetFillColor(shadedBlue)
#        graphSysts.SetFillStyle(3001)
        graphSysts.Draw("2same0")


            
        def cleanList(l, tolerance):
            sl = sorted(list(set(l)))
            lc=[]
#            for i in sl[:-1]:
            for i,j in zip(sl[:-1],sl[1:]):
                    
                print 'difference'
                print i
                print j
                print abs(i - j)
                if abs(i - j) < tolerance:
                    
                    print 'is below threshold'
                    print tolerance
                    print 'not appending'
                    continue
                print 'appending'
                lc.append(i)

            lc.append(sl[-1])
            return list(set(lc))
            
        def staggerPoints(graph,pitch,offset):
            print 'staggerPoints'
            print 'offset'
            print offset
            print graph.GetN()
            size = graph.GetN()*pitch
            xl = map(lambda x, edx: x-edx, getListFromBuffer(graph.GetX(),graph.GetN()),getListFromBuffer(graph.GetEXlow(),graph.GetN()))
            xr =( map(lambda x, edx: x+edx, getListFromBuffer(graph.GetX(),graph.GetN()),getListFromBuffer(graph.GetEXhigh(),graph.GetN())))
            print xl
            print xr
            x = sorted(list(set().union(xl,xr)))
            new_xb = []
            for ix in range(len(x)-1):
                step = float((x[ix+1]-x[ix])/float(pitch))
                for i in np.arange(x[ix],x[ix+1],step):
                    print 'adding to new_xb'
                    print i
#                    new_xb.append(float("{0:.5f}".format(i))) ##rounding to the 5th decimal digit to avoid numerical problems in splitting bins
                    new_xb.append(float(i)) 
                new_xb.append(x[ix+1])    
###            new_xb = sorted(list(set(new_xb)))
            print 'new_xb before cleaning'
            print sorted(list(set(new_xb)))
            new_xb = sorted(list(set( cleanList(new_xb,1e-05) )))
            print 'new_xb after cleaning'
            print new_xb
            new_x = map(lambda x1, x2 : x1+ (x2-x1)/2, new_xb[:-1], new_xb[1:])
            new_xEup = map(lambda xc, xb : xb-xc, new_x, new_xb[1:])
            new_xEdown = map(lambda xc, xb : abs(xb-xc), new_x, new_xb[:-1])
            print 'Eup'
            print new_xEup
            print len(new_xEup)
            print 'Edown'
            print new_xEdown
            print len(new_xEdown)
            y = getListFromBuffer(graph.GetY(),graph.GetN())
            yEup = getListFromBuffer(graph.GetEYhigh(),graph.GetN())
            yEdown= getListFromBuffer(graph.GetEYlow(),graph.GetN())
            new_y = [0.]*graph.GetN()*pitch
            new_yEup = [0.]*graph.GetN()*pitch
            new_yEdown= [0.]*graph.GetN()*pitch
            for iy in range(len(y)):
#                print iy
#                print pitch
#                print int(pitch/2)
#                print offset
                ##fill middle point (i.e. old position) + offset
                new_y[iy*pitch + int(pitch/2) + offset] = y[iy]
                new_yEup[iy*pitch + int(pitch/2) + offset] = yEup[iy]
                new_yEdown[iy*pitch + int(pitch/2) + offset] = yEdown[iy]
            print new_y
            print new_yEup
            print new_yEdown
            print new_x
            print new_xb
            print 'len(new_xb)'
            print len(new_xb)
            print 'len(new_x)'
            print len(new_x)
            print 'len(new_y)'
            print len(new_y)
            new_graph = TGraphAsymmErrors(graph.GetN()*pitch, array('d',new_x), array('d',new_y), array('d',new_xEdown), array('d',new_xEup), array('d',new_yEdown), array('d',new_yEup))
            return new_graph

        graphsTheory={}
        graphsTheory_noEY={}
        graphsTheory_staggered={}
        print 'graphsTheory'
        myKeys = ["HX", "ggHpowheg", "ggHamcatnlo", "ggHamcatnloNNLOPS", "SM"]
        if options.skipPOWHEG:
            myKeys = ["HX", "ggHamcatnlo", "ggHamcatnloNNLOPS", "SM"]
#        for pr in predictions.keys():        
        for pr in myKeys:        
###            graphLine = TGraphAsymmErrors(graph.GetN(), graph.GetX(), array('d', [a['y%s'%pr] for a in sortedData]), graph.GetEXlow(), graph.GetEXhigh(), array('d',[0]*graph.GetN()), array('d',[0]*graph.GetN()))
            graphLine = TGraphAsymmErrors(graph.GetN(), graph.GetX(), array('d', [a['y%s'%pr] for a in sortedData]), graph.GetEXlow(), graph.GetEXhigh(), array('d', [a['y%s'%pr]*a['y%sUnc'%pr] for a in sortedData]), array('d', [a['y%s'%pr]*a['y%sUnc'%pr] for a in sortedData]) ) 
            graphLine.SetName(graph.GetName()+"_line_%s" %pr)
            print pr
            if pr not in ['HX','SM']:
                graphLine = AddGraphs(graphLine, graphsTheory['HX']) 
            graphLine_staggered = staggerPoints(graphLine, 11, predictions[pr]['shift']).Clone()
            graphLine_staggered.SetName(graph.GetName()+"_line_staggered_%s" %pr)            
            print 'N new'
            print graphLine_staggered.GetN()
            graphLine_staggered.SetMarkerColor(predictions[pr]['color'])
            graphLine_staggered.SetLineColor(predictions[pr]['color'])
#            shadeColor = TColor.GetColorTransparent(predictions[pr]['color'], 0.9)
#            graphLine_staggered.SetFillColor(shadeColor)
#            graphLine_staggered.SetMarkerSize(0)
            graphLine_staggered.SetFillStyle(predictions[pr]['fill'])
#            graphLine_staggered.SetMarkerStyle(predictions[pr]['marker'])
            graphLine.Print("all")
            graphLine_staggered.Print("all")
            graphLine.SetLineColor(predictions[pr]['color'])
            graphLine_staggered.SetFillColor(predictions[pr]['color'])
            graphLine.SetLineWidth(1)
            graphLine.SetMarkerSize(0)
        
    #    gStyle.SetErrorY(0.)
            graphLine_noEY = setYerrorsToZero(graphLine)
            print 'after yerrors to zero'
            graphLine.Print('all')
            graphLine_noEY.Print('all')
            graphLine_noEY.SetLineColor(predictions[pr]['color'])
            graphLine_noEY.SetLineWidth(1)
            graphLine_noEY.SetMarkerSize(0)
            graphLine_noEY.SetMarkerColor(predictions[pr]['color'])

            if pr is not 'SM':
                if pr is not 'HX':
                    graphLine_noEY.Draw("psameZ")
#                graphLine.Draw("psameZ")
                    graphLine_staggered.Draw("same2")
#                graphLine_staggered.Draw("psame2")
                    graphsTheory_noEY[pr]=graphLine_noEY
                else:
#                    graphLine_noEY.Draw("psameZ")
#                    hHX = graphLine_noEY.GetHistogram()
                    hHX = convertToHist(graphLine_noEY)
                    hHX.SetFillColor(0)
                    hHX.SetFillStyle(0)
                    graphsTheory_noEY[pr]=hHX
                    hHX.Draw('samehist')
#                    graphLine_staggered.Draw("same2")
                    
            graphsTheory[pr]=graphLine

            graphsTheory_staggered[pr]=graphLine_staggered
        

##        print 'LATEXTAB   \multirow{%d}{*}{\textbf{%s}} & $%f$-$%f$  & $%f^{+%.2g}_{-%.2g}$  &  $%f^{+%.2g}_{-%.2g}$  \\\\' % (graph.GetN(), observables[obs]['xlabel'], sortedData[0]['x'] - sortedData[0]['errx'], sortedData[0]['x'] + sortedData[0]['errx'], sortedData[0]['y'], sortedData[0]['erryup'], sortedData[0]['errydown'],  graphsTheory['SM'].GetY()[0], graphsTheory['SM'].GetEYhigh()[0], graphsTheory['SM'].GetEYlow()[0])
##        for ip in range(1,graph.GetN()):
##            print 'LATEXTAB    & $%f$-$%f$  & $%f^{+%.2g}_{-%.2g}$   &  $%f^{+%.2g}_{-%.2g}$  \\\\' % ( sortedData[ip]['x'] - sortedData[ip]['errx'], sortedData[ip]['x'] + sortedData[ip]['errx'], sortedData[ip]['y'], sortedData[ip]['erryup'], sortedData[ip]['errydown'],  graphsTheory['SM'].GetY()[ip], graphsTheory['SM'].GetEYhigh()[ip], graphsTheory['SM'].GetEYlow()[ip])


        print 'LATEXTAB   \multirow{%d}{*}{\textbf{%s}} & $%f$-$%f$  & $%f^{+%f}_{-%f}$  &  $%f^{+%f}_{-%f}$  \\\\' % (graph.GetN(), observables[obs]['xlabel'], sortedData[0]['x'] - sortedData[0]['errx'], sortedData[0]['x'] + sortedData[0]['errx'], sortedData[0]['y'], sortedData[0]['erryup'], sortedData[0]['errydown'],  graphsTheory['SM'].GetY()[0], graphsTheory['SM'].GetEYhigh()[0], graphsTheory['SM'].GetEYlow()[0])
        for ip in range(1,graph.GetN()):
            print 'LATEXTAB    & $%f$-$%f$  & $%f^{+%f}_{-%f}$   &  $%f^{+%f}_{-%f}$  \\\\' % ( sortedData[ip]['x'] - sortedData[ip]['errx'], sortedData[ip]['x'] + sortedData[ip]['errx'], sortedData[ip]['y'], sortedData[ip]['erryup'], sortedData[ip]['errydown'],  graphsTheory['SM'].GetY()[ip], graphsTheory['SM'].GetEYhigh()[ip], graphsTheory['SM'].GetEYlow()[ip])
        

        print 'NNLOPS: powheg / amcatnlo'
       ## for a in sortedData:
       ##     print str(a['yggHpowhegNNLOPS'])+'	'+str(a['yggHamcatnloNNLOPS'])+'	'+str(a['yggHpowhegNNLOPS']/a['yggHamcatnloNNLOPS'])+'	'+str( (a['yggHpowhegNNLOPS']+a['yHX'])/(a['yggHamcatnloNNLOPS']+a['yHX']))
        graph.Draw("psame0")

    #    histo = graph.GetHistogram()
    #    histo.SetLineColor(myColorA0)
    #    histo.Print("all")
    #    histo.Draw('histsame')
        legData = TLegend(0.18,0.66,0.48,0.79)
        legData.SetFillStyle(0)
        legData.SetTextFont(43)
        legData.SetBorderSize(0)
#        legData.AddEntry(graphLine, 'aMC@NLO + NNLOPS', 'l')
        legData.AddEntry(graph, 'Data, stat #oplus syst unc.', 'lpe')
        legData.AddEntry(graphSysts, 'Systematic uncertainty', 'f')

        legRight = TLegend(0.50,0.60,0.85,0.83)
        legRight.SetFillStyle(0)
        legRight.SetTextFont(43)
        legRight.SetBorderSize(0)
#        legRight.AddEntry(graphLine, 'aMC@NLO + NNLOPS', 'l')
        legRight.AddEntry(graphsTheory_noEY['HX'], predictions['HX']['leg'] , 'l')
        for pr in predictions.keys():
 #           legRight.AddEntry(graphsTheory[pr], pr, 'l')
            if pr is not 'SM' and pr is not 'HX':
                if options.skipPOWHEG and pr is 'ggHpowheg': 
                    continue
                legRight.AddEntry(graphsTheory_staggered[pr], predictions[pr]['leg'] , 'lf')
 
        


        legLeft = TLegend(0.22,0.45,0.52,0.82)
        legLeft.SetFillStyle(0)
        legLeft.SetBorderSize(0)
        legLeft.AddEntry(graph, 'Best fit, stat #oplus syst unc.', 'lpe')
        legLeft.AddEntry(graphSysts, 'Syst unc.', 'f')
        for pr in predictions.keys():
 #           legLeft.AddEntry(graphsTheory[pr], pr, 'l')
            if pr is not 'SM':
                if options.skipPOWHEG and pr is 'ggHpowheg': 
                    continue
                legLeft.AddEntry(graphsTheory_staggered[pr], predictions[pr]['leg'] , 'lf')
#        legLeft.AddEntry(graphLine, 'aMC@NLO + NNLOPS', 'l')
#        legLeft.AddEntry(graph, 'Best fit, stat #oplus syst unc.', 'l')
#        legLeft.AddEntry(graphSysts, 'Syst unc.', 'f')

        tex_m=TLatex()
        tex_m.SetNDC()
        tex_m.SetTextAlign(12)
        tex_m.SetTextFont(42)
        tex_m.SetTextSize(0.035)
    #        tex_m.SetLineWidth(2)



        if (obs in observables.keys()):
            if 'legend' in observables[obs].keys():
                print "variable "+str(obs)+" has legend option"
                if observables[obs]['legend'] == 'r':
                    legRight.Draw()
                    tex_m.DrawLatex(0.515,0.57,"#sigma_{SM}(#scale[0.7]{H #rightarrow #gamma#gamma}) from DOI:#scale[0.8]{10.23731/CYRM-2017-002}")
                else:
                    legRight.Draw()
                    tex_m.DrawLatex(0.51,0.845,"LHC HXSWG YR4, m_{H}=125.09")

        else:
            legRight.Draw()
            tex_m.DrawLatex(0.56,0.845,"LHC HXSWG YR4, m_{H}=125.09")
        legData.Draw()



##extra text
        if (obs in observables.keys()):
            if observables[obs]['rescaleLastToSecondLast']:
                tex_m=TLatex()
                tex_m.SetTextAlign(22)
                tex_m.SetTextFont(42)
                tex_m.SetTextSize(0.023)
                tex_m.SetLineWidth(2)
                tex_m.SetTextAngle(0)
                lastBinText = "#sigma_{fid} ("+str(observables[obs]['name'])+">"+(str(lastBinBoundary).replace('.0',''))+str(observables[obs]['unit'])+")/"+(str(secondLastBinWidth).replace('.0',''))
                tex_m.DrawLatex(float(lastX),float(lastMaxY), str(lastBinText) )
            if 'extratext1' in observables[obs].keys():
                print "variable "+str(obs)+" has extratext1"
                tex_m=TLatex()
                tex_m.SetNDC()
                tex_m.SetTextAlign(12)
                tex_m.SetTextFont(62)
                tex_m.SetTextSize(0.03)
                tex_m.SetLineWidth(2)
                tex_m.DrawLatex(0.2,0.62, str(observables[obs]['extratext1']))
            if 'extratext2' in observables[obs].keys():
                print "variable "+str(obs)+" has extratext2"
                tex_m=TLatex()
                tex_m.SetNDC()
                tex_m.SetTextAlign(12)
                tex_m.SetTextFont(62)
                tex_m.SetTextSize(0.03)
                tex_m.SetLineWidth(2)
                tex_m.DrawLatex(0.2,0.56, str(observables[obs]['extratext2']))




        tex_m=TLatex()
        tex_m.SetNDC()
        tex_m.SetTextAlign(12)
        tex_m.SetTextFont(42)
        tex_m.SetTextSize(0.075)
        tex_m.SetLineWidth(2)
#        tex_m.DrawLatex(0.18,0.93,"#bf{CMS}, #it{Preliminary})"
        tex_m.DrawLatex(0.15,0.94,"#bf{CMS}")
            
        tex_m=TLatex()
        tex_m.SetNDC()
        tex_m.SetTextAlign(12)
        tex_m.SetTextFont(42)
        tex_m.SetTextSize(0.065)
    #        tex_m.SetLineWidth(2)
        tex_m.DrawLatex(0.60,0.94,"35.9 fb^{-1} (13 TeV)")

        tex_m=TLatex()
        tex_m.SetNDC()
        tex_m.SetTextAlign(12)
        tex_m.SetTextFont(42)
        tex_m.SetTextSize(0.065)
    #        tex_m.SetLineWidth(2)
        tex_m.DrawLatex(0.19,0.85,"H #rightarrow #gamma#gamma")



        if scaleLastBin != 1.:
            tex_m=TLatex()
            tex_m.SetNDC()
            tex_m.SetTextAlign(12)
            tex_m.SetTextFont(42)
            tex_m.SetTextSize(0.03)
    #        tex_m.SetLineWidth(2)
            tex_m.DrawLatex(0.85,0.85,"SF =%s"%scaleLastBin)
            
##        if 'customLabels' in observables[obs].keys():
##            hlabels.Draw('AXISsame')
##        else:
##            graph.Draw('Asame')



        pads[1].cd()
        pads[1].SetGridy()
        ratio = TGraphAsymmErrors(graph.GetN(), graph.GetX(), array('d', [a['y']/a['ySM'] for a in sortedData]), graph.GetEXhigh(), graph.GetEXlow(), array('d',[a['errydown']/a['ySM'] for a in sortedData]), array('d', [a['erryup']/a['ySM'] for a in sortedData]) )
        y_rgR = getListFromBuffer(ratio.GetY(), ratio.GetN())
        yEu_rgR = getListFromBuffer(ratio.GetEYhigh(), ratio.GetN())
        yEd_rgR = getListFromBuffer(ratio.GetEYlow(), ratio.GetN())
        yMinusE = map(lambda y, ey: y-ey, y_rgR,yEd_rgR)
        print 'yMinusE'
        print yMinusE
        yPlusE  = map(lambda y, ey: y+ey, y_rgR,yEu_rgR)
        print 'yPlusE'
        print yPlusE
        ratio.GetYaxis().SetRangeUser(-0.5,2.5)
        yMinusEmin = float(min(yMinusE))
        yPlusEMax  = float(max(yPlusE))
        if( yMinusEmin<-0.5 or yPlusEMax>2.5):
            print 'ymins below min or yplus above'
            if( yMinusEmin<-0.5 and yPlusEMax>2.5):
                print 'both'
                print 'setting range user to %f, %f' %(yMinusEmin - abs(0.2*yMinusEmin), yPlusEMax + abs(0.2*yPlusEMax))
                ratio.GetYaxis().SetRangeUser(yMinusEmin - abs(0.2*yMinusEmin), yPlusEMax + abs(0.2*yPlusEMax))
            elif yMinusEmin<-0.5:
                print 'ymins below min'
                ratio.GetYaxis().SetRangeUser(yMinusEmin - abs(0.2*yMinusEmin), 2.5)
            elif yPlusEMax>2.5:
                print 'ymins below min or yplus above'
                ratio.GetYaxis().SetRangeUser(-0.5, yPlusEMax + abs(0.2*yPlusEMax))
        ratio.SetMarkerStyle(20)
        ratio.SetMarkerColor(kBlack)
        ratio.SetMarkerSize(0.6)
        print "ratio graph"
        ratio.Print("all")
        ratio.Draw("ap")



        if (obs in observables.keys()):
            if 'xlabel' in observables[obs].keys():
                print "variable "+str(obs)+" has a label for xaxis"
                print "setting xlabel to "+str(observables[obs]['xlabel'])
                ratio.GetXaxis().SetTitle(observables[obs]['xlabel'])
    ###        graph.GetXaxis().SetTitleSize(1.2)
        ratio.GetYaxis().SetTitle("#splitline{Ratio to HX + ggH}{aMC@NLO, #scale[0.9]{NNLOPS}}")
        ratio.GetYaxis().CenterTitle(True)
        ratio.GetXaxis().SetTitleFont(43)
        ratio.GetYaxis().SetTitleFont(43)
        ratio.GetXaxis().SetTitleSize(20)
        ratio.GetYaxis().SetTitleSize(12)
        ratio.GetYaxis().SetTitleOffset(2.6)
        ratio.GetXaxis().SetTitleOffset(3.4)
        ratio.GetYaxis().SetLabelFont(43)
        ratio.GetXaxis().SetLabelFont(43)
#        ratio.GetYaxis().SetLabelSize(0.10)
#        ratio.GetXaxis().SetLabelSize(0.107)
        ratio.GetYaxis().SetLabelSize(16)
        ratio.GetXaxis().SetLabelSize(16)
        ratio.GetYaxis().SetLabelOffset(0.007)
        ratio.GetXaxis().SetLabelOffset(0.007)

        
        x_rgR = getListFromBuffer(ratio.GetX(), ratio.GetN())
        xEu_rgR = getListFromBuffer(ratio.GetEXhigh(), ratio.GetN())
        xEd_rgR = getListFromBuffer(ratio.GetEXlow(), ratio.GetN())
        print 'xaxis range ratio'
        print float(x_rgR[0]) - float(xEd_rgR[0])
        print float(x_rgR[-1]) + float(xEu_rgR[-1])
        ratio.GetXaxis().SetRangeUser(float(x_rgR[0]) - float(xEd_rgR[0]), float(x_rgR[-1]) + float(xEu_rgR[-1]) )

        ratio.GetYaxis().SetNdivisions(507)

        ratio.Draw("ap")
        pads[1].Update()

        if 'customLabels' in observables[obs].keys():
            hlabels_ratio = convertToHist(ratio)
            hlabels_ratio.Reset()
            hlabels_ratio.GetXaxis().SetTitle( ratio.GetXaxis().GetTitle() )
            hlabels_ratio.GetYaxis().SetTitle( ratio.GetYaxis().GetTitle() )

            hlabels_ratio.GetXaxis().SetTitleFont( ratio.GetXaxis().GetTitleFont() )
            hlabels_ratio.GetYaxis().SetTitleFont( ratio.GetYaxis().GetTitleFont() )

            hlabels_ratio.GetXaxis().SetTitleSize( ratio.GetXaxis().GetTitleSize() )
            hlabels_ratio.GetYaxis().SetTitleSize( ratio.GetYaxis().GetTitleSize() )

            hlabels_ratio.GetXaxis().SetTitleOffset( ratio.GetXaxis().GetTitleOffset() )
            hlabels_ratio.GetYaxis().SetTitleOffset( ratio.GetYaxis().GetTitleOffset() )

            hlabels_ratio.GetXaxis().SetLabelFont( ratio.GetXaxis().GetLabelFont() )
            hlabels_ratio.GetYaxis().SetLabelFont( ratio.GetYaxis().GetLabelFont() )

            hlabels_ratio.GetXaxis().SetLabelSize( ratio.GetXaxis().GetLabelSize() )
            hlabels_ratio.GetYaxis().SetLabelSize( ratio.GetYaxis().GetLabelSize() )

            hlabels_ratio.GetXaxis().SetLabelOffset( ratio.GetXaxis().GetLabelOffset() )
            hlabels_ratio.GetYaxis().SetLabelOffset( ratio.GetYaxis().GetLabelOffset() )

            hlabels_ratio.GetYaxis().SetNdivisions( 507 )

            print 'min max hist'
            print ratio.GetYaxis().GetXmin(),ratio.GetYaxis().GetXmax()
            
            hlabels_ratio.GetYaxis().SetRangeUser( ratio.GetYaxis().GetXmin() , ratio.GetYaxis().GetXmax()  )
            
            for (idx,label) in enumerate(observables[obs]['customLabels']):
                print 'setting label ',label
                print 'to bin ',idx
                hlabels_ratio.GetXaxis().SetBinLabel(idx+1, str(label))
            hlabels_ratio.Draw('AXIS')
            hlabels_ratio.Draw('AXIGsame')
            ratio.Draw('psame')
        
#        if 'customLabels' in observables[obs].keys():
#            hlabels_ratio.Draw('AXISsame')
#        else:
#            ratio.Draw('Asame')
    
        for fmt in savefmts:
            savename = str(options.outdir)+"/expectedPrecisionNNLOPS_ub_"+str(obs) 
            if options.resizeLast != -1:
                savename = str(savename)+"_lastBinResizedTo"+str(options.resizeLast)
            if options.resizeFirst != -1000:
                savename = str(savename)+"_firstBinResizedTo"+str(options.resizeFirst)
            if options.hideFirstBin:
                savename = str(savename)+"_hideFirstBin"
            c1.SaveAs(str(savename)+str(fmt))    
        SymmError = [0.5*(a+b) for a,b in zip(ErrUp,ErrDown)]
        MeanError = np.mean(SymmError)
        StdDevError = np.std(SymmError)
        numberOfBins = len(SymmError)
        MedianError = np.median(SymmError)
        Quant25Error = np.percentile(SymmError,25)
        Quant75Error = np.percentile(SymmError,75)
        if(options.skipFirstInMean):
            MeanError = np.mean(SymmError[1:])
            StdDevError = np.std(SymmError[1:])
            numberOfBins = len(SymmError[1:])
            MedianError = np.median(SymmError[1:])
            Quant25Error = np.percentile(SymmError[1:],25)
            Quant75Error = np.percentile(SymmError[1:],75)
    
        out_file = open(str(options.outdir)+"/meanPrecision.txt", "a")
        out_file.write(str(obs)+"	"+str(numberOfBins)+"	"+str(MeanError)+"	"+str(StdDevError)+"	"+str(MedianError)+"	"+str(Quant25Error)+"	"+str(Quant75Error)+"\n")
        out_file.close
##    infile = TFile(options.files, "READ")
##    ws = infile.Get(options.wsname)
##    alldata = ws.allData()
##    genBoundaries=[]
##    recoBoundaries=[]
##    for d in alldata:
##        genVar,recoVar = getVarsName(d.GetName())
##        d.Print()
##        genBs,recoBs = getBinBoundaries(d.GetName())
##        genBoundaries.append(genBs[0])
##        genBoundaries.append(genBs[1])
##        recoBoundaries.append(recoBs[0])
##        recoBoundaries.append(recoBs[1])
##    genBoundaries = sorted(list(set(genBoundaries)))
##    recoBoundaries = sorted(list(set(recoBoundaries)))
##    if options.hideFirstBin:
##        if -1000 in genBoundaries:
##            genBoundaries.remove(-1000)
##        if -1000 in recoBoundaries:
##            recoBoundaries.remove(-1000)
##    lastBinB = genBoundaries[-1]
##    firstBinB = genBoundaries[0]
##    if options.resizeLast != -1:
##        genBoundaries[-1]=options.resizeLast
##        recoBoundaries[-1]=options.resizeLast
##    if options.resizeFirst != -1000:
##        genBoundaries[0]=options.resizeFirst
##        recoBoundaries[0]=options.resizeFirst
##    print array('d', genBoundaries)
##    print array('d', recoBoundaries)
##    
##    label=""
##    respMs={}
##    for cat in options.categories.split(","):
##        respMs[cat] = TH2D("resp_matrix_"+str(label)+"_"+str(cat), "resp_matrix_"+str(label)+"_"+str(cat), len(genBoundaries)-1, array('d',genBoundaries), len(recoBoundaries)-1, array('d',recoBoundaries))
##        respMs[cat].GetXaxis().SetTitle(genVar)
##        respMs[cat].GetYaxis().SetTitle(recoVar)
###        respM.GetXaxis().Print()
###        respM.GetYaxis().Print()
##        for d in alldata:
##            if cat in d.GetName():
##                genBs,recoBs = getBinBoundaries(d.GetName())
##                if genBs[1] == lastBinB and options.resizeLast != -1:
##                    genBs[1]=options.resizeLast
##                if recoBs[1] == lastBinB and options.resizeLast != -1:
##                    recoBs[1]=options.resizeLast
##
##                if genBs[0] == firstBinB and options.resizeFirst != -1000:
##                    genBs[0]=options.resizeFirst
##                if recoBs[0] == firstBinB and options.resizeFirst != -1000:
##5B                    recoBs[0]=options.resizeFirst
##                print "dname: "+str(d.GetName())+", sumW: "+str(d.sumEntries())
##                respMs[cat].SetBinContent(respMs[cat].GetBin( respMs[cat].GetXaxis().FindBin(0.5*(genBs[0]+genBs[1])), respMs[cat].GetYaxis().FindBin(0.5*(recoBs[0]+recoBs[1]) ) ), respMs[cat].GetBinContent( respMs[cat].GetBin( respMs[cat].GetXaxis().FindBin(0.5*(genBs[0]+genBs[1])), respMs[cat].GetYaxis().FindBin(0.5*(recoBs[0]+recoBs[1]) ) ) )  +  d.sumEntries() )
##    cat="all"
##    respMAll = TH2D("resp_matrix_"+str(label)+"_"+str(cat), "resp_matrix_"+str(label)+"_"+str(cat), len(genBoundaries)-1, array('d',genBoundaries), len(recoBoundaries)-1, array('d',recoBoundaries))
##    respMAll.GetXaxis().SetTitle(genVar)
##    respMAll.GetYaxis().SetTitle(recoVar)
##    for rm in respMs.values():
##        respMAll.Add(rm)
##    respMs[cat]=respMAll
##
##    c1=TCanvas()
##    for respM in respMs.keys():
##        respMs[respM].Draw("colz")
##        if options.logx:
##            c1.SetLogx()
##        if options.logy:
##            c1.SetLogy()
##        if options.logz:
##            c1.SetLogz()
##        respMs[respM].Print("ALL")
##        for fmt in savefmts:
##            savename = str(options.outdir)+"/respMatrix_"+str(respM) 
##            if options.resizeLast != -1:
##                savename = str(savename)+"_lastBinResizedTo"+str(options.resizeLast)
##            if options.resizeFirst != -1000:
##                savename = str(savename)+"_firstBinResizedTo"+str(options.resizeFirst)
##            if options.hideFirstBin:
##                savename = str(savename)+"_hideFirstBin"
##            c1.SaveAs(str(savename)+str(fmt))
    

def splitCanvas(c):
    pads=[]
    pads.append(TPad("pad1","pad1",0,0.35,1.,1.))
    pads.append(TPad("pad2","pad2",0,0.0,1.,0.35))
    pads[0].SetBottomMargin(0.06)

    pads[1].SetBottomMargin(0.40)
    pads[1].SetTopMargin(0.038)

    pads[0].SetLeftMargin(0.15)
    pads[1].SetLeftMargin(0.15)
    for pad in pads:
        pad.Draw()
    return pads        
    
    
## ------------------------------------------------------------------------------------------------------------------------------------------------------    
if __name__ == "__main__":
    parser = OptionParser(option_list=[
            make_option("-i", "--indir",
                        action="store", type="string", dest="indir",
                        default="./",
                        help="input directory", metavar="DIR"
                        ),
            make_option("-f", "--files",
                        action="store", type="string", dest="files",
                        default="allSig125IA.root",
                        help="pattern of files to be read", metavar="PATTERN"
                        ), 
            make_option("-F", "--filesFreezeNuis",
                        action="store", type="string", dest="filesFreezeNuis",
                        default="allSig125IA.root",
                        help="pattern of files to be read", metavar="PATTERN"
                        ), 
            make_option("-w", "--wsname",
                        action="store", type="string", dest="wsname",
                        default="cms_hgg_13TeV",
                        help="name of ws to be read", metavar="PATTERN"
                        ), 
            make_option("-t", "--treeName",
                        action="store", type="string", dest="treename",
                        default="TestTree",
                        help="TTree name", metavar="TREENAME"
                        ),
            make_option("-o", "--outfile",
                        action="store", type="string", dest="outfile",
                        default="reduced.root",
                        help="outputfile", metavar="FILE"
                        ),
            make_option("-l", "--label",
                        action="store", type="string", dest="label",
                        default="",
                        help="label", metavar="LABEL"
                        ),

            make_option("-V", "--variable",
                        action="store", dest="variable", type="string",
                        default="PtNJets2p5",
                        help="list of variable"
                        ),
            make_option("-T", "--tmvaSettings",
                        action="store", dest="tmvaSettings", type="string",
                        default="dipho.json",
                        help="settings for the TMVA training"
                        ),
            make_option("-v", "--verbose",
                        action="store_true", dest="verbose",
                        default=False,
                        ),
            make_option("-O", "--optimize",
                        action="store_true", dest="optimize",
                        default=False,
                        ),
            make_option("-D", "--outdir",
                        action="store", type="string", dest="outdir",
                        default="plots_mul",
                        ),
            make_option("-L", "--logz",
                        action="store_true", dest="logz",
                        default=False,
                        ),
            make_option("-x", "--logx",
                        action="store_true", dest="logx",
                        default=False,
                        ),
            make_option("-y", "--logy",
                        action="store_true", dest="logy",
                        default=False,
                        ),
            make_option("-H", "--hideFirstBin",
                        action="store_true", dest="hideFirstBin",
                        default=False,
                        ),
            make_option("-N", "--maxEntries",
                        action="store", type="int", dest="maxEntries",
                        default=-1,
                        ),
            make_option("-R", "--resizeLast",
                        action="store", type="float", dest="resizeLast",
                        default=-1,
                        ),
            make_option("--resizeFirst",
                        action="store", type="float", dest="resizeFirst",
                        default=-1000,
                        ),
            make_option("-C", "--categories",
                        action="store", type="string", dest="categories",
                        default="SigmaMpTTag_0,SigmaMpTTag_1,SigmaMpTTag_2",
                        ),
            make_option("-S", "--skipFirstInMean",
                        action="store_true", dest="skipFirstInMean",
                        default=False,
                        ),
            make_option("-s", "--spectrum",
                        action="store_true", dest="spectrum",
                        default=False,
                        ),

            make_option("-W", "--skipPOWHEG",
                        action="store_true", dest="skipPOWHEG",
                        default=False,
                        ),

            ])

    (options, args) = parser.parse_args()

    sys.argv.append("-b")
    main(options, args)






