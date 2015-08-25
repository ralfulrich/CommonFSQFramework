#!/usr/bin/env python
import CommonFSQFramework.Core.ExampleProofReader
from BadChannels2013 import badChannelsSecMod

import sys, os, time
sys.path.append(os.path.dirname(__file__))
from os import listdir
from os.path import isfile, join

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import edm
from math import sqrt, log10
from array import *
import copy

from outsource_analzye_muon import *

from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

outfolder = os.environ["HaloMuonOutput"]

class Step2_Selection_2(CommonFSQFramework.Core.ExampleProofReader.ExampleProofReader):
    def init( self, maxEvents = None):
        #####
        self.flag_use_merjin_electronic_channel_noise = False # mine RMS
        #####


#        self.maxFileNo = slaveParameters["maxFileNo"]
        firstRun = (self.maxFileNo == -1)

        self.maxEvents = maxEvents
        self.hist = {}
        self.nEventsRandom = 0
        #self.tree = {}
        self.hist["EventCount"] = ROOT.TH1D("EventCount","EventCount",10,0,20)
       
        #MC
        #if not self.isData:

        #Getting sector RMS and Mean
        if firstRun:
            inputFile = ROOT.TFile(join(outfolder,"mean_rms_2.root"))
            inputFile_noise = ROOT.TFile(join(outfolder,"Histograms_StatFitInformationPedNoise.root"))
            


        else:
            inputFile = ROOT.TFile(join(outfolder,"plotsMuonselectioncuts_2_{n:04d}.root".format(n=self.maxFileNo)))
            # check naming if rerunning step1 again
        
        
        # create the branches and assign the fill-variables to them

        # ftree = ROOT.TFile("tree.root", "recreate")
        # tree = ROOT.TTree("nt", "nt")
        # np= np.zeros(1, dtype=float)
        # tree.Branch('ch_energy', ch_energy, 'ch_energy')

        hist_sec_Mean = inputFile.Get("data_PPMinBias_Run2013/hist_sec_Mean")
        hist_sec_RMS = inputFile.Get("data_PPMinBias_Run2013/hist_sec_RMS")

        #Getting channel RMS and Mean
        hist_ch_Mean = inputFile.Get("data_PPMinBias_Run2013/hist_ch_Mean")
        hist_ch_RMS =  inputFile.Get("data_PPMinBias_Run2013/hist_ch_RMS")

        #Getting electronic noise from in the Merijns file
        hist_ch_noise_RMS= inputFile_noise.Get("SumHistoNoiseFit")

        #hist_sec_noise_RMS= inputFile_noise.Get("SumHistoNoiseFit")
        self.hist["2DMuonCountMap"] =  ROOT.TH2D("2DMuonCountMap","2DMuonCountMap", 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist["2DMuonCountFor3Ch"] =  ROOT.TH2D("2DMuonCountFor3Ch","2DMuonCountFor3Ch", 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist["2DMuonactiveregion"] =  ROOT.TH2D("2DMuonactiveregion","2DMuonactiveregion", 14, 0.5, 14.5, 16, 0.5, 16.5)
        #self.hist["2DMuonCountMap_allch"] =  ROOT.TH2D("2DMuonCountMap_allch","2DMuonCountMap_allch", 14, 0.5, 14.5, 16, 0.5, 16.5)

        self.hist["2DcountChannelAboveTwoSigma_AllEvt"] =  ROOT.TH2D("2DcountChannelAboveTwoSigma_AllEvt","2DcountChannelAboveTwoSigma_AllEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist["2DcountChannelAboveTwoSigma_RndEvt"] =  ROOT.TH2D("2DcountChannelAboveTwoSigma_RndEvt","2DcountChannelAboveTwoSigma_RndEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist["2DcountChannelAboveTwoSigmaFor3SigmaSectors_AllEvt"] =  ROOT.TH2D("2DcountChannelAboveTwoSigmaFor3SigmaSectors_AllEvt","2DcountChannelAboveTwoSigmaFor3SigmaSectors_AllEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist["2DcountChannelAboveTwoSigmaFor3SigmaSectors_RndEvt"] =  ROOT.TH2D("2DcountChannelAboveTwoSigmaFor3SigmaSectors_RndEvt","2DcountChannelAboveTwoSigmaFor3SigmaSectors_RndEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)
        # self.hist["2DcountChannelAboveTwoSigmaForMuonCandidateSector_AllEvt"] =  ROOT.TH2D("2DcountChannelAboveTwoSigmaForMuonCandidateSector_AllEvt","2DcountChannelAboveTwoSigmaForMuonCandidateSector_AllEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)
        # self.hist["2DcountChannelAboveTwoSigmaForMuonCandidateSector_RndEvt"] =  ROOT.TH2D("2DcountChannelAboveTwoSigmaForMuonCandidateSector_RndEvt","2DcountChannelAboveTwoSigmaForMuonCandidateSector_RndEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)

        self.hist["2DcountChannelsAboveNoiseForAllSectors_AllEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseForAllSectors_AllEvt","2DcountChannelsAboveNoiseForAllSectors_AllEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
        self.hist["2DcountChannelsAboveNoiseForAllSectors_RndEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseForAllSectors_RndEvt","2DcountChannelsAboveNoiseForAllSectors_RndEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
        self.hist["2DcountChannelsAboveNoiseFor3SigmaSectors_AllEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseFor3SigmaSectors_AllEvt","2DcountChannelsAboveNoiseFor3SigmaSectors_AllEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
        self.hist["2DcountChannelsAboveNoiseFor3SigmaSectors_RndEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseFor3SigmaSectors_RndEvt","2DcountChannelsAboveNoiseFor3SigmaSectors_RndEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
        # self.hist["2DcountChannelsAboveNoiseForMuonCandidateSector_AllEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseForMuonCandidateSector_AllEvt","2DcountChannelsAboveNoiseForMuonCandidateSector_AllEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
        # self.hist["2DcountChannelsAboveNoiseForMuonCandidateSector_RndEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseForMuonCandidateSector_RndEvt","2DcountChannelsAboveNoiseForMuonCandidateSector_RndEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)

        self.hist["2DSecRMSHot_Vs_SecRMSSecondHot_AllEvt"] = ROOT.TH2D("2DSecRMSHot_Vs_SecRMSSecondHot_AllEvt","2DSecRMSHot_Vs_SecRMSSecondHot_AllEvt",88, -2, 20, 88, -2, 20)
        self.hist["2DSecRMSHot_Vs_SecRMSSecondHot_RndEvt"] = ROOT.TH2D("2DSecRMSHot_Vs_SecRMSSecondHot_RndEvt","2DSecRMSHot_Vs_SecRMSSecondHot_RndEvt",88, -2, 20, 88, -2, 20)
        self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_AllEvt"] = ROOT.TH2D("2DDeltaSigma_Vs_SecRMSSecondHot_AllEvt","2DDeltaSigma_Vs_SecRMSSecondHot_AllEvt",88,-2,20, 88, -2, 20)
        self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_RndEvt"] = ROOT.TH2D("2DDeltaSigma_Vs_SecRMSSecondHot_RndEvt","2DDeltaSigma_Vs_SecRMSSecondHot_RndEvt",88, -2, 20, 88,-2, 20)

        # self.hist["2DsectorChMeanForAllSectors"] = ROOT.TH2D("2DsectorChMeanForAllSectors","2DsectorChMeanForAllSectors", 14, -0.5, 13.5, 16, -0.5, 15.5)
        # self.hist["2DsectorChMeanForAllSectors_RndEvt"] = ROOT.TH2D("2DsectorChMeanForAllSectors_RndEvt","2DsectorChMeanForAllSectors_RndEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
        # self.hist["2DsectorChMeanFor3sigmaSectors"] = ROOT.TH2D("2DsectorChMeanFor3sigmaSectors","2DsectorChMeanFor3sigmaSectors", 14, -0.5, 13.5, 16, 0.5, 16.5)
        # self.hist["2DsectorChMeanFor3sigmaSectors_RndEvt"] = ROOT.TH2D("2DsectorChMeanFor3sigmaSectors_RndEvt","2DsectorChMeanFor3sigmaSectors_RndEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
        # self.hist["2DsectorChMeanForMuonCandidateSectors"] = ROOT.TH2D("2DsectorChMeanForMuonCandidateSectors","2DsectorChMeanForMuonCandidateSectors", 14, -0.5, 13.5, 16, 0.5, 16.5)
        # self.hist["2DsectorChMeanForMuonCandidateSectors_RndEvt"] = ROOT.TH2D("2DsectorChMeanForMuonCandidateSectors_RndEvt","2DsectorChMeanForMuonCandidateSectors_RndEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)

        # self.hist["2DsectorChRMSForAllSectors"] = ROOT.TH2D("2DsectorChRMSForAllSectors","2DsectorChRMSForAllSectors", 14, -0.5, 13.5, 16, 0.5, 16.5)
        # self.hist["2DsectorChRMSForAllSectors_RndEvt"] = ROOT.TH2D("2DsectorChRMSForAllSectors_RndEvt","2DsectorChRMSForAllSectors_RndEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
        # self.hist["2DsectorChRMSFor3sigmaSectors"] = ROOT.TH2D("2DsectorChRMSFor3sigmaSectors","2DsectorChRMSFor3sigmaSectors", 14, -0.5, 13.5, 16, 0.5, 16.5)
        # self.hist["2DsectorChRMSFor3sigmaSectors_RndEvt"] = ROOT.TH2D("2DsectorChRMSFor3sigmaSectors_RndEvt","2DsectorChRMSFor3sigmaSectors_RndEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
        # self.hist["2DsectorChRMSForMuonCandidateSectors"] = ROOT.TH2D("2DsectorChRMSForMuonCandidateSectors","2DsectorChRMSForMuonCandidateSectors", 14, -0.5, 13.5, 16, 0.5, 16.5)
        # self.hist["2DsectorChRMSForMuonCandidateSectors_RndEvt"] = ROOT.TH2D("2DsectorChRMSForMuonCandidateSectors_RndEvt","2DsectorChRMSForMuonCandidateSectors_RndEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
        # self.hist["2DsectorChRMS_Vs_ChAboveNoise"] = ROOT.TH2D("2DsectorChRMS_Vs_ChAboveNoise","2DsectorChRMS_Vs_ChAboveNoise", 14, -0.5, 13.5, 14, -0.5, 13.5)
        # self.hist["2DsectorChRMS_Vs_ChAboveNoise_RndEvt"] = ROOT.TH2D("2DsectorChRMS_Vs_ChAboveNoise_RndEvt","2DsectorChRMS_Vs_ChAboveNoise_RndEvt", 14, -0.5, 13.5, 14, -0.5, 13.5)
        self.hist["2DcountSectorRMS_AllEvt"] =  ROOT.TH2D("2DcountSectorRMS_AllEvt","2DcountSectorRMS_AllEvt", 88, -2, 20, 16, 0.5, 16.5)
        self.hist["2DcountSectorRMS_RndEvt"] =  ROOT.TH2D("2DcountSectorRMS_RndEvt","2DcountSectorRMS_RndEvt", 88, -2, 20, 16, 0.5, 16.5)

        # self.hist["2DMuonNoTriggerCountMap"] =  ROOT.TH2D("2DMuonNoTriggerCountMap","2DMuonNoTriggerCountMap", 14, -0.5, 13.5, 16, -0.5, 15.5)
        self.hist["GoodMuonCountPerSec"] =  ROOT.TH1D("GoodMuonCountPerSec","GoodMuonCountPerSec", 16, 0.5, 16.5)
        self.hist["RunsWithGoodMuons"] =  ROOT.TH1D("RunsWithGoodMuons","RunsWithGoodMuons", 10000, 247000-0.5, 257000-0.5)
        self.hist["RunsAllTrigger"] =  ROOT.TH1D("RunsAllTrigger","RunsAllTrigger", 10000, 247000-0.5, 257000-0.5)
        self.hist["DeltaSigma_AllEvt"] = ROOT.TH1D("DeltaSigma_AllEvt","DeltaSigma_AllEvt",100, 0, 50)
        self.hist["DeltaSigma_RndEvt"] = ROOT.TH1D("DeltaSigma_RndEvt","DeltaSigma_RndEvt",100, 0, 50)

        histcalibrationname = '2DMuonSignalMap'
        histcalibrationKatname = '2DMuonSignalMap_katerina'
        histcalibrationIgorname = '2DMuonSignalMap_igor'
        histcalibrationRatioIgorname = "2DMuonSignalMapIgor_Ratio"
        histcalibrationRatioIgor_Katname = "2DMuonSignalMapIgor_KatRatio"
        histcalibrationRationame = "2DMuonSignalMap_Ratio"
        histcalibrationname_rms = '2DMuonSignalMap_rms'
        histcalibration_notdividedRefname = '2DMuonSignalMap_notdividedRef'
        self.hist[histcalibration_notdividedRefname] =  ROOT.TH2D(histcalibration_notdividedRefname,histcalibration_notdividedRefname, 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist[histcalibrationKatname] =  ROOT.TH2D(histcalibrationKatname,histcalibrationKatname, 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist[histcalibrationIgorname] =  ROOT.TH2D(histcalibrationIgorname,histcalibrationIgorname, 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist[histcalibrationRationame] =  ROOT.TH2D(histcalibrationRationame,histcalibrationRationame, 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist[histcalibrationRatioIgorname] =  ROOT.TH2D(histcalibrationRatioIgorname,histcalibrationRatioIgorname, 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist[histcalibrationRatioIgor_Katname] =  ROOT.TH2D(histcalibrationRatioIgor_Katname,histcalibrationRatioIgor_Katname, 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist[histcalibrationname_rms] =  ROOT.TH2D(histcalibrationname_rms,histcalibrationname_rms, 14, 0.5, 14.5, 16, 0.5, 16.5)
        if not firstRun:
            self.hist[histcalibrationname] = inputFile.Get("data_MinimumBias_Run2015A/2DMuonSignalMap")
           # print "Extracted histogram from file. Checking entries:",  self.hist[histcalibrationname].GetEntries()
        else: #first time running analyser the calibration constants are all set to 1
            self.hist[histcalibrationname] =  ROOT.TH2D(histcalibrationname,histcalibrationname, 14, 0.5, 14.5, 16, 0.5, 16.5)
           
            for imod in xrange(0,14):
                for isec in xrange(0,16):
                    self.hist[histcalibrationname].SetBinContent(self.hist[histcalibrationname].FindBin(imod+1,isec+1), 1.) #all factors set to 1
                     
                    print "Created new 2d calibration map. Checking entries:",  self.hist[histcalibrationname].GetEntries(), ". maxFileNo =", self.maxFileNo, firstRun
        #histCalibration.SetBit(ROOT.TH1.kIsAverage) # does not seem to work. use TProfile2D instead

        #make new histograms
        for isec in xrange(0,16):
            for imod in xrange(0,14):
                hname = 'MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                hwithouttrigger = 'MuonSignalSecCh_withoutrigger_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                hdifferenmuonselect = 'MuonSignalSecCh_differentSelection_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                hnoise_neighbor = 'CastorNoise_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                hnoise_electronic_RMS = 'CastorNoise_electronic_RMS_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                hnoise_randomtrg = 'CastorNoise_randomtrg_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
               
                #hempty = 'CastorEmptysectors_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                self.hist[hname] = ROOT.TH1D( hname, hname, 50, -100, 400)
                self.hist[hwithouttrigger] = ROOT.TH1D(hwithouttrigger, hwithouttrigger, 50, -100, 400)
                self.hist[hdifferenmuonselect] = ROOT.TH1D(hdifferenmuonselect, hdifferenmuonselect, 50, -100, 400)
                self.hist[hnoise_neighbor] = ROOT.TH1D(hnoise_neighbor, hnoise_neighbor, 50, -100, 400)
                self.hist[hnoise_randomtrg] = ROOT.TH1D(hnoise_randomtrg, hnoise_randomtrg, 50, -100, 400)
                self.hist[hnoise_electronic_RMS] = ROOT.TH1D(hnoise_electronic_RMS,hnoise_electronic_RMS, 224,-0.5,223.5)
                
            henergy= 'MuonSignalSec_energy_sec_{sec}'.format(sec=str(isec+1))
            self.hist[henergy] = ROOT.TH1D( henergy, henergy, 50, 0, 400)
        
        histMuonSignalMean = "1DMuonsignalMean"
        self.hist[histMuonSignalMean] = ROOT.TH1D(histMuonSignalMean,histMuonSignalMean,224,-0.5,223.5)
        histMuonSignalKatMean = "1DMuonsignalKaterinaMean"
        self.hist[histMuonSignalKatMean] = ROOT.TH1D(histMuonSignalKatMean,histMuonSignalKatMean,224,-0.5,223.5)
        histMuonSignalRatioMean = "1DMuonsignalRatioMean"
        self.hist[histMuonSignalRatioMean] = ROOT.TH1D(histMuonSignalRatioMean, histMuonSignalRatioMean,224,-0.5,223.5)
        histMuonSignalRatioMeanPull = "1DMuonsignalRatioMeanPull"
        self.hist[histMuonSignalRatioMeanPull] = ROOT.TH1D(histMuonSignalRatioMeanPull, histMuonSignalRatioMeanPull,224,-0.5,223.5)
        histMuonSignalRatioLog = "1DMuonsignalRatioLog"
        self.hist[histMuonSignalRatioLog] = ROOT.TH1D(histMuonSignalRatioLog, histMuonSignalRatioLog,80,-2,2)
        
        histMuonSignalPull_Meli_Kat = "1DMuonsignalPull_Meli_Kat"
        self.hist[histMuonSignalPull_Meli_Kat] = ROOT.TH1D(histMuonSignalPull_Meli_Kat, histMuonSignalPull_Meli_Kat,40,-4,4)
        histMuonSignalPull_Meli_Igor = "1DMuonsignalPull_Meli_Igor"
        self.hist[histMuonSignalPull_Meli_Igor] = ROOT.TH1D(histMuonSignalPull_Meli_Igor, histMuonSignalPull_Meli_Igor,40,-4,4)
        histMuonSignalPull_Kat_Igor = "1DMuonsignalPull_Kat_Igor"
        self.hist[histMuonSignalPull_Kat_Igor] = ROOT.TH1D(histMuonSignalPull_Kat_Igor, histMuonSignalPull_Kat_Igor,40,-4,4)
        
        hnameAllsec= 'MuonSignalAllSec_energy'
        self.hist[hnameAllsec] = ROOT.TH1D(hnameAllsec, hnameAllsec,50, 0, 400)
       


               #self.hist[hempty] = ROOT. TH1D( hempty, hempty, 500, -100, 400)

         
       

        #get channel energies from input file
        self.ch_mean = [[0 for _ in xrange(14)] for _ in xrange(16)]
        self.ch_RMS = [[0 for _ in xrange(14)] for _ in xrange(16)]
        for i in xrange(0, 224):
            #how I get Channel information
            isec = i//14
            imod = i%14
            #print sec, mod
            #i = sec*14+module   -> sec0 mod0 -> i=0   ; sec9 mod 4
            self.ch_mean[isec][imod] = hist_ch_Mean.GetBinContent(i+1) #+1 because of underflow
            if self.flag_use_merjin_electronic_channel_noise:
                self.ch_RMS[isec][imod] = hist_ch_noise_RMS.GetBinContent(i+1) #+1 because of underflow
            else:
                self.ch_RMS[isec][imod] = hist_ch_RMS.GetBinContent(i+1) #+1 because of underflow

            #print i, "sec,mod", sec, mod, hist_ch_RMS.GetBinContent(i+1)


        #get mean and rms from sector histogram
        self.sec_mean = [0] * 16
        self.sec_RMS = [0] * 16
        for i in xrange(0,16):
            self.sec_mean[i] = hist_sec_Mean.GetBinContent(i+1) #+1 because of underflow
            self.sec_RMS[i] = hist_sec_RMS.GetBinContent(i+1) #+1 because of underflow

        #these ones are needed to update mean and RMS for the new calibration constants
        self.new_ch_mean = [[0 for _ in xrange(14)] for _ in xrange(16)]
        self.new_ch_RMS = [[0 for _ in xrange(14)] for _ in xrange(16)]
        self.new_sec_mean = [0] * 16
        self.new_sec_RMS = [0] * 16

        #TProfile for storing means and RMS. and when jobs are merged, the averge is taken
        self.hist['hist_sec_Mean'] = ROOT.TProfile('hist_sec_Mean','hist_sec_Mean',16,0.5,16.5)
        self.hist['hist_sec_RMS'] = ROOT.TProfile('hist_sec_RMS','hist_sec_RMS',16,0.5,16.5) #change later to hist_sec_mean
        self.hist['hist_ch_Mean'] = ROOT.TProfile('hist_ch_Mean','hist_ch_Mean',224,0.5,224.5)
        self.hist['hist_ch_RMS'] = ROOT.TProfile('hist_ch_RMS','hist_ch_RMS',224,0.5,224.5)
        

        for h in self.hist:
            self.hist[h].Sumw2()
            self.GetOutputList().Add(self.hist[h])


    def getListSigmaSector(self, energy_sec):
        listSigmaSector = [0] * 16

        for i in xrange(0,16):
            if self.sec_RMS[i]:
                sigma = (energy_sec[i] - self.sec_mean[i]) / self.sec_RMS[i]
            else:
                print "Warning: RMS of sec", i, "is zero. Assume triggered"
                sigma = np.sign(energy_sec[i] - self.sec_mean[i]) * 1e9
            listSigmaSector[i] = sigma

        return listSigmaSector

    def getListSigmaChannel(self, energy_ch):
        listSigmaChannel = [[0 for _ in xrange(14)] for _ in xrange(16)]
        for isec in xrange(16):
            for imod in xrange(14):
                energy = energy_ch[isec][imod]
                # print "RMS mod {mod} sec {sec}: {rms}".format(mod=str(imod),sec=str(isec),rms=self.ch_RMS[isec][imod])
                if [isec+1,imod+1] in badChannelsSecMod:
                    #print "skipping channel", imod, isec
                    continue

                sigma = 0
                if self.ch_RMS[isec][imod] != 0:
                    sigma = (energy - self.ch_mean[isec][imod]) / self.ch_RMS[isec][imod]

                listSigmaChannel[isec][imod] = sigma

        return listSigmaChannel


    def filterListSector(self, l, sigma_th):
        listSectorsAboveThr = []
        for isec in xrange(16):
            if l[isec] > sigma_th: listSectorsAboveThr.append(isec)
        return listSectorsAboveThr

    def filterListChannel(self, l, sigma_th):
        listChannelsAboveThr = []
        for isec in xrange(16):
            for imod in xrange(14):
                if l[isec][imod] > sigma_th: listChannelsAboveThr.append([isec,imod])
        return listChannelsAboveThr


    def filterListSigma(self, l, sigma_th):
        return [[ii,isigma] for ii, isigma in l if isigma > sigma_th]

    def getlistSectorMeanRMSInZ(self, ch_energy):
        zmean = [0] * 16
        zrms = [0] * 16
        esecsum = [0] * 16
        for isec in xrange(0,16):
            for imod in xrange(0,14):
                if [isec+1, imod+1] in badChannelsSecMod:
                    #print "skipping channel", imod, isec
                    continue
                esecsum[isec] += ch_energy[isec][imod]
                zmean[isec] += ch_energy[isec][imod] * imod
                zrms[isec] += ch_energy[isec][imod] * imod*imod
            zmean[isec] /= esecsum[isec]
            zrms[isec] /= esecsum[isec]

            tmpvalue = zrms[isec]-(zmean[isec]*zmean[isec])
            if tmpvalue <= 0:
                zrms[isec] = -100
            else:
                zrms[isec] = sqrt(tmpvalue)

        return zmean, zrms


    def analyze(self):
        weight = 1
        num = 0
        goodMuonEvent = False
        goodMuonEventWithoutAnyTriggerSelection = False
        goodMuonEventDifferentSelection = False


        # genTracks
        #num = self.fChain.genTracks.size()
        #print num
        #print self.maxEta # see slaveParams below
        
       

        self.hist["EventCount"].Fill("all",1)

        isBptxminus = True
            #if not (self.fChain.trgl1L1GTTech[1] and self.fChain.trgl1L1GTTech[2]): #not (bptx+ and bptx-)
        if (self.fChain.trgl1L1GTTech[1]) and not (self.fChain.trgl1L1GTTech[2]):
           isBptxminus = False


        
        self.hist["EventCount"].Fill("bptx +",1)

        # if self.flag_use_merjin_electronic_channel_noise:
        for isec in xrange(16):
            for imod in xrange(14):
                hnoise_electronic_RMS = 'CastorNoise_electronic_RMS_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                self.hist[hnoise_electronic_RMS].Fill(self.ch_RMS[isec][imod])
                
        histcalibrationname = '2DMuonSignalMap'
        histCalibration = self.hist[histcalibrationname]
        histcalibrationname_rms = '2DMuonSignalMap_rms'
       
       
       
        ch_energy = [[0.0 for _ in xrange(14)] for _ in xrange(16)]
        sec_energy = [0.0] * 16
        sec_front_energy = [0.0] * 16

        if self.fChain.CastorRecHitEnergy.size() != 224:
            return 0

        self.hist["EventCount"].Fill("size rh",1)

        isRandom = False
        if self.fChain.trgRandom:
            if self.fChain.trgl1L1GTTech[7]: #not bptx+ or bptx-
                self.nEventsRandom += 1
                isRandom = 1

        if isRandom: self.hist["EventCount"].Fill("rnd trg",1)


        for i in xrange(224):
            isec = self.fChain.CastorRecHitSector.at(i)-1
            imod = self.fChain.CastorRecHitModule.at(i)-1
            rh_energy = self.fChain.CastorRecHitEnergy.at(i)

            binnumber = histCalibration.FindBin(imod+1,isec+1)
            #binnumber_notdividedRef = histcalibration_notdividedRef.FindBin(imod+1,isec+1)
            calibration = histCalibration.GetBinContent(binnumber)

            ich_energy = calibration * rh_energy
            ch_energy[isec][imod] = ich_energy

            if [isec+1, imod+1] not in badChannelsSecMod:
                sec_energy[isec] += ich_energy
                if imod < 5: sec_front_energy[isec] += ich_energy

            if isRandom:
                hnoise_randomtrg = 'CastorNoise_randomtrg_mod_{imod}_sec_{isec}'.format(imod=str(imod+1), isec=str(isec+1))
                self.hist[hnoise_randomtrg].Fill(ch_energy[isec][imod])
                self.new_ch_mean[isec][imod] += ich_energy
                self.new_ch_RMS[isec][imod] += ich_energy**2
                if [isec+1, imod+1] not in badChannelsSecMod:
                    self.new_sec_mean[isec] += ich_energy
                    self.new_sec_RMS[isec] += ich_energy**2






        ######################
        # Setup muon trigger #
        ######################
        hasMuonTrigger = (self.fChain.trgCastorHaloMuon or self.fChain.trgl1L1GTAlgo[102])
        # if not hasMuonTrigger: return 1
        if hasMuonTrigger: self.hist["EventCount"].Fill("muon trg",1)


            #  print 'sec', sec, 'mod', mod, 'e_ch=', ch_energy, "e_sec", sec_energy[sec]
            #  print "all sector energies:  ", sec_energy


        # zmean, zrms = self.getlistSectorMeanRMSInZ(ch_energy)

        listSigmaSector = self.getListSigmaSector(sec_energy)
        listSigmaChannel = self.getListSigmaChannel(ch_energy)

       

        SigmaHottestSector = None
        SigmaSecHottestSector = None
        HottestSector = None
        SecondHottestSector = None
        for isec in xrange(16):
            sigma = listSigmaSector[isec]

            self.hist["2DcountSectorRMS_AllEvt"].Fill(sigma,isec+1)
            if isRandom:
                self.hist["2DcountSectorRMS_RndEvt"].Fill(sigma,isec+1)

            if sigma > SigmaHottestSector:
                SigmaSecHottestSector = SigmaHottestSector
                SecondHottestSector = HottestSector

                SigmaHottestSector = sigma
                HottestSector = isec

            elif sigma > SigmaSecHottestSector:
                SigmaSecHottestSector = sigma
                SecondHottestSector = isec

        DeltaSigma= (SigmaHottestSector - SigmaSecHottestSector)
        self.hist["2DSecRMSHot_Vs_SecRMSSecondHot_AllEvt"].Fill(SigmaHottestSector,SigmaSecHottestSector)
        self.hist["DeltaSigma_AllEvt"].Fill(DeltaSigma)
        self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_AllEvt"].Fill(DeltaSigma,SigmaSecHottestSector)
        if isRandom:
            self.hist["2DSecRMSHot_Vs_SecRMSSecondHot_RndEvt"].Fill(SigmaHottestSector,SigmaSecHottestSector)
            self.hist["DeltaSigma_RndEvt"].Fill(DeltaSigma)
            self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_RndEvt"].Fill(DeltaSigma,SigmaSecHottestSector)
        # selection of all channels above noise
        listAllChannelsAboveNoise = [[] for _ in xrange(16)]


        for isec in xrange(16):
            for imod in xrange(14):
                sigma = listSigmaChannel[isec][imod]
                if sigma > 2 and [isec+1,imod+1] not in badChannelsSecMod:
                    listAllChannelsAboveNoise[isec].append(imod)
                    self.hist["2DcountChannelAboveTwoSigma_AllEvt"].Fill(imod+1,isec+1)
                    if isRandom:
                        self.hist["2DcountChannelAboveTwoSigma_RndEvt"].Fill(imod+1,isec+1)

        for isec in xrange(0,16):
            self.hist["2DcountChannelsAboveNoiseForAllSectors_AllEvt"].Fill(len(listAllChannelsAboveNoise[isec]),isec+1)
            # self.hist["2DsectorChMeanForAllSectors"].Fill(zmean[isec],isec+1)
            # self.hist["2DsectorChRMSForAllSectors"].Fill(zrms[isec],isec+1)
            # self.hist["2DsectorChRMS_Vs_ChAboveNoise"].Fill(len(listAllChannelsAboveNoise[isec]),zrms[isec])
            if isRandom:
                self.hist["2DcountChannelsAboveNoiseForAllSectors_RndEvt"].Fill(len(listAllChannelsAboveNoise[isec]),isec+1)
                # self.hist["2DsectorChMeanForAllSectors_RndEvt"].Fill(zmean[isec],isec+1)
                # self.hist["2DsectorChRMSForAllSectors_RndEvt"].Fill(zrms[isec],isec+1)
                # self.hist["2DsectorChRMS_Vs_ChAboveNoise_RndEvt"].Fill(len(listAllChannelsAboveNoise[isec]),zrms[isec])

        ################################
        # cut on muon candidate sector #
        ################################
        if SigmaSecHottestSector > 2.5:
            return 0

        muonSec = HottestSector

        for imod in xrange(14):
            sigma = listSigmaChannel[muonSec][imod]
            if [muonSec+1,imod+1] in badChannelsSecMod:
                #print "skipping channel", mod, sec
                continue
            if sigma > 2: 
                # listAllChannelsAboveNoise[muonSec].append(imod)
                self.hist["2DcountChannelAboveTwoSigmaFor3SigmaSectors_AllEvt"].Fill(imod+1,muonSec+1)
                if isRandom:
                    self.hist["2DcountChannelAboveTwoSigmaFor3SigmaSectors_RndEvt"].Fill(imod+1,muonSec+1)
                    

        self.hist["2DcountChannelsAboveNoiseFor3SigmaSectors_AllEvt"].Fill(len(listAllChannelsAboveNoise[muonSec]),muonSec+1)
        # self.hist["2DsectorChMeanFor3sigmaSectors"].Fill(zmean[muonSec],muonSec+1)
        # self.hist["2DsectorChRMSFor3sigmaSectors"].Fill(zrms[muonSec],muonSec+1)
        if isRandom:
            self.hist["2DcountChannelsAboveNoiseFor3SigmaSectors_RndEvt"].Fill(len(listAllChannelsAboveNoise[muonSec]),muonSec+1)
            self.hist["EventCount"].Fill("Rntrgsigma sec cut",1)
            # self.hist["2DsectorChMeanFor3sigmaSectors_RndEvt"].Fill(zmean[muonSec],muonSec+1)
            # self.hist["2DsectorChRMSFor3sigmaSectors_RndEvt"].Fill(zrms[muonSec],muonSec+1)
        if hasMuonTrigger:
           self.hist["EventCount"].Fill("MuTrgsigma sec cut",1)

        if not isBptxminus:
           self.hist["EventCount"].Fill("Bptx(+)sigma sec cut",1)



        self.hist["EventCount"].Fill("sigma sec cut",1)

        ################################
        # cut oon noise veto sector #
        ################################
        # if DeltaSigma < 1:
        #     return 0

        # for imod in xrange(14):
        #     sigma = listSigmaChannel[isec][imod]
        #     if [isec+1,imod+1] in badChannelsSecMod:
        #         #print "skipping channel", imod, isec
        #         continue
        #     if sigma > 2:
        #         # listAllChannelsAboveNoise[isec].append(imod)
        #         self.hist["2DcountChannelAboveTwoSigmaForMuonCandidateSector_AllEvt"].Fill(imod+1,muonSec+1)
        #         if isRandom:
        #             self.hist["2DcountChannelAboveTwoSigmaForMuonCandidateSector_RndEvt"].Fill(imod+1,muonSec+1)

        # self.hist["2DcountChannelsAboveNoiseForMuonCandidateSector_AllEvt"].Fill(len(listAllChannelsAboveNoise[muonSec]),muonSec+1)
        # self.hist["2DsectorChMeanForMuonCandidateSectors"].Fill(zmean[muonSec],muonSec+1)
        # self.hist["2DsectorChRMSForMuonCandidateSectors"].Fill(zrms[muonSec],muonSec+1)
        # if isRandom:
        #     self.hist["2DcountChannelsAboveNoiseForMuonCandidateSector_RndEvt"].Fill(len(listAllChannelsAboveNoise[muonSec]),muonSec+1)
        #     self.hist["2DsectorChMeanForMuonCandidateSectors_RndEvt"].Fill(zmean[muonSec],muonSec+1)
        #     self.hist["2DsectorChRMSForMuonCandidateSectors_RndEvt"].Fill(zrms[muonSec],muonSec+1)

        self.hist["EventCount"].Fill("delta sigma cut",1)


        #a selected sector have a signal above the noise threshold and remaining sectors are below a threshold.excluded the active sector + 2 neigbours + "the mirrored-sectors"
        isec_pn = ((muonSec + 1) + 16)%16
        isec_mn = ((muonSec - 1) + 16)%16

        isec_pns = ((muonSec + 1) + 8)%16
        isec_mns = ((muonSec- 1) + 8)%16
        #print "neighbors plus and minus:" , isec_pn, isec_mn, "mirrored neighbors plus and minus:" , isec_pns, isec_mns
        cond1 =  muonSec or isec_pn or isec_mn
        cond2= ((muonSec+8)%16 or isec_pns or isec_mns)
        for isec in xrange(0,16):
            cond1 =  isec == muonSec or isec == isec_pn or isec == isec_mn
            cond2 = (isec == (muonSec+8)%16 or isec == isec_pns or isec == isec_mns)
            for imod in xrange(0,14):
                hnoise_neighbor = 'CastorNoise_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                if not cond1 or not cond2:
                    self.hist[hnoise_neighbor].Fill(ch_energy[isec][imod])




        #selection on channels above noise
        # listChannelsAboveNoise = []
        # for i in xrange(0,224):
        #     isec = self.fChain.CastorRecHitSector.at(i)-1
        #     imod = self.fChain.CastorRecHitModule.at(i)-1hanme
        #     if [isec+1,imod+1] in badChannelsSecMod:
        #         #print "skipping channel", imod, isec
        #         continue

        #     if ch_energy[isec][imod]> self.ch_mean[isec][imod] + 2*self.ch_RMS[isec][imod]:
        #         self.hist["2DMuonCountMap_allch"].Fill(imod+1,isec+1)
        #         if isec == muonSec:
        #             listChannelsAboveNoise.append(i)
        #     # hchannel = 'Channelenergy_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
        #     # self.hist[hchannel].Fill(ch_energy[isec][imod])
        countChannelsAboveNoise = len(listAllChannelsAboveNoise[muonSec])
        #print "Channels above noise:", listChannelsAboveNoise


        NchannelNoiseCut = countChannelsAboveNoise > 3 # before we used the countChannelsAboveNoise > 3
        # if muonSec in [6,7,10,11,12,13]:
        #     if countChannelsAboveNoise > 4:
        #         NchannelNoiseCut = True
        
        if isRandom:
           self.hist["EventCount"].Fill("RndTrg ch cut",1)
        if hasMuonTrigger:
           self.hist["EventCount"].Fill("MuTrg ch cut",1)
        if not isBptxminus:
            self.hist["EventCount"].Fill("Bptx(+) ch cut",1)




        if countChannelsAboveNoise > 6:
            goodMuonEventDifferentSelection = True

        Front_Module = False
        Mid_Module   = False
        Rear_Module  = False

        if NchannelNoiseCut:
            self.hist["EventCount"].Fill("N_ch noise cut",1)
            
            for imod in xrange(0,14):
                self.hist["2DMuonCountFor3Ch"].Fill(imod+1,muonSec+1)

            for imod in listAllChannelsAboveNoise[muonSec]:
                # imod = self.fChain.CastorRecHitModule.at(i)-1
                if imod <= 3:# [0,1,2,3]
                   Front_Module = True
                   self.hist["2DMuonactiveregion"].Fill(1,muonSec+1)  
                elif imod <= 8: #[4,5,6,7,8]
                     Mid_Module = True
                     self.hist["2DMuonactiveregion"].Fill(2,muonSec+1)  
                else: # [9,10,11,12,13]:
                      Rear_Module = True
                      self.hist["2DMuonactiveregion"].Fill(3,muonSec+1) 

            if (Front_Module + Mid_Module + Rear_Module) >= 3: #maybe change to 3
                goodMuonEventWithoutAnyTriggerSelection = True
                if hasMuonTrigger:
                   goodMuonEvent = True

            if muonSec == 7 or muonSec == 6:
                if (Front_Module + Mid_Module + Rear_Module) >= 2: #maybe change to 3
                    goodMuonEventWithoutAnyTriggerSelection = True
                    if hasMuonTrigger:
                        goodMuonEvent = True
               
        #found an interesting event. now fill histograms for channels above noise
        if goodMuonEvent:
            # print "Good event in (sec,mod)", sec, mod, "Front,Mid,Back", Front_Module, Mid_Module, Rear_Module
            self.hist["EventCount"].Fill("good muon evt",1)
            if isRandom:
               self.hist["EventCount"].Fill("RndtrggoodMuonEvent cut",1)
            if hasMuonTrigger:
               self.hist["EventCount"].Fill("MuTrg goodMuonEvent cut",1)

            if not isBptxminus:
               self.hist["EventCount"].Fill("Bptx(+) goodMuonEvent cut",1)



         
            #Energy of muon
             
            energy_secsum = [0.0] * 16
            for isec in xrange(16):
                for imod in xrange(14):
                    if [isec+1,imod+1] in badChannelsSecMod:
                        continue

                    energy_secsum[isec] += ch_energy[isec][imod]

            henergy = 'MuonSignalSec_energy_sec_{sec}'.format(sec=str(muonSec+1))
            hnameAllsec ='MuonSignalAllSec_energy'
            self.hist[henergy].Fill(energy_secsum[muonSec])

            self.hist[hnameAllsec].Fill(energy_secsum[muonSec])



            self.hist["RunsWithGoodMuons"].Fill(self.fChain.run)


            for imod in xrange(0,14):
                hname = 'MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(muonSec+1))
                self.hist[hname].Fill(ch_energy[muonSec][imod])
            self.hist["GoodMuonCountPerSec"].Fill(muonSec+1)
            for imod in listAllChannelsAboveNoise[muonSec]:
                self.hist["2DMuonCountMap"].Fill(imod+1,muonSec+1)

            #self.tree["nt"].Fill()
        else:
            self.hist["RunsAllTrigger"].Fill(self.fChain.run)


        if goodMuonEventWithoutAnyTriggerSelection:
            self.hist["EventCount"].Fill("good muon NO trg",1)

            # print "Good event in (sec,mod)", sec, mod, "Front,Mid,Back", Front_Module, Mid_Module, Rear_Module
            for imod in xrange(0,14):
                hname = 'MuonSignalSecCh_withoutrigger_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(muonSec+1))
                self.hist[hname].Fill(ch_energy[muonSec][imod])
           
        if goodMuonEventDifferentSelection:
            for imod in xrange(0,14):
                hname = 'MuonSignalSecCh_differentSelection_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(muonSec+1))
                self.hist[hname].Fill(ch_energy[muonSec][imod])




        #plot all the good muons
        # if goodMuonEvent:
        #     new_dir = join(outfolder,"goodMuonEvents_{n:04d}/".format(n=self.maxFileNo+1))
        #     if not os.path.exists(new_dir):
        #         os.makedirs(new_dir)

        #     trgEvtFileName  = "Muon_run_" + str(self.fChain.run) + "_event_" + str(self.fChain.event) + ".pdf"
        #     fig = plt.figure()
        #     ax = fig.add_subplot(111, projection='3d')
        #     ch_energy_wo_bad_ch = np.asarray(ch_energy)
        #     for [isec, imod] in badChannelsSecMod:
        #         ch_energy_wo_bad_ch[isec-1][imod-1] = 0
        #     for isec in xrange(16):
        #         xs = np.arange(14)
        #         ys = ch_energy_wo_bad_ch[isec]
        #         ax.bar(xs, ys, zs=isec, zdir='y', color="g" if isec == muonSec else "r", alpha=0.8)
        #     ax.set_xlabel('Module')
        #     ax.set_ylabel('Sector')
        #     ax.set_zlabel('Channel Energy')
        #     ax.set_title('Sec: {sec} ({sig:.1f}sigma)   Trigger: {trg}   Mod[F/M/R]: {f}/{m}/{r} \n Run: {run}   Event: {evt}'.format(sec=muonSec, sig=listSigmaSector[muonSec], trg="Yes" if hasMuonTrigger else "no", f=Front_Module, m=Mid_Module, r=Rear_Module, run=self.fChain.run, evt=self.fChain.event), multialignment='center')
        #     fig.savefig(join(join(outfolder,new_dir),trgEvtFileName))
        #     plt.close("all")


                # hwithouttrigger = 'MuonSignalSecCh_withoutrigger_mod_{mod}_sec_{sec}'.format(mod=str(mod), sec=str(sec))
                #     self.hist[hwithouttrigger].Fill(ch_energy[sec][mod])
                #     self.hist["2DMuonNoTriggerCountMap"].Fill(mod,sec)


        return 1

    def finalize(self):
        print "Finalize:"
        # normFactor = self.getNormalizationFactor()
        # print "  applying norm", normFactor
        # for h in self.hist:
        #     self.hist[h].Scale(normFactor)

        
                # bin=  hist_merge_MuonSignal_Noise.FindBin(imod, isec)
                # Muon_Noise_signal = hname + hNoise
                #hist_merge_MuonSignal_Noise.SetBinContent(bin,Muon_Noise_signal)

                # From 16 channels in one module, only 10 are filled into histogram. Scale to get per channel noise...
                # self.hist[hnoise_neighbor].Scale(1./10.) # hnoise_neighbor.Scale(1./10.)AttributeError: 'str' object has no attribute 'Scale'

        #self.hist[hname].Scale(1./self.hist[hname].Integral());

        # self.hist[hname].Scale(1./self.hist[hname].GetBinWidth(i));


        # if self.nEventsRandom > 0:
        #     for i in xrange(0,16):
        #             print "Bin ", i, ": ", self.new_sec_mean[i], self.new_sec_RMS[i] / float(self.nEventsRandom)
        #             self.new_sec_mean[i] /= float(self.nEventsRandom)
        #             self.new_sec_RMS[i] = sqrt(self.new_sec_RMS[i]/float(self.nEventsRandom) - self.new_sec_mean[i]**2)
        #             self.hist['hist_sec_Mean'].Fill(i, self.new_sec_mean[i] )
        #             self.hist['hist_sec_RMS'].Fill(i, self.new_sec_RMS[i] )

        #     for isec in xrange(0,16):
        #         for imod in xrange(0,14):
        #             if self.nEventsRandom>0:
        #                 i = isec * 14 + imod
        #                 self.new_ch_mean[isec][imod] /= float(self.nEventsRandom)
        #                 self.new_ch_RMS[isec][imod] = sqrt(self.new_ch_RMS[isec][imod]/float(self.nEventsRandom) - self.new_ch_mean[isec][imod]**2)
        #                 self.hist['hist_ch_Mean'].Fill(i, self.new_ch_mean[isec][imod] )
        #                 self.hist['hist_ch_RMS'].Fill(i, self.new_ch_RMS[isec][imod] )
               #ftree.Write()            
    def finalizeWhenMerged(self):
        olist = self.GetOutputList()
        histos = {}
        for o in olist:
            if not "TH1" in o.ClassName():
                if not "TH2" in o.ClassName():
                    continue
            histos[o.GetName()] = o
            # print " TH1/2 histogram in output: ", o.GetName()

        for isec in xrange(0,16):
            for imod in xrange(0,14):
                hnoise_neighbor = 'CastorNoise_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                hnoise_randomtrg = 'CastorNoise_randomtrg_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))

                binnumber = histos["EventCount"].GetXaxis().FindBin("rnd trg")
                if not histos["EventCount"].GetBinContent(binnumber) == 0:
                    histos[hnoise_randomtrg].Scale( 1./histos["EventCount"].GetBinContent(binnumber) )

                #  From 16 channels in one module, only 10 are filled into histogram
                if not histos[hnoise_neighbor].GetEntries() == 0:
                    histos[hnoise_neighbor].Scale( 1./histos[hnoise_neighbor].GetEntries() )

                hname = 'MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                if not histos["GoodMuonCountPerSec"].GetBinContent(isec+1) == 0:
                    histos[hname].Scale( 1./histos["GoodMuonCountPerSec"].GetBinContent(isec+1) )

        inputFile = ROOT.TFile(join(outfolder,"mean_rms_2.root"))
        hist_ch_Mean = inputFile.Get("data_PPMinBias_Run2013/hist_ch_Mean")
        histcalibrationname = '2DMuonSignalMap'
        histcalibrationKatname = '2DMuonSignalMap_katerina'
        histcalibrationIgorname = '2DMuonSignalMap_igor'
        histcalibrationRationame = "2DMuonSignalMap_Ratio"
        histcalibrationRatioIgorname = "2DMuonSignalMapIgor_Ratio"
        histcalibrationRatioIgor_Katname = "2DMuonSignalMapIgor_KatRatio"
        histcalibrationname_rms = '2DMuonSignalMap_rms'
        histcalibration_notdividedRefname = '2DMuonSignalMap_notdividedRef'
        histcalibration = histos[histcalibrationname]
        histCalibrationKat = histos[histcalibrationKatname]
        histCalibrationIgor = histos[histcalibrationIgorname]
        histcalibrationRatio = histos[histcalibrationRationame]
        histcalibrationRatioIgor_Kat = histos[histcalibrationRatioIgor_Katname]
        histcalibrationRatioIgor = histos[histcalibrationRatioIgorname]
        histcalibration_rms = histos[histcalibrationname_rms]
        histcalibration_notdividedRef = histos[histcalibration_notdividedRefname]
        hreferencename = 'MuonSignalSecCh_mod_4_sec_9' #reference channel is (counting from one) mod=4 sec=9
        referenceMean= histos[hreferencename].GetMean()
        referenceRMS = histos[hreferencename].GetRMS()
        referenceNmuons = histos[hreferencename].GetEntries()
        meanReferenceNoise = hist_ch_Mean.GetBinContent(8*14 + 3 +1)
        # referenceMean -= meanReferenceNoise #switch on at some point or da a fit

        def getKaterinacalibrationnumbers():
            calib=[[  0.7510,    0.8700,    2.7370,    0.0000,    0.3630,    0.6430,    0.0000,    0.3100,    0.2120,    0.2740,    0.3030,    0.1690,    0.2650,    0.1550], 
                  [  0.6190,    0.6160,    1.8130,    0.8690,    0.1820,    0.6280,    0.0000,    0.5070,    0.1680,    0.2910,    0.3380,    0.1460,    0.2490,    0.1250], 
                  [  1.0700,    0.6510,    1.4250,    0.7660,    0.3040,    0.1930,    8.2170,   13.2900,    0.4650,    0.2350,    0.0000,    0.2950,    0.3430,    0.3510], 
                  [  0.5310,    0.3300,    0.8910,    0.8260,    0.1170,    0.3300,    0.0000,    0.0000,    0.0000,    0.6390,    0.0000,    0.3170,    0.0000,    0.4580], 
                  [  0.6120,    0.0000,    1.3410,    0.7020,    0.1560,    0.5690,    0.8360,    0.0000,    0.0000,    0.5230,    0.2360,    0.3290,    0.3990,    0.3420], 
                  [  1.3130,    0.4870,    1.4000,    0.6320,    0.1990,    0.7950,    1.2090,    0.0000,    0.5100,    0.7060,    0.2330,    0.2800,    0.4830,    0.4410], 
                  [  0.4160,    0.2820,    1.0430,    0.3130,    0.1140,    0.0860,  250.6690,    0.1950,    0.4200,    6.9160,    3.4790,    1.5110,    4.8590,    3.5340], 
                  [  0.3420,    0.2950,    1.1980,    1.4030,    0.2130,    1.0730,    0.0000,    0.2060,    0.6350,   27.2690,    9.4210,    3.3400,    3.4880,    1.0100], 
                  [  0.3030,    0.8460,    1.4120,    1.0000,    0.2180,    0.8830,    0.0000,    0.1320,    0.1950,    0.2490,    0.2250,    0.2270,    0.2990,    0.2780], 
                  [  0.9040,    1.4030,    2.6580,    1.1900,    0.2350,    1.5570,    0.0000,    0.3160,    0.1990,    0.3100,    0.1790,    0.2510,    0.2510,    0.2520], 
                  [  1.0160,    0.9930,    1.6950,    0.8870,    0.2850,    0.6230,    0.0000,   10.0790,    0.3730,    0.2440,    9.6350,    0.5240,    0.6990,    0.3790], 
                  [  1.1690,    1.1300,    2.1400,    1.3920,    0.2630,    1.2470,    0.0000,    0.0000,    0.5670,    0.3030,   99.3510,    0.3510,    0.1980,    0.3560], 
                  [  0.9160,    1.2700,    1.6430,    0.8070,    0.2310,    2.3020,    0.0000,    0.0000,    0.3230,    0.2910,    0.0000,    0.3430,    0.1280,    0.3080], 
                  [  0.6010,    0.9840,    2.1400,    0.8210,    0.1770,    1.0970,    0.0000,    0.0000,    0.2030,    0.2920,   16.6350,    0.3020,    0.3510,    0.3680], 
                  [  0.7590,    1.3650,    2.9620,    1.1740,    0.3800,    2.3370,    0.0000,  517.2540,    0.2690,    0.0000,    0.1940,    0.2740,    0.2800,    0.4100], 
                  [  0.7420,    0.9720,    2.4600,    0.9240,    0.2200,    0.1630,    3.9070,    0.1970,    0.2700,    0.2580,    0.1510,    0.1340,    0.2790,    0.2620]]
            # calib2 = list(reversed(calib))
            # fig = plt.figure(figsize=(15, 8))
            # ax = fig.add_subplot(111)
            # map2d = plt.imshow(calib2,interpolation="none")
            # ax.set_aspect('equal')
            # cbar = plt.colorbar(orientation='vertical',norm=mpl.colors.Normalize(vmin=0, vmax=8))
            # map2d.set_clim(0, 8.0)
            
            for sec,line in enumerate(calib):
                for mod,z in enumerate(line):
                    #print "hist->SetBinContent({x},{y},{z});".format(x=mod+1,y=sec+1,z=z)
                    histCalibrationKat.SetBinContent(mod+1,sec+1,z);

        def getIgorcalibrationnumbers():
            calib_igor= [[1.000425,  0.931195,  1.722751,  0.0,        0.249538,    0.394333,   0.0,   0.155133,    0.123307,   0.264286,    0.217681,   0.125848,   0.234006,  0.075464], 
                        [1.28372,    1.08159,   1.974339,  0.851123,   0.352999,    0.377796,   0.0,   0.246074,    0.191197,   0.263945,    0.289805,   0.120486,   0.298206,  0.076201],
                        [1.252213,   0.925163,  1.754341,  0.803017,   0.209004,    0.424808,   0.0,   0.0, 0.0,    0.221043,   0.0,         0.249512,   0.269058,   0.15739], 
                        [1.100644,   0.966124,  1.63636,   1.019983,   0.264205,    0.492828,   0.0,   0.0, 0.0,    0.240881,   0.0,         0.168351,   0.0,        0.128583], 
                        [1.000826,   0.0,       1.293415,  0.86196,    0.295809,    0.475899,   0.0,   0.0, 0.0,    0.314342,   0.20608,     0.164808,   0.271375,   0.156374],
                        [0.722059,   0.780461,  1.995963,  0.661209,   0.259085,    0.487429,   0.0,   0.0, 0.0,    0.377144,   0.26507,     0.206327,   0.336299,   0.199795], 
                        [1.096988,   1.218608,  1.745835,  0.888802,   0.225927,    0.45836,    0.0,   0.020096,    0.142001,   0.385257,    0.256415,   0.189381,   0.376582,   0.299507],
                        [1.091334,   1.025778,  2.106153,  1.199047,   0.337598,    0.58613,    0.0,   0.0,         0.154015,   0.424145,    0.292679,   0.204659,   0.404868,   0.113908], 
                        [1.054333,   1.474295,  1.124023,  1.0,        0.315433,    1.090999,   0.0,   0.166044,    0.154629,   0.215142,    0.154745,   0.230574,   0.0,        0.188459], 
                        [1.282433,   1.126571,  2.271284,  1.148607,   0.315623,    1.124756,   0.0,   0.076184,    0.150419,   0.244338,    0.149832,   0.202471,   0.158535,   0.221797], 
                        [1.207038,   1.494127,  1.223583,  0.99518,    0.282358,    0.342161,   0.0,   0.0,         0.158282,   0.160647,    0.317039,   0.414463,   0.45517,    0.183079],
                        [1.634185,   1.708066,  2.259616,  1.518867,   0.310965,    0.90365,    0.0,   0.0,         0.104076,   0.191959,    0.281285,   0.261916,   0.120645,   0.158001], 
                        [1.528451,   1.548834,  2.212892,  1.048122,   0.276433,    0.0,        0.0,   0.0,         0.072358,   0.188554,    0.254051,   0.237827,   0.082658,   0.146762], 
                        [1.335749,   1.727109,  2.040581,  1.309911,   0.245352,    1.458605,   0.0,   0.0,         0.10903,    0.278077,    0.327009,   0.203752,   0.221086,   0.162093], 
                        [1.215644,   1.620683,  2.06946,   1.029579,   0.269562,    1.893122,   0.0,   0.0,         0.176633,   0.0,         0.146783,   0.183976,   0.177801,   0.265278],
                        [1.282735,   1.439203,  1.96242,   0.842085,   0.197566,    0.0,        0.0,   0.117696,    0.18644,    0.191907,    0.127682,   0.134317,   0.216069,   0.14071]]

            for isec,line in enumerate(calib_igor):
                for imod,z in enumerate(line):
                    #print "hist->SetBinContent({x},{y},{z});".format(x=mod+1,y=sec+1,z=z)
                    histCalibrationIgor.SetBinContent(imod+1,isec+1,z);
                  

        getKaterinacalibrationnumbers()
        getIgorcalibrationnumbers()
        
        sigma_refCh = 0
        if referenceNmuons > 0:
            sigma_refCh = referenceRMS/sqrt(referenceNmuons)

        if referenceMean != 0:
            print "Warning reference channel (Sec 9 Mod 4) empty. Not dividing"
        for isec in xrange(0,16):
            for imod in xrange(0,14):
                hname = 'MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                hnoise_neighbor = 'CastorNoise_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                mean= histos[hname].GetMean()
                rms = histos[hname].GetRMS()
                

                Nmuons = histos[hname].GetEntries()
                sigma_mean = 0
                if Nmuons > 0:
                    sigma_mean = rms/sqrt(Nmuons)


                i = isec * 14 + imod
                meanNoise = hist_ch_Mean.GetBinContent(i+1) #noise is not estimated well. do iterative procedure from ralf? or random trigger?
                binnumber = histcalibration.FindBin(imod+1, isec+1)
                noiseSubtractedMean = (mean)# - meanNoise)
                histcalibration_notdividedRef.SetBinContent(binnumber, noiseSubtractedMean)
                
                mean_Kat = histCalibrationKat.GetBinContent(binnumber)
                mean_Igor = histCalibrationIgor.GetBinContent(binnumber)
                # histMuonSignalMean.SetBinContent(i+1,mean)
                # histMuonSignalMean.SetBinError(i+1,rms)
                
                sigma_mean_diff_refCh = 0
                if referenceMean != 0:
                    noiseSubtractedMean /= referenceMean #divide by reference channel
                    if mean > 0 and referenceMean > 0:
                        sigma_mean_diff_refCh = sqrt(noiseSubtractedMean**2 * ((sigma_mean/mean)**2 + (sigma_refCh/referenceMean)**2))
                   

                #if [isec+1,imod+1] in badChannelsSecMod:
                #   noiseSubtractedMean = 1 #set ba channels always to 1
                
                inv_mean = 0
                sigma_inv_mean = 0
                if noiseSubtractedMean>0:
                    inv_mean = 1./noiseSubtractedMean #invert to get a "correction factor"
                    sigma_inv_mean = sigma_mean_diff_refCh/(noiseSubtractedMean**2)
                if imod>=0 and imod<2:
                    inv_mean /= 2 #em modules are half the size
                    sigma_inv_mean /= 2
                
                meanRatio = 0
                sigma_meanRatio = 0
                meanRatio_igor=0
                meanRatioIgor_Kat=0
                if mean_Kat>0 and inv_mean>0 and mean_Igor>0 :
                    meanRatio = inv_mean/mean_Kat
                    meanRatio_igor = inv_mean/mean_Igor
                    meanRatioIgor_Kat= mean_Kat/mean_Igor
                    sigma_meanRatio = sqrt( meanRatio**2 * ( (sigma_inv_mean/inv_mean)**2 + 0.2**2 ) )
                
                    pull_Meli_Kat= (inv_mean -mean_Kat)/sqrt((mean_Kat*0.2)**2 + (sigma_inv_mean)**2)
                    pull_Meli_Igor= (inv_mean -mean_Igor)/sqrt((mean_Igor*0.2)**2 + (sigma_inv_mean)**2)
                    pull_Kat_Igor= (mean_Kat -mean_Igor)/sqrt((mean_Igor*0.2)**2 + (mean_Kat*0.2)**2)
                
                    histcalibration.SetBinContent(binnumber, inv_mean)
                    histcalibration_rms.SetBinContent(binnumber,sigma_inv_mean)
                    histcalibrationRatio.SetBinContent(binnumber, meanRatio)
                    histcalibrationRatioIgor.SetBinContent(binnumber, meanRatio_igor)
                    histcalibrationRatioIgor_Kat.SetBinContent(binnumber, meanRatio_igor)
                    histMuonSignalMean = histos["1DMuonsignalMean"] 
                    histMuonSignalKatMean= histos["1DMuonsignalKaterinaMean"]
                    histMuonSignalRatioMean = histos["1DMuonsignalRatioMean"] 
                    histMuonSignalRatioLog = histos["1DMuonsignalRatioLog"] 
                    histMuonSignalRatioMeanPull = histos["1DMuonsignalRatioMeanPull"] 
                    histMuonSignalPull_Meli_Kat = histos["1DMuonsignalPull_Meli_Kat"] 
                    histMuonSignalPull_Meli_Igor = histos["1DMuonsignalPull_Meli_Igor"] 
                    histMuonSignalPull_Kat_Igor = histos["1DMuonsignalPull_Kat_Igor"] 
                    histMuonSignalKatMean.SetBinContent(i+1,mean_Kat)
                    histMuonSignalKatMean.SetBinError(i+1,mean_Kat*0.2)
                    histMuonSignalMean.SetBinContent(i+1,inv_mean)
                    histMuonSignalMean.SetBinError(i+1,sigma_inv_mean)
                    histMuonSignalRatioMean.SetBinContent(i+1,meanRatio)
                    histMuonSignalRatioMean.SetBinError(i+1,sigma_meanRatio)
                    histMuonSignalRatioMeanPull.SetBinContent(i+1,pull_Meli_Kat)
                    histMuonSignalPull_Meli_Kat.Fill(pull_Meli_Kat)
                    histMuonSignalPull_Meli_Igor.Fill(pull_Meli_Igor)
                    histMuonSignalPull_Kat_Igor.Fill(pull_Kat_Igor)
                if meanRatio > 0:
                    histMuonSignalRatioLog.Fill(log10(meanRatio))
                else:
                    print "Warning mean ratio for channel ID", i, "is <= 0. Cannot calculate log10"
                #histMuonSignalRatioMean.SetBinError(i+1,sigma_meanRatio)

                print "checking means for muons", imod, isec, mean, histos[hname].GetEntries()
         
        

if __name__ == "__main__":
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.AutoLibraryLoader.enable()

    # sampleList = None # run through all

    # debug config:
    # Run printTTree.py alone to get the samples list
    sampleList = []
    sampleList.append("data_PPMinBias_Run2013")


    slaveParams = {}

    #check which output files already exist
    maxFileNo = -1
    filenames = [ f for f in listdir(outfolder) if isfile(join(outfolder,f)) ]
    listOfOutputFiles = []
    for ifilename in filenames:
        if not "plotsMuonselectioncuts" in ifilename:
            continue
        else:
            print "Found previous output file with name: ", ifilename
            maxFileNo = max(maxFileNo,int(ifilename.split("_")[1].strip(".root")))
            listOfOutputFiles.append(join(outfolder,ifilename))

    #prompt for deletion of previous outputfiles
    if maxFileNo != -1:
        print "\nDo you want to delete previous output filenames and start from beginning? (y/n)"
        case = ""
        while not (case == "y" or case == "n"):
            case = raw_input(">> ").lower()
            if case == "y":
                maxFileNo = -1
                for ifilename in listOfOutputFiles:
                    if os.path.exists(ifilename):
                        print "Deleting file:", ifilename
                        os.remove(ifilename)
                break
            if case == "n":
                break

    #decide on output filename
    

    outFileName = join(outfolder,"plotsMuonselectioncuts_2_{n:04d}.root".format(n=maxFileNo+1))
    if maxFileNo == -1:
        print "No previous output file found"
    else:
        print "Found prevous output file(s). Setting new output file name to:", outFileName

    # use printTTree.py <sampleName> to see what trees are avaliable inside the skim file

    slaveParams["maxFileNo"] = maxFileNo
    Step2_Selection_2.runAll(treeName="MuonCastorVTwo",
                           slaveParameters = slaveParams,
                           sampleList = sampleList,
                           maxFilesMC = None,
                           maxFilesData =None,
                           nWorkers =8,
                           outFile = outFileName,
                           verbosity = 2)
