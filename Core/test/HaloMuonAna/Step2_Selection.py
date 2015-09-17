#!/usr/bin/env python
import CommonFSQFramework.Core.ExampleProofReader
from BadChannels2015 import badChannelsSecMod

import sys, os, time
sys.path.append(os.path.dirname(__file__))
from os import listdir
from os.path import isfile, join

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import edm
from math import sqrt
from array import *
import copy

from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

outfolder = os.environ["HaloMuonOutput"]

class Step2_Selection(CommonFSQFramework.Core.ExampleProofReader.ExampleProofReader):
    def init( self, maxEvents = None):
        #####
        self.flag_use_merjin_electronic_channel_noise = True
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
            inputFile = ROOT.TFile(join(outfolder,"mean_rms.root"))
            inputFile_noise = ROOT.TFile(join(outfolder,"Histograms_StatFitInformationPedNoise.root"))
        else:
            inputFile = ROOT.TFile(join(outfolder,"plotsMuonselectioncuts_{n:04d}.root".format(n=self.maxFileNo)))
            # check naming if rerunning step1 again
        
        
       # create the branches and assign the fill-variables to them

        # ftree = ROOT.TFile("tree.root", "recreate")
        # tree = ROOT.TTree("nt", "nt")
        # np= np.zeros(1, dtype=float)
        # tree.Branch('ch_energy', ch_energy, 'ch_energy')

        hist_sec_Mean = inputFile.Get("data_MinimumBias_Run2015A/hist_sec_Mean")
        hist_sec_RMS = inputFile.Get("data_MinimumBias_Run2015A/hist_sec_RMS")

        #Getting channel RMS and Mean
        hist_ch_Mean = inputFile.Get("data_MinimumBias_Run2015A/hist_ch_Mean")
        hist_ch_RMS =  inputFile.Get("data_MinimumBias_Run2015A/hist_ch_RMS")

        #Getting electronic noise from in the Merijns file
        hist_ch_noise_RMS= inputFile_noise.Get("SumHistoNoiseFit")

        #hist_sec_noise_RMS= inputFile_noise.Get("SumHistoNoiseFit")
        self.hist["2DMuonCountMap"] =  ROOT.TH2D("2DMuonCountMap","2DMuonCountMap", 14, 0.5, 14.5, 16, 0.5, 16.5)
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
        histcalibration_notdividedRefname = '2DMuonSignalMap_notdividedRef'
        histratioLEDname = '2DMuonSignalMap_ratio_Mel_LED'

        self.hist[histcalibration_notdividedRefname] =  ROOT.TH2D(histcalibration_notdividedRefname,histcalibration_notdividedRefname, 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist[histratioLEDname] =  ROOT.TH2D(histratioLEDname,histratioLEDname, 14, 0.5, 14.5, 16, 0.5, 16.5)
        histLED1800Vname = '2DMuonSignalMapLED_1800V2D'
        histLEDerror1800Vname = '2DMuonSignalMapLEDerror_1800V2D'
        self.hist[histLED1800Vname] =  ROOT.TH2D(histLED1800Vname,histLED1800Vname, 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist[histLEDerror1800Vname] =  ROOT.TH2D(histLEDerror1800Vname,histLEDerror1800Vname, 14, 0.5, 14.5, 16, 0.5, 16.5)
        histGainPPHVname = 'GainPPHV_to_MuonHV'
        histGainErrorPPHVname = 'GainErrorPPHV_to_MuonHV'
        self.hist[histGainPPHVname] =  ROOT.TH2D(histGainPPHVname,histGainPPHVname, 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist[histGainErrorPPHVname] =  ROOT.TH2D(histGainErrorPPHVname,histGainErrorPPHVname, 14, 0.5, 14.5, 16, 0.5, 16.5)

        histMuonSignalPull_Meli_Gain= "1DMuonsignalPull_Meli_Gain"
        self.hist[histMuonSignalPull_Meli_Gain] = ROOT.TH1D(histMuonSignalPull_Meli_Gain, histMuonSignalPull_Meli_Gain,40,-4,4)
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

        #histcalibration_notdividedRefname = '2DMuonSignalMap_notdividedRef'
        #histcalibration_notdividedRef = self.hist[histcalibration_notdividedRefname]

        ch_energy = [[0.0 for _ in xrange(14)] for _ in xrange(16)]
        sec_energy = [0.0] * 16
        sec_front_energy = [0.0] * 16

        if self.fChain.CastorRecHitEnergy.size() != 224:
            return 0

        self.hist["EventCount"].Fill("size rh",1)

        isRandom = False
        if self.fChain.trgRandom:
            if not (self.fChain.trgl1L1GTTech[1] or self.fChain.trgl1L1GTTech[2]): #not bptx+ or bptx-
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


        NchannelNoiseCut = countChannelsAboveNoise > 3 # before we used the countChannelsAboveNoise > 5
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

            for imod in listAllChannelsAboveNoise[muonSec]:
                # imod = self.fChain.CastorRecHitModule.at(i)-1
                if imod <= 3:# [0,1,2,3]
                    Front_Module = True
                elif imod <= 8: #[4,5,6,7,8]
                    Mid_Module = True
                else: # [9,10,11,12,13]:
                    Rear_Module = True

            if (Front_Module + Mid_Module + Rear_Module) >= 3: #maybe change to 3
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
        if goodMuonEvent:
            new_dir = join(outfolder,"goodMuonEvents_{n:04d}/".format(n=self.maxFileNo+1))
            if not os.path.exists(new_dir):
                os.makedirs(new_dir)

            trgEvtFileName  = "Muon_run_" + str(self.fChain.run) + "_event_" + str(self.fChain.event) + ".pdf"
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ch_energy_wo_bad_ch = np.asarray(ch_energy)
            for [isec, imod] in badChannelsSecMod:
                ch_energy_wo_bad_ch[isec-1][imod-1] = 0
            for isec in xrange(16):
                xs = np.arange(14)
                ys = ch_energy_wo_bad_ch[isec]
                ax.bar(xs, ys, zs=isec, zdir='y', color="g" if isec == muonSec else "r", alpha=0.8)
            ax.set_xlabel('Module')
            ax.set_ylabel('Sector')
            ax.set_zlabel('Channel Energy')
            ax.set_title('Sec: {sec} ({sig:.1f}sigma)   Trigger: {trg}   Mod[F/M/R]: {f}/{m}/{r} \n Run: {run}   Event: {evt}'.format(sec=muonSec, sig=listSigmaSector[muonSec], trg="Yes" if hasMuonTrigger else "no", f=Front_Module, m=Mid_Module, r=Rear_Module, run=self.fChain.run, evt=self.fChain.event), multialignment='center')
            fig.savefig(join(join(outfolder,new_dir),trgEvtFileName))
            plt.close("all")


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


        if self.nEventsRandom > 0:
            for i in xrange(0,16):
                    print "Bin ", i, ": ", self.new_sec_mean[i], self.new_sec_mean[i] / self.nEventsRandom
                    self.new_sec_mean[i] /= float(self.nEventsRandom)
                    self.new_sec_RMS[i] = sqrt(self.new_sec_RMS[i]/float(self.nEventsRandom) - self.new_sec_mean[i]**2)
                    self.hist['hist_sec_Mean'].Fill(i, self.new_sec_mean[i] )
                    self.hist['hist_sec_RMS'].Fill(i, self.new_sec_RMS[i] )

            for isec in xrange(0,16):
                for imod in xrange(0,14):
                    if self.nEventsRandom>0:
                        i = isec * 14 + imod
                        self.new_ch_mean[isec][imod] /= float(self.nEventsRandom)
                        self.new_ch_RMS[isec][imod] = sqrt(self.new_ch_RMS[isec][imod]/float(self.nEventsRandom) - self.new_ch_mean[isec][imod]**2)
                        self.hist['hist_ch_Mean'].Fill(i, self.new_ch_mean[isec][imod] )
                        self.hist['hist_ch_RMS'].Fill(i, self.new_ch_RMS[isec][imod] )
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

        inputFile = ROOT.TFile(join(outfolder,"mean_rms.root"))
        hist_ch_Mean = inputFile.Get("data_MinimumBias_Run2015A/hist_ch_Mean")
        histcalibrationname = '2DMuonSignalMap'
      
        histcalibration_notdividedRefname = '2DMuonSignalMap_notdividedRef'
        histratioLEDname = '2DMuonSignalMap_ratio_Mel_LED'
        histLED1800Vname = '2DMuonSignalMapLED_1800V2D'
        histLEDerror1800Vname = '2DMuonSignalMapLEDerror_1800V2D'
        histGainPPHVname = 'GainPPHV_to_MuonHV'
        histGainErrorPPHVname = 'GainErrorPPHV_to_MuonHV'
        histcalibration = histos[histcalibrationname]
        histcalibration_notdividedRef = histos[histcalibration_notdividedRefname]
        histcalibrationLED=histos[histratioLEDname]
        hreferencename = 'MuonSignalSecCh_mod_4_sec_9' #reference channel is (counting from one) mod=4 sec=9
        referenceMean= histos[hreferencename].GetMean()
        referenceNmuons = histos[hreferencename].GetEntries()
        referenceRMS = histos[hreferencename].GetRMS()
        histLED1800V = histos[histLED1800Vname]
        histLEDerror1800V = histos[histLEDerror1800Vname]
        histGainPPHV= histos[histGainPPHVname]
        histGainErrorPPHV= histos[histGainErrorPPHVname]
        
        meanReferenceNoise = hist_ch_Mean.GetBinContent(8*14 + 3 +1)
        


        # referenceMean -= meanReferenceNoise #switch on at some point or da a fit
        
        def getGainValues_LED_1800V():
            gain= [[75500.0, 1300000.0, 568000.0, 79100.0, 0.0, 362000.0, 330000.0, 512000.0, 2130000.0, 121000.0, 133000.0, 1760000.0, 104000.0, 1560000.0],
                  [72200.0, 530000.0, 344000.0, 0.0, 793000.0, 270000.0, 47400.0, 467000.0, 1140000.0, 0.0, 55200.0, 2460000.0, 611000.0, 525000.0], 
                  [362000.0, 557000.0, 423000.0, 843000.0, 320000.0, 274000.0, 382000.0, 547000.0, 457000.0, 271000.0, 683000.0, 515000.0, 369000.0, 462000.0], 
                  [579000.0, 522000.0, 377000.0, 746000.0, 293000.0, 248000.0, 297000.0, 463000.0, 1080000.0, 477000.0, 542000.0, 647000.0, 788000.0, 561000.0],
                  [354000.0, 94100.0, 326000.0, 0.0, 0.0, 307000.0, 486000.0, 456000.0, 325000.0, 0.0, 1570000.0, 393000.0, 365000.0, 638000.0],
                  [644000.0, 0.0, 0.0, 0.0, 310000.0, 282000.0, 226000.0, 486000.0, 463000.0, 420000.0, 0.0, 160000.0, 346000.0, 624000.0], 
                  [535000.0, 486000.0, 414000.0, 703000.0, 285000.0, 83900.0, 397000.0, 661000.0, 832000.0, 311000.0, 423000.0, 748000.0, 309000.0, 320000.0], 
                  [645000.0, 509000.0, 406000.0, 501000.0, 335000.0, 279000.0, 483000.0, 676000.0, 757000.0, 278000.0, 503000.0, 647000.0, 374000.0, 857000.0], 
                  [0.0, 172000.0, 310000.0, 543000.0, 213000.0, 159000.0, 448000.0, 569000.0, 974000.0, 710000.0, 734000.0, 296000.0, 0.0, 403000.0], 
                  [551000.0, 550000.0, 298000.0, 583000.0, 263000.0, 142000.0, 568000.0, 581000.0, 846000.0, 613000.0, 279000.0, 320000.0, 291000.0, 385000.0], 
                  [327000.0, 515000.0, 303000.0, 506000.0, 160000.0, 140000.0, 506000.0, 532000.0, 756000.0, 920000.0, 392000.0, 324000.0, 288000.0, 304000.0], 
                  [356000.0, 214000.0, 359000.0, 624000.0, 214000.0, 247000.0, 524000.0, 48600.0, 717000.0, 732000.0, 533000.0, 464000.0, 856000.0, 664000.0], 
                  [483000.0, 564000.0, 0.0, 469000.0, 349000.0, 164000.0, 522000.0, 558000.0, 355000.0, 596000.0, 442000.0, 472000.0, 377000.0, 789000.0],
                  [573000.0, 431000.0, 293000.0, 401000.0, 346000.0, 180000.0, 435000.0, 0.0, 473000.0, 556000.0, 440000.0, 573000.0, 453000.0, 0.0],
                  [441000.0, 426000.0, 363000.0, 0.0, 345000.0, 122000.0, 459000.0, 436000.0, 219000.0, 663000.0, 314000.0, 631000.0, 517000.0, 363000.0],
                  [536000.0, 487000.0, 354000.0, 0.0, 375000.0, 212000.0, 274000.0, 483000.0, 504000.0, 806000.0, 844000.0, 825000.0, 506000.0, 552000.0]]


            for isec,line in enumerate(gain):
                for imod,z in enumerate(line):
                    #print "hist->SetBinContent({x},{y},{z});".format(x=mod+1,y=sec+1,z=z)
                    histLED1800V.SetBinContent(imod+1,isec+1,z);
                  

        
        getGainValues_LED_1800V()
        


        def getgainError_LED_1800V():
            # gain_error = [[1,2,3],[4,5,6]]
            gain_error = [[0.297684, 0.303952, 0.154135, 0.143911, 0.00635945, 0.00696305, 0.013834, 0.0321267, 0.0131778, 0.00933203, 0.050539, 0.202872, 0.050539, 0.050539], 
                          [0.273256, 0.160228, 0.139015, 0.313545, 0.00771067, 0.00761597, 0.0624354, 0.0248926, 0.0202406, 0.00548391, 0.00676534, 0.234323, 0.050539, 0.050539], 
                          [0.208304, 0.209137, 0.20524, 0.190448, 0.0232386, 0.00471184, 0.0277715, 0.0154876, 0.008732, 0.0134014, 0.030463, 0.0378628, 0.0626422, 0.0364017], 
                          [0.542134, 0.393696, 0.132817, 0.242715, 0.0507952, 0.13058, 0.0520935, 0.0290168, 0.0121874, 0.0140919, 0.0130064, 0.0335709, 0.026466, 0.0332873], 
                          [0.27, 0.160646, 0.16222, 0.195976, 0.00751524, 0.0067034, 0.0115408, 0.0185129, 0.00929156, 0.0129714, 0.0114801, 0.00506028, 1.0658, 0.606506], 
                          [0.154156, 0.153823, 0.119356, 0.202521, 0.00716295, 0.0070185, 0.0120852, 0.0464483, 0.00874689, 0.0116634, 0.050539, 0.151539, 0.377385, 0.410735], 
                          [0.184748, 0.160822, 0.393225, 0.178964, 0.00593741, 0.00470421, 0.017487, 0.0232038, 0.00917938, 0.0107265, 0.0160861, 0.356823, 0.0939621, 0.0386872], 
                          [0.157523, 0.175195, 0.16801, 0.105009, 0.0290363, 0.00498382, 0.0212775, 0.022056, 0.0157854, 0.0104235, 0.194334, 0.206459, 0.261539, 0.0521384], 
                          [0.151579, 0.143932, 0.135856, 0.171023, 0.253775, 0.0120092, 0.0237112, 0.0268974, 0.0142285, 0.0211197, 0.067731, 0.050539, 0.050539, 0.336312], 
                          [0.163748, 0.161818, 0.181251, 0.253084, 0.00782319, 0.00626064, 0.0318485, 0.0249402, 0.00646081, 0.0111822, 0.050539, 0.189923, 0.050539, 0.653284], 
                          [0.238548, 0.188285, 0.12227, 0.172699, 0.00682779, 0.00786426, 0.0209686, 0.0204322, 0.0101018, 0.0140415, 0.221228, 0.0962536, 0.365371, 0.050539], 
                          [0.166892, 0.243934, 0.206652, 0.140395, 0.0051994, 0.00605701, 0.0256196, 0.0165193, 0.00994767, 0.0158936, 0.132698, 0.125007, 0.374697, 0.385424], 
                          [0.115767, 0.173322, 0.126637, 0.12721, 0.00871158, 0.00421985, 0.0232185, 0.0220267, 0.00926294, 0.0144162, 0.0723647, 0.171277, 0.0925839, 0.208289], 
                          [0.168528, 0.173265, 0.146743, 0.202867, 0.00571115, 0.00731292, 0.0177893, 0.0232495, 0.00932107, 0.00953935, 0.22871, 0.108933, 0.174756, 0.050539], 
                          [0.183391, 0.137414, 0.124799, 0.168779, 0.0485851, 0.00782375, 0.0183694, 0.0220257, 0.00824439, 0.00902181, 0.0328126, 0.138288, 0.112013, 0.173584], 
                          [0.257366, 0.167887, 0.272659, 0.184186, 0.00944811, 0.0108569, 0.0476646, 0.0224154, 0.00838952, 0.0083331, 0.13013, 0.0684296, 0.156254, 0.13299]]

            for isec, line in enumerate(gain_error):
                for imod, z in enumerate(line):
                    #print "hist->SetBinContent({x},{y},{z});".format(x=mod+1,y=sec+1,z=z)
                    histLEDerror1800V.SetBinContent(imod+1,isec+1,z);
                  



        getgainError_LED_1800V()
        def getgain_PPHV_to_MuonHV():
            gain=   [[17.8879, 13.2482, 14.0743, 13.3278, 2.12091, 2.12543, 3.21527, 2.98481, 2.11992, 2.06513, 2.10157, 4.18052, 2.10157, 2.10157],
                    [15.5813, 14.3688, 13.9962, 15.8174, 2.14671, 2.16158, 2.94851, 3.03287, 2.07367, 2.09709, 2.02886, 5.6682, 2.10157, 2.10157],
                    [17.1485, 14.2758, 13.6352, 13.7265, 2.19081, 2.16161, 3.12191, 3.12444, 2.11671, 2.14295, 1.68227, 2.62755, 2.90551, 2.25366], 
                    [25.4937, 15.0353, 13.5396, 14.5396, 2.00951, 2.01052, 2.94567, 3.10935, 2.10436, 2.10581, 1.9966, 2.14238, 1.67304, 2.10112], 
                    [22.6664, 14.0919, 14.2022, 13.6725, 2.14257, 2.13443, 3.07885, 3.15034, 2.12721, 2.111, 2.12661, 2.04572, 2.1345, 2.23077], 
                    [13.8189, 13.0366, 12.7392, 14.4154, 2.14965, 2.13784, 2.99889, 3.00984, 2.09463, 2.09944, 2.10157, 0.507937, 2.24675, 2.14433],
                    [14.1375, 14.3402, 12.6179, 14.1134, 2.14431, 2.11088, 2.99728, 3.23413, 2.10955, 2.10675, 2.0835, 6.92593, 1.99355, 2.31884], 
                    [13.8083, 13.8735, 13.0352, 13.8626, 2.00958, 2.13991, 3.18717, 3.15437, 2.09151, 2.11572, 2.21586, 2.13531, 2.51007, 2.03081],
                    [15.0498, 14.8564, 12.6001, 13.3532, 2.00998, 2.16025, 3.17426, 3.14376, 2.0985, 2.09938, 2.12139, 2.10157, 2.10157, 1.72961],
                    [13.4029, 14.6153, 13.4451, 13.2316, 2.12537, 2.11262, 3.21792, 3.06312, 2.16952, 2.08848, 2.10157, 1.15108, 2.10157, 1.96429], 
                    [13.8487, 14.6547, 13.9845, 14.004, 2.16241, 1.98043, 3.01777, 3.02263, 2.12504, 2.16387, 1.9798, 2.0, 2.21538, 2.10157], 
                    [13.5454, 14.1516, 11.5071, 14.1275, 2.08736, 2.11662, 3.24957, 3.10665, 2.13867, 2.04457, 2.22083, 2.04405, 1.98608, 1.98802], 
                    [13.3621, 14.021, 13.4167, 11.0612, 2.15244, 2.09452, 3.13056, 3.04276, 2.04811, 2.10842, 2.125, 3.80645, 1.75349, 5.33108],
                    [13.4709, 14.0436, 13.4366, 14.2121, 2.12454, 2.12266, 2.99075, 3.13996, 2.06635, 2.07508, 2.31579, 2.23828, 2.27638, 2.10157],
                    [13.2229, 13.5576, 13.6674, 14.794, 2.00953, 2.02293, 3.12361, 3.20957, 2.08386, 2.10041, 0.850949, 2.31136, 4.57522, 2.21341],
                    [13.4984, 13.0597, 13.3461, 13.7108, 2.1543, 2.14047, 2.89888, 3.05443, 2.07617, 2.12138, 2.25668, 2.10459, 2.3318, 2.07519]]
        # calib2 = list(reversed(calib))
        # fig = plt.figure(figsize=(15, 8)) 
        # ax = fig.add_subplot(111)
        # map2d = plt.imshow(calib2,interpolation="none")
        # ax.set_aspect('equal')
        # cbar = plt.colorbar(orientation='vertical',norm=mpl.colors.Normalize(vmin=0, vmax=8))
        # map2d.set_clim(0, 8.0)
        
            for sec,line in enumerate(gain):
                for mod,z in enumerate(line):
                    #print "hist->SetBinContent({x},{y},{z});".format(x=mod+1,y=sec+1,z=z)
                    histGainPPHV.SetBinContent(mod+1,sec+1,z);
             

        getgain_PPHV_to_MuonHV()


        def getgain_PPHVError_to_MuonHV():
            gain_error =    [[0.297684, 0.303952, 0.154135, 0.143911, 0.00635945, 0.00696305, 0.013834, 0.0321267, 0.0131778, 0.00933203, 0.050539, 0.202872, 0.050539, 0.050539],
                            [0.273256, 0.160228, 0.139015, 0.313545, 0.00771067, 0.00761597, 0.0624354, 0.0248926, 0.0202406, 0.00548391, 0.00676534, 0.234323, 0.050539, 0.050539],
                            [0.208304, 0.209137, 0.20524, 0.190448, 0.0232386, 0.00471184, 0.0277715, 0.0154876, 0.008732, 0.0134014, 0.030463, 0.0378628, 0.0626422, 0.0364017],
                            [0.542134, 0.393696, 0.132817, 0.242715, 0.0507952, 0.13058, 0.0520935, 0.0290168, 0.0121874, 0.0140919, 0.0130064, 0.0335709, 0.026466, 0.0332873],
                            [0.27, 0.160646, 0.16222, 0.195976, 0.00751524, 0.0067034, 0.0115408, 0.0185129, 0.00929156, 0.0129714, 0.0114801, 0.00506028, 1.0658, 0.606506],
                            [0.154156, 0.153823, 0.119356, 0.202521, 0.00716295, 0.0070185, 0.0120852, 0.0464483, 0.00874689, 0.0116634, 0.050539, 0.151539, 0.377385, 0.410735],
                            [0.184748, 0.160822, 0.393225, 0.178964, 0.00593741, 0.00470421, 0.017487, 0.0232038, 0.00917938, 0.0107265, 0.0160861, 0.356823, 0.0939621, 0.0386872],
                            [0.157523, 0.175195, 0.16801, 0.105009, 0.0290363, 0.00498382, 0.0212775, 0.022056, 0.0157854, 0.0104235, 0.194334, 0.206459, 0.261539, 0.0521384],
                            [0.151579, 0.143932, 0.135856, 0.171023, 0.253775, 0.0120092, 0.0237112, 0.0268974, 0.0142285, 0.0211197, 0.067731, 0.050539, 0.050539, 0.336312],
                            [0.163748, 0.161818, 0.181251, 0.253084, 0.00782319, 0.00626064, 0.0318485, 0.0249402, 0.00646081, 0.0111822, 0.050539, 0.189923, 0.050539, 0.653284],
                            [0.238548, 0.188285, 0.12227, 0.172699, 0.00682779, 0.00786426, 0.0209686, 0.0204322, 0.0101018, 0.0140415, 0.221228, 0.0962536, 0.365371, 0.050539],
                            [0.166892, 0.243934, 0.206652, 0.140395, 0.0051994, 0.00605701, 0.0256196, 0.0165193, 0.00994767, 0.0158936, 0.132698, 0.125007, 0.374697, 0.385424],
                            [0.115767, 0.173322, 0.126637, 0.12721, 0.00871158, 0.00421985, 0.0232185, 0.0220267, 0.00926294, 0.0144162, 0.0723647, 0.171277, 0.0925839, 0.208289],
                            [0.168528, 0.173265, 0.146743, 0.202867, 0.00571115, 0.00731292, 0.0177893, 0.0232495, 0.00932107, 0.00953935, 0.22871, 0.108933, 0.174756, 0.050539],
                            [0.183391, 0.137414, 0.124799, 0.168779, 0.0485851, 0.00782375, 0.0183694, 0.0220257, 0.00824439, 0.00902181, 0.0328126, 0.138288, 0.112013, 0.173584],
                            [0.257366, 0.167887, 0.272659, 0.184186, 0.00944811, 0.0108569, 0.0476646, 0.0224154, 0.00838952, 0.0083331, 0.13013, 0.0684296, 0.156254, 0.13299]]


            for sec,line in enumerate(gain_error):
                for mod,z in enumerate(line):
                    #print "hist->SetBinContent({x},{y},{z});".format(x=mod+1,y=sec+1,z=z)
                    histGainErrorPPHV.SetBinContent(mod+1,sec+1,z);
             

        getgain_PPHVError_to_MuonHV()

       
        sigma_refCh = 0
        if referenceNmuons > 0:
            sigma_refCh = (referenceRMS/sqrt(referenceNmuons) )*mean_gainerror_PPHV
            
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
                #binnumber_notdividedRef = histcalibration_notdividedRef.FindBin(imod+1, isec+1)
                meandividedRefmean = (mean)# - meanNoise)
                histcalibration_notdividedRef.SetBinContent(binnumber, meandividedRefmean)
                mean_gain_LED1800V = histLED1800V.GetBinContent(binnumber)
                mean_gainerror_LED1800V = histLEDerror1800V.GetBinContent(binnumber)
                mean_gain_PPHV = histGainPPHV.GetBinContent(binnumber)
                mean_gainerror_PPHV = histGainErrorPPHV.GetBinContent(binnumber)
                
                sigma_mean_PPHV = 0
                mean_gain =0
                if  referenceMean != 0:
                    
                    mean_gain *= mean_gain_PPHV
                    sigma_mean_PPHV = sigma_mean* mean_gainerror_PPHV
                    
                    referenceMean*= mean_gain_PPHV 
                    sigma_mean_gaincorrection = sqrt((mean_gain)**2 *(sigma_mean)**2 + (mean)**2 *(sigma_mean_PPHV)**2)
                    
                    meandividedRefmean = mean_gain /referenceMean
       
                   
                    if mean > 0 and referenceMean > 0:
                       sigma_mean_diff_refCh = sqrt((meandividedRefmean)**2 * (sigma_refCh)**2 + (referenceMean)**2*(sigma_mean_gaincorrection)**2)
                       
                       

                       # if imod>=0 and imod<2:
                       # inv_mean /= 2 #em modules are half the size
                       # sigma_inv_mean /= 2

                

                inv_mean = 0
                sigma_inv_mean = 0
                sigma_meanerror_PPHV= 0
                if  meandividedRefmean>0:
                    
                    inv_mean = 1./meandividedRefmean #invert to get a "correction factor"
                    
                    sigma_inv_mean = sigma_mean_diff_refCh/(meandividedRefmean**2)

                if  imod>=0 and imod<2:
                    inv_mean /= 2        #em modules are half the size
                    sigma_inv_mean /= 2
                
                
                
               
                if inv_mean>0 and mean_gain_LED1800V>0 :
                   
                    meanRatio_gain = inv_mean/mean_gain_LED1800V
                   
                    # sigma_meanRatio = sqrt( meanRatio**2 * ( (sigma_inv_mean/inv_mean)**2 + 0.2**2 ) )
                
                   
                    pull_Meli_Gain= (inv_mean- mean_gain_LED1800V)/sqrt((mean_gainerror_LED1800V)**2 + (sigma_inv_mean)**2)
                    
                    histcalibration.SetBinContent(binnumber, inv_mean)
                    histcalibrationLED.SetBinContent(binnumber, meanRatio_gain)
                
                
                    histcalibration.SetBinContent(binnumber, meandividedRefmean)
                    histMuonSignalPull_Meli_Gain = histos["1DMuonsignalPull_Meli_Gain"] 
 
               
                    histMuonSignalPull_Meli_Gain.Fill(pull_Meli_Gain)
               
                    print "checking means for muons", imod, isec, mean, histos[hname].GetEntries()
         
        

if __name__ == "__main__":
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.AutoLibraryLoader.enable()

    sampleList = None # run through all

    # debug config:
    # Run printTTree.py alone to get the samples list
    #sampleList = []
    #sampleList.append("QCD_Pt-15to3000_TuneZ2star_Flat_HFshowerLibrary_7TeV_pythia6")


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
    

    outFileName = join(outfolder,"plotsMuonselectioncuts_{n:04d}.root".format(n=maxFileNo+1))
    if maxFileNo == -1:
        print "No previous output file found"
    else:
        print "Found prevous output file(s). Setting new output file name to:", outFileName

    # use printTTree.py <sampleName> to see what trees are avaliable inside the skim file
    #max_events = 100
    slaveParams["maxFileNo"] = maxFileNo
    Step2_Selection.runAll(treeName="MuonCastorVTwo",
                           slaveParameters = slaveParams,
                           sampleList = sampleList,
                           maxFilesMC = None,
                           maxFilesData =None,
                           #maxNevents = max_events,
                           nWorkers =8,
                           outFile = outFileName,
                           verbosity = 2)
