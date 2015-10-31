#!/usr/bin/env python
DATASOURCE = "DATA" #options are "TOYMC" "DATA". first generates events, the second uses the sample specified
NBRHIST =  418 # the number of hist for toyMc
NBREVTPERHIST = 10000 #2387


import CommonFSQFramework.Core.ExampleProofReader
from BadChannels2015 import badChannelsSecMod

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
# import time

#from outsource_analzye_muon import *


from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt


if DATASOURCE == "TOYMC":
    from ToyMC import toyMCClass
    theToyMC = toyMCClass()


outfolder = os.environ["HaloMuonOutput"]




class Step2_Selection_2(CommonFSQFramework.Core.ExampleProofReader.ExampleProofReader):
    def init( self, maxEvents = None):
        #####
#        self.flag_use_merjin_electronic_channel_noise = False # mine RMS
        self.flag_use_merjin_electronic_channel_noise = True # RU CHECK: is this Gaussian fit ?? 
        #####


        # self.maxFileNo = slaveParameters["maxFileNo"]
        firstRun = (self.maxFileNo == -1)

        self.maxEvents = maxEvents
        self.hist = {}
        self.nEventsRandom = 0
        #self.tree = {}
        self.hist["EventCount"] = ROOT.TH1D("EventCount","EventCount", 10, 0, 20)

#        self.hist["RunsAllTrigger"] =  ROOT.TH1D("RunsAllTrigger", "RunsAllTrigger", 10000, 247000-0.5, 257000-0.5)
        
        #if not self.isData:

        #Getting sector RMS and Mean
        if firstRun:
            inputFile = ROOT.TFile(join(outfolder,"mean_rms.root"))
            inputFile_noise = ROOT.TFile(join(outfolder,"Histograms_StatFitInformationPedNoise.root"))
            


        else:
            inputFile = ROOT.TFile(join(outfolder,"plotsMuonselectioncuts_2_{n:04d}.root".format(n=self.maxFileNo)))
            # check naming if rerunning step1 again
        #
        



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

#        print "Pointer Address of Noise level hist:", hist_sec_Mean

        # we need a set of histogramms for different exclusivity selections around the muon
        for iDelta in xrange(self.iDeltaStart,self.iDeltaEnd):
            
            histcalibrationname = '2DMuonSignalMap_Excl' + str(iDelta)
            histcalibrationname_sigma = '2DMuonSignalMap_Excl' + str(iDelta)+'_sigma'
            histcalibrationname_RMS = '2DMuonSignalMap_Excl'+ str(iDelta)+'_RMS'
            histcalibration_notdividedRefname = '2DMuonSignalMap_Excl' + str(iDelta)+'_notdividedRef'
            
            self.hist[histcalibration_notdividedRefname] =  ROOT.TH2D(histcalibration_notdividedRefname,histcalibration_notdividedRefname, 14, 0.5, 14.5, 16, 0.5, 16.5)
            self.hist[histcalibrationname_RMS] =  ROOT.TH2D(histcalibrationname_RMS,histcalibrationname_RMS, 14, 0.5, 14.5, 16, 0.5, 16.5)
            self.hist[histcalibrationname_sigma] =  ROOT.TH2D(histcalibrationname_sigma,histcalibrationname_sigma, 14, 0.5, 14.5, 16, 0.5, 16.5)


            if not firstRun:
                self.hist[histcalibrationname] = inputFile.Get("data_MinimumBias_Run2015A/2DMuonSignalMap")
                print "Extracted histogram from file. Checking entries:",  self.hist[histcalibrationname].GetEntries()

            else: #first time running analyser the calibration constants are all set to 1
                self.hist[histcalibrationname] =  ROOT.TH2D(histcalibrationname,histcalibrationname, 14, 0.5, 14.5, 16, 0.5, 16.5)

                for imod in xrange(0,14):
                    for isec in xrange(0,16):
                        self.hist[histcalibrationname].SetBinContent(self.hist[histcalibrationname].FindBin(imod+1,isec+1), 1.) #all factors set to 1



            self.hist["EventCount_Excl" + str(iDelta)] = ROOT.TH1D("EventCount_Excl" + str(iDelta),"EventCount",10,0,20)

            self.hist["2DMuonCountMap_Excl" + str(iDelta)] =  ROOT.TH2D("2DMuonCountMap_Excl" + str(iDelta),"2DMuonCountMap", 14, 0.5, 14.5, 16, 0.5, 16.5)
            self.hist["2DMuonCountFor3Ch_Excl" + str(iDelta)] =  ROOT.TH2D("2DMuonCountFor3Ch_Excl" + str(iDelta),"2DMuonCountFor3Ch", 14, 0.5, 14.5, 16, 0.5, 16.5)
            self.hist["2DMuonactiveregion_Excl" + str(iDelta)] =  ROOT.TH2D("2DMuonactiveregion_Excl" + str(iDelta),"2DMuonactiveregion", 14, 0.5, 14.5, 16, 0.5, 16.5)

#            self.hist["2DcountChannelAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_AllEvt"] =  ROOT.TH2D("2DcountChannelAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_AllEvt","2DcountChannelAboveNoiseFor3SigmaSectors_AllEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)
            self.hist["2DcountChannelAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_MuoEvt"] =  ROOT.TH2D("2DcountChannelAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_MuoEvt","2DcountChannelAboveNoiseFor3SigmaSectors_MuoEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)
            self.hist["2DcountChannelAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_RndEvt"] =  ROOT.TH2D("2DcountChannelAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_RndEvt","2DcountChannelAboveNoiseFor3SigmaSectors_RndEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)

            self.hist["2DcountChannelsAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_MuoEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_MuoEvt","2DcountChannelsAboveNoiseFor3SigmaSectors_MuoEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
#            self.hist["2DcountChannelsAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_AllEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_AllEvt","2DcountChannelsAboveNoiseFor3SigmaSectors_AllEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
            self.hist["2DcountChannelsAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_RndEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_RndEvt","2DcountChannelsAboveNoiseFor3SigmaSectors_RndEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)

#            self.hist["2DSecRMSHot_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_AllEvt"] = ROOT.TH2D("2DSecRMSHot_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_AllEvt","2DSecRMSHot_Vs_SecRMSSecondHot_AllEvt",88, -2, 20, 88, -2, 20)
#            self.hist["2DSecRMSHot_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_MuoEvt"] = ROOT.TH2D("2DSecRMSHot_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_MuoEvt","2DSecRMSHot_Vs_SecRMSSecondHot_MuoEvt",88, -2, 20, 88, -2, 20)
#            self.hist["2DSecRMSHot_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_RndEvt"] = ROOT.TH2D("2DSecRMSHot_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_RndEvt","2DSecRMSHot_Vs_SecRMSSecondHot_RndEvt",88, -2, 20, 88, -2, 20)

            self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_MuoEvt"] = ROOT.TH2D("2DDeltaSigma_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_MuoEvt","2DDeltaSigma_Vs_SecRMSSecondHot_MuoEvt",88,-2,20, 88, -2, 20)
#            self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_AllEvt"] = ROOT.TH2D("2DDeltaSigma_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_AllEvt","2DDeltaSigma_Vs_SecRMSSecondHot_AllEvt",88,-2,20, 88, -2, 20)
            self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_RndEvt"] = ROOT.TH2D("2DDeltaSigma_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_RndEvt","2DDeltaSigma_Vs_SecRMSSecondHot_RndEvt",88, -2, 20, 88,-2, 20)

        # self.hist["2DMuonNoTriggerCountMap"] =  ROOT.TH2D("2DMuonNoTriggerCountMap","2DMuonNoTriggerCountMap", 14, -0.5, 13.5, 16, -0.5, 15.5)
            self.hist["GoodMuonCountPerSec_Excl"+str(iDelta)] =  ROOT.TH1D("GoodMuonCountPerSec_Excl"+str(iDelta),"GoodMuonCountPerSec", 16, 0.5, 16.5)
#            self.hist["RunsWithGoodMuons_Excl"+str(iDelta)] =  ROOT.TH1D("RunsWithGoodMuons_Excl"+str(iDelta),"RunsWithGoodMuons", 10000, 247000-0.5, 257000-0.5)
            self.hist["DeltaSigma_Excl"+str(iDelta)+"_MuoEvt"] = ROOT.TH1D("DeltaSigma_Excl"+str(iDelta)+"_MuoEvt","DeltaSigma_MuoEvt",100, 0, 50)
#            self.hist["DeltaSigma_Excl"+str(iDelta)+"_AllEvt"] = ROOT.TH1D("DeltaSigma_Excl"+str(iDelta)+"_AllEvt","DeltaSigma_AllEvt",100, 0, 50)
            self.hist["DeltaSigma_Excl"+str(iDelta)+"_RndEvt"] = ROOT.TH1D("DeltaSigma_Excl"+str(iDelta)+"_RndEvt","DeltaSigma_RndEvt",100, 0, 50)

            for isec in xrange(0,16):

                henergy= 'MuonSignalSec_energy_Excl'+str(iDelta)+'_sec_{sec}'.format(sec=str(isec+1))
                self.hist[henergy] = ROOT.TH1D( henergy, henergy, 50, 0, 400)      

                for imod in xrange(0,14):
                    hname = 'MuonSignalSecCh_Excl'+str(iDelta)+'_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))               
#                    hwithouttrigger = 'MuonSignalSecCh_Excl'+str(iDelta)+'_withoutrigger_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                    self.hist[hname] = ROOT.TH1D( hname, hname, 50, -100, 400)
 #                   self.hist[hwithouttrigger] = ROOT.TH1D(hwithouttrigger, hwithouttrigger, 50, -100, 400)

                    if DATASOURCE == "TOYMC":
                        hMC1D = '1DSignalMCCh_Excl'+str(iDelta)+'_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                        self.hist[hMC1D] = ROOT.TH1D(hMC1D,hMC1D,100, -10, 90)

                        # I have 10000 events and 10000/1000 i need to have 10 histograms
                        for iRep in xrange(NBRHIST): 
                            hname_MC = 'MuonSignalMCCh_Excl'+str(iDelta)+'_mod_{mod}_sec_{sec}_number_{number}'.format(mod=str(imod+1), sec=str(isec+1), number=str(iRep))
                            self.hist[hname_MC] = ROOT.TH1D( hname_MC, hname_MC, 50, -100, 400)   


            hnameAllsec= 'MuonSignalAllSec_Excl'+str(iDelta)+'_energy'
            self.hist[hnameAllsec] = ROOT.TH1D(hnameAllsec, hnameAllsec,50, 0, 400)


            histMuonSignalMean = "1DMuonsignalMean_Excl"+str(iDelta)
            self.hist[histMuonSignalMean] = ROOT.TH1D(histMuonSignalMean,histMuonSignalMean,224,-0.5,223.5)






        
        self.hist["2DcountChannelAboveNoise_MuoEvt"] =  ROOT.TH2D("2DcountChannelAboveNoise_MuoEvt","2DcountChannelAboveNoise_MuoEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)
#        self.hist["2DcountChannelAboveNoise_AllEvt"] =  ROOT.TH2D("2DcountChannelAboveNoise_AllEvt","2DcountChannelAboveNoise_AllEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist["2DcountChannelAboveNoise_RndEvt"] =  ROOT.TH2D("2DcountChannelAboveNoise_RndEvt","2DcountChannelAboveNoise_RndEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)

        self.hist["2DcountChannelsAboveNoiseForAllSectors_MuoEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseForAllSectors_MuoEvt","2DcountChannelsAboveNoiseForAllSectors_MuoEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
#        self.hist["2DcountChannelsAboveNoiseForAllSectors_AllEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseForAllSectors_AllEvt","2DcountChannelsAboveNoiseForAllSectors_AllEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
        self.hist["2DcountChannelsAboveNoiseForAllSectors_RndEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseForAllSectors_RndEvt","2DcountChannelsAboveNoiseForAllSectors_RndEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)


        self.hist["2DcountSectorRMS_MuoEvt"] =  ROOT.TH2D("2DcountSectorRMS_MuoEvt","2DcountSectorRMS_MuoEvt", 88, -2, 20, 16, 0.5, 16.5)
#        self.hist["2DcountSectorRMS_AllEvt"] =  ROOT.TH2D("2DcountSectorRMS_AllEvt","2DcountSectorRMS_AllEvt", 88, -2, 20, 16, 0.5, 16.5)
        self.hist["2DcountSectorRMS_RndEvt"] =  ROOT.TH2D("2DcountSectorRMS_RndEvt","2DcountSectorRMS_RndEvt", 88, -2, 20, 16, 0.5, 16.5)
                            


        # print "Created new 2d calibration map. Checking entries:",  self.hist[histcalibrationname].GetEntries(), ". maxFileNo =", self.maxFileNo, firstRun
        #histCalibration.SetBit(ROOT.TH1.kIsAverage) # does not seem to work. use TProfile2D instead

        # make new histograms for noise
#
 #       for isec in xrange(0,16):
  #          for imod in xrange(0,14):
   #             hnoise_electronic_RMS = 'CastorNoise_electronic_RMS_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
    #            hnoise_randomtrg = 'CastorNoise_randomtrg_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
     #           #hempty = 'CastorEmptysectors_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
      #          self.hist[hnoise_electronic_RMS] = ROOT.TH1D(hnoise_electronic_RMS,hnoise_electronic_RMS, 224,-0.5,223.5)
       #         self.hist[hnoise_randomtrg] = ROOT.TH1D(hnoise_randomtrg, hnoise_randomtrg, 50, -100, 400)


                

                
        # Getting electronic noise from in the Merijns file
        hist_ch_noise_RMS= inputFile_noise.Get("SumHistoNoiseFit")

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
#        self.new_ch_mean = [[0 for _ in xrange(14)] for _ in xrange(16)]
#        self.new_ch_RMS = [[0 for _ in xrange(14)] for _ in xrange(16)]
#        self.new_sec_mean = [0] * 16
#        self.new_sec_RMS = [0] * 16

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



    def getListSigmaSectorFromChannels(self, listSigmaChannel):
        listSigmaSector = [0] * 16

        for isec in xrange(0,16):
            sigma_this_sector = 0

            for imod in xrange(0,14):
                if listSigmaChannel[isec][imod] is not None:
                    sigma_this_sector += listSigmaChannel[isec][imod]

                listSigmaSector[isec] = sigma_this_sector

        return listSigmaSector



    # find all channel sigma values, and 
    # count channels above noise-cut in sectors
    def getListSigmaChannel(self, energy_ch, sigma_cut, badChannelsSecMod):
        listSigmaChannel = [[0 for _ in xrange(14)] for _ in xrange(16)]
        listAllChannelsAboveNoise = [[] for _ in xrange(16)]
        for isec in xrange(16):
            for imod in xrange(14):
                energy = energy_ch[isec][imod]
                # print "RMS mod {mod} sec {sec}: {rms}".format(mod=str(imod),sec=str(isec),rms=self.ch_RMS[isec][imod])
                if [isec+1,imod+1] in badChannelsSecMod:
                    listSigmaChannel[isec][imod] = None
                    #print "skipping channel", imod, isec
                    continue
                
                sigma = 0
                if self.ch_RMS[isec][imod] != 0:
                    sigma = (energy - self.ch_mean[isec][imod]) / self.ch_RMS[isec][imod]
                
                listSigmaChannel[isec][imod] = sigma
                
                if sigma > sigma_cut:
                    listAllChannelsAboveNoise[isec].append(imod)

        return listSigmaChannel, listAllChannelsAboveNoise






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
        weight = 1 # Ralf: what is this used for ? CFF internal?
        num = 0 # Ralf: what is this used for ? CFF internal?
                

        # genTracks
        #num = self.fChain.genTracks.size()
        #print num
        #print self.maxEta # see slaveParams below
        
        self.hist["EventCount"].Fill("all",1)
        #print "EventCount 1 :", self.hist["EventCount"].GetEntries()

#        self.hist["RunsAllTrigger"].Fill(self.fChain.run)


        if DATASOURCE == "DATA":
            trgl1L1GTTech1    = self.fChain.trgl1L1GTTech[1]
            trgl1L1GTTech2    = self.fChain.trgl1L1GTTech[2]
            trgl1L1GTTech7    = self.fChain.trgl1L1GTTech[7]
            trgRandom         = self.fChain.trgRandom
            trgCastorHaloMuon = self.fChain.trgCastorHaloMuon
            trgl1L1GTAlgo102  = self.fChain.trgl1L1GTAlgo[102]
#            histcalibrationname = '2DMuonSignalMap'
            histcalibrationname = '2DMuonSignalMap_Excl' + str(self.iDeltaEnd-1)

        elif DATASOURCE == "TOYMC":
            theToyMC.GenerateEvent()
            trgl1L1GTTech1    = True
            trgl1L1GTTech2    = True
            trgl1L1GTTech7    = False
            trgRandom         = False
            trgCastorHaloMuon = True
            trgl1L1GTAlgo102  = True
            
            histcalibrationname = '2DMuonSignalMap'
            #EventCount.GetBin("all") // 10000 = number for the plot 2DSignalMap_number
            #histcalibrationname = '2DMuonSignalMap_number'
            #if histos.has_key(2DSignalMap_number):
            #    pass
            #else
            #    create new 2DSign

            #if EventCount.GetBin("all") % 10000 == 0:
            #    for channels:
            #        histogram for signal.Reset()
           
            
            Nevent = 0
            #10000/1000 10000 events are produced and divided by 1000
            Nevent = int((self.hist["EventCount"].GetBinContent(self.hist["EventCount"].GetXaxis().FindBin("all"))-1)/NBREVTPERHIST) # subtracted by -1 because it is need to start from 0
            # print  "EventCount (bin \'all\') :" ,self.hist["EventCount"].GetBinContent(self.hist["EventCount"].GetXaxis().FindBin("all"))
            # print  "EventCount 2 :" ,self.hist["EventCount"].GetXaxis().FindBin("all")
            # print "Nevent" ,Nevent
            
            if Nevent >= NBRHIST: Nevent = NBRHIST-1
           
                   

        isBptxplus = False
        isBptxminus = False
        # if not (trgl1L1GTTech1 and trgl1L1GTTech2): #not (bptx+ and bptx-)
        if (trgl1L1GTTech1): # bptx PLUS
            isBptxplus = True
        if (trgl1L1GTTech2): # bptx MINUS
            isBptxminus = True
                

        isZeroBias = isBptxminus and isBptxplus

        # if self.flag_use_merjin_electronic_channel_noise:
#        for isec in xrange(16):
 #           for imod in xrange(14):
  #              hnoise_electronic_RMS = 'CastorNoise_electronic_RMS_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
   #             self.hist[hnoise_electronic_RMS].Fill(self.ch_RMS[isec][imod])

                
        histCalibration = self.hist[histcalibrationname]
        #histcalibrationname_sigma = '2DMuonSignalMap_sigma'


        ch_energy = [[0.0 for _ in xrange(14)] for _ in xrange(16)]
#        sec_energy = [0.0] * 16
#        sec_front_energy = [0.0] * 16

        if self.fChain.CastorRecHitEnergy.size() != 224:
            return 0

        self.hist["EventCount"].Fill("size rh",1)
        
        
        if isZeroBias:
            self.hist["EventCount"].Fill("ZB",1)
                     
       
        isRandom = False
        if trgRandom:
            if trgl1L1GTTech7: # not bptx+ or bptx- (bptx quit)
                self.nEventsRandom += 1
                isRandom = 1

        if isRandom:
            self.hist["EventCount"].Fill("rnd trg",1)

        
        for i in xrange(224):
            isec = i//14
            imod = i%14
            if DATASOURCE == "DATA":
                rh_energy = self.fChain.CastorRecHitEnergy.at(i)
            elif DATASOURCE == "TOYMC":
                rh_energy = theToyMC.CastorRecHitEnergy(i)

            binnumber = histCalibration.FindBin(imod+1,isec+1)
            calibration = histCalibration.GetBinContent(binnumber)

            ich_energy = calibration * rh_energy
            ch_energy[isec][imod] = ich_energy

#            if [isec+1, imod+1] not in badChannelsSecMod:
#                sec_energy[isec] += ich_energy
#                if imod < 5: sec_front_energy[isec] += ich_energy

#            if isRandom:
#                hnoise_randomtrg = 'CastorNoise_randomtrg_mod_{imod}_sec_{isec}'.format(imod=str(imod+1), isec=str(isec+1))
 #               self.hist[hnoise_randomtrg].Fill(ch_energy[isec][imod])
#                self.new_ch_mean[isec][imod] += ich_energy
#                self.new_ch_RMS[isec][imod] += ich_energy**2
#                if [isec+1, imod+1] not in badChannelsSecMod:
 #                   self.new_sec_mean[isec] += ich_energy
  #                  self.new_sec_RMS[isec] += ich_energy**2
            




        ######################
        # Setup muon trigger #
        ######################
        hasMuonTrigger = (trgCastorHaloMuon or trgl1L1GTAlgo102) # RU CHECK: why is this not the same ?????
        # if not hasMuonTrigger: return 1
        if hasMuonTrigger: self.hist["EventCount"].Fill("muon trg",1)

        if hasMuonTrigger and isBptxplus and not isBptxminus:
            self.hist["EventCount"].Fill("muon+bptx",1)
       

            #  print 'sec', sec, 'mod', mod, 'e_ch=', ch_energy, "e_sec", sec_energy[sec]
            #  print "all sector energies:  ", sec_energy


        # zmean, zrms = self.getlistSectorMeanRMSInZ(ch_energy)
        

        sigma_channel_cut = 2 # WARNING: CUT VALUE ON CHANNEL-SIGMA !!!!!!!!!!!!! (Melike: 2)
        listSigmaChannel,listAllChannelsAboveNoise = self.getListSigmaChannel(ch_energy, sigma_channel_cut, badChannelsSecMod)

        # plot channels above noise
        if len(listAllChannelsAboveNoise) is not 16:
            print ("ERROR: listAllChannelsAboveNoise has not length 16 but " + str(len(listAllChannelsAboveNoise)))
            return 0
        
        for isec in xrange(16):
            for imod in listAllChannelsAboveNoise[isec]:                
        #        self.hist["2DcountChannelAboveNoise_AllEvt"].Fill(imod+1,isec+1)
                if isRandom:
                    self.hist["2DcountChannelAboveNoise_RndEvt"].Fill(imod+1,isec+1)
                if hasMuonTrigger:
                    self.hist["2DcountChannelAboveNoise_MuoEvt"].Fill(imod+1,isec+1)

        
#        listSigmaSector = self.getListSigmaSector(sec_energy)                  # Melike
        listSigmaSector = self.getListSigmaSectorFromChannels(listSigmaChannel) # ->  Ralf



        # print sector noise distributions
        for isec in xrange(16):
            sigma = listSigmaSector[isec]

      #      self.hist["2DcountSectorRMS_AllEvt"].Fill(sigma,isec+1)
            if hasMuonTrigger:
                self.hist["2DcountSectorRMS_MuoEvt"].Fill(sigma,isec+1)
            if isRandom:
                self.hist["2DcountSectorRMS_RndEvt"].Fill(sigma,isec+1)



        # printout of number of channels above noise
        for isec in xrange(0,16):
#            self.hist["2DcountChannelsAboveNoiseForAllSectors_AllEvt"].Fill(len(listAllChannelsAboveNoise[isec]),isec+1)
            # self.hist["2DsectorChMeanForAllSectors"].Fill(zmean[isec],isec+1)
            # self.hist["2DsectorChRMSForAllSectors"].Fill(zrms[isec],isec+1)
            # self.hist["2DsectorChRMS_Vs_ChAboveNoise"].Fill(len(listAllChannelsAboveNoise[isec]),zrms[isec])
            if hasMuonTrigger:
                self.hist["2DcountChannelsAboveNoiseForAllSectors_MuoEvt"].Fill(len(listAllChannelsAboveNoise[isec]),isec+1)
            if isRandom:
                self.hist["2DcountChannelsAboveNoiseForAllSectors_RndEvt"].Fill(len(listAllChannelsAboveNoise[isec]),isec+1)
                # self.hist["2DsectorChMeanForAllSectors_RndEvt"].Fill(zmean[isec],isec+1)
                # self.hist["2DsectorChRMSForAllSectors_RndEvt"].Fill(zrms[isec],isec+1)
                # self.hist["2DsectorChRMS_Vs_ChAboveNoise_RndEvt"].Fill(len(listAllChannelsAboveNoise[isec]),zrms[isec])




        ###############################################
        # the main muon selection loop starts here
        
        for iSectorTest in xrange(16):
            for iDeltaSector in xrange(self.iDeltaStart,self.iDeltaEnd): # exclusivity sectors next to muon: 1 ... 8
                
                HottestSector = iSectorTest
                SigmaHottestSector = listSigmaSector[iSectorTest]
                SecondHottestSector = None
                SigmaSecHottestSector = None
                
                for iSectorTestNeighbors in xrange(1, iDeltaSector+1):
                    iSectorTestNeighborsLeft = (iSectorTest+iSectorTestNeighbors) % 16
                    iSectorTestNeighborsRight = (iSectorTest-iSectorTestNeighbors) % 16
                    
                    if SigmaHottestSector < listSigmaSector[iSectorTestNeighborsRight] or SigmaHottestSector < listSigmaSector[iSectorTestNeighborsLeft] :
                        HottestSector = None
                        break # this sector cannot be a muon candidate for iDeltaSector (or larger)
                    
                    if SigmaSecHottestSector is None: # initialize
                        if listSigmaSector[iSectorTestNeighborsLeft] > listSigmaSector[iSectorTestNeighborsRight]:
                            SigmaSecHottestSector = listSigmaSector[iSectorTestNeighborsLeft]
                            SecondHottestSector = iSectorTestNeighborsLeft
                        else :
                            SigmaSecHottestSector = listSigmaSector[iSectorTestNeighborsRight]
                            SecondHottestSector = iSectorTestNeighborsRight
                    
                    else :
                        if listSigmaSector[iSectorTestNeighborsLeft] > SigmaSecHottestSector:
                            SigmaSecHottestSector = listSigmaSector[iSectorTestNeighborsLeft]
                            SecondHottestSector = iSectorTestNeighborsLeft

                        if listSigmaSector[iSectorTestNeighborsRight] > SigmaSecHottestSector:
                            SigmaSecHottestSector = listSigmaSector[iSectorTestNeighborsRight]
                            SecondHottestSector = iSectorTestNeighborsRight
                                        
                if HottestSector is None:
                    break # this sector cannot be a muon candidate for iDeltaSector (or larger) -> next sector


                ################################################################################
                # we found a muon candidate in iSectorTest with exclusivity iDeltaSector !!!
                                
                DeltaSigma = (SigmaHottestSector - SigmaSecHottestSector)
#                self.hist["2DSecRMSHot_Vs_SecRMSSecondHot_Excl" + str(iDeltaSector) + "_AllEvt"].Fill(SigmaHottestSector,SigmaSecHottestSector)
 #               self.hist["DeltaSigma_Excl" + str(iDeltaSector) + "_AllEvt"].Fill(DeltaSigma)
  #              self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_Excl" + str(iDeltaSector) + "_AllEvt"].Fill(DeltaSigma,SigmaSecHottestSector)
                if hasMuonTrigger:
#                    self.hist["2DSecRMSHot_Vs_SecRMSSecondHot_Excl" + str(iDeltaSector) + "_MuoEvt"].Fill(SigmaHottestSector,SigmaSecHottestSector)
                    self.hist["DeltaSigma_Excl" + str(iDeltaSector) + "_MuoEvt"].Fill(DeltaSigma)
                    self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_Excl" + str(iDeltaSector) + "_MuoEvt"].Fill(DeltaSigma,SigmaSecHottestSector)
                if isRandom:
#                    self.hist["2DSecRMSHot_Vs_SecRMSSecondHot_Excl" + str(iDeltaSector) + "_RndEvt"].Fill(SigmaHottestSector,SigmaSecHottestSector)
                    self.hist["DeltaSigma_Excl" + str(iDeltaSector) + "_RndEvt"].Fill(DeltaSigma)
                    self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_Excl" + str(iDeltaSector) + "_RndEvt"].Fill(DeltaSigma,SigmaSecHottestSector)



#                goodMuonEvent = False
                goodMuonEventWithoutAnyTriggerSelection = False

                ################################
                # cut on muon candidate sector #
                ################################
                if SigmaSecHottestSector > 2.5: # WARNING: CUT VALUE ON NOISE HARDCODED HERE !!!!
#                    return 0
                    continue

                muonSec = HottestSector
        
                self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("sec-cut (all)",1)


                # for muon candidate sectors, make some plots
                for imod in xrange(14):
                    sigma = listSigmaChannel[muonSec][imod]
                    if sigma is not None: # check for bad channel
                        #print "skipping channel", mod, sec
                        continue
                    if sigma > sigma_channel_cut:  
                        # listAllChannelsAboveNoise[muonSec].append(imod)
#                        self.hist["2DcountChannelAboveNoiseFor3SigmaSectors_Excl" + str(iDeltaSector) + "_AllEvt"].Fill(imod+1,muonSec+1)
                        if hasMuonTrigger:
                            self.hist["2DcountChannelAboveNoiseFor3SigmaSectors_Excl" + str(iDeltaSector) + "_MuoEvt"].Fill(imod+1,muonSec+1)
                        if isRandom:
                            self.hist["2DcountChannelAboveNoiseFor3SigmaSectors_Excl" + str(iDeltaSector) + "_RndEvt"].Fill(imod+1,muonSec+1)


#                self.hist["2DcountChannelsAboveNoiseFor3SigmaSectors_Excl" + str(iDeltaSector) + "_AllEvt"].Fill(len(listAllChannelsAboveNoise[muonSec]),muonSec+1)
                # self.hist["2DsectorChMeanFor3sigmaSectors"].Fill(zmean[muonSec],muonSec+1)
                # self.hist["2DsectorChRMSFor3sigmaSectors"].Fill(zrms[muonSec],muonSec+1)
                if hasMuonTrigger:
                    self.hist["2DcountChannelsAboveNoiseFor3SigmaSectors_Excl" + str(iDeltaSector) + "_MuoEvt"].Fill(len(listAllChannelsAboveNoise[muonSec]),muonSec+1)
                    self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("sec-cut (muo)",1)
                if isRandom:
                    self.hist["2DcountChannelsAboveNoiseFor3SigmaSectors_Excl" + str(iDeltaSector) + "_RndEvt"].Fill(len(listAllChannelsAboveNoise[muonSec]),muonSec+1)
                    self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("sec-cut (rnd)",1)
                    # self.hist["2DsectorChMeanFor3sigmaSectors_RndEvt"].Fill(zmean[muonSec],muonSec+1)
                    # self.hist["2DsectorChRMSFor3sigmaSectors_RndEvt"].Fill(zrms[muonSec],muonSec+1)

                if hasMuonTrigger and isBptxplus and not isBptxminus:
                   self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("sec-cut (muo+bptx)",1)


        
#                self.hist["EventCount"].Fill("delta sigma cut",1)


                countChannelsAboveNoise = len(listAllChannelsAboveNoise[muonSec])
                # print "Channels above noise:", listChannelsAboveNoise
        
        
                NchannelNoiseCut = countChannelsAboveNoise > 3 # before we used the countChannelsAboveNoise > 3
                # if muonSec in [6,7,10,11,12,13]:
                #     if countChannelsAboveNoise > 4:
                #         NchannelNoiseCut = True
                    
            
        
            
                Front_Module = 0
                Mid_Module   = 0
                Rear_Module  = 0

                if NchannelNoiseCut:

                    self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("N_ch-cut (all)",1)
    
                    if isRandom:
                        self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("N_ch-cut (rnd)",1)
                    if hasMuonTrigger:
                        self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("N_ch-cut (muo)",1)
                    if hasMuonTrigger and isBptxplus and not isBptxminus:
                        self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("N_ch-cut (muo+bptx)",1)


                    for imod in xrange(0,14):
                        self.hist["2DMuonCountFor3Ch_Excl" + str(iDeltaSector)].Fill(imod+1,muonSec+1)

                    for imod in listAllChannelsAboveNoise[muonSec]:
                        # imod = self.fChain.CastorRecHitModule.at(i)-1
                        if imod <= 3:# [0,1,2,3]
                           Front_Module = 1
                           self.hist["2DMuonactiveregion_Excl" + str(iDeltaSector)].Fill(1,muonSec+1)  
                        elif imod <= 8: #[4,5,6,7,8]
                             Mid_Module = 1
                             self.hist["2DMuonactiveregion_Excl" + str(iDeltaSector)].Fill(2,muonSec+1)  
                        else: # [9,10,11,12,13]:
                              Rear_Module = 1
                              self.hist["2DMuonactiveregion_Excl" + str(iDeltaSector)].Fill(3,muonSec+1) 


                    if (Front_Module + Mid_Module + Rear_Module) >= 3: #maybe change to 3
                        goodMuonEventWithoutAnyTriggerSelection = True
#                        if hasMuonTrigger:
#                           goodMuonEvent = True

                    if muonSec == 7 or muonSec == 6:
                        if (Front_Module + Mid_Module + Rear_Module) >= 2: #maybe change to 3
                            goodMuonEventWithoutAnyTriggerSelection = True
#                            if hasMuonTrigger:
#                                goodMuonEvent = True

                #found an interesting event. now fill histograms for channels above noise
                if goodMuonEventWithoutAnyTriggerSelection:
# RU was: if goodMuonEvent:                    
                    # print "Good event in (sec,mod)", sec, mod, "Front,Mid,Back", Front_Module, Mid_Module, Rear_Module


#                    self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("good muon (all)",1)

                    if isRandom:
                       self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("good muon (rnd)",1)

                    if hasMuonTrigger:
                       self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("good muon (muo)",1)

                    if hasMuonTrigger and isBptxplus and not isBptxminus:
                       self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("good muon (muo+bptx)",1)

                    if isZeroBias:
                       self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("good muon (ZB)",1)                        

                    if isZeroBias and hasMuonTrigger:
                       self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("good muon (ZB+muon)",1)                        




                    # Energy of muons

#                    for imod in xrange(0,14):
 #                       hname = 'MuonSignalSecCh_Excl' + str(iDeltaSector) + '_withoutrigger_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(muonSec+1))
  #                      self.hist[hname].Fill(ch_energy[muonSec][imod])



                    if hasMuonTrigger:
                       
                        energy_secsum = [0.0] * 16
                        for isec in xrange(16):
                            for imod in xrange(14):
                                if [isec+1,imod+1] in badChannelsSecMod:
                                    continue
                                energy_secsum[isec] += ch_energy[isec][imod]
                                
                                
                        henergy = 'MuonSignalSec_energy_Excl' + str(iDeltaSector) + '_sec_{sec}'.format(sec=str(muonSec+1))
                        self.hist[henergy].Fill(energy_secsum[muonSec])


                        if DATASOURCE == "TOYMC":
                            for imod in xrange(14):
                                hname_MC = 'MuonSignalMCCh_Excl' + str(iDeltaSector) + '_mod_{mod}_sec_{sec}_number_{number}'.format(mod=str(imod+1), sec=str(muonSec+1), number=str(Nevent))   
                                self.hist[hname_MC].Fill(ch_energy[muonSec][imod])
                                # mean_mc= self.hist[hname_MC].GetMean()


                        hnameAllsec ='MuonSignalAllSec_Excl' + str(iDeltaSector) + '_energy'
                        self.hist[hnameAllsec].Fill(energy_secsum[muonSec])


#                        self.hist["RunsWithGoodMuons_Excl" + str(iDeltaSector)].Fill(self.fChain.run)

                        for imod in xrange(0,14):
                            hname = 'MuonSignalSecCh_Excl' + str(iDeltaSector) + '_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(muonSec+1))
                            self.hist[hname].Fill(ch_energy[muonSec][imod])

                        self.hist["GoodMuonCountPerSec_Excl" + str(iDeltaSector)].Fill(muonSec+1)

                        for imod in listAllChannelsAboveNoise[muonSec]:
                            self.hist["2DMuonCountMap_Excl" + str(iDeltaSector)].Fill(imod+1,muonSec+1)

                        #self.tree["nt"].Fill()


        return 1





    def finalize(self):
        print "Finalize:"




    def finalizeWhenMerged(self):
        olist = self.GetOutputList()
        histos = {}
        for o in olist:
            if not "TH1" in o.ClassName():
                if not "TH2" in o.ClassName():
                    continue
            histos[o.GetName()] = o

        #     print " TH1/2 histogram in output: ", o.GetName()
        # print "finsihed printing olist of", self
          

        for isec in xrange(0,16):
            for imod in xrange(0,14):

#                hnoise_randomtrg = 'CastorNoise_randomtrg_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))  
 #               binnumber = histos["EventCount"].GetXaxis().FindBin("rnd trg")
 #               if not histos["EventCount"].GetBinContent(binnumber) == 0:
  #                  histos[hnoise_randomtrg].Scale( 1./histos["EventCount"].GetBinContent(binnumber) )

                # scale muon energy distributions by number of muons
                for iDeltaSector in xrange(self.iDeltaStart,self.iDeltaEnd):
                    hname = 'MuonSignalSecCh_Excl' + str(iDeltaSector) + '_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                
                    if not histos["GoodMuonCountPerSec_Excl" + str(iDeltaSector)].GetBinContent(isec+1) == 0:
                        histos[hname].Scale( 1./histos["GoodMuonCountPerSec_Excl"+str(iDeltaSector)].GetBinContent(isec+1) )

                



#       inputFile = ROOT.TFile(join(outfolder,"mean_rms.root"))
#       hist_ch_Mean = inputFile.Get("data_MinimumBias_Run2015A/hist_ch_Mean")

        for iDeltaSector in xrange(self.iDeltaStart,self.iDeltaEnd):

            histcalibration = histos['2DMuonSignalMap_Excl'+str(iDeltaSector)]        
            histcalibration_sigma = histos['2DMuonSignalMap_Excl'+str(iDeltaSector)+'_sigma']
            histcalibration_RMS= histos['2DMuonSignalMap_Excl'+str(iDeltaSector)+'_RMS']
            histcalibration_notdividedRef = histos['2DMuonSignalMap_Excl'+str(iDeltaSector)+'_notdividedRef']
            
            for isec in xrange(0,16):
                for imod in xrange(0,14):
                    i = isec * 14 + imod
                    #meanNoise = hist_ch_Mean.GetBinContent(i+1) #noise is not estimated well. do iterative procedure from ralf? or random trigger?
                    binnumber = histcalibration.FindBin(imod+1, isec+1)

                    if  DATASOURCE == "DATA":
                        hname = 'MuonSignalSecCh_Excl'+str(iDeltaSector)+'_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                        Nmuons = histos[hname].GetEntries()
                        mean= histos[hname].GetMean()
                        rms = histos[hname].GetRMS()
                        sigma_mean = 0
                        if Nmuons > 0 :
                           sigma_mean = rms/sqrt(Nmuons)

                        hreferencename = 'MuonSignalSecCh_Excl'+str(iDeltaSector)+'_mod_4_sec_9' #reference channel is (counting from one) mod=4 sec=9



                    elif DATASOURCE == "TOYMC":
                          for ihist in xrange(NBRHIST):
                              hname_MC = 'MuonSignalMCCh_Excl'+str(iDeltaSector)+'_mod_{mod}_sec_{sec}_number_{number}'.format(mod=str(imod+1), sec=str(isec+1), number=str(ihist))    
                              hMC1D = '1DSignalMCCh_Excl'+str(iDeltaSector)+'_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                              mean_mc= histos[hname_MC].GetMean()
                              histos[hMC1D].Fill(mean_mc)
                              Nmuons = histos[hMC1D].GetEntries()
                              mean= histos[hMC1D].GetMean()       
                              rms= histos[hMC1D].GetRMS()
                              
                              hreferencename = '1DSignalMCCh_Excl'+str(iDeltaSector)+'_mod_4_sec_9' #reference channel is (counting from one) mod=4 sec=9
                              sigma_mean = 0
                              if Nmuons > 0 :
                                 sigma_mean = rms/sqrt(Nmuons)


                    referenceMean= histos[hreferencename].GetMean()
                    referenceRMS = histos[hreferencename].GetRMS()
                    referenceNmuons = histos[hreferencename].GetEntries()

                    print  "1:referenceRMS, referenceMean", referenceMean,referenceRMS
                    if referenceNmuons > 0 :
                        sigma_refCh = referenceRMS/sqrt(referenceNmuons)
                    if  referenceMean == 0:
                        print "Warning reference channel (Sec 9 Mod 4) empty. Not dividing"         



                    # and sigma_mean",  sigma_refCh, sigma_inv_mean
                    # print "Muon Signal mean Ch_{sec}_{mod} = {mval}".format(sec=str(isec+1),mod=str(imod+1),mval=str(mean))
                    print  "2 : referenceRMS, referenceMean", referenceRMS, referenceMean                 

                    histcalibration_notdividedRef.SetBinContent(binnumber,mean)
                    histcalibration_RMS.SetBinContent(binnumber,rms)

                    MeanDividedRefMean = mean
                    sigma_mean_diff_refCh = 0


                    if  referenceMean != 0 :
                        MeanDividedRefMean /= referenceMean #divide by reference channel
                        print "sigma_refCh", sigma_refCh
                        if mean > 0 and referenceMean > 0 :
                           sigma_mean_diff_refCh = sqrt(MeanDividedRefMean**2 * ((sigma_mean/mean)**2 + (sigma_refCh/referenceMean)**2))

                    inv_mean = 0
                    sigma_inv_mean = 0
                    if  MeanDividedRefMean>0 :
                        inv_mean = 1./MeanDividedRefMean #invert to get a "correction factor"
                        sigma_inv_mean = sigma_mean_diff_refCh/(MeanDividedRefMean**2)


                    if  imod>=0 and imod<2:
                        inv_mean /= 2 #em modules are half the size
                        sigma_inv_mean /= 2


                    print "inv_mean and sigma_mean", inv_mean, sigma_inv_mean

                    histcalibration.SetBinContent(binnumber, inv_mean)
                    histcalibration_sigma.SetBinContent(binnumber,sigma_inv_mean)

                    histMuonSignalMean = histos["1DMuonsignalMean_Excl"+str(iDeltaSector)] 
                    histMuonSignalMean.SetBinContent(i+1,inv_mean)
                    histMuonSignalMean.SetBinError(i+1,sigma_inv_mean)


                    # you use the value in the following at next step which is pull distribution 
                    # if  inv_TOYMCmean>0 and sigma_inv_TOYMCmean>0:    
                    #     pull_TOYMC_Melike= (inv_mean -inv_TOYMCmean)/sqrt((sigma_inv_TOYMCmean)**2 + (sigma_inv_mean)**2)

                    # if  meanRatio > 0:
                    #     histMuonSignalRatioLog.Fill(log10(meanRatio))
               


        

if __name__ == "__main__":

    if len(sys.argv) != 2 :
        print ("ERROR: SPECIFY iDeltaStart as command line argument !!!! Try again.  ")
        sys.exit(1)

    iDeltaStart = int(sys.argv[1])    # min 1
    iDeltaEnd = iDeltaStart+1         # max 9


    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.AutoLibraryLoader.enable()
    
    # sampleList = None # run through all
    histos_from_main = {}
    # debug config:
    # Run printTTree.py alone to get the samples list
    sampleList = []
    sampleList.append("data_MinimumBias_Run2015A")
    
    

    slaveParams = {}
    
    slaveParams['iDeltaStart'] = iDeltaStart
    slaveParams['iDeltaEnd'] = iDeltaEnd
    

    if DATASOURCE == "DATA":
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

    
        # decide on output filename
        slaveParams["maxFileNo"] = maxFileNo
#        outFileName = join(outfolder,"plotsMuonselectioncuts_2_{n:04d}.root".format(n=maxFileNo+1))
        outFileName = "muons_" + str(iDeltaStart) + "_" + str(iDeltaEnd) + ".root"
        
        if maxFileNo == -1:
            print "No previous output file found"
        else:
            print "Found prevous output file(s). Setting new output file name to:", outFileName

        # use printTTree.py <sampleName> to see what trees are avaliable inside the skim file

        print ("Writing output to file: " + outFileName)

        Step2_Selection_2.runAll(treeName="MuonCastorVTwo",
                               slaveParameters = slaveParams,
                               sampleList = sampleList,
                               maxFilesMC = None,
                               maxFilesData = None,
                               nWorkers =8,
                               outFile = outFileName,
                               verbosity = 2)
#                               maxNevents = 10000,




    elif DATASOURCE == "TOYMC":
        print "Im in TOYMC!!!!!!!!!!!!!!!!!"
        # repetitions = 2
        # for iRep in xrange(repetitions):
        #max_events = 100 #if command out when you want to run whole files 
        max_events = int(NBRHIST*NBREVTPERHIST-1)
        slaveParams["maxFileNo"] =-1
        Step2_Selection_2.runAll(treeName="MuonCastorVTwo",
                               slaveParameters = slaveParams,
                               sampleList = sampleList,
                               maxFilesMC = None,
                               maxFilesData = None,
                               maxNevents = max_events,
                               nWorkers = 1,
                               outFile = join(outfolder,"output_for_toyMC.root"),
                               verbosity = 2,
                               donotvalidate = True)
       
