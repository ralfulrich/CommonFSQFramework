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
import json
# import time

# from outsource_analzye_muon import *


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
#        self.flag_use_merijn_electronic_channel_noise = False # mine RMS
        self.flag_use_merijn_electronic_channel_noise = True 
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

        # Getting sector RMS and Mean
        if firstRun:
            inputFile = ROOT.TFile(join(outfolder,"mean_rms.root"))
            inputFile_noise = ROOT.TFile(join(outfolder,"Histograms_StatFitInformationPedNoise.root"))
        else:
            inputFile = ROOT.TFile(join(outfolder,"plotsMuonselectioncuts_2_{n:04d}.root".format(n=self.maxFileNo)))
            # check naming if rerunning step1 again
        #
        


        #############################################################
        # also write a ttree with the final muon events!            #
        # in same format at CFF, but reduced content                #
        # create the branches and assign the fill-variables to them #
        #############################################################
        
        outTree = ROOT.TTree("muons", "selected muon events")
        
        self.OUTtrgl1L1GTTech = array( 'i', 103 * [0] )
        self.OUTtrgRandom = array( 'i', [0] )
        self.OUTtrgCastorHaloMuon = array( 'i', [0] )
        self.OUTCastorRecHitEnergy = array( 'd', 224 * [0.0])
        self.OUTrun = array( 'i', [0] )
        self.OUTlumi = array( 'i', [0] )

        outTree.Branch('trgl1L1GTTech', self.OUTtrgl1L1GTTech, 'trgl1L1GTTech[103]/I' )
        outTree.Branch('trgRandom', self.OUTtrgRandom, 'trgRandom/I')
        outTree.Branch('trgCastorHaloMuon', self.OUTtrgCastorHaloMuon, 'trgCastorHaloMuon/I')
        outTree.Branch('CastorRecHitEnergy', self.OUTCastorRecHitEnergy, 'CastorRecHitEnergy[224]/D')
        outTree.Branch('run', self.OUTrun, 'run/I')
        outTree.Branch('lumi', self.OUTlumi, 'lumi/I')

        setattr(self, "outTree", outTree)
        self.addToOutput(self.outTree)
        
        
        
        
        # #Getting channel RMS and Mean
        # hist_ch_Mean = inputFile.Get("data_MinimumBias_Run2015A/hist_ch_Mean")
        # hist_ch_RMS =  inputFile.Get("data_MinimumBias_Run2015A/hist_ch_RMS")

        
        hist_sec_Mean = inputFile.Get("data_MinimumBias_Run2015A/hist_sec_Mean")
        hist_sec_RMS = inputFile.Get("data_MinimumBias_Run2015A/hist_sec_RMS")
        #Getting channel RMS and Mean
        hist_ch_Mean = inputFile.Get("data_MinimumBias_Run2015A/hist_ch_Mean")
        hist_ch_RMS =  inputFile.Get("data_MinimumBias_Run2015A/hist_ch_RMS")

        # print "Pointer Address of Noise level hist:", hist_sec_Mean

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
                self.hist[histcalibrationname] = inputFile.Get("data_Cosmics_MuonHLTSkim_2015E_4T/2DMuonSignalMap")
                print "Extracted histogram from file. Checking entries:",  self.hist[histcalibrationname].GetEntries()

            else: #first time running analyser the calibration constants are all set to 1
                self.hist[histcalibrationname] =  ROOT.TH2D(histcalibrationname,histcalibrationname, 14, 0.5, 14.5, 16, 0.5, 16.5)

                for imod in xrange(0,14):
                    for isec in xrange(0,16):
                        self.hist[histcalibrationname].SetBinContent(self.hist[histcalibrationname].FindBin(imod+1,isec+1), 1.) #all factors set to 1



            self.hist["EventCount_Excl" + str(iDelta)] = ROOT.TH1D("EventCount_Excl" + str(iDelta),"EventCount",10,0,20)

            self.hist["MuonCount_Muo_Excl" + str(iDelta)] = ROOT.TH1D("MuonCount_Muo_Excl"+str(iDelta),"MuonCount (Muo);N_{muon}/evt;;", 17, -0.5, 16.5)
            self.hist["MuonCount_Rnd_Excl" + str(iDelta)] = ROOT.TH1D("MuonCount_Rnd_Excl"+str(iDelta),"MuonCount (Rnd);N_{muon}/evt;;", 17, -0.5, 16.5)

            self.hist["2DMuonCountMap_Excl" + str(iDelta)] =  ROOT.TH2D("2DMuonCountMap_Excl" + str(iDelta),"2DMuonCountMap;module;sector;", 14, 0.5, 14.5, 16, 0.5, 16.5)
            self.hist["2DMuonCountFor3Ch_Excl" + str(iDelta)] =  ROOT.TH2D("2DMuonCountFor3Ch_Excl" + str(iDelta),"2DMuonCountFor3Ch;module;sector;", 14, 0.5, 14.5, 16, 0.5, 16.5)
            self.hist["2DMuonactiveregion_Excl" + str(iDelta)] =  ROOT.TH2D("2DMuonactiveregion_Excl" + str(iDelta),"2DMuonactiveregion;module;sector", 14, 0.5, 14.5, 16, 0.5, 16.5)

            self.hist["2DcountChannelAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_MuoEvt"] =  ROOT.TH2D("2DcountChannelAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_MuoEvt","2DcountChannelAboveNoiseFor3SigmaSectors_MuoEvt;module;sector;", 14, 0.5, 14.5, 16, 0.5, 16.5)
            self.hist["2DcountChannelAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_RndEvt"] =  ROOT.TH2D("2DcountChannelAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_RndEvt","2DcountChannelAboveNoiseFor3SigmaSectors_RndEvt;module;sector;", 14, 0.5, 14.5, 16, 0.5, 16.5)

            self.hist["2DcountChannelsAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_MuoEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_MuoEvt","2DcountChannelsAboveNoiseFor3SigmaSectors_MuoEvt;module;sector;", 14, -0.5, 13.5, 16, 0.5, 16.5)
            self.hist["2DcountChannelsAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_RndEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseFor3SigmaSectors_Excl"+str(iDelta)+"_RndEvt","2DcountChannelsAboveNoiseFor3SigmaSectors_RndEvt;module;sector;", 14, -0.5, 13.5, 16, 0.5, 16.5)

            self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_MuoEvt"] = ROOT.TH2D("2DDeltaSigma_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_MuoEvt","2DDeltaSigma_Vs_SecRMSSecondHot_MuoEvt;#Delta#sigma;#sigma(second);",88,-2,20, 88, -2, 20)
            self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_RndEvt"] = ROOT.TH2D("2DDeltaSigma_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_RndEvt","2DDeltaSigma_Vs_SecRMSSecondHot_RndEvt;#Delta#sigma;#sigma(second);",88, -2, 20, 88,-2, 20)
            self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_MuoCand"] = ROOT.TH2D("2DDeltaSigma_Vs_SecRMSSecondHot_Excl"+str(iDelta)+"_MuoCand","2DDeltaSigma_Vs_SecRMSSecondHot_MuoCand;#Delta#sigma;#sigma(second);",88,-2,20, 88, -2, 20)

            # self.hist["2DMuonNoTriggerCountMap"] =  ROOT.TH2D("2DMuonNoTriggerCountMap","2DMuonNoTriggerCountMap", 14, -0.5, 13.5, 16, -0.5, 15.5)
            self.hist["GoodMuonCountPerSec_Excl"+str(iDelta)] =  ROOT.TH1D("GoodMuonCountPerSec_Excl"+str(iDelta),"GoodMuonCountPerSec;sector;N_{muon};", 16, 0.5, 16.5)
            self.hist["RunsWithGoodMuons_Excl"+str(iDelta)] =  ROOT.TH1D("RunsWithGoodMuons_Excl"+str(iDelta), "RunsWithGoodMuons", 1, 1, 0)
            self.hist["DeltaSigma_Excl"+str(iDelta)+"_MuoEvt"] = ROOT.TH1D("DeltaSigma_Excl"+str(iDelta)+"_MuoEvt","DeltaSigma_MuoEvt;#Delta#sigma;;",100, 0, 50)
            self.hist["DeltaSigma_Excl"+str(iDelta)+"_RndEvt"] = ROOT.TH1D("DeltaSigma_Excl"+str(iDelta)+"_RndEvt","DeltaSigma_RndEvt;Delta#sigma;;",100, 0, 50)
            self.hist["DeltaSigma_Excl"+str(iDelta)+"_MuoCand"] = ROOT.TH1D("DeltaSigma_Excl"+str(iDelta)+"_MuoCand","DeltaSigma_MuoCand;#Delta#sigma;;",100, 0, 50)


            for isec in xrange(0,16):

                henergy= 'MuonSignalSec_energy_Excl'+str(iDelta)+'_sec_{sec}'.format(sec=str(isec+1))
                self.hist[henergy] = ROOT.TH1D( henergy, henergy+";Energy;;", 100, 0, 1000)
               
                for imod in xrange(0,14):
                    hname = 'MuonSignalSecCh_Excl'+str(iDelta)+'_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))               
                    self.hist[hname] = ROOT.TH1D( hname, hname+";Energy;;", 50, -100, 400)

                    if DATASOURCE == "TOYMC":
                        hMC1D = '1DSignalMCCh_Excl'+str(iDelta)+'_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                        self.hist[hMC1D] = ROOT.TH1D(hMC1D,hMC1D+";Energy;;",100, -10, 90)
                        
                        # I have 10000 events and 10000/1000 i need to have 10 histograms
                        for iRep in xrange(NBRHIST): 
                            hname_MC = 'MuonSignalMCCh_Excl'+str(iDelta)+'_mod_{mod}_sec_{sec}_number_{number}'.format(mod=str(imod+1), sec=str(isec+1), number=str(iRep))
                            self.hist[hname_MC] = ROOT.TH1D( hname_MC, hname_MC, 50, -100, 400)   


            hnameAllsec= 'MuonSignalAllSec_Excl'+str(iDelta)+'_energy'
            self.hist[hnameAllsec] = ROOT.TH1D(hnameAllsec, hnameAllsec+";Energy;;", 150, 0, 2500)
            
            henergyGeV='MuonSignalAllSec_GeV_Excl'+str(iDelta)+'_energy'
            self.hist[henergyGeV] = ROOT.TH1D(henergyGeV, henergyGeV+";Energy;;", 500, 0, 75)

            histMuonSignalMean = "1DMuonsignalMean_Excl"+str(iDelta)
            self.hist[histMuonSignalMean] = ROOT.TH1D(histMuonSignalMean,histMuonSignalMean+";channel;<Energy>;",224,-0.5,223.5)

            histMuonFluct = "MuonFluctuations_Excl"+str(iDelta)
            self.hist[histMuonFluct] = ROOT.TH2D(histMuonFluct, histMuonFluct+";E_{sector};E_{channel}^{max};", 50,0,1000,70,0,400)







        
        self.hist["2DcountChannelAboveNoise_MuoEvt"] =  ROOT.TH2D("2DcountChannelAboveNoise_MuoEvt","2DcountChannelAboveNoise_MuoEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)
        self.hist["2DcountChannelAboveNoise_RndEvt"] =  ROOT.TH2D("2DcountChannelAboveNoise_RndEvt","2DcountChannelAboveNoise_RndEvt", 14, 0.5, 14.5, 16, 0.5, 16.5)

        self.hist["2DcountChannelsAboveNoiseForAllSectors_MuoEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseForAllSectors_MuoEvt","2DcountChannelsAboveNoiseForAllSectors_MuoEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)
        self.hist["2DcountChannelsAboveNoiseForAllSectors_RndEvt"] = ROOT.TH2D("2DcountChannelsAboveNoiseForAllSectors_RndEvt","2DcountChannelsAboveNoiseForAllSectors_RndEvt", 14, -0.5, 13.5, 16, 0.5, 16.5)


        self.hist["2DcountSectorRMS_MuoEvt"] =  ROOT.TH2D("2DcountSectorRMS_MuoEvt","2DcountSectorRMS_MuoEvt", 88, -2, 20, 16, 0.5, 16.5)
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



        # TProfile for storing (means) and RMS
        self.hist['hist_sec_Mean'] = ROOT.TProfile('hist_sec_Mean','hist_sec_Mean',16,0.5,16.5)
        self.hist['hist_sec_RMS'] = ROOT.TProfile('hist_sec_RMS','hist_sec_RMS',16,0.5,16.5) #change later to hist_sec_mean
        self.hist['hist_ch_RMS'] = ROOT.TProfile('hist_ch_RMS','hist_ch_RMS',224,0.5,224.5)
        self.hist['hist_ch_good'] = ROOT.TProfile('hist_ch_good','hist_ch_good',224,0.5,224.5)
                
                
        # Getting electronic noise from in the Merijns file
        hist_ch_noise_RMS = inputFile_noise.Get("SumHistoNoiseFit")  # this is the Gaussian Fit !!
        
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
            if self.flag_use_merijn_electronic_channel_noise:
                self.ch_RMS[isec][imod] = hist_ch_noise_RMS.GetBinContent(i+1) #+1 because of underflow
            else:
                self.ch_RMS[isec][imod] = hist_ch_RMS.GetBinContent(i+1) #+1 because of underflow

            # and also remember for later reference
            self.hist['hist_ch_RMS'].Fill(i+1, self.ch_RMS[isec][imod])
            if [isec+1,imod+1] not in badChannelsSecMod:
                self.hist['hist_ch_good'].Fill(i+1, 1)


            #print i, "sec,mod", sec, mod, hist_ch_RMS.GetBinContent(i+1)


        #get mean and rms from sector histogram
        self.sec_mean = [0] * 16
        self.sec_RMS = [0] * 16
        for i in xrange(0,16):
            self.sec_mean[i] = hist_sec_Mean.GetBinContent(i+1) #+1 because of underflow
            self.sec_RMS[i] = hist_sec_RMS.GetBinContent(i+1) #+1 because of underflow

            # and also remember for later reference
            self.hist['hist_sec_RMS'].Fill(i+1, self.sec_RMS[i])
            self.hist['hist_sec_Mean'].Fill(i+1, self.sec_mean[i])

        #these ones are needed to update mean and RMS for the new calibration constants
#        self.new_ch_mean = [[0 for _ in xrange(14)] for _ in xrange(16)]
#        self.new_ch_RMS = [[0 for _ in xrange(14)] for _ in xrange(16)]
#        self.new_sec_mean = [0] * 16
#        self.new_sec_RMS = [0] * 16

        

        
        
        for h in self.hist:
            self.hist[h].Sumw2()
            self.GetOutputList().Add(self.hist[h])




        # from Merijn using the characterization data
        self.TranspositionFactor = [[17.8879, 13.2482, 14.0743, 13.3278, 2.12091, 2.12543, 3.21527, 2.98481, 2.11992, 2.06513, 2.10157, 4.18052, 2.10157, 2.10157],
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


        # from Melike using everything, also muon intercalib results
        self.Absolutecalibration2015 = [[0.016874, 0.016531, 0.02974, 0.015972, 0.004336, 0.005311, 0.005042, 0.005669, 0.00246, 0.003017, 0.00382, 0.006415, 0.008673, 0.004469], 
                    [0.019528, 0.013637, 0.039159, 0.017504, 0.005906, 0.004969, 0.0, 0.004563, 0.002357, 0.0, 0.002439, 0.005647, 0.00515, 0.004313], 
                    [0.018033, 0.015976, 0.032389, 0.016185, 0.003734, 0.00437, 0.005185, 0.0, 0.001718, 0.002401, 0.0, 0.004594, 0.008272, 0.005952],
                    [0.029058, 0.014194, 0.033104, 0.016516, 0.003732, 0.005391, 0.006832, 0.006688, 0.002633, 0.003109, 0.003426, 0.003281, 0.003259, 0.008384], 
                    [0.019439, 0.01289, 0.029525,  0.0, 0.002122, 0.002438, 0.003519, 0.004777, 0.001497, 0.002614, 0.00131, 0.001461, 0.003673, 0.002588], 
                    [0.020009, 0.014693, 0.032381, 0.015669, 0.002941, 0.002977, 0.005145, 0.005011, 0.0021, 0.003598, 0.001715, 0.0005, 0.006963, 0.003981],
                    [0.011084, 0.016251, 0.020354, 0.017343, 0.0, 0.013322, 0.007125, 0.004629, 0.005993, 0.021396, 0.013318, 0.023119, 0.024841, 0.057863], 
                    [0.012645, 0.011933, 0.025745, 0.016747, 0.003168, 0.006409, 0.005631, 0.0, 0.00398, 0.020326, 0.012471, 0.00766, 0.02901, 0.008533], 
                    [0.017433, 0.02335, 0.022966, 0.0272, 0.002624, 0.009199, 0.00517, 0.003675, 0.001882, 0.002163, 0.001921, 0.002456, 0.003948, 0.003556], 
                    [0.016628, 0.021822, 0.038431, 0.02047, 0.003518, 0.011943, 0.004984, 0.004376, 0.002686, 0.002735, 0.002116, 0.001831, 0.00444, 0.005347], 
                    [0.022451, 0.019186, 0.029669, 0.024445, 0.005458, 0.002917, 0.007484, 0.008021, 0.003736, 0.002226, 0.008541, 0.007556, 0.0187, 0.01019], 
                    [0.024777, 0.019928, 0.028925, 0.032663, 0.004808, 0.005331, 0.006216, 0.008464, 0.003521, 0.002246, 0.005665, 0.0, 0.005555, 0.007251], 
                    [0.01675, 0.017174, 0.039164, 0.017898, 0.004786, 0.024493, 0.007055, 0.006875, 0.003386, 0.003246, 0.011542, 0.014695, 0.015738, 0.031568],
                    [0.013626, 0.026727, 0.034807, 0.026609, 0.004608, 0.015818, 0.006961, 0.0, 0.002699, 0.002787, 0.009467, 0.005394, 0.013983, 0.121789],
                    [0.014901, 0.017273, 0.035037, 0.071443, 0.002436, 0.008082, 0.005174, 0.004774, 0.002453, 0.001865, 0.000842, 0.003679, 0.011773, 0.00654], 
                    [0.014908, 0.01537, 0.033796, 0.078791, 0.002585, 0.014259, 0.005098, 0.005297, 0.002803, 0.002401, 0.002467, 0.002542, 0.006237, 0.004821]]



        # eventually read additional json file for LS selection
        json_file = open(join(outfolder,"../muon_2015_3.8T.json"))
        self.JSON = json.loads(json_file.read())
        json_file.close()






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


    ##############################################
    # find all channel sigma values, and         #
    # count channels above noise-cut in sectors  #
    ##############################################
    def getListSigmaChannel(self, energy_ch, sigma_cut, badChannelsSecMod):
        listSigmaChannel = [[0 for _ in xrange(14)] for _ in xrange(16)]
        listAllChannelsAboveNoise = [[] for _ in xrange(16)]
        for isec in xrange(16):
            for imod in xrange(14):
                energy = energy_ch[isec][imod]
                # print "RMS mod {mod} sec {sec}: {rms}".format(mod=str(imod),sec=str(isec),rms=self.ch_RMS[isec][imod])
                if [isec+1,imod+1] in badChannelsSecMod:
                    listSigmaChannel[isec][imod] = None
                    # print "skipping channel", imod, isec
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
        ##################################################   
        ##################################################
        ##################################################
        ##################################################
        # if self.fChain.run <= 247483: return 1
        ##################################################
        ##################################################
        ##################################################
        ##################################################
        weight = 1 # Ralf: what is this used for ? CFF internal?
        num = 0 # Ralf: what is this used for ? CFF internal?
                
        
        self.hist["EventCount"].Fill("all",1)
        
        
        if self.JSON != None:
            goodLS = False
            if str(self.fChain.run) in self.JSON: 
                for LS in self.JSON[str(self.fChain.run)]:
                    if self.fChain.lumi>=LS[0] and self.fChain.lumi<=LS[1]:
                        goodLS = True
            if not goodLS:
                return 0
            self.hist["EventCount"].Fill("goodLS",1)
            


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

        # if self.flag_use_merijn_electronic_channel_noise:
#        for isec in xrange(16):
 #           for imod in xrange(14):
  #              hnoise_electronic_RMS = 'CastorNoise_electronic_RMS_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
   #             self.hist[hnoise_electronic_RMS].Fill(self.ch_RMS[isec][imod])

                
        histCalibration = self.hist[histcalibrationname]
        #histcalibrationname_sigma = '2DMuonSignalMap_sigma'


        ch_energy = [[0.0 for _ in xrange(14)] for _ in xrange(16)]
        sec_energy = [0.0] * 16

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

            if [isec+1, imod+1] not in badChannelsSecMod:
                sec_energy[isec] += ich_energy
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

        
        listSigmaSector = self.getListSigmaSector(sec_energy)                  # Melike
#        listSigmaSector = self.getListSigmaSectorFromChannels(listSigmaChannel) # ->  Ralf



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
        # the main muon selection loop starts here    #
        ###############################################

        event_with_muon = False
        
        num_muons = [0] * (self.iDeltaEnd-self.iDeltaStart)
        num_muons_rnd = [0] * (self.iDeltaEnd-self.iDeltaStart)
        
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



                ######################################################################################
                # we found a POSSIBLE muon candidate in iSectorTest with exclusivity iDeltaSector    #
                ######################################################################################
                                
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
                if SigmaSecHottestSector > 3: # WARNING: CUT VALUE ON NOISE HARDCODED HERE !!!! Melike: 2.5
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



                ################################################################################
                # found an interesting event: MUON CANDIDATE                                   #
                ################################################################################

                if goodMuonEventWithoutAnyTriggerSelection:

                    event_with_muon = True

                    # RU was: if goodMuonEvent:                    
                    # print "Good event in (sec,mod)", sec, mod, "Front,Mid,Back", Front_Module, Mid_Module, Rear_Module


#                    self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("good muon (all)",1)

                    if isRandom:
                       self.hist["EventCount_Excl" + str(iDeltaSector)].Fill("good muon (rnd)",1)
                       num_muons_rnd[iDeltaSector-self.iDeltaStart] += 1

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
                        
                        self.hist["DeltaSigma_Excl" + str(iDeltaSector) + "_MuoCand"].Fill(DeltaSigma)
                        self.hist["2DDeltaSigma_Vs_SecRMSSecondHot_Excl" + str(iDeltaSector) + "_MuoCand"].Fill(DeltaSigma,SigmaSecHottestSector)

                        num_muons[iDeltaSector-self.iDeltaStart] += 1
                       
                        energy_secsum = [0.0] * 16
                        energyGeV_secsum = [0.0] * 16
                        max_channel = [0.0] * 16
                        for isec in xrange(16):
                            for imod in xrange(14):
                                if [isec+1,imod+1] in badChannelsSecMod:
                                    continue
                                energy_secsum[isec] += ch_energy[isec][imod]
                                trpvalues = self.TranspositionFactor[isec][imod]
                                absintercalibvalues = self.Absolutecalibration2015[isec][imod]
                                energyGeV_secsum [isec] += ch_energy[isec][imod] * absintercalibvalues / trpvalues
                                if (ch_energy[isec][imod] > max_channel[isec]):
                                    max_channel[isec] = ch_energy[isec][imod]
                                
                                
                        henergy = 'MuonSignalSec_energy_Excl' + str(iDeltaSector) + '_sec_{sec}'.format(sec=str(muonSec+1))
                        self.hist[henergy].Fill(energy_secsum[muonSec])
                        

                        if DATASOURCE == "TOYMC":
                            for imod in xrange(14):
                                hname_MC = 'MuonSignalMCCh_Excl' + str(iDeltaSector) + '_mod_{mod}_sec_{sec}_number_{number}'.format(mod=str(imod+1), sec=str(muonSec+1), number=str(Nevent))   
                                self.hist[hname_MC].Fill(ch_energy[muonSec][imod])
                                # mean_mc= self.hist[hname_MC].GetMean()


                        self.hist['MuonSignalAllSec_Excl' + str(iDeltaSector) + '_energy'].Fill(energy_secsum[muonSec])
                        self.hist['MuonSignalAllSec_GeV_Excl' + str(iDeltaSector) + '_energy'].Fill(energyGeV_secsum[muonSec])
                        
                        self.hist["MuonFluctuations_Excl"+str(iDeltaSector)].Fill(energy_secsum[muonSec], max_channel[muonSec])

                        self.hist["RunsWithGoodMuons_Excl" + str(iDeltaSector)].Fill(str(self.fChain.run), 1)

                        for imod in xrange(0,14):
                            hname = 'MuonSignalSecCh_Excl' + str(iDeltaSector) + '_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(muonSec+1))
                            self.hist[hname].Fill(ch_energy[muonSec][imod])

                        self.hist["GoodMuonCountPerSec_Excl" + str(iDeltaSector)].Fill(muonSec+1)

                        for imod in listAllChannelsAboveNoise[muonSec]:
                            self.hist["2DMuonCountMap_Excl" + str(iDeltaSector)].Fill(imod+1,muonSec+1)

                        #self.tree["nt"].Fill()


                            
        for iDeltaSector in xrange(self.iDeltaStart,self.iDeltaEnd): 
            self.hist["MuonCount_Muo_Excl" + str(iDeltaSector)].Fill(num_muons[iDeltaSector-self.iDeltaStart])
            self.hist["MuonCount_Rnd_Excl" + str(iDeltaSector)].Fill(num_muons_rnd[iDeltaSector-self.iDeltaStart])
        
        
        if event_with_muon:
            
            for i in xrange(0, len(self.OUTtrgl1L1GTTech)):
                if i<len(self.fChain.trgl1L1GTTech):
                    self.OUTtrgl1L1GTTech[i] = self.fChain.trgl1L1GTTech[i]
                else:
                    self.OUTtrgl1L1GTTech[i] = 0
            self.OUTtrgRandom = self.fChain.trgRandom
            self.OUTtrgCastorHaloMuon = self.fChain.trgCastorHaloMuon
            for i in xrange(0, len(self.OUTCastorRecHitEnergy)):
                if i<len(self.fChain.CastorRecHitEnergy):
                    self.OUTCastorRecHitEnergy[i] = self.fChain.CastorRecHitEnergy[i]
                else:
                    self.OUTCastorRecHitEnergy[i] = 0
            self.OUTrun = self.fChain.run
            self.OUTlumi = self.fChain.lumi
            self.outTree.Fill()
            
        return 1




    def finalize(self):
        
        if hasattr(self, 'outTree'):
            self.outTree.AutoSave()
        
        print "Finalize:"




    def finalizeWhenMerged(self):
        # print "Finalize when merged"
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

               # hnoise_randomtrg = 'CastorNoise_randomtrg_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))  
               # binnumber = histos["EventCount"].GetXaxis().FindBin("rnd trg")
               # if not histos["EventCount"].GetBinContent(binnumber) == 0:
               #     histos[hnoise_randomtrg].Scale( 1./histos["EventCount"].GetBinContent(binnumber) )

                # scale muon energy distributions by number of muons
                for iDeltaSector in xrange(self.iDeltaStart,self.iDeltaEnd):
                    hname = 'MuonSignalSecCh_Excl' + str(iDeltaSector) + '_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                
                    if not histos["GoodMuonCountPerSec_Excl" + str(iDeltaSector)].GetBinContent(isec+1) == 0:
                        histos[hname].Scale( 1./histos["GoodMuonCountPerSec_Excl"+str(iDeltaSector)].GetBinContent(isec+1) )

                



        # inputFile = ROOT.TFile(join(outfolder,"mean_rms.root"))
        # hist_ch_Mean = inputFile.Get("data_MinimumBias_Run2015A/hist_ch_Mean")

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
        pass


        

if __name__ == "__main__":

    if len(sys.argv) != 4 :
        print ("ERROR: SPECIFY iDeltaStart, iDeltaEnd, dataset.  ")
        print ("       iDeltaStart: minimal muon exclusivity (number of empty sectors next to muon)")
        print ("       iDeltaEnd: maximal muon exclusivity around muon")
        print ("       dataset: data_Cosmics_MuonHLTSkim_2015E_4T or data_MinimumBias_Run2015A")
        print ("Try again!!")
        sys.exit(1)
    
    iDeltaStart = int(sys.argv[1])    # min 1
#    iDeltaEnd = iDeltaStart+1         # max 9
    
#    if len(sys.argv) == 3:
    iDeltaEnd = int(sys.argv[2]) + 1
    
    if iDeltaStart<1 or iDeltaStart>9 or iDeltaEnd<1 or iDeltaEnd>9 or iDeltaEnd<iDeltaStart+1:
        print ("start and end must be between 1 and 8, and start<end")
        sys.exit(1)
        
    datasetname = str(sys.argv[3])
    if datasetname != "data_Cosmics_MuonHLTSkim_2015E_4T" and datasetname != "data_MinimumBias_Run2015A" :
        print ("dataset must be either: data_Cosmics_MuonHLTSkim_2015E_4T or data_MinimumBias_Run2015A")
        sys.exit(1)


    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.AutoLibraryLoader.enable()
    
    # sampleList = None # run through all
    histos_from_main = {}
    # debug config:
    # Run printTTree.py alone to get the samples list
    sampleList = []
#    sampleList.append("data_MinimumBias_Run2015A") # all triggers, 18MEvents
#    sampleList.append("data_Cosmics_MuonHLTSkim_2015E_4T") # only the muon triggers 
    sampleList.append(datasetname) 

    

#    print json.__class__.__name__
#    sys.exit(1)

    slaveParams = {}
    
    slaveParams['iDeltaStart'] = iDeltaStart
    slaveParams['iDeltaEnd'] = iDeltaEnd
    # slaveParams['json'] = json


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
       
