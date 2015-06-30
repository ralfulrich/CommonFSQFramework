#!/usr/bin/env python
import CommonFSQFramework.Core.ExampleProofReader
from BadChannels2015 import badChannelsModSec

import sys, os, time
sys.path.append(os.path.dirname(__file__))
from os import listdir
from os.path import isfile, join

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import edm
from math import sqrt
from array import *

from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

outfolder = os.environ["HaloMuonOutput"]

class Step2_Selection(CommonFSQFramework.Core.ExampleProofReader.ExampleProofReader):
    def init( self, maxEvents = None):
#        self.maxFileNo = slaveParameters["maxFileNo"]
        firstRun = (self.maxFileNo == -1)

        self.maxEvents = maxEvents
        self.hist = {}
        self.nEventsRandom = 0

        self.hist["EventCount"] = ROOT.TH1D("EventCount","EventCount",10,0.5,10.5)

        #Getting sector RMS and Mean
        if firstRun:
            inputFile = ROOT.TFile(join(outfolder,"mean_rms.root"))
        else:
            inputFile = ROOT.TFile(join(outfolder,"plotsMuonselectioncuts_{n:04d}.root".format(n=self.maxFileNo)))
            # check naming if rerunning step1 again
        hist_sec_Mean = inputFile.Get("data_MinimumBias_Run2015A/hist_sec_Mean")
        hist_sec_RMS = inputFile.Get("data_MinimumBias_Run2015A/hist_sec_RMS")

        #Getting channel RMS and Mean
        hist_ch_Mean = inputFile.Get("data_MinimumBias_Run2015A/hist_ch_Mean")
        hist_ch_RMS =  inputFile.Get("data_MinimumBias_Run2015A/hist_ch_RMS")

        self.hist["2DMuonCountMap"] =  ROOT.TH2D("2DMuonCountMap","2DMuonCountMap", 14, -0.5, 13.5, 16, -0.5, 15.5)
        self.hist["2DMuonCountMap_allch"] =  ROOT.TH2D("2DMuonCountMap_allch","2DMuonCountMap_allch", 14, -0.5, 13.5, 16, -0.5, 15.5)
        self.hist["2DCountMap_allch"] =  ROOT.TH2D("2DCountMap_allch","2DCountMap_allch", 14, -0.5, 13.5, 16, -0.5, 15.5)
        self.hist["2DcountChannelsAboveNoiseInSectorAbove3Sigma"] =  ROOT.TH2D("2DcountChannelsAboveNoiseInSectorAbove3Sigma","2DcountChannelsAboveNoiseInSectorAbove3Sigma", 14, -0.5, 13.5, 16, -0.5, 15.5)
        self.hist["2DcountChannelsAboveNoiseForAllSectors"] = ROOT.TH2D("2DcountChannelsAboveNoiseForAllSectors","2DcountChannelsAboveNoiseForAllSectors", 14, -0.5, 13.5, 16, -0.5, 15.5)
        self.hist["2DcountSectorAboveRMS"] =  ROOT.TH2D("2DcountSectorAboveRMS","2DcountSectorAboveRMS", 500, -50, 450, 16, -0.5, 15.5)
        self.hist["2DMuonNoTriggerCountMap"] =  ROOT.TH2D("2DMuonNoTriggerCountMap","2DMuonNoTriggerCountMap", 14, -0.5, 13.5, 16, -0.5, 15.5)
        self.hist["GoodMuonCountPerSec"] =  ROOT.TH1D("GoodMuonCountPerSec","GoodMuonCountPerSec", 16, -0.5, 15.5)
        self.hist["RunsWithGoodMuons"] =  ROOT.TH1D("RunsWithGoodMuons","RunsWithGoodMuons", 10000, 247000-0.5, 257000-0.5)
        self.hist["RunsAllTrigger"] =  ROOT.TH1D("RunsAllTrigger","RunsAllTrigger", 10000, 247000-0.5, 257000-0.5)

        histcalibrationname = '2DMuonSignalMap'
        if not firstRun:
            self.hist[histcalibrationname] = inputFile.Get("data_MinimumBias_Run2015A/2DMuonSignalMap")
           # print "Extracted histogram from file. Checking entries:",  self.hist[histcalibrationname].GetEntries()
        else: #first time running analyser the calibration constants are all set to 1
            self.hist[histcalibrationname] =  ROOT.TH2D(histcalibrationname,histcalibrationname, 14, -0.5, 13.5, 16, -0.5, 15.5)
            for imod in xrange(0,14):
                for isec in xrange(0,16):
                    self.hist[histcalibrationname].SetBinContent(self.hist[histcalibrationname].FindBin(imod,isec), 1.) #all factors set to 1
                    print "Created new 2d calibration map. Checking entries:",  self.hist[histcalibrationname].GetEntries(), ". maxFileNo =", self.maxFileNo, firstRun
        histCalibration= self.hist[histcalibrationname]
        #histCalibration.SetBit(ROOT.TH1.kIsAverage) # does not seem to work. use TProfile2D instead

        #make new histograms
        for isec in xrange(0,16):
            for imod in xrange(0,14):
                hname = 'MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                hwithouttrigger = 'MuonSignalSecCh_withoutrigger_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                hnoise_neighbor = 'CastorNoise_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))

                hnoise_randomtrg = 'CastorNoise_randomtrg_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                #hempty = 'CastorEmptysectors_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                self.hist[hname] = ROOT.TH1D( hname, hname, 50, -100, 400)
                self.hist[hwithouttrigger] = ROOT.TH1D(hwithouttrigger, hwithouttrigger, 50, -100, 400)
                self.hist[hnoise_neighbor] = ROOT.TH1D(hnoise_neighbor, hnoise_neighbor, 50, -100, 400)
                self.hist[hnoise_randomtrg] = ROOT.TH1D(hnoise_randomtrg, hnoise_randomtrg, 50, -100, 400)
                #self.hist[hempty] = ROOT. TH1D( hempty, hempty, 500, -100, 400)

        #get channel energies from input file
        self.ch_mean = [[0 for _ in range(14)] for _ in range(16)]
        self.ch_RMS = [[0 for _ in range(14)] for _ in range(16)]
        for i in xrange(0, 224):
            #how I get Channel information
            sec = i//14
            mod = i%14
            #print sec, mod
            #i = sec*14+module   -> sec0 mod0 -> i=0   ; sec9 mod 4
            self.ch_mean[sec][mod] = hist_ch_Mean.GetBinContent(i+1) #+1 because of underflow
            self.ch_RMS[sec][mod] = hist_ch_RMS.GetBinContent(i+1) #+1 because of underflow

            #print i, "sec,mod", sec, mod, hist_ch_RMS.GetBinContent(i+1)

        #get mean and rms from sector histogram
        self.sec_mean = [0] * 16
        self.sec_RMS = [0] * 16
        for i in xrange(0,16):
            self.sec_mean[i] = hist_sec_Mean.GetBinContent(i+1) #+1 because of underflow
            self.sec_RMS[i] = hist_sec_RMS.GetBinContent(i+1) #+1 because of underflow

        #these ones are needed to update mean and RMS for the new calibration constants
        self.new_ch_mean = [[0 for _ in range(14)] for _ in range(16)]
        self.new_ch_RMS = [[0 for _ in range(14)] for _ in range(16)]
        self.new_sec_mean = [0] * 16
        self.new_sec_RMS = [0] * 16

        #TProfile for storing means and RMS. and when jobs are merged, the averge is taken
        self.hist['hist_sec_Mean'] = ROOT.TProfile('hist_sec_Mean','hist_sec_Mean',16,-0.5,15.5)
        self.hist['hist_sec_RMS'] = ROOT.TProfile('hist_sec_RMS','hist_sec_RMS',16,-0.5,15.5) #change later to hist_sec_mean
        self.hist['hist_ch_Mean'] = ROOT.TProfile('hist_ch_Mean','hist_ch_Mean',224,-0.5,223.5)
        self.hist['hist_ch_RMS'] = ROOT.TProfile('hist_ch_RMS','hist_ch_RMS',224,-0.5,223.5)


        for h in self.hist:
            self.hist[h].Sumw2()
            self.GetOutputList().Add(self.hist[h])

    def analyze(self):
        weight = 1
        num = 0
        goodMuonEvent = False
        # genTracks
        #num = self.fChain.genTracks.size()
        #print num
        #print self.maxEta # see slaveParams below

        self.hist["EventCount"].Fill(1)

        histcalibrationname = '2DMuonSignalMap'
        histCalibration= self.hist[histcalibrationname]

        ch_energy = [[0 for _ in range(14)] for _ in range(16)]
        sec_energy =[0] * 16
        sec_front_energy = [0] * 16

        if self.fChain.CastorRecHitEnergy.size() != 224:
            return 0

        self.hist["EventCount"].Fill(2)

        isRandom = False
        if self.fChain.trgRandom:
            if not (self.fChain.trgl1L1GTTech[1] or self.fChain.trgl1L1GTTech[2]): #not bptx+ or bptx-
                self.nEventsRandom += 1
                isRandom = 1

        if isRandom: self.hist["EventCount"].Fill(3)

        for i in xrange(0, 224):
            sec = self.fChain.CastorRecHitSector.at(i)-1
            mod = self.fChain.CastorRecHitModule.at(i)-1
            binnumber = histCalibration.FindBin(mod,sec)
            calibration =  histCalibration.GetBinContent(binnumber)
            ich_energy = calibration * self.fChain.CastorRecHitEnergy.at(i)
            ch_energy[sec][mod] = ich_energy
            hnoise_randomtrg = 'CastorNoise_randomtrg_mod_{mod}_sec_{sec}'.format(mod=str(mod), sec=str(sec))
            if [mod,sec] not in badChannelsModSec:
                sec_energy[sec] += ich_energy
                if mod < 5: sec_front_energy[sec] += ich_energy
            if isRandom:
                self.hist[hnoise_randomtrg].Fill(ch_energy[sec][mod])
                self.new_ch_mean[sec][mod] += ich_energy
                self.new_ch_RMS[sec][mod] += ich_energy**2
                if [mod,sec] not in badChannelsModSec:
                    self.new_sec_mean[sec] += ich_energy
                    self.new_sec_RMS[sec] += ich_energy**2




            #  print 'sec', sec, 'mod', mod, 'e_ch=', ch_energy, "e_sec", sec_energy[sec]
            #  print "all sector energies:  ", sec_energy

        listSectorsAboveNoise = []
        for i in xrange(0,16):
            if self.sec_RMS[i]:
                sigma = (sec_energy[i] - self.sec_mean[i]) / self.sec_RMS[i]
            else:
                print "Warning: RMS of sec", i, "is zero. Assume triggered"
                sigma = np.sign(sec_energy[i] - self.sec_mean[i]) * 1e9
            listSectorsAboveNoise.append([i,sigma])
            self.hist["2DcountSectorAboveRMS"].Fill(sigma,i)


        # selection of all channels above noise
        listAllChannelsAboveNoise = [[] for _ in range(16)]

        for i in xrange(0,224):
            sec = self.fChain.CastorRecHitSector.at(i)-1
            mod = self.fChain.CastorRecHitModule.at(i)-1
            if [mod,sec] in badChannelsModSec:
                #print "skipping channel", mod, sec
                continue

            if ch_energy[sec][mod]> self.ch_mean[sec][mod] + 2*self.ch_RMS[sec][mod]:
                # self.hist["2DMuonCountMap_allch"].Fill(mod,sec)
                listAllChannelsAboveNoise[sec].append(i)
                self.hist["2DCountMap_allch"].Fill(mod,sec)

        for isec in xrange(0,16):
            self.hist["2DcountChannelsAboveNoiseForAllSectors"].Fill(len(listAllChannelsAboveNoise[isec]),isec)


        def filterListSigma(l,sigma_th):
            return [[ii,isigma] for ii, isigma in l if isigma > sigma_th]

        #selection on sectors (bad channels already excluded)
        if not (len(filterListSigma(listSectorsAboveNoise,2)) == len(filterListSigma(listSectorsAboveNoise,3)) == 1): #only 1 sector above 3 Sigma and all others below 2 Sigma
            return 0
        # for i in xrange(0,16):
        #     if self.sec_RMS[i]:
        #         sigma = (sec_front_energy[i] - self.sec_mean[i]) / self.sec_RMS[i]
        #         if sigma >= 5:
        #             return 0
        muonSec = (filterListSigma(listSectorsAboveNoise,3))[0] #selected sector with muon. (secId,sigmaAboveNoise)
            #if listSectorsAboveNoise[0] != 8:
            #   return

        self.hist["EventCount"].Fill(4)

        #a selected sector have a signal above the noise threshold and remaining sectors are below a threshold.excluded the active sector + 2 neigbours + "the mirrored-sectors"
        isec_pn = ((muonSec[0] + 1) + 16)%16
        isec_mn = ((muonSec[0] - 1) + 16)%16

        isec_pns = ((muonSec[0] + 1) + 8)%16
        isec_mns = ((muonSec[0]- 1) + 8)%16
        #print "neighbors plus and minus:" , isec_pn, isec_mn, "mirrored neighbors plus and minus:" , isec_pns, isec_mns
        cond1 =  muonSec[0] or isec_pn or isec_mn
        cond2= ((muonSec[0]+8)%16 or isec_pns or isec_mns)
        for isec in xrange(0,16):
            cond1 =  isec == muonSec[0] or isec == isec_pn or isec == isec_mn
            cond2 = (isec == (muonSec[0]+8)%16 or isec == isec_pns or isec == isec_mns)
            for imod in xrange(0,14):
                hnoise_neighbor = 'CastorNoise_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                if not cond1 or not cond2:
                    self.hist[hnoise_neighbor].Fill(ch_energy[isec][imod])

                    
        #selection on channels above noise
        listChannelsAboveNoise = []
        for i in xrange(0,224):
            sec = self.fChain.CastorRecHitSector.at(i)-1
            mod = self.fChain.CastorRecHitModule.at(i)-1
            if [mod,sec] in badChannelsModSec:
                #print "skipping channel", mod, sec
                continue

            if ch_energy[sec][mod]> self.ch_mean[sec][mod] + 2*self.ch_RMS[sec][mod]:
                self.hist["2DMuonCountMap_allch"].Fill(mod,sec)
                if sec == muonSec[0]:
                    listChannelsAboveNoise.append(i)
            # hchannel = 'Channelenergy_mod_{mod}_sec_{sec}'.format(mod=str(mod), sec=str(sec))
            # self.hist[hchannel].Fill(ch_energy[sec][mod])
        countChannelsAboveNoise = len(listChannelsAboveNoise)
        #print "Channels above noise:", listChannelsAboveNoise

        self.hist["2DcountChannelsAboveNoiseInSectorAbove3Sigma"].Fill(countChannelsAboveNoise,muonSec[0])

        if countChannelsAboveNoise > 5:
            self.hist["EventCount"].Fill(5)

            hasMuonTrigger = (self.fChain.trgCastorHaloMuon or self.fChain.trgl1L1GTAlgo[102])

            Front_Module = False
            Mid_Module= False
            Rear_Module =False

            for i in listChannelsAboveNoise:
                sec = self.fChain.CastorRecHitSector.at(i)-1
                mod = self.fChain.CastorRecHitModule.at(i)-1
                if mod in [0,1,2,3,4]:
                    Front_Module =True

                if mod in [5,6,7,8,9]:
                    Mid_Module= True

                if mod in [10,11,12,13]:
                    Rear_Module= True

            if (Front_Module + Mid_Module + Rear_Module) >= 3: #maybe change to 3
                if hasMuonTrigger:
                    goodMuonEvent = True
        

        #found an interesting event. now fill histograms for channels above noise
        if goodMuonEvent:
            self.hist["EventCount"].Fill(6)

            print "Good event in (sec,mod)", sec, mod, "Front,Mid,Back", Front_Module, Mid_Module, Rear_Module
            self.hist["RunsWithGoodMuons"].Fill(self.fChain.run)
            for imod in xrange(0,14):
                hname = 'MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(muonSec[0]))
                self.hist[hname].Fill(ch_energy[muonSec[0]][imod])




            self.hist["GoodMuonCountPerSec"].Fill(muonSec[0])
            for ich in listChannelsAboveNoise:
                sec = self.fChain.CastorRecHitSector.at(ich)-1
                mod = self.fChain.CastorRecHitModule.at(ich)-1
                self.hist["2DMuonCountMap"].Fill(mod,sec)
        else:
            self.hist["RunsAllTrigger"].Fill(self.fChain.run)

        #plot all the good muons
        if goodMuonEvent:
            new_dir = join(outfolder,"goodMuonEvents_{n:04d}/".format(n=self.maxFileNo+1))
            if not os.path.exists(new_dir):
                os.makedirs(new_dir)

            trgEvtFileName  = "Muon_run_" + str(self.fChain.run) + "_event_" + str(self.fChain.event) + ".pdf"
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ch_energy_wo_bad_ch = np.asarray(ch_energy)
            for imod,isec in badChannelsModSec:
                ch_energy_wo_bad_ch[isec][imod] = 0
            for isec in np.arange(16):
                xs = np.arange(14)
                ys = ch_energy_wo_bad_ch[isec]
                ax.bar(xs, ys, zs=isec, zdir='y', color="g" if isec == muonSec[0] else "r", alpha=0.8)
            ax.set_xlabel('Module')
            ax.set_ylabel('Sector')
            ax.set_zlabel('Channel Energy')
            ax.set_title('Sec: {sec} ({sig:.1f}sigma)   Trigger: {trg}   Mod[F/M/R]: {f}/{m}/{r} \n Run: {run}   Event: {evt}'.format(sec=muonSec[0], sig=muonSec[1], trg="Yes" if hasMuonTrigger else "no", f=Front_Module, m=Mid_Module, r=Rear_Module, run=self.fChain.run, evt=self.fChain.event), multialignment='center')
            fig.savefig(join(join(outfolder,new_dir),trgEvtFileName))


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

        inputFile = ROOT.TFile(join(outfolder,"mean_rms.root"))
        hist_ch_Mean = inputFile.Get("data_MinimumBias_Run2015A/hist_ch_Mean")
        histcalibrationname = '2DMuonSignalMap'
        histcalibration = self.hist[histcalibrationname]
        hreferencename = 'MuonSignalSecCh_mod_3_sec_8' #reference channel is (counting from one) mod=4 sec=9

        referenceMean= self.hist[hreferencename].GetMean()
        meanReferenceNoise = hist_ch_Mean.GetBinContent(8*14 + 3 +1)
#        referenceMean -= meanReferenceNoise #switch on at some point or da a fit


        if referenceMean != 0:
            print "Warning reference channel (Sec 3 Mod 8) empty. Not dividing"
        for isec in xrange(0,16):
            for imod in xrange(0,14):
                hname = 'MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                hnoise_neighbor = 'CastorNoise_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                mean= self.hist[hname].GetMean()
                i = isec * 14 + imod
                meanNoise = hist_ch_Mean.GetBinContent(i+1) #noise is not estimated well. do iterative procedure from ralf? or random trigger?
                binnumber = histcalibration.FindBin(imod, isec)
                noiseSubtractedMean = (mean)# - meanNoise)
                if referenceMean != 0:
                    noiseSubtractedMean /= referenceMean
                if [imod,isec] in badChannelsModSec:
                    noiseSubtractedMean = 1 #set bad channels always to 1
                histcalibration.SetBinContent(binnumber, noiseSubtractedMean)
                print "checking means for muons", imod, isec, mean, self.hist[hname].GetEntries()

                # bin=  hist_merge_MuonSignal_Noise.FindBin(imod, isec)
                # Muon_Noise_signal = hname + hNoise
                #hist_merge_MuonSignal_Noise.SetBinContent(bin,Muon_Noise_signal)

                #From 16 channels in one module, only 10 are filled into histogram. Scale to get per channel noise...
                self.hist[hnoise_neighbor].Scale(1./10.) # hnoise_neighbor.Scale(1./10.)AttributeError: 'str' object has no attribute 'Scale'
        
        #self.hist[hname].Scale(1./self.hist[hname].Integral());                          

        # self.hist[hname].Scale(1./self.hist[hname].GetBinWidth(i)); 
                                                        
       
        assert self.nEventsRandom > 0

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

    def finalizeWhenMerged(self):
        olist = self.GetOutputList()
        histos = {}
        for o in olist:
            if not "TH1" in o.ClassName():
                if not "TH2" in o.ClassName():
                    continue
            histos[o.GetName()] = o
            print " TH1/2 histogram in output: ", o.GetName()

        for isec in xrange(0,16):
            for imod in xrange(0,14):
                hnoise_neighbor = 'CastorNoise_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                hnoise_randomtrg = 'CastorNoise_randomtrg_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))

                histos[hnoise_neighbor].Scale( 1./histos["EventCount"].GetBinContent(4) )
                histos[hnoise_randomtrg].Scale( 1./histos["EventCount"].GetBinContent(3) )



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

    slaveParams["maxFileNo"] = maxFileNo
    Step2_Selection.runAll(treeName="MuonCastorVTwo",
                           slaveParameters = slaveParams,
                           sampleList = sampleList,
                           maxFilesMC = None,
                           maxFilesData =None,  
                           nWorkers = 8,
                           outFile = outFileName,
                           verbosity = 2)