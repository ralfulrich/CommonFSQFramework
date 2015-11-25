#!/usr/bin/env python
import CommonFSQFramework.Core.ExampleProofReader
from BadChannels2015 import badChannelsSecMod

import sys, os, time
sys.path.append(os.path.dirname(__file__))

import ROOT

ROOT.gROOT.SetBatch(True)
from ROOT import edm
from math import sqrt


from array import *

class Step1_Mean_RMS(CommonFSQFramework.Core.ExampleProofReader.ExampleProofReader):
    def init( self, maxEvents= None): #for all events put maxEvents to None

        self.maxEvents = maxEvents

        self.hist = {}
        self.nEvents = float(0) #float for division
        self.nEventsRandom = float(0)
        
        for isec in xrange(0,16):
            for imod in xrange(0,14):
                hname ='ChannelNoise_sec_{sec}_mod_{mod}'.format(sec=str(isec+1),mod=str(imod+1))
                self.hist[hname] = ROOT. TH1D( hname, hname, 500, -50,200)

            hname = 'SectorNoise_sec_{sec}'.format(sec=str(isec+1))
            self.hist[hname] = ROOT.TH1D( hname, hname, 500, -50, 200)
            
        self.new_ch_mean = [[0 for _ in range(14)] for _ in range(16)]
        self.new_ch_RMS = [[0 for _ in range(14)] for _ in range(16)]
        self.new_sec_mean = [0] * 16
        self.new_sec_RMS = [0] * 16

        #TProfile will only be filled in finalise. and when jobs are merged, the averge is taken
        self.hist['hist_sec_Mean'] = ROOT.TProfile('hist_sec_Mean','hist_sec_Mean',16,-0.5,15.5)
        self.hist['hist_sec_RMS'] = ROOT.TProfile('hist_sec_RMS','hist_sec_RMS',16,-0.5,15.5)
        self.hist['hist_ch_Mean'] = ROOT.TProfile('hist_ch_Mean','hist_ch_Mean',224,-0.5,223.5)
        self.hist['hist_ch_RMS'] = ROOT.TProfile('hist_ch_RMS','hist_ch_RMS',224,-0.5,223.5)

        for h in self.hist:
            self.hist[h].Sumw2()
            self.GetOutputList().Add(self.hist[h])


    #called for every event
    def analyze(self): # for each event
        self.nEvents += 1.
        if self.nEvents%1000 == 0:
            print "analyze(): Event:", self.nEvents

        if self.fChain.trgl1L1GTTech[1] or self.fChain.trgl1L1GTTech[2]:
            return 0
        if not self.fChain.trgRandom:
            return 0


        Runs = self.fChain.run
        self.nEventsRandom += 1

        weight = 1
        num = 0
        # genTracks
        #num = self.fChain.genTracks.size()
        #print num
        #print self.maxEta # see slaveParams below
        #self.hist["numGenTracks"].Fill(num, weight)
        sec_sum = [0] * 16
        for i in xrange(0, self.fChain.CastorRecHitEnergy.size()):
            # how I get Channel information
            sec = self.fChain.CastorRecHitSector.at(i)-1
            mod = self.fChain.CastorRecHitModule.at(i)-1
            
            ich_energy = self.fChain.CastorRecHitEnergy.at(i)

            self.new_ch_mean[sec][mod] += ich_energy
            self.new_ch_RMS[sec][mod] += ich_energy**2

            hname = 'ChannelNoise_sec_{sec}_mod_{mod}'.format(sec=str(sec+1),mod=str(mod+1))
            self.hist[hname].Fill(ich_energy)

            if [sec+1,mod+1] not in badChannelsSecMod:
                sec_sum[sec] += ich_energy

                self.new_sec_mean[sec] += ich_energy
                self.new_sec_RMS[sec] += ich_energy**2
             
        for isec in xrange(0,16):
            hname = 'SectorNoise_sec_{sec}'.format(sec=str(isec+1))
            self.hist[hname].Fill(sec_sum[isec])

        return 1

    #called for each worker unit but after event loop
    def finalize(self):
        print "finalize(): called"
        # normFactor = self.getNormalizationFactor()
        # print "  applying norm", normFactor
        # for h in self.hist:
        #     self.hist[h].Scale(normFactor)

        if self.nEventsRandom <= 0:
            print "Error: Worker did not process any events"

        for isec in xrange(0,16):
            if self.nEventsRandom > 0:
                print "Bin ", isec, ": ", self.new_sec_mean[isec], self.new_sec_mean[isec] / self.nEventsRandom
                self.new_sec_mean[isec] /= float(self.nEventsRandom)
                self.new_sec_RMS[isec] = sqrt(self.new_sec_RMS[isec]/float(self.nEventsRandom) - self.new_sec_mean[isec]**2)
                self.hist['hist_sec_Mean'].Fill(isec, self.new_sec_mean[isec] )
                self.hist['hist_sec_RMS'].Fill(isec, self.new_sec_RMS[isec] )


        for isec in xrange(0,16):
            for imod in xrange(0,14):
                if self.nEventsRandom>0:
                    i = isec * 14 + imod
                    self.new_ch_mean[isec][imod] /= float(self.nEventsRandom)
                    self.new_ch_RMS[isec][imod] = sqrt(self.new_ch_RMS[isec][imod]/float(self.nEventsRandom) - self.new_ch_mean[isec][imod]**2)
                    self.hist['hist_ch_Mean'].Fill(i, self.new_ch_mean[isec][imod] )
                    self.hist['hist_ch_RMS'].Fill(i, self.new_ch_RMS[isec][imod] )






    #called after combination of all worker unit
#    def finalizeWhenMerged(self):
#        print "nothing to do here"



if __name__ == "__main__":
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.AutoLibraryLoader.enable()
    sampleList = []
    sampleList.append("data_MinimumBias_MuonHLTSkim_Run2015A")
    #sampleList = None # run through all
    maxFilesMC = None # run through all ffiles found
    maxFilesData = None # same
    nWorkers = 8 # Use all cpu cores

    # debug config:
    # Run printTTree.py alone to get the samples list
    #sampleList = []
    #sampleList.append("QCD_Pt-15to3000_TuneZ2star_Flat_HFshowerLibrary_7TeV_pythia6")
    #maxFilesMC = 1
    #maxFilesData = 1
    #maxFilesData = 1
   # nWorkers = 1


    slaveParams = {}
    #slaveParams["maxEta"] = 2.


    # use printTTree.py <sampleName> to see what trees are avaliable inside the skim file
    Step1_Mean_RMS.runAll(treeName="MuonCastorVTwo",
           slaveParameters=slaveParams,
           sampleList=sampleList,
           maxFilesMC = maxFilesMC,
           maxFilesData = maxFilesData,
           nWorkers=nWorkers,
           outFile = "output/mean_rms.root" ) #2013 mean_RMS_2, for 2015 mean_RMs
