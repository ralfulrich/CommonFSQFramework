#!/usr/bin/env python
import CommonFSQFramework.Core.ExampleProofReader
from BadChannels2015 import badChannelsModSec

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

        for isec in xrange(0,16):
            hname ='SectorsNoise_sec_{sec}'.format(sec=str(isec))
            self.hist[hname] = ROOT. TH1D( hname, hname, 200, 1,0)

        self.new_ch_mean = [[0 for _ in range(14)] for _ in range(16)]
        self.new_ch_RMS = [[0 for _ in range(14)] for _ in range(16)]
        self.new_sec_mean = [0] * 16
        self.new_sec_RMS = [0] * 16

        for h in self.hist:
            self.hist[h].Sumw2()
            self.GetOutputList().Add(self.hist[h])


    #called for every event
    def analyze(self): # for each event

#        Trigger= self.fChain.trgHalomuon
#        if Trigger!=1 :
#            return 0

        weight = 1
        num = 0
        self.nEvents += 1.
        if self:nEvents%1000 == 0:
            print "analyze(): Event:", self.nEvents
        # genTracks
        #num = self.fChain.genTracks.size()
        #print num
        #print self.maxEta # see slaveParams below
        #self.hist["numGenTracks"].Fill(num, weight)

        for i in xrange(0, self.fChain.CastorRecHitEnergy.size()):
            # how I get Channel information
            sec = self.fChain.CastorRecHitSector.at(i)-1
            mod = self.fChain.CastorRecHitModule.at(i)-1
            ich_energy = self.fChain.CastorRecHitEnergy.at(i)
            hname = 'SectorsNoise_sec_{sec}'.format(sec=str(sec))
            self.hist[hname].Fill(ich_energy)

            self.new_ch_mean[sec][mod] += ich_energy
            self.new_ch_RMS[sec][mod] += ich_energy**2
            if [mod,sec] not in badChannelsModSec:
                self.new_sec_mean[sec] += ich_energy
                self.new_sec_RMS[sec] += ich_energy**2


        return 1

    #called for each worker unit but after event loop
    def finalize(self):
        print "finalize(): called"
        normFactor = self.getNormalizationFactor()
        print "  applying norm", normFactor
        for h in self.hist:
            self.hist[h].Scale(normFactor)

        #TProfile will only be filled in finalise. and when jobs are merged, the averge is taken
        self.hist['hist_sec_Mean'] = ROOT.TProfile('hist_sec_Mean','hist_sec_Mean',16,-0.5,15.5)
        self.hist['hist_sec_RMS'] = ROOT.TProfile('hist_sec_RMS','hist_sec_RMS',16,-0.5,15.5)
        self.hist['hist_ch_Mean'] = ROOT.TProfile('hist_ch_Mean','hist_ch_Mean',224,-0.5,223.5)
        self.hist['hist_ch_RMS'] = ROOT.TProfile('hist_ch_RMS','hist_ch_RMS',224,-0.5,223.5)


        if self.nEvents <= 0:
            print "Error: Worker did not process any events"

        for i in xrange(0,16):
            if self.nEvents > 0:
                print "Bin ", i, ": ", self.new_sec_mean[i], self.new_sec_mean[i] / self.nEvents
                self.new_sec_mean[i] /= float(self.nEvents)
                self.new_sec_RMS[i] = sqrt(self.new_sec_RMS[i]/float(self.nEvents) - self.new_sec_mean[i]**2)
                self.hist['hist_sec_Mean'].Fill(i, self.new_sec_mean[i] )
                self.hist['hist_sec_RMS'].Fill(i, self.new_sec_RMS[i] )


        for imod in xrange(0,16):
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

    sampleList = None # run through all
    maxFilesMC = None # run through all ffiles found
    maxFilesData = None # same
    nWorkers = None # Use all cpu cores

    # debug config:
    # Run printTTree.py alone to get the samples list
    #sampleList = []
    #sampleList.append("QCD_Pt-15to3000_TuneZ2star_Flat_HFshowerLibrary_7TeV_pythia6")
    #maxFilesMC = 1
    #maxFilesData = 1
   # maxFilesData = 1
    nWorkers = 8


    slaveParams = {}
    #slaveParams["maxEta"] = 2.


    # use printTTree.py <sampleName> to see what trees are avaliable inside the skim file
    Step1_Mean_RMS.runAll(treeName="CastorTree",
           slaveParameters=slaveParams,
           sampleList=sampleList,
           maxFilesMC = maxFilesMC,
           maxFilesData = maxFilesData,
           nWorkers=nWorkers,
           outFile = "output/mean_rms.root" )
