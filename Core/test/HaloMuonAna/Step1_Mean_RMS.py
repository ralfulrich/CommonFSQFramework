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
        self.sect_mean = []
        self.sect_RMS = []
        self.ch_mean = []
        self.ch_RMS = []
        self.nEvents = float(0) #float for division

        #TProfile will only be filled in finalise. and when jobs are merged, the averge is taken
        self.hist['hSector_Mean'] = ROOT.TProfile('hSector_Mean','hSector_Mean',16,-0.5,15.5)
        self.hist['hSector_RMS'] = ROOT.TProfile('hSector_RMS','hSector_RMS',16,-0.5,15.5)
        self.hist['hist_ch_Mean'] = ROOT.TProfile('hist_ch_Mean','hist_ch_Mean',224,-0.5,223.5)
        self.hist['hist_ch_RMS'] = ROOT.TProfile('hist_ch_RMS','hist_ch_RMS',224,-0.5,223.5)

        for isec in xrange(0,16):
            self.sect_mean.append(0)
            self.sect_RMS.append(0)
            for imod in xrange(0,14):
                hname ='SectorsNoise_sec_{sec}'.format(sec=str(isec))
                self.hist[hname] = ROOT. TH1D( hname, hname, 200, 1,0)


                self.ch_mean.append(0)
                self.ch_RMS.append(0)

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
            if sec in [isec for imod,isec in badChannelsModSec]:
                continue
            channel_energy = self.fChain.CastorRecHitEnergy.at(i)
            hname = 'SectorsNoise_sec_{sec}'.format(sec=str(sec))
            self.hist[hname].Fill(channel_energy)


            self.ch_mean[i] += channel_energy
            self.ch_RMS[i] += channel_energy**2
            self.sect_mean[sec] += channel_energy
            self.sect_RMS[sec] += channel_energy**2


        return 1

    #called for each worker unit but after event loop
    def finalize(self):
        print "finalize(): called"
        normFactor = self.getNormalizationFactor()
        print "  applying norm", normFactor
        for h in self.hist:
            self.hist[h].Scale(normFactor)

#        self.hist['hSector_Mean_RMS'].Fill( sector_mean, sector_rms )

        for i in xrange(0,16):
            if self.nEvents > 0:
                print "Bin ", i, ": ", self.sect_mean[i], self.sect_mean[i] / self.nEvents
                self.sect_mean[i] /= float(self.nEvents)
                self.sect_RMS[i] = sqrt(self.sect_RMS[i]/float(self.nEvents) - self.sect_mean[i]**2)
                self.hist['hSector_Mean'].Fill(i, self.sect_mean[i] )
                self.hist['hSector_RMS'].Fill(i, self.sect_RMS[i] )


        for i in xrange(0,224):
            if self.nEvents>0:
               self.ch_mean[i] /= float(self.nEvents)
               self.ch_RMS[i] = sqrt(self.ch_RMS[i]/float(self.nEvents) - self.ch_mean[i]**2)
               self.hist['hist_ch_Mean'].Fill(i, self.ch_mean[i] )
               self.hist['hist_ch_RMS'].Fill(i, self.ch_RMS[i] )
            else:
                print "Worker did not process any events"






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
    nWorkers = 1


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
