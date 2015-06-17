#!/usr/bin/env python
import CommonFSQFramework.Core.ExampleProofReader

import sys, os, time
sys.path.append(os.path.dirname(__file__))

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import edm

from array import *

class RunEventsLumi(CommonFSQFramework.Core.ExampleProofReader.ExampleProofReader):
    def init( self, maxEvents= None): #for all events put maxEvents to None
    

        self.maxEvents = maxEvents
        self.hist = {}

      #  self.hist["numEvents"] =  ROOT.TH1F("numEvents", "numEvents",  100, 0, 1)
        self.hist["Runs"] =  ROOT.TH1F("Runs", "Runs",  2000, 246000-0.5, 278000-0.5)
        self.hist["Lumi"] =  ROOT.TH1F("Lumi", "Lumi",  1600, 0-0.5, 1600-0.5)
        self.hist["BC"] =  ROOT.TH1F("BC", "BC",  3600, 0-0.5, 3600-0.5)



        for h in self.hist:
            self.hist[h].Sumw2()
            self.GetOutputList().Add(self.hist[h])

    def analyze(self):
        weight = 1
        num = 0
      
        Run_number = self.fChain.run
        Events = self.fChain.event
        LS = self.fChain.lumi
        Bunch_crossing= self.fChain.bx


        print Events,  Run_number, LS, Bunch_crossing
        #print self.maxEta # see slaveParams below
    
       # self.hist["numEvents"].Fill( Events)
        self.hist["Runs"].Fill(Run_number)
        self.hist["Lumi"].Fill(LS)
        self.hist["BC"].Fill(self.fChain.bx)

        return 1

    def finalize(self):
        print "Finalize:"
        normFactor = self.getNormalizationFactor()
        print "  applying norm", normFactor
        for h in self.hist:
            self.hist[h].Scale(normFactor)

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
    #nWorkers = 1


    slaveParams = {}
    #slaveParams["maxEta"] = 2.


    # use printTTree.py <sampleName> to see what trees are avaliable inside the skim file
    RunEventsLumi.runAll(treeName="CastorTree",
           slaveParameters=slaveParams,
           sampleList=sampleList,
           maxFilesMC = maxFilesMC,
           maxFilesData = maxFilesData,
           nWorkers=nWorkers,
           outFile = "plotsRunEventsLumi.root" )
