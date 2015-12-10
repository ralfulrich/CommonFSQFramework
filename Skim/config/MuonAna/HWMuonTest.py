#!/usr/bin/env python
import CommonFSQFramework.Core.ExampleProofReader

import sys, os, time
sys.path.append(os.path.dirname(__file__))

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import edm

import numpy as np
import gc
import matplotlib.pyplot as plt

from array import *

class HWMuonTest(CommonFSQFramework.Core.ExampleProofReader.ExampleProofReader):
    def init( self):

        self.hist = {}
        self.hist["heventbx"] =  ROOT.TH1F("heventbx",   "heventbx",  3600, 0, 3600)
        self.hist["hrun"] = ROOT.TH1F("hrun", "hrun", 3000, 0, 3000)
        self.hist["hcasenergy"] = ROOT.TH1F("hcasenergy", "hcasenergy", 10000, 0, 10000)

        # for i in range(100):
        #     hname = "h2_test_mem_{num}".format(num=str(i))
        #     self.hist[hname] = ROOT.TH2F(hname,hname, 1000, 0, 100, 1000, 0, 100)


        for h in self.hist:
            self.hist[h].Sumw2()
            self.GetOutputList().Add(self.hist[h])

        self.rhs = None
        self.rhm = None
        self.rhe = None

    def analyze(self):
        weight = 1
        bx = self.fChain.bx
        run = self.fChain.run

        if self.fChain.CastorRecHitEnergy.size() != 224:
            return -1

        # self.getCastorEnergy()

        self.rhe = np.array( map(lambda x: x,self.fChain.CastorRecHitEnergy) )
        # if self.rhs == None:
        #     self.rhs = np.array( map(lambda x: x,self.fChain.CastorRecHitSector) )
        #     self.rhm = np.array( map(lambda x: x,self.fChain.CastorRecHitSector) )
        sumecas = self.rhe.sum()

        # genTracks
        #num = self.fChain.genTracks.size()
        #print num
        #print self.maxEta # see slaveParams below
        self.hist["heventbx"].Fill(bx, weight)
        self.hist["hrun"].Fill(run-261000, weight)
        self.hist["hcasenergy"].Fill(sumecas, weight)

        # gc.collect()

        # fig = plt.figure()
        # plt.close('all')

        return 1

    def finalize(self):
        print "Finalize:"
        # normFactor = self.getNormalizationFactor()
        # print "  applying norm", normFactor
        # for h in self.hist:
        #     self.hist[h].Scale(normFactor)

    def finalizeWhenMerged(self):

        pass

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
    sampleList = []
    # sampleList.append("data_Cosmics_MuonHLTSkim_2015E_4T")
    sampleList.append("data_Cosmics_MuonHLTSkim_HI2015E_4T")
    # maxFilesMC = 1
    # maxFilesData = 1
    # maxFilesData = 1
    nWorkers = 8


    # slaveParams = {}
    # slaveParams["maxEta"] = 2.


    # use printTTree.py <sampleName> to see what trees are avaliable inside the skim file
    HWMuonTest.runAll(treeName="MuonCastorVTwo",
           # slaveParameters=slaveParams,
           sampleList=sampleList,
           maxFilesMC = maxFilesMC,
           maxFilesData = maxFilesData,
           nWorkers=nWorkers,
           outFile = "plotsHWMuonTest.root" )
