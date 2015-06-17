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

maxFileNo = -1

class Step2_Selection(CommonFSQFramework.Core.ExampleProofReader.ExampleProofReader):
    def init( self, maxEvents = None):
        global maxFileNo

        self.maxEvents = maxEvents
        self.hist = {}

#Sector RMS and Mean
        inputFile = ROOT.TFile("/afs/cern.ch/work/m/makbiyik/public/HALOMUON/CMSSW_7_3_5/src/CommonFSQFramework/Core/test/HalomuonAnalyzer/d15Jun/mean_rms.root")
        hSectorMean = inputFile.Get("data_MinimumBias/hSector_Mean")
        hSectorRMS = inputFile.Get("data_MinimumBias/hSector_RMS")

#Channel RMS and Mean

        hChMean = inputFile.Get("data_MinimumBias/hist_ch_Mean")
        hChRMS =  inputFile.Get("data_MinimumBias/hist_ch_RMS")

        self.hist["2DMuonCountMap"] =  ROOT.TH2D("2DMuonCountMap","2DMuonCountMap", 14, -0.5, 13.5, 16, -0.5, 15.5)
        self.hist["2DMuonNoTriggerCountMap"] =  ROOT.TH2D("2DMuonNoTriggerCountMap","2DMuonNoTriggerCountMap", 14, -0.5, 13.5, 16, -0.5, 15.5)
        self.hist["SecAboveNoiseCount"] =  ROOT.TH1D("SecAboveNoiseCount","SecAboveNoiseCount", 16, -0.5, 15.5)
        histcalibrationname = '2DMuonSignalMap'
        if maxFileNo != -1: #if -1 -> does not exist
            inFileName = "output/plotsMuonselectioncuts_{n:04d}.root".format(n=maxFileNo)
            inputFile2 = ROOT.TFile(inFileName)
            self.hist[histcalibrationname] = inputFile2.Get("data_MinimumBias/2DMuonSignalMap")
            print "Extracted histogram from file. Checking entries:",  self.hist[histcalibrationname].GetEntries()
        else: #first time running analyser
            self.hist[histcalibrationname] =  ROOT.TH2D(histcalibrationname,histcalibrationname, 14, -0.5, 13.5, 16, -0.5, 15.5)
            for binx in xrange(0,16):
                for biny in xrange(0,14):
                    self.hist[histcalibrationname].SetBinContent(binx,biny, 1.) #all factors set to 1
            print "Created new 2d calibration map. Checking entries:",  self.hist[histcalibrationname].GetEntries()
        histCalibration= self.hist[histcalibrationname]

        #get sector energies from input file
        for isec in xrange(0,16):
            for imod in xrange(0,14):
                hname = 'MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                hwithouttrigger = 'MuonSignalSecCh_withtoutrigger_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                hempty = 'CastorEmptysectors_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                self.hist[hname] = ROOT.TH1D( hname, hname, 100, 1, 0)
                self.hist[ hwithouttrigger] = ROOT.TH1D( hwithouttrigger, hwithouttrigger, 100, 1,0)
                self.hist[hempty] = ROOT. TH1D( hempty, hempty, 100, 1,0)

        for h in self.hist:
            self.hist[h].Sumw2()
            self.GetOutputList().Add(self.hist[h])


        #get channel energies from input file
        self.ch_mean = [[0 for _ in range(14)] for _ in range(16)]
        self.ch_RMS = [[0 for _ in range(14)] for _ in range(16)]
        for i in xrange(0, 224):
           #how I get Channel information
            sec = i/14
            mod = i%14
       #     print sec, mod
            #i = sec*14+module   -> sec0 mod0 -> i=0   ; sec9 mod 4
            self.ch_mean[sec][mod] = hChMean.GetBinContent(i+1) #+1 because of underflow
            self.ch_RMS[sec][mod] = hChRMS.GetBinContent(i+1) #+1 because of underflow
            binnumber = histCalibration.FindBin(sec,mod)
            calibration =  histCalibration.GetBinContent(binnumber)
            self.ch_mean[sec][mod] = calibration* self.ch_mean[sec][mod]
            self.ch_RMS[sec][mod] = calibration* self.ch_RMS[sec][mod]


         #get mean and rms from sector histogram
        self.sector_mean = [0] * 16
        self.sector_RMS = [0] * 16
        for i in xrange(0,16):
            self.sector_mean[i] = hSectorMean.GetBinContent(i+1) #+1 because of underflow
            self.sector_RMS[i] = hSectorRMS.GetBinContent(i+1) #+1 because of underflow
            binnumber = histCalibration.FindBin(sec)
            calibration =  histCalibration.GetBinContent(binnumber)
            self.sector_mean[i] = calibration* self.sector_mean[i]
            self.sector_RMS[i] = calibration* self.sector_RMS[i]
#do same for rms
        #print "Sector means are: ", self.sector_mean
       # print "Sector RMS are: ", self.sector_RMS
       # inputFile.Close()
        # if inputFile2.IsOpen():
        #     inputFile.Close()

    def analyze(self):
        weight = 1
        num = 0
        # genTracks
        #num = self.fChain.genTracks.size()
        #print num
        #print self.maxEta # see slaveParams below

        histcalibrationname = '2DMuonSignalMap'
        histCalibration= self.hist[histcalibrationname]

        channel_energy = [[0 for _ in range(14)] for _ in range(16)]
        sector_energy =[0] * 16
        for i in xrange(0, self.fChain.CastorRecHitEnergy.size()):
           # how I get Channel information
            sec = self.fChain.CastorRecHitSector.at(i)-1
            mod = self.fChain.CastorRecHitModule.at(i)-1
            binnumber = histCalibration.FindBin(sec,mod)
            calibration =  histCalibration.GetBinContent(binnumber)
            channel_energy[sec][mod] += calibration*self.fChain.CastorRecHitEnergy.at(i)
            sector_energy[sec] += calibration*self.fChain.CastorRecHitEnergy.at(i)

          #  print 'sec', sec, 'mod', mod, 'e_ch=', channel_energy, "e_sec", sector_energy[sec]

      #  print "all sector energies:  ", sector_energy

        listSectorsAboveNoise = []
        for i in xrange(0,16):
            if sector_energy[i] > self.sector_mean[i] + 3*self.sector_RMS[i]:
                listSectorsAboveNoise.append(i)
                self.hist["SecAboveNoiseCount"].Fill(i)
              # print i,  sector_energy[i]


        #selection on sectors (bad channels alread excluded)
        countSectorsAboveNoise = len(listSectorsAboveNoise)
        if countSectorsAboveNoise != 1:
            return
        #if listSectorsAboveNoise[0] != 8:
        #   return

        isec_pn = ((listSectorsAboveNoise[0] + 1) + 16)%16
        isec_mn = ((listSectorsAboveNoise[0] - 1) + 16)%16

        isec_pns = ((listSectorsAboveNoise[0] + 1) + 8)%16
        isec_mns = ((listSectorsAboveNoise[0] - 1) + 8)%16

        cond1 = listSectorsAboveNoise[0] or isec_pn or isec_mn
        cond2= ((listSectorsAboveNoise[0]+8)%16 or isec_pns or isec_mns)

        for isec in xrange(0,16):
            for imod in xrange(0,14):
                hempty = 'CastorEmptysectors_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                if not cond1 or not cond2:
                    self.hist[hempty].Fill(channel_energy[imod][isec])


        #selection on channels above noise
        listChannelsAboveNoise = []
        for i in xrange(0,224):
            sec = self.fChain.CastorRecHitSector.at(i)-1
            mod = self.fChain.CastorRecHitModule.at(i)-1
            if [mod,sec] in badChannelsModSec:
                print "skipping channel", mod, sec
                continue


            if sec == listSectorsAboveNoise[0]:
              # channel_energy[sec][mod] += self.fChain.CastorRecHitEnergy.at(i)
               if channel_energy[sec][mod]> self.ch_mean[sec][mod] + 3*self.ch_RMS[sec][mod]:
                  listChannelsAboveNoise.append(i)
              # hchannel = 'Channelenergy_mod_{mod}_sec_{sec}'.format(mod=str(mod), sec=str(sec))
              # self.hist[hchannel].Fill(channel_energy[sec][mod])


        countChannelsAboveNoise = len(listChannelsAboveNoise)

        if countChannelsAboveNoise >3 :
           Trigger= self.fChain.trgmuon
           Front_Module = False
           Mid_Module= False
           Back_Module =False

           for i in listChannelsAboveNoise:
               sec = self.fChain.CastorRecHitSector.at(i)-1
               mod = self.fChain.CastorRecHitModule.at(i)-1
               if mod in [0,1,2,3]:
                  Front_Module =True

               if mod in [4,5,6,7]:
                  Mid_Module= True

               if mod in [8,9,10,11,12]:
                  Back_Module= True

           for i in listChannelsAboveNoise:
               sec = self.fChain.CastorRecHitSector.at(i)-1
               mod = self.fChain.CastorRecHitModule.at(i)-1
               hname = 'MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(mod), sec=str(sec))
               hwithouttrigger = 'MuonSignalSecCh_withtoutrigger_mod_{mod}_sec_{sec}'.format(mod=str(mod), sec=str(sec))
               if Front_Module + Mid_Module + Back_Module >= 2: #maybe change to 3
                   print "events in ", sec, mod, Trigger, "Front,Mid,Back", Front_Module, Mid_Module, Back_Module
                   if Trigger ==1 :
                     self.hist[hname].Fill(channel_energy[sec][mod])
                     self.hist["2DMuonCountMap"].Fill(mod,sec)
                   else:
                     self.hist[hwithouttrigger].Fill(channel_energy[sec][mod])
                     self.hist["2DMuonNoTriggerCountMap"].Fill(mod,sec)


        return 1

    def finalize(self):
        print "Finalize:"
        normFactor = self.getNormalizationFactor()
        print "  applying norm", normFactor
        for h in self.hist:
            self.hist[h].Scale(normFactor)

        inputFile = ROOT.TFile("/afs/cern.ch/work/m/makbiyik/public/HALOMUON/CMSSW_7_3_5/src/CommonFSQFramework/Core/test/HalomuonAnalyzer/output/mean_rms.root")
        hChMean = inputFile.Get("data_MinimumBias/hist_ch_Mean")
        histcalibrationname = '2DMuonSignalMap'
        histcalibration = self.hist[histcalibrationname]
        hreferencename = 'MuonSignalSecCh_mod_3_sec_8' #refernce channel is (counting from one) mod=4 sec=9
        referenceMean= self.hist[hreferencename].GetMean()
        meanReferenceNoise = hChMean.GetBinContent(8*14 + 3 +1) #noise is not estimated well. do iterative procedure from ralf? or random trigger?
        referenceMean -= meanReferenceNoise

        for isec in xrange(0,16):
            for imod in xrange(0,14):
                hname = 'MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                mean= self.hist[hname].GetMean()
                i = isec * 14 + imod
                meanNoise = hChMean.GetBinContent(i+1) #noise is not estimated well. do iterative procedure from ralf? or random trigger?
                binnumber = histcalibration.FindBin(imod, isec)
                noiseSubtractedMean = (mean - meanNoise)
                if referenceMean != 0:
                    noiseSubtractedMean /= referenceMean
                else:
                    print "Warning reference channel empty"
                histcalibration.SetBinContent(binnumber, noiseSubtractedMean)
                print "checking means for muons", imod, isec, mean, self.hist[hname].GetEntries()

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
   # maxFilesData = 10
   # nWorkers = 8


    slaveParams = {}
    slaveParams["maxEta"] = 2.

    #check which output files already exist
    path = "./"
    filenames = [ f for f in listdir(path) if isfile(join(path,f)) ]
    for ifilename in filenames:
        if not "plotsMuonselectioncuts" in ifilename:
            continue
        else:
            print "Found previous output file with name: ", ifilename
            maxFileNo = max(maxFileNo,int(ifilename.split("_")[1].strip(".root")))

    #decide on output filename
    outFileName = "output/plotsMuonselectioncuts_{n:04d}.root".format(n=maxFileNo+1)
    if maxFileNo == -1:
        print "No previous output file found"
    else:
        print "Found prevous output file(s). Setting new output file name to:", outFileName

    # use printTTree.py <sampleName> to see what trees are avaliable inside the skim file
    Step2_Selection.runAll(treeName="CastorTree",
           slaveParameters=slaveParams,
           sampleList=sampleList,
           maxFilesMC = maxFilesMC,
           maxFilesData = maxFilesData,
           nWorkers=nWorkers,
           outFile = outFileName)
