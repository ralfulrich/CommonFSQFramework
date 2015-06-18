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

maxFileNo = -1
outfolder = "/afs/cern.ch/user/c/cbaus/pp13TeV2015/HaloMuons/CMSSW_7_4_4_patch4/src/CommonFSQFramework/Core/test/HaloMuonAna/output/"

class Step2_Selection(CommonFSQFramework.Core.ExampleProofReader.ExampleProofReader):
    def init( self, maxEvents = None):
        global maxFileNo
        firstRun = (maxFileNo == -1)

        self.maxEvents = maxEvents
        self.hist = {}

        #Getting sector RMS and Mean
        if firstRun:
            inputFile = ROOT.TFile(join(outfolder,"mean_rms.root"))
        else:
            inputFile = join(outfolder,"plotsMuonselectioncuts_{n:04d}.root".format(n=maxFileNo))
        hSectorMean = inputFile.Get("data_MinimumBias/hSector_Mean")
        hSectorRMS = inputFile.Get("data_MinimumBias/hSector_RMS")

        #Getting channel RMS and Mean
        hChMean = inputFile.Get("data_MinimumBias/hist_ch_Mean")
        hChRMS =  inputFile.Get("data_MinimumBias/hist_ch_RMS")

        self.hist["2DMuonCountMap"] =  ROOT.TH2D("2DMuonCountMap","2DMuonCountMap", 14, -0.5, 13.5, 16, -0.5, 15.5)
        self.hist["2DMuonNoTriggerCountMap"] =  ROOT.TH2D("2DMuonNoTriggerCountMap","2DMuonNoTriggerCountMap", 14, -0.5, 13.5, 16, -0.5, 15.5)
        self.hist["SecAboveNoiseCount"] =  ROOT.TH1D("SecAboveNoiseCount","SecAboveNoiseCount", 16, -0.5, 15.5)
        self.hist["RunsWithGoodMuons"] =  ROOT.TH1D("SecAboveNoiseCount","SecAboveNoiseCount", 10000, 247000-0.5, 248000-0.5)
        self.hist["Runs"] =  ROOT.TH1D("SecAboveNoiseCount","SecAboveNoiseCount", 10000, 247000-0.5, 248000-0.5)

        histcalibrationname = '2DMuonSignalMap'
        if not firstRun:
            self.hist[histcalibrationname] = inputFile.Get("data_MinimumBias/2DMuonSignalMap")
            print "Extracted histogram from file. Checking entries:",  self.hist[histcalibrationname].GetEntries()
        else: #first time running analyser the calibration constants are all set to 1
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
                #hempty = 'CastorEmptysectors_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                self.hist[hname] = ROOT.TH1D( hname, hname, 100, 1, 0)
                self.hist[ hwithouttrigger] = ROOT.TH1D( hwithouttrigger, hwithouttrigger, 100, 1,0)
                #self.hist[hempty] = ROOT. TH1D( hempty, hempty, 100, 1,0)

        for h in self.hist:
            self.hist[h].Sumw2()
            self.GetOutputList().Add(self.hist[h])


        #get channel energies from input file
        self.ch_mean = [[0 for _ in range(14)] for _ in range(16)]
        self.ch_RMS = [[0 for _ in range(14)] for _ in range(16)]
        self.new_ch_mean = [[0 for _ in range(14)] for _ in range(16)]
        self.new_ch_RMS = [[0 for _ in range(14)] for _ in range(16)]
        for i in xrange(0, 224):
            #how I get Channel information
            sec = i/14
            mod = i%14
            #print sec, mod
            #i = sec*14+module   -> sec0 mod0 -> i=0   ; sec9 mod 4
            self.ch_mean[sec][mod] = hChMean.GetBinContent(i+1) #+1 because of underflow
            self.ch_RMS[sec][mod] = hChRMS.GetBinContent(i+1) #+1 because of underflow
            binnumber = histCalibration.FindBin(sec,mod)
            calibration =  histCalibration.GetBinContent(binnumber)
            self.ch_mean[sec][mod] = calibration* self.ch_mean[sec][mod]
            self.ch_RMS[sec][mod] = calibration* self.ch_RMS[sec][mod]


        #get mean and rms from sector histogram
        self.sec_mean = [0] * 16
        self.sec_RMS = [0] * 16
        self.new_sec_mean = [0] * 16
        self.new_sec_RMS = [0] * 16
        for i in xrange(0,16):
            self.sec_mean[i] = hSectorMean.GetBinContent(i+1) #+1 because of underflow
            self.sec_RMS[i] = hSectorRMS.GetBinContent(i+1) #+1 because of underflow
            binnumber = histCalibration.FindBin(sec)
            calibration =  histCalibration.GetBinContent(binnumber)
            self.sec_mean[i] = calibration* self.sec_mean[i]
            self.sec_RMS[i] = calibration* self.sec_RMS[i]
#do same for rms
        #print "Sector means are: ", self.sec_mean
       # print "Sector RMS are: ", self.sec_RMS
       # inputFile.Close()
       # if inputFile2.IsOpen():
       #     inputFile.Close()

    def analyze(self):
        weight = 1
        num = 0
        goodMuonEvent = False
        # genTracks
        #num = self.fChain.genTracks.size()
        #print num
        #print self.maxEta # see slaveParams below

        histcalibrationname = '2DMuonSignalMap'
        histCalibration= self.hist[histcalibrationname]

        ch_energy = [[0 for _ in range(14)] for _ in range(16)]
        sec_energy =[0] * 16
        for i in xrange(0, self.fChain.CastorRecHitEnergy.size()):
            sec = self.fChain.CastorRecHitSector.at(i)-1
            mod = self.fChain.CastorRecHitModule.at(i)-1
            binnumber = histCalibration.FindBin(sec,mod)
            calibration =  histCalibration.GetBinContent(binnumber)
            ch_energy[sec][mod] += calibration*self.fChain.CastorRecHitEnergy.at(i)
            if [mod,sec] not in badChannelsModSec:
                sec_energy[sec] += calibration*self.fChain.CastorRecHitEnergy.at(i)

            #  print 'sec', sec, 'mod', mod, 'e_ch=', ch_energy, "e_sec", sec_energy[sec]
            #  print "all sector energies:  ", sec_energy

        listSectorsAboveNoise = []
        for i in xrange(0,16):
            sigma = (sec_energy[i] - self.sec_mean[i]) / self.sec_RMS[i]
            listSectorsAboveNoise.append([i,sigma])
        def filterListSigma(l,sigma_th):
            return [[ii,isigma] for ii, isigma in l if isigma > sigma_th]

        #selection on sectors (bad channels already excluded)
        if len(filterListSigma(listSectorsAboveNoise,2)) != 1 or len(filterListSigma(listSectorsAboveNoise,3)) != 1: #only 1 sector above 3 Sigma and all others below 2 Sigma
            return
        muonSec = (filterListSigma(listSectorsAboveNoise,3))[0]
        #print "Sectors above noise", listSectorsAboveNoise
            #if listSectorsAboveNoise[0] != 8:
            #   return

        #what is this in the following lines? please comment code properly
        # isec_pn = ((listSectorsAboveNoise3Sigma[0] + 1) + 16)%16
        # isec_mn = ((listSectorsAboveNoise3Sigma[0] - 1) + 16)%16

        # isec_pns = ((listSectorsAboveNoise3Sigma[0] + 1) + 8)%16
        # isec_mns = ((listSectorsAboveNoise3Sigma[0] - 1) + 8)%16

        # cond1 = listSectorsAboveNoise3Sigma[0] or isec_pn or isec_mn
        # cond2= ((listSectorsAboveNoise3Sigma[0]+8)%16 or isec_pns or isec_mns)

        # for isec in xrange(0,16):
        #     for imod in xrange(0,14):
        #         hempty = 'CastorEmptysectors_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
        #         if not cond1 or not cond2:
        #             self.hist[hempty].Fill(ch_energy[imod][isec])


        #selection on channels above noise
        listChannelsAboveNoise = []
        for i in xrange(0,224):
            sec = self.fChain.CastorRecHitSector.at(i)-1
            mod = self.fChain.CastorRecHitModule.at(i)-1
            if [mod,sec] in badChannelsModSec:
                print "skipping channel", mod, sec
                continue

            if sec == muonSec[0]:
              # ch_energy[sec][mod] += self.fChain.CastorRecHitEnergy.at(i)
               if ch_energy[sec][mod]> self.ch_mean[sec][mod] + 2*self.ch_RMS[sec][mod]:
                  listChannelsAboveNoise.append(i)
              # hchannel = 'Channelenergy_mod_{mod}_sec_{sec}'.format(mod=str(mod), sec=str(sec))
              # self.hist[hchannel].Fill(ch_energy[sec][mod])


        countChannelsAboveNoise = len(listChannelsAboveNoise)
        #print "Channels above noise:", listChannelsAboveNoise

        if countChannelsAboveNoise > 4:
            hasMuonTrigger= self.fChain.trgmuon
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

            if Front_Module + Mid_Module + Rear_Module >= 2: #maybe change to 3
                if hasMuonTrigger:
                    goodMuonEvent = True

        #found an interesting event. now fill histograms for channels above noise
        if goodMuonEvent:
            print "Good event in (sec,mod)", sec, mod, "Front,Mid,Back", Front_Module, Mid_Module, Rear_Module
            self.hist["RunsWithGoodMuons"].Fill(self.fChain.run)
            for i in xrange(0,224):
                sec = self.fChain.CastorRecHitSector.at(i)-1
                mod = self.fChain.CastorRecHitModule.at(i)-1
                hname = 'MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(mod), sec=str(sec))
                self.hist[hname].Fill(ch_energy[sec][mod])
                self.hist["2DMuonCountMap"].Fill(mod,sec)
        else:
            self.hist["Runs"].Fill(self.fChain.run)

        #plot all the good muons
        if goodMuonEvent:
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
            ax.set_title('Sec: {sec} ({sig:.1f}sigma)   Trigger: {trg}   Mod[F/M/R]: {f}/{m}/{r} \n Run: {run}   Event: {evt}'.format(sec=sec, sig=muonSec[1], trg="Yes" if hasMuonTrigger else "no", f=Front_Module, m=Mid_Module, r=Rear_Module, run=self.fChain.run, evt=self.fChain.event), multialignment='center')
            fig.savefig(join(outfolder,trgEvtFileName))


                # hwithouttrigger = 'MuonSignalSecCh_withtoutrigger_mod_{mod}_sec_{sec}'.format(mod=str(mod), sec=str(sec))
                #     self.hist[hwithouttrigger].Fill(ch_energy[sec][mod])
                #     self.hist["2DMuonNoTriggerCountMap"].Fill(mod,sec)

        return 1

    def finalize(self):
        print "Finalize:"
        normFactor = self.getNormalizationFactor()
        print "  applying norm", normFactor
        for h in self.hist:
            self.hist[h].Scale(normFactor)

        inputFile = ROOT.TFile(join(outfolder,"mean_rms.root"))
        hChMean = inputFile.Get("data_MinimumBias/hist_ch_Mean")
        histcalibrationname = '2DMuonSignalMap'
        histcalibration = self.hist[histcalibrationname]
        hreferencename = 'MuonSignalSecCh_mod_3_sec_8' #refernce channel is (counting from one) mod=4 sec=9
        referenceMean= self.hist[hreferencename].GetMean()
        meanReferenceNoise = hChMean.GetBinContent(8*14 + 3 +1) #noise is not estimated well. do iterative procedure from ralf? or random trigger?
#        referenceMean -= meanReferenceNoise

        if referenceMean != 0:
            print "Warning reference channel empty. Not dividing"
        for isec in xrange(0,16):
            for imod in xrange(0,14):
                hname = 'MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod), sec=str(isec))
                mean= self.hist[hname].GetMean()
                i = isec * 14 + imod
                meanNoise = hChMean.GetBinContent(i+1) #noise is not estimated well. do iterative procedure from ralf? or random trigger?
                binnumber = histcalibration.FindBin(imod, isec)
                noiseSubtractedMean = (mean)# - meanNoise)
                if referenceMean != 0:
                    noiseSubtractedMean /= referenceMean
                histcalibration.SetBinContent(binnumber, noiseSubtractedMean)
                print "checking means for muons", imod, isec, mean, self.hist[hname].GetEntries()

if __name__ == "__main__":
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.AutoLibraryLoader.enable()

    sampleList = None # run through all
    maxFilesMC = None # run through all ffiles found

    # debug config:
    # Run printTTree.py alone to get the samples list
    #sampleList = []
    #sampleList.append("QCD_Pt-15to3000_TuneZ2star_Flat_HFshowerLibrary_7TeV_pythia6")
    #maxFilesMC = 1


    slaveParams = {}
    slaveParams["maxEta"] = 2.

    #check which output files already exist
    filenames = [ f for f in listdir(outfolder) if isfile(join(outfolder,f)) ]
    for ifilename in filenames:
        if not "plotsMuonselectioncuts" in ifilename:
            continue
        else:
            print "Found previous output file with name: ", ifilename
            maxFileNo = max(maxFileNo,int(ifilename.split("_")[1].strip(".root")))

    #decide on output filename
    outFileName = join(outfolder,"plotsMuonselectioncuts_{n:04d}.root".format(n=maxFileNo+1))
    if maxFileNo == -1:
        print "No previous output file found"
    else:
        print "Found prevous output file(s). Setting new output file name to:", outFileName

    # use printTTree.py <sampleName> to see what trees are avaliable inside the skim file

    Step2_Selection.runAll(treeName="CastorTree",
           slaveParameters=slaveParams,
           sampleList=sampleList,
           maxFilesMC = maxFilesMC,
           maxFilesData = None,
           nWorkers=8,
           outFile = outFileName,
                           verbosity=2)
