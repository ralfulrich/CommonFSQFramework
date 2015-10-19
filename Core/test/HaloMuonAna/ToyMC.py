import ROOT
ROOT.gROOT.SetBatch(True)
from math import sqrt, log10
from array import *
import copy
from EmptyChannels2015 import emptyChannelsSecMod
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

import os
from os.path import isfile, join

class toyMCClass:
    def __init__(self):
        outfolder = os.environ["HaloMuonOutput"]
        self.inputFile = ROOT.TFile(join(outfolder,"input_for_toyMC.root"))
        self.outputFile = ROOT.TFile(join(outfolder,"Eventgenerator.root"),"RECREATE")

        self.hist = {}

        self.hNoise= [[0 for _ in xrange(14)] for _ in xrange(16)]
        self.hMuon = [[0 for _ in xrange(14)] for _ in xrange(16)]
        
        for isec in xrange(0,16):
            for imod in xrange(0,14):
                hnoisename = 'data_MinimumBias_Run2015A/CastorNoise_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                #hist[hnoisename] = ROOT.TH1D(hnoisename, hnoisename, 50, -100, 400)
                self.hNoise[isec][imod] = self.inputFile.Get(hnoisename)
                hname = 'data_MinimumBias_Run2015A/MuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))   
                #self.hist[hname] = ROOT.TH1D(hname, hname, 50,-100, 400)   
                self.hMuon[isec][imod] = self.inputFile.Get(hname)
                h_noisesignal= 'HNoiseSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod+1),sec=str(isec+1))
                self.hist[h_noisesignal]= ROOT.TH1D(h_noisesignal,h_noisesignal,50,-100,400)
                h_muonsignal= 'hMuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
                self.hist[h_muonsignal]= ROOT.TH1D(h_muonsignal,h_muonsignal,50,-100,400)
         
        self.CastorEnergy = [[0 for _ in xrange(14)] for _ in xrange(16)]
 
    def GenerateEvent(self):
        muonsec = np.random.randint(0,16,1)[0]

        # print "Muonsec =",muonsec
        for i in xrange(224):
            isec = i//14
            imod = i%14
            #print "isec=",isec, h_noisesignal
            
            # cas_e_noise = 0
            # cas_e_muon  = 0
            if [isec+1,imod+1] in emptyChannelsSecMod:
               continue
            cas_e_noise = self.hNoise[isec][imod].GetRandom()
            cas_e_muon = self.hMuon[isec][imod].GetRandom()
            

            # print "cas_e_muon && cas_e_noise ", cas_e_muon, cas_e_noise 
            cas_e = cas_e_muon if muonsec == isec else cas_e_noise 

            self.CastorEnergy[isec][imod] = cas_e

            h_noisesignal= 'HNoiseSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod+1),sec=str(isec+1))
            self.hist[h_noisesignal].Fill(cas_e_noise)
            h_muonsignal= 'hMuonSignalSecCh_mod_{mod}_sec_{sec}'.format(mod=str(imod+1), sec=str(isec+1))
            self.hist[h_muonsignal].Fill(cas_e_muon)

    def CastorRecHitEnergy(self,i):
        isec = i//14
        imod = i%14 
        return self.CastorEnergy[isec][imod]
            

           

   