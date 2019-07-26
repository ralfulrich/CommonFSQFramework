import FWCore.ParameterSet.Config as cms
import CommonFSQFramework.Core.Util
import os
import sys

useCustomCond = False
useRAW = False # True: RAW, False: AOD

# this is important to get the right trigger setup
from Configuration.StandardSequences.Eras import eras
process = cms.Process("Treemaker", eras.Run2_2018_pp_on_AA) # eras.Run2_2018)

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
if (useRAW):
    #process.load("Configuration.StandardSequences.Reconstruction_cff")
    process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
    process.load("Configuration.StandardSequences.RawToDigi_cff")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_Express_v2', '')

from os import walk

#if "RU_JOBINPUT" not in os.environ:
#    print ("you need to specif RU_JOBINPUT")
#    sys.exit(0)

# Source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        'file:@RU_JOBINPUT@'
        )
    )             
        
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )

# HLT path filter 
import HLTrigger.HLTfilters.hltHighLevel_cfi 
process.TriggerFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
    #HLTPaths = ['HLT_L1CastorMuon_v*','HLT_L1Tech59_CASTORHaloMuon_v*'], #  # provide list of HLT paths (or patterns) you want
    #HLTPaths = ['HLT_HIRandom_v*'],
    #HLTPaths = ['HLT_HIZeroBias_v*'],
    HLTPaths = ['*'],
    #andOr = cms.bool(True),   # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = True  ## NNED THis FOR 2015 June since the HLT trigger was renamed! 
    ) 
    
# Here starts the CFF specific part
import CommonFSQFramework.Core.customizePAT
process = CommonFSQFramework.Core.customizePAT.customize(process)
process = CommonFSQFramework.Core.customizePAT.customizeGT(process)

# define treeproducer
process.CFFTree = cms.EDAnalyzer("CFFTreeProducer")

#import CommonFSQFramework.Core.CaloRecHitViewsConfigs
import CommonFSQFramework.Core.CaloTowerViewsConfigs
import CommonFSQFramework.Core.CastorViewsConfigs
#import CommonFSQFramework.Core.CastorDigiViewConfigs
#import CommonFSQFramework.Core.PFObjectsViewsConfigs
import CommonFSQFramework.Core.TriggerResultsViewsConfigs
import CommonFSQFramework.Core.RecoTrackViewsConfigs

process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CastorViewsConfigs.get(["CastorRecHitViewFull"]))
#process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CastorViewsConfigs.get(["CastorRecHitViewFull","CastorTowerView"]))
#process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CastorDigiViewConfigs.get(["CastorDigiView"]))
process.CFFTree._Parameterizable__setParameters({"TriggerResultsViewConfig" : cms.PSet(
            miniView = cms.string("TriggerResultsView"),
            branchPrefix = cms.untracked.string("HLT"),
            #process = cms.string("HLT"), 
            storePrescales = cms.bool(False),
            useGlobalReadoutRecord = cms.bool(False),
            triggers = cms.vstring("ZB", "Random", "CastorMediumJet", "CastorHighJet", "NotBptxOR", "BptxMinusOnly", "BptxPlusOnly"),
            ZB = cms.vstring("HLT_HIZeroBias_v*"),
            Random = cms.vstring("HLT_HIRandom_v*"),
            CastorMediumJet = cms.vstring("HLT_HICastor_MediumJet_BptxAND_v*"),
            CastorHighJet = cms.vstring("HLT_HICastor_HighJet_BptxAND_v*"), 
            NotBptxOR = cms.vstring("HLT_HIL1NotBptxOR_v*"),
            BptxMinusOnly = cms.vstring("HLT_HIL1UnpairedBunchBptxMinus_v*"),
            BptxPlusOnly = cms.vstring("HLT_HIL1UnpairedBunchBptxPlus_v*"),
            )}
    )
#process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CaloRecHitViewsConfigs.get(["HFRecHitView"]))
process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CaloTowerViewsConfigs.get(["CaloTowerView"]))
process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.RecoTrackViewsConfigs.get(["RecoTrackView"]))
#process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.PFObjectsViewsConfigs.get(["PFCandidateView","ecalPFClusterView","hcalPFClusterView","hfPFClusterView"]))

process.dump = cms.EDAnalyzer('EventContentAnalyzer')

if (useRAW) :
    process.raw2digi_step = cms.Path(process.RawToDigi)
    process.reconstruction_step = cms.Path(process.reconstructionHeavyIons)
    #process.reconstruction_step = cms.Path(process.reconstruction)

#process.CFF = cms.Path(process.dump * process.TriggerFilter * process.CFFTree)
process.CFF = cms.Path(process.TriggerFilter * process.CFFTree)
#process.castorreco_step = cms.Path(process.castorreco) # *process.CastorFullReco)

if (useRAW) :
    process = CommonFSQFramework.Core.customizePAT.addPath(process, process.raw2digi_step)
    process = CommonFSQFramework.Core.customizePAT.addPath(process, process.reconstruction_step)
process = CommonFSQFramework.Core.customizePAT.addPath(process, process.CFF)
