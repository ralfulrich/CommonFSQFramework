import FWCore.ParameterSet.Config as cms
import CommonFSQFramework.Core.Util
import os
import sys

useCustomCond = False
rawOrStreamer = True # True: RAW, False: streamer

# this is important to get the right trigger setup
from Configuration.StandardSequences.Eras import eras
process = cms.Process("Treemaker", eras.Run2_2018_pp_on_AA) # eras.Run2_2018)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_Express_v2', '')

from os import walk

#if "RU_JOBINPUT" not in os.environ:
#    print ("you need to specif RU_JOBINPUT")
#    sys.exit(0)

# Source
if not rawOrStreamer:
    process.source = cms.Source(
        "NewEventStreamFileReader",
        fileNames = cms.untracked.vstring(
            #        'file:' + os.environ['RU_JOBINPUT'] 
            'file:@RU_JOBINPUT@'
            )
        )
else:
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
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )

# HLT path filter 
import HLTrigger.HLTfilters.hltHighLevel_cfi 
if not rawOrStreamer:
    process.TriggerFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
        TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
        #HLTPaths = ['HLT_L1CastorMuon_v*','HLT_L1Tech59_CASTORHaloMuon_v*'], #  # provide list of HLT paths (or patterns) you want
        #HLTPaths = ['HLT_HIRandom_v*'],
        #HLTPaths = ['HLT_HIZeroBias_v*'],
        HLTPaths = ['*'],
        #andOr = cms.bool(True),   # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
        throw = True  ## NNED THis FOR 2015 June since the HLT trigger was renamed! 
        ) 
else:
    process.TriggerFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
        TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
        #HLTPaths = ['HLT_L1CastorMuon_v*','HLT_L1Tech59_CASTORHaloMuon_v*'], #  # provide list of HLT paths (or patterns) you want
        #HLTPaths = ['HLT_HIRandom_v*'],
        #HLTPaths = ['HLT_HIZeroBias_v*'],
        HLTPaths = ['*'],
        #andOr = cms.bool(True),   # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
        throw = True  ## NNED THis FOR 2015 June since the HLT trigger was renamed! 
        ) 
    
    

# get custom CASTOR conditions to mark/remove bad channels
if useCustomCond:
    process.load("CondCore.DBCommon.CondDBSetup_cfi")
    process.CastorDbProducer = cms.ESProducer("CastorDbProducer")

    process.es_ascii = cms.ESSource("CastorTextCalibrations",
                                    input = cms.VPSet(
            cms.PSet(
                object = cms.string('ChannelQuality'),
                file = cms.FileInPath('data/customcond/castor/BadChannels2015.txt')
                ),
            )
                                    )

    process.es_prefer_castor = cms.ESPrefer('CastorTextCalibrations','es_ascii')


# Here starts the CFF specific part
import CommonFSQFramework.Core.customizePAT
process = CommonFSQFramework.Core.customizePAT.customize(process)
process = CommonFSQFramework.Core.customizePAT.customizeGT(process)

# define treeproducer
process.CFFTree = cms.EDAnalyzer("CFFTreeProducer")


#import CommonFSQFramework.Core.CaloRecHitViewsConfigs
import CommonFSQFramework.Core.CaloTowerViewsConfigs
import CommonFSQFramework.Core.CastorViewsConfigs
import CommonFSQFramework.Core.CastorDigiViewConfigs
#import CommonFSQFramework.Core.PFObjectsViewsConfigs
import CommonFSQFramework.Core.TriggerResultsViewsConfigs


process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CastorViewsConfigs.get(["CastorRecHitViewFull"]))
#process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CastorViewsConfigs.get(["CastorRecHitViewFull","CastorTowerView"]))
process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CastorDigiViewConfigs.get(["CastorDigiView"]))
process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.TriggerResultsViewsConfigs.get(["L1GTriggerResultsView"]))
#process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CaloRecHitViewsConfigs.get(["HFRecHitView"]))
process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CaloTowerViewsConfigs.get(["CaloTowerView"]))
#process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.PFObjectsViewsConfigs.get(["PFCandidateView","ecalPFClusterView","hcalPFClusterView","hfPFClusterView"]))


#from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
#MassReplaceInputTag(process, new="rawDataRepacker", old="rawDataCollector")

#if not rawOrStreamer:
process.castorDigis.InputLabel           = cms.InputTag("rawDataRepacker")
process.csctfDigis.producer              = cms.InputTag("rawDataRepacker")
process.dttfDigis.DTTF_FED_Source        = cms.InputTag("rawDataRepacker")
process.ecalDigis.InputLabel             = cms.InputTag("rawDataRepacker")
process.ecalPreshowerDigis.sourceTag     = cms.InputTag("rawDataRepacker")
process.gctDigis.inputLabel              = cms.InputTag("rawDataRepacker")
process.gtDigis.DaqGtInputTag            = cms.InputTag("rawDataRepacker")
process.hcalDigis.InputLabel             = cms.InputTag("rawDataRepacker")
process.muonCSCDigis.InputObjects        = cms.InputTag("rawDataRepacker")
process.muonDTDigis.inputLabel           = cms.InputTag("rawDataRepacker")
process.muonRPCDigis.InputLabel          = cms.InputTag("rawDataRepacker")
process.scalersRawToDigi.scalersInputTag = cms.InputTag("rawDataRepacker")
process.siPixelDigis.InputLabel          = cms.InputTag("rawDataRepacker")
process.siStripDigis.ProductLabel        = cms.InputTag("rawDataRepacker")
process.twinMuxStage2Digis.DTTM7_FED_Source = cms.InputTag("rawDataRepacker")
#process.l1tStage2Fed.rawTag = cms.InputTag('rawDataRepacker')
process.omtfStage2Digis.inputLabel   = cms.InputTag("rawDataRepacker")
process.emtfStage2Digis.InputLabel   = cms.InputTag("rawDataRepacker")
process.gtStage2Digis.InputLabel     = cms.InputTag("rawDataRepacker")
process.gmtStage2Digis.InputLabel    = cms.InputTag("rawDataRepacker")
process.caloStage2Digis.InputLabel   = cms.InputTag("rawDataRepacker")
process.RPCTwinMuxRawToDigi.inputTag = cms.InputTag("rawDataRepacker")
process.caloLayer1Digis.InputLabel   = cms.InputTag("rawDataRepacker") 
process.bmtfDigis.InputLabel         = cms.InputTag("rawDataRepacker")

#from EventFilter.RPCRawToDigi.rpcTwinMuxRawToDigi_cfi import rpcTwinMuxRawToDigi
#process.rpcCPPFRawToDigi.inputTag = 'rawDataRepacker'
#from EventFilter.L1TRawToDigi.gtStage2Digis_cfi import gtStage2Digis


#from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
#MassReplaceInputTag(process, new="rawDataRepacker", old="rawDataCollector")


process.dump=cms.EDAnalyzer('EventContentAnalyzer')


#process.CFF = cms.Path(process.dump * process.TriggerFilter * process.CFFTree)
process.CFF = cms.Path(process.TriggerFilter * process.CFFTree)

##process.raw2digi_custom_step = cms.Path(process.RawToDigi) # L1TRawToDigi*process.castorDigis)
#process.raw2digi_custom_step = cms.Path(process.L1TRawToDigi*process.castorDigis)
#process.raw2digi_custom_step = cms.Path(process.L1TRawToDigi)
process.raw2digi_custom_step = cms.Path(process.gtStage2Digis*process.castorDigis)
# #process.raw2digi_step = cms.Path(process.castorDigis)
#process.castorreco_step = cms.Path(process.reconstruction)
process.castorreco_step = cms.Path(process.castorreco) # *process.CastorFullReco)

process = CommonFSQFramework.Core.customizePAT.addPath(process, process.raw2digi_custom_step)
process = CommonFSQFramework.Core.customizePAT.addPath(process, process.castorreco_step)
process = CommonFSQFramework.Core.customizePAT.addPath(process, process.CFF)
