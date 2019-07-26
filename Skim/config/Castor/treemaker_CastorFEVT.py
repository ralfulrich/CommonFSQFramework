import FWCore.ParameterSet.Config as cms
import CommonFSQFramework.Core.Util
import os

useCustomCond = False

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

files = []
mypath='/eos/cms/store/t0streamer/Data/HIExpress/000/326/722'
for (dirpath, dirnames, filenames) in walk(mypath):
    for filename in filenames:
        files.append("file:" + dirpath + "/" + filename)
print files

# Source
process.source = cms.Source(
    "NewEventStreamFileReader",
    fileNames = cms.untracked.vstring(
        files
#        'file:/eos/cms/store/t0streamer/Data/HIExpress/000/326/722/run326722_ls0033_streamHIExpress_StorageManager.dat',
#        
        )
    )

# process.source = cms.Source("PoolSource",

#    fileNames = cms.untracked.vstring( # EXPRESS DATA
# #   'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/7F667FFC-BF8C-7445-9EE2-A86D3B4755A0.root',
# #   'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/22585A1A-6D5C-F64F-87EE-B9D0493CDE02.root'

# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/22585A1A-6D5C-F64F-87EE-B9D0493CDE02.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/24DE0AC3-9367-934E-A7F3-EFB1DF42B4D5.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/38576913-F59C-4042-A2D8-0C83F9927DE5.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/3F8B6F5D-A2BE-014D-874F-A7B997C343BC.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/4C3D07C1-E8BC-EC4B-9F8D-32DF82FC4C09.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/53B422B0-BE43-8F48-A1DA-929728C02977.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/5BA9EB73-46F6-704F-96E5-061CF98D1312.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/5FC5EDD5-0329-A04C-A1E1-64B1B98DA741.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/6C41E335-288C-F240-AAB4-AE7E3A3AC078.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/7F667FFC-BF8C-7445-9EE2-A86D3B4755A0.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/8FF2A315-68FB-2345-B5B7-278DD4614256.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/94C684E7-5993-1E41-9B1A-0B2105F39F3C.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/A798E38A-7832-D14D-B8A8-6E18A91D223F.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/A917CD68-42C0-F340-94D8-5167FE74FB2B.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/B50BFCD8-10BA-0141-88EA-7E17934DDD00.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/BC01935D-6BB8-724E-AB49-DFD0E1118008.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/C1BCD67F-ECFD-0344-A61A-56459069F6B9.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/C6C7CB7F-02FE-A54D-9E72-D6A73D8D23A2.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/C7662A3C-96EC-9743-A56A-57E9DF207EF3.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/C99FBAC6-901F-ED40-BFE9-639D8E52C288.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/CD0A99A3-5180-4B49-9B82-3FBB7508CA93.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/CDCF84DA-21DA-9D44-AB47-448613F80C34.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/E30EBBCC-E23F-8141-9D7A-C7CFE9A373C2.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/E840777D-278A-EA46-BACC-B47509BD564D.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/F3FC3B5B-9AE2-344D-8480-CC9955E0C130.root',
# 'file:/eos/cms/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/499/00000/FA991A29-CCFD-D847-B1D5-35EDF5B1E093.root',

#         )

# )                            


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )

 # HLT path filter 
import HLTrigger.HLTfilters.hltHighLevel_cfi 
process.TriggerFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
     TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
     #HLTPaths = ['HLT_L1CastorMuon_v*','HLT_L1Tech59_CASTORHaloMuon_v*'], #  # provide list of HLT paths (or patterns) you want
     HLTPaths = ['HLT_HIRandom_v*'],
     #HLTPaths = ['*'],
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
#import CommonFSQFramework.Core.CaloTowerViewsConfigs
import CommonFSQFramework.Core.CastorViewsConfigs
#import CommonFSQFramework.Core.PFObjectsViewsConfigs
import CommonFSQFramework.Core.TriggerResultsViewsConfigs


process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CastorViewsConfigs.get(["CastorRecHitViewFull"]))
#process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CastorViewsConfigs.get(["CastorRecHitViewFull","CastorTowerView"]))
process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.TriggerResultsViewsConfigs.get(["L1GTriggerResultsView"]))
#process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CaloRecHitViewsConfigs.get(["HFRecHitView"]))
#process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.CaloTowerViewsConfigs.get(["CaloTowerView"]))
#process.CFFTree._Parameterizable__setParameters(CommonFSQFramework.Core.PFObjectsViewsConfigs.get(["PFCandidateView","ecalPFClusterView","hcalPFClusterView","hfPFClusterView"]))


#from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
#MassReplaceInputTag(process, new="rawDataRepacker", old="rawDataCollector")

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
process.omtfStage2Digis.inputLabel = cms.InputTag("rawDataRepacker")
process.emtfStage2Digis.InputLabel = cms.InputTag("rawDataRepacker")
process.gtStage2Digis.InputLabel = cms.InputTag("rawDataRepacker")
process.gmtStage2Digis.InputLabel = cms.InputTag("rawDataRepacker")
process.caloStage2Digis.InputLabel = cms.InputTag("rawDataRepacker")
process.RPCTwinMuxRawToDigi.inputTag = cms.InputTag("rawDataRepacker")
process.caloLayer1Digis.InputLabel = cms.InputTag("rawDataRepacker") 
process.bmtfDigis.InputLabel = cms.InputTag("rawDataRepacker")

#from EventFilter.RPCRawToDigi.rpcTwinMuxRawToDigi_cfi import rpcTwinMuxRawToDigi
#process.rpcCPPFRawToDigi.inputTag = 'rawDataRepacker'
#from EventFilter.L1TRawToDigi.gtStage2Digis_cfi import gtStage2Digis


#from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
#MassReplaceInputTag(process, new="rawDataRepacker", old="rawDataCollector")


process.dump=cms.EDAnalyzer('EventContentAnalyzer')

#process.FiltererdTree = cms.Path(process.dump * process.TriggerFilter * process.CFFTree)
process.FiltererdTree = cms.Path(process.TriggerFilter * process.CFFTree)

##process.raw2digi_custom_step = cms.Path(process.RawToDigi) # L1TRawToDigi*process.castorDigis)
#process.raw2digi_custom_step = cms.Path(process.L1TRawToDigi*process.castorDigis)
#process.raw2digi_custom_step = cms.Path(process.L1TRawToDigi)
process.raw2digi_custom_step = cms.Path(process.gtStage2Digis*process.castorDigis)
# #process.raw2digi_step = cms.Path(process.castorDigis)
#process.castorreco_step = cms.Path(process.reconstruction)
process.castorreco_step = cms.Path(process.castorreco) # *process.CastorFullReco)

process = CommonFSQFramework.Core.customizePAT.addPath(process, process.raw2digi_custom_step)
process = CommonFSQFramework.Core.customizePAT.addPath(process, process.castorreco_step)
process = CommonFSQFramework.Core.customizePAT.addPath(process, process.FiltererdTree)
