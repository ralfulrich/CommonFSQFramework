import FWCore.ParameterSet.Config as cms
import CommonFSQFramework.Core.Util
import os

isData = True

if "TMFSampleName" not in os.environ:
    print "TMFSampleName not found, assuming we are running on MC"
else:
    s = os.environ["TMFSampleName"]
    sampleList=CommonFSQFramework.Core.Util.getAnaDefinition("sam")
    isData =  sampleList[s]["isData"]
    if isData: print "Disabling MC-specific features for sample",s

# for test purpose
# isData = False

process = cms.Process("Treemaker")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Source
process.source = cms.Source("PoolSource",
    # fileNames = cms.untracked.vstring('/store/user/hvanhaev/ZeroBias1/Run2015A-v1_RERECO_Run247324_GR_P_V54_withCustomCond-v1/150608_213851/0000/output_data_rereco_1.root')
    #fileNames = cms.untracked.vstring('/store/user/hvanhaev/MinBias_TuneMonash13_13TeV-pythia8/RunIISpring15DR74-NoPU0T_MCRUN2_740TV0_step2-v2/150610_055012/0000/step2_RAW2DIGI_L1Reco_RECO_1.root')
    #fileNames = cms.untracked.vstring('/store/mc/RunIISpring15DR74/ReggeGribovPartonMC_13TeV-EPOS/GEN-SIM-RECO/NoPURawReco_castor_MCRUN2_74_V8B-v1/10000/BC62D29E-7707-E511-A6D9-AC853D9F5344.root')
    #fileNames = cms.untracked.vstring('/store/data/Run2015A/ZeroBias/RECO/PromptReco-v1/000/247/607/00000/52EA626D-9210-E511-843F-02163E01451D.root')
    fileNames = cms.untracked.vstring('/store/data/Run2015A/MinimumBias/RECO/PromptReco-v1/000/247/403/00000/061A4C94-3B0F-E511-B492-02163E0140E0.root')
    # fileNames = cms.untracked.vstring('/store/hidata/HIRun2013/PAMinBiasUPC/RECO/PromptReco-v1/000/209/948/00000/184CCE65-8C60-E211-B2D4-003048D2BE12.root')
)

# from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValProdTTbarGENSIMRECO
# process.source.fileNames = filesRelValProdTTbarGENSIMRECO

# Geometry and Detector Conditions
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
if isData: process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
if not isData: process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.load("Configuration.StandardSequences.MagneticField_cff")

# for test purpose
# process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_75_V5', '')

print process.GlobalTag.globaltag
# process.load("Configuration.StandardSequences.Reconstruction_cff")

# # get custom CASTOR conditions to mark/remove bad channels
# process.load("CondCore.DBCommon.CondDBSetup_cfi")
# process.CastorDbProducer = cms.ESProducer("CastorDbProducer")

# process.es_ascii = cms.ESSource("CastorTextCalibrations",
#    input = cms.VPSet(
#        cms.PSet(
#            object = cms.string('ChannelQuality'),
#            file = cms.FileInPath('data/customcond/castor/quality__2015.txt')
#        ),
#    )
# )

# process.es_prefer_castor = cms.ESPrefer('CastorTextCalibrations','es_ascii')

# # for MC reproduce the CastorTowers and CastorJets to remove the bad channels there
# if not isData:
#     process.load('RecoLocalCalo.Castor.Castor_cff')
#     process.CastorReReco = cms.Path(process.CastorFullReco)

# produce HF PFClusters
# process.PFClustersHF = cms.Path(process.particleFlowRecHitHF*process.particleFlowClusterHF)




###############################################################################
###############################################################################
## switch to uncheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
#addMETCollection(process, labelName='patMETCalo', metSource='met')
addMETCollection(process, labelName='patMETPF', metSource='pfMetT1')
#addMETCollection(process, labelName='patMETTC', metSource='tcMet') # FIXME: removed from RECO/AOD; needs functionality to add to processing

## uncomment the following line to add different jet collections
## to the event content
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection

# uncomment the following lines to switch to ak4CaloJets in your PAT output
labelAK4Calo = 'AK4Calo'
addJetCollection(
   process,
   labelName = labelAK4Calo,
   jetSource = cms.InputTag('ak4CaloJets'),
   jetCorrections = ('AK7Calo', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-1'), # FIXME: Use proper JECs, as soon as available
   btagDiscriminators = [
       'pfJetBProbabilityBJetTags'
     , 'pfJetProbabilityBJetTags'
     , 'pfTrackCountingHighPurBJetTags'
     , 'pfTrackCountingHighEffBJetTags'
     , 'pfSimpleSecondaryVertexHighEffBJetTags'
     , 'pfSimpleSecondaryVertexHighPurBJetTags'
     , 'pfCombinedInclusiveSecondaryVertexV2BJetTags'
     ]
   )
#process.out.outputCommands.append( 'drop *_selectedPatJets%s_pfCandidates_*'%( labelAK4Calo ) )
## JetID works only with RECO input for the CaloTowers (s. below for 'process.source.fileNames')
#process.patJets.addJetID=True
#process.load("RecoJets.JetProducers.ak4JetID_cfi")
#process.patJets.jetIDMap="ak4JetID"
process.patJetsAK4Calo.useLegacyJetMCFlavour=True # Need to use legacy flavour since the new flavour requires jet constituents which are dropped for CaloJets from AOD
###############################################################################
###############################################################################

###############################################################################
# HLT path filter
process.hltcastormuon = cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     HLTPaths = cms.vstring('HLT_L1CastorMuon_v1','HLT_L1Tech59_CASTORHaloMuon_v1'), # provide list of HLT paths (or patterns) you want
     eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key #HLT_MinBiasBSC # HLT_L1Tech_BSC_minBias
     andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
     throw = cms.bool(False)    # throw exception on unknown path names
) 
###############################################################################
###############################################################################

# Here starts the CFF specific part
import CommonFSQFramework.Core.customizePAT
process = CommonFSQFramework.Core.customizePAT.customize(process)

# GT customization
process = CommonFSQFramework.Core.customizePAT.customizeGT(process)

# define treeproducer
process.MuonCastorVTwo = cms.EDAnalyzer("CFFTreeProducer")

process.FiltererdTree = cms.Path(process.hltcastormuon*process.MuonCastorVTwo)

import CommonFSQFramework.Core.VerticesViewsConfigs
import CommonFSQFramework.Core.CaloRecHitViewsConfigs
import CommonFSQFramework.Core.CaloTowerViewsConfigs
import CommonFSQFramework.Core.CastorViewsConfigs
import CommonFSQFramework.Core.PFObjectsViewsConfigs
import CommonFSQFramework.Core.TriggerResultsViewsConfigs
import CommonFSQFramework.Core.JetViewsConfigs

if not isData:
    import CommonFSQFramework.Core.GenLevelViewsConfigs
    

# process.MuonCastorVTwo._Parameterizable__setParameters(CommonFSQFramework.Core.VerticesViewsConfigs.get(["VerticesView"]))
# process.MuonCastorVTwo._Parameterizable__setParameters(CommonFSQFramework.Core.CaloRecHitViewsConfigs.get(["HBHERecHitView","HFRecHitView"]))
# process.MuonCastorVTwo._Parameterizable__setParameters(CommonFSQFramework.Core.CaloTowerViewsConfigs.get(["CaloTowerView"]))
# process.MuonCastorVTwo._Parameterizable__setParameters(CommonFSQFramework.Core.CastorViewsConfigs.get(["ak5CastorJetView"]))
process.MuonCastorVTwo._Parameterizable__setParameters(CommonFSQFramework.Core.CastorViewsConfigs.get(["CastorRecHitViewFull"]))
# process.MuonCastorVTwo._Parameterizable__setParameters(CommonFSQFramework.Core.PFObjectsViewsConfigs.get(["PFCandidateView","ecalPFClusterView","hcalPFClusterView","hfPFClusterView"]))
process.MuonCastorVTwo._Parameterizable__setParameters(CommonFSQFramework.Core.TriggerResultsViewsConfigs.get(["CastorSpecialMuonTriggerResultsView","L1GTriggerResultsView"]))
#process.MuonCastorVTwo._Parameterizable__setParameters(CommonFSQFramework.Core.JetViewsConfigs.get(["JetViewAK4Calo"]))

if not isData:
    process.MuonCastorVTwo._Parameterizable__setParameters(CommonFSQFramework.Core.GenLevelViewsConfigs.get(["GenPartView"]))
    process.MuonCastorVTwo._Parameterizable__setParameters(CommonFSQFramework.Core.GenLevelViewsConfigs.get(["ak4GenJetView"]))

# # add paths
# if not isData:
#     process = CommonFSQFramework.Core.customizePAT.addPath(process, process.CastorReReco)

# process = CommonFSQFramework.Core.customizePAT.addPath(process, process.PFClustersHF)
# process = CommonFSQFramework.Core.customizePAT.addTreeProducer(process, process.MuonCastorVTwo)
process = CommonFSQFramework.Core.customizePAT.addPath(process, process.FiltererdTree)
