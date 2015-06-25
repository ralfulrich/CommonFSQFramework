from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'treemaker_CastorMuon.py'

config.section_("Data")
config.Data.inputDataset = '/A/B/C'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased' # alt: LumiBased
config.Data.unitsPerJob = 100000
#config.Data.totalUnits = 1000 # havent worked last time, use lumi mask?
#config.Data.totalUnits = 10000 # havent worked last time, use lumi mask?
config.Data.lumiMask = "CommonFSQFramework/Skim/lumi/MinBias_CastorMuonRuns_v2.json"
#config.Data.dbsUrl = "global"


config.Data.publication = False
config.Data.publishDataName = 'HW_CastorMuon'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/cwohrman/CFF/CastorMuon'
config.General.workArea='MuonAna_0T_20150624'
config.General.requestName='MuonAna_0T_20150624_data_MinimumBias_Run2015A'
config.Data.publishDataName='MuonAna_0T_20150624_data_MinimumBias_Run2015A'
config.Data.inputDataset='/MinimumBias/Run2015A-PromptReco-v1/RECO'
config.Data.splitting='LumiBased'
config.Data.unitsPerJob=10
config.Data.lumiMask='/afs/cern.ch/work/c/cwohrman/CastorTriggerL1EmulationTest/ForJetReco/CMSSW_7_5_0_pre5/src/CommonFSQFramework/Skim/lumi/MinBias_CastorMuonRuns_v2.json'
