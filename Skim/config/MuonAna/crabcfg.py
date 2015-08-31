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
config.Data.unitsPerJob = 200000
#config.Data.totalUnits = 1000 # havent worked last time, use lumi mask?
#config.Data.totalUnits = 10000 # havent worked last time, use lumi mask?
config.Data.lumiMask = "CommonFSQFramework/Skim/lumi/MinBias_CastorMuonRuns_v2.json"
#config.Data.dbsUrl = "global"


config.Data.publication = False
config.Data.publishDataName = 'HW_CastorMuon'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/cwohrman/CFF/CastorMuon'
