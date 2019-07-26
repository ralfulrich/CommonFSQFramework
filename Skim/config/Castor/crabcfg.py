from WMCore.Configuration import Configuration 
config = Configuration() 
config.section_("General") 
config.General.workArea = 'crab_projects' 
config.section_("User") 
config.User.voGroup = 'dcms' 
config.section_("JobType") 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'treemaker_CastorRAW.py' 
config.section_("Data") 
config.Data.inputDataset = '/A/B/C' 
config.Data.inputDBS = 'global' 
config.Data.splitting = 'EventAwareLumiBased' 
config.Data.unitsPerJob = 200000 
#config.Data.totalUnits = 1000 
# havent worked last time, use lumi mask? 
#config.Data.totalUnits = 10000 
# havent worked last time, use lumi mask? 
config.Data.lumiMask = "circulating.json" 
#config.Data.dbsUrl = "global" 
config.Data.ignoreLocality = True 
config.Data.publication = False 
# config.Data.outLFNDirBase = 'srm://dgridsrm-fzk.gridka.de:8443/srm/managerv2?SFN=/pnfs/gridka.de/dcms/disk-only/store/user/makbiyik/CastorMuon' 
config.Data.outLFNDirBase = '/store/user/rulrich/HIRun2018A/circulating' 
# config.Data.publishDataName = 'HW_CastorMuon' 
config.section_("Site") 
config.Site.whitelist = ['T2_CH_*','T2_DE_*'] 
# config.Site.storageSite = 'T2_CH_CERN' 
# config.Data.outLFNDirBase = '/store/group/phys_heavyions/cwohrman/CFF/CastorMuon' 
# config.Site.storageSite = 'T2_DE_RWTH' 
config.Site.storageSite = 'T2_DE_DESY'

