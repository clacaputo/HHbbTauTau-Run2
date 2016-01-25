from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'TTbar_ext3_SyncTree'
config.General.workArea = 'TTbar'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/miniAOD_skim_Sync.py'

config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2_ext3-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 40000
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/FirstProduction/' # or '/store/group/<subdir>'
config.Data.publication = True

config.Site.storageSite = 'T2_IT_Bari'
