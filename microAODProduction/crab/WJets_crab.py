from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'WJetsToLNu'
config.General.workArea = 'WJetsToLNu'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/miniAOD_skim_EleID.py'

config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 50000
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/FirstProduction/' # or '/store/group/<subdir>'
config.Data.publication = True

config.Site.storageSite = 'T2_IT_Bari'
