from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'RunD_SyncTree'
config.General.workArea = '2015Data'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/syncTreeProducer_cfg.py'
config.JobType.pyCfgParams = ['isData=True','runOnCrab=True','sampleType=Run2015D']

config.Data.inputDataset = '/SingleMuon/ccaputo-crab_RunD-ce4e02672c88fc064906a22e024aadae/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 8000
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/RunOnPublishedDataset/' # or '/store/group/<subdir>'
config.Data.publication = False

config.Site.storageSite = 'T2_IT_Bari'
