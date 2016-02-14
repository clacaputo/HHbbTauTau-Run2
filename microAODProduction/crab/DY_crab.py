from CRABClient.UserUtilities import config
config = config()

config.General.workArea = 'DYJetsToLL_HT'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/miniAOD_skim_Sync.py'
config.JobType.pyCfgParams = ['sampleType=Spring15MC','computeHT=True']

config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/user/ccaputo/HHbbtautau/Run2/ThirdProduction/HTBinning/' # or '/store/group/<subdir>'
config.Data.publication = True

config.Site.storageSite = 'T2_IT_Bari'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException


    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    ###################
    ##     DYJets    ##
    ###################
    config.General.requestName = 'DYJetsToLL_M-50'
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    config.Data.unitsPerJob = 15000
    submit(config)

    config.General.requestName = 'DYJetsToLL_M-50_HT-100to200'
    config.Data.inputDataset = '/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    config.Data.unitsPerJob = 6000
    submit(config)

    config.General.requestName = 'DYJetsToLL_M-50_HT-200to400'
    config.Data.inputDataset = '/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    config.Data.unitsPerJob = 2000
    submit(config)

    config.General.requestName = 'DYJetsToLL_M-50_HT-400to600'
    config.Data.inputDataset = '/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2/MINIAODSIM'
    config.Data.unitsPerJob = 3000
    submit(config)

    config.General.requestName = 'DYJetsToLL_M-50_HT-600toInf'
    config.Data.inputDataset = '/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    config.Data.unitsPerJob = 2000
    submit(config)
