from CRABClient.UserUtilities import config
config = config()

config.General.workArea = 'WJetsToLNu_HT'

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
    ##     WJets     ##
    ###################
    config.General.requestName = 'WJetsToLNu'
    config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    config.Data.unitsPerJob = 50000
    submit(config)

    config.General.requestName = 'WJetsToLNu_HT-100to200'
    config.Data.inputDataset = '/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    config.Data.unitsPerJob = 10000
    submit(config)

    config.General.requestName = 'WJetsToLNu_HT-200to400'
    config.Data.inputDataset = '/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    config.Data.unitsPerJob = 8000
    submit(config)

    config.General.requestName = 'WJetsToLNu_HT-400to600'
    config.Data.inputDataset = '/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    config.Data.unitsPerJob = 8000
    submit(config)

    config.General.requestName = 'WJetsToLNu_HT-600toInf'
    config.Data.inputDataset = '/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    config.Data.unitsPerJob = 8000
    submit(config)
