import FWCore.ParameterSet.Config as cms

process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for PHYS14 scenario PU4bx50 : global tag is ???
#    for PHYS14 scenario PU20bx25: global tag is PHYS14_25_V1
#  as a rule, find the global tag in the DAS under the Configs for given dataset
#process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

inputFilesAOD = cms.untracked.vstring(
    # AOD test files from /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/AODSIM
       '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/002F7FDD-BA13-E511-AA63-0026189437F5.root',
       '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00610EE7-C213-E511-842C-00304833529A.root',
       '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00623DCC-A813-E511-A302-0025905B85A2.root',
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM
    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root',
    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00C0BECF-6F14-E511-96F8-0025904B739A.root',
    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0260F225-7614-E511-A79F-00A0D1EE8EB4.root',
    )

personalFileMiniAOD = cms.untracked.vstring('file:allBranches_eleIDs_signal_v2.root',)

# Set up input/output depending on the format
# You can list here either AOD or miniAOD files, but not both types mixed
#
useAOD = False

if useAOD == True :
    inputFiles = inputFilesAOD
    outputFile = "electron_ntuple.root"
    print("AOD input files are used")
else :
    inputFiles = personalFileMiniAOD
    outputFile = "tauID_ntuple.root"
    print("MiniAOD input files are used")
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             

#
# Set up electron ID (VID framework)
#

#from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
## turn on VID producer, indicate data format  to be
## DataFormat.AOD or DataFormat.MiniAOD, as appropriate
#if useAOD == True :
#    dataFormat = DataFormat.AOD
#else :
#    dataFormat = DataFormat.MiniAOD

#switchOnVIDElectronIdProducer(process, dataFormat)

## define which IDs we want to produce
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']

##add them to the VID producer
#for idmod in my_id_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#
# Configure an example module for user analysis with electrons
#
print "="*50
print "Setting the ntupler"
print "="*50

process.synctupler = cms.EDAnalyzer('SyncTreeProducer',

                                 genParticles = cms.InputTag("genParticles"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #

                                 tauSrc  = cms.InputTag("slimmedTaus"),
                                 muonSrc    = cms.InputTag("slimmedMuons"),
                                 vtxSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 bits = cms.InputTag("TriggerResults","","HLT"),
                                 prescales = cms.InputTag("patTrigger"),
                                 objects = cms.InputTag("selectedPatTrigger"),
                                )

process.ntupler = cms.EDAnalyzer('ElectronsIDAnalyzer',
                                 # The module automatically detects AOD vs miniAOD, so we configure both
                                 #
                                 # Common to all formats objects
                                 #
                                 # ... none ...
                                 #
                                 electronBool = cms.untracked.bool(False), ##Switch btw electron and tau ntupler
                                                                           ## TO BE IMPROVED
                                 # Objects specific to AOD format
                                 #
                                 electrons    = cms.InputTag("gedGsfElectrons"),
                                 genParticles = cms.InputTag("genParticles"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                 tauSrc  = cms.InputTag("slimmedTaus"),
                                 #
                                 # ID decisions (common to all formats)
                                 eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                                 eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                                 #
                                 # ValueMaps with MVA results
                                 #
                                 mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                                 mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories")
                                )



process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

print "="*50
print "Starting the ntupler"
print "="*50

# Make sure to add the ID sequence upstream from the user analysis module
#process.p = cms.Path(process.egmGsfElectronIDSequence * process.ntupler)
process.p = cms.Path(process.synctupler)
