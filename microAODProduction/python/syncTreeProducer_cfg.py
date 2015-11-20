import FWCore.ParameterSet.Config as cms
process = cms.Process("HTauTauSyncTree2015")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500


from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('globalTag',
                  'MCRUN2_74_V9::All',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Global Tag to use. Default: MCRUN2_74_V9::All")
options.register ('fileList',
                  'fileList.txt',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "List of root files to process.")
options.register ('fileNamePrefix',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Prefix to add to input file names.")
options.register ('includeSim',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Include Sim. Default: False")

options.parseArguments()

#--------------------------------------
# Event Source & # of Events to process
#---------------------------------------

from HHbbTauTau.RunTools.readFileList import *

process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring())
readFileList(process.source.fileNames, options.fileList, options.fileNamePrefix)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

#-----------------------------
# Geometry
#-----------------------------
process.load("Configuration.StandardSequences.Geometry_cff")
#-----------------------------
# Magnetic Field
#-----------------------------
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#-------------
# Global Tag
#-------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for PHYS14 scenario PU4bx50 : global tag is ???
#    for PHYS14 scenario PU20bx25: global tag is PHYS14_25_V1
#  as a rule, find the global tag in the DAS under the Configs for given dataset
#process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
process.GlobalTag.globaltag = options.globalTag

#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile) )

#-------------
# SyncTree Producer
#-------------
process.synctupler = cms.EDAnalyzer('SyncTreeProducer',

                                 genParticles = cms.InputTag("genParticles"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #

                                 tauSrc    = cms.InputTag("slimmedTaus"),
                                 muonSrc   = cms.InputTag("slimmedMuons"),
                                 vtxSrc    = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 pfMETSrc  = cms.InputTag("slimmedMETs"),
                                 bits      = cms.InputTag("TriggerResults","","HLT"),
                                 prescales = cms.InputTag("patTrigger"),
                                 objects   = cms.InputTag("selectedPatTrigger"),
                                )

process.p = cms.Path(process.synctupler)
