import FWCore.ParameterSet.Config as cms
process = cms.Process("HTauTauTree")
#------------------------
# Message Logger Settings
#------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10


from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('globaltag',
                  'START53_V21::All',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Global Tag to use. Default: START53_V21::All")
options.register ('fileList',
                  'fileList.py',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "List of root files to process.")

options.parseArguments()

#--------------------------------------
# Event Source & # of Events to process
#---------------------------------------
fileList = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = fileList)
execfile(options.fileList)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

#-----------------------------
# Geometry
#-----------------------------
process.load("Configuration.Geometry.GeometryIdeal_cff")
#-----------------------------
# Magnetic Field
#-----------------------------
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#-------------
# Global Tag
#-------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = options.globaltag

#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile) )
#--------------------------------------------------
# VHTauTau Tree Specific
#--------------------------------------------------
process.load("HHbbTauTau.TreeProduction.TreeCreator_cfi")
process.load("HHbbTauTau.TreeProduction.TreeWriter_cfi")
process.load("HHbbTauTau.TreeProduction.TreeContentConfig_cff")

#-------------------------------------------------------
# PAT 
#------------------------------------------------------
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")

import PhysicsTools.PatAlgos.tools.tauTools as tauTools
#tauTools.switchToPFTauHPS(process) # For HPS Taus

## --
## Switch on PAT trigger
## --
import PhysicsTools.PatAlgos.tools.trigTools as trigTools
trigTools.switchOnTrigger( process, outputModule='' ) # This is optional and can be omitted.

process.p = cms.Path(
  process.treeCreator +
  process.treeContentSequence +
  process.treeWriter
)
