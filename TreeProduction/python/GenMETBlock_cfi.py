import FWCore.ParameterSet.Config as cms

genMETBlock = cms.EDAnalyzer("GenMETBlock",
    verbosity = cms.int32(0),
    genMETSrc = cms.InputTag('genMetTrue')
)
