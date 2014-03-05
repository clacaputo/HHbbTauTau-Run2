import FWCore.ParameterSet.Config as cms

treeWriter = cms.EDAnalyzer("TreeMakerModule",
  verbosity  = cms.int32(0),
  createTree = cms.bool(False)
)
