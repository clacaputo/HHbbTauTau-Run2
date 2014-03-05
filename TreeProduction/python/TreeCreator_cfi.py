import FWCore.ParameterSet.Config as cms

treeCreator = cms.EDAnalyzer("TreeMakerModule",
  verbosity  = cms.int32(0),
  createTree = cms.bool(True)
)
