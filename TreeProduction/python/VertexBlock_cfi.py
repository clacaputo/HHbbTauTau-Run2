import FWCore.ParameterSet.Config as cms

vertexBlock = cms.EDAnalyzer("VertexBlock",
  verbosity = cms.int32(1),
  vertexSrc = cms.InputTag('offlinePrimaryVerticesWithBS')
)
