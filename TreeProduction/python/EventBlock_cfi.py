import FWCore.ParameterSet.Config as cms

eventBlock = cms.EDAnalyzer("EventBlock",
    verbosity = cms.int32(0),
    l1InputTag  = cms.InputTag('gtDigis'),
    vertexInputTag = cms.InputTag('offlinePrimaryVerticesWithBS'),
    vertexMinimumNDOF = cms.uint32(4),
    vertexMaxAbsZ = cms.double(24.),
    vertexMaxd0 = cms.double(2.),
    hpTrackThreshold = cms.double(0.25)
)
