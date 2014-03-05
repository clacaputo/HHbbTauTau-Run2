import FWCore.ParameterSet.Config as cms

electronBlock = cms.EDAnalyzer("ElectronBlock",
  verbosity       = cms.int32(0),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
  trackSrc        = cms.InputTag('generalTracks'),
  vertexSrc       = cms.InputTag('offlinePrimaryVertices'),
  electronSrc     = cms.InputTag('electronVariables')
)
