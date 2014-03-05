import FWCore.ParameterSet.Config as cms

muonBlock = cms.EDAnalyzer("MuonBlock",
  verbosity       = cms.int32(0),
  muonSrc         = cms.InputTag('muonVariables'),
  vertexSrc       = cms.InputTag('offlinePrimaryVertices'),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
  beamSpotCorr    = cms.bool(True),
  muonID          = cms.string('GlobalMuonPromptTight')
)
