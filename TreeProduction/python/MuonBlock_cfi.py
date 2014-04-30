import FWCore.ParameterSet.Config as cms

muonBlock = cms.EDAnalyzer("MuonBlock",
  muonSrc         = cms.InputTag('patMuonsWithEmbeddedVariables'),
  vertexSrc       = cms.InputTag('patVertices'),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
  beamSpotCorr    = cms.bool(True),
  muonID          = cms.string('GlobalMuonPromptTight')
)
