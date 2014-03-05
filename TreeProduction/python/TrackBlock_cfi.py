import FWCore.ParameterSet.Config as cms

trackBlock = cms.EDAnalyzer("TrackBlock",
  verbosity       = cms.int32(0),
  trackSrc        = cms.InputTag('generalTracks'),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot')
)
