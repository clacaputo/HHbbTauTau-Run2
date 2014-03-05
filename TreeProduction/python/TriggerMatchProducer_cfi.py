import FWCore.ParameterSet.Config as cms
## ---
## PAT trigger matching
## --
matchHLTMuons = cms.EDProducer(
   # matching in DeltaR, sorting by best DeltaR
   "PATTriggerMatcherDRLessByR"
   # matcher input collections
   , src     = cms.InputTag( 'selectedPatMuons' )
   , matched = cms.InputTag( 'patTrigger' )
   # selections of trigger objects
   , matchedCuts = cms.string( 'type( "TriggerMuon" ) && ( path( "HLT_IsoMu24_v*" ) || path( "HLT_IsoMu24" ) )' )
   # selection of matches
   , maxDPtRel   = cms.double( 0.5 ) # no effect here
   , maxDeltaR   = cms.double( 0.5 )
   , maxDeltaEta = cms.double( 0.2 ) # no effect here
   # definition of matcher output
   , resolveAmbiguities    = cms.bool( True )
   , resolveByMatchQuality = cms.bool( True )
  )
