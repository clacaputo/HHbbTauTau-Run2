import FWCore.ParameterSet.Config as cms

metBlock = cms.EDAnalyzer("METBlock",
  verbosity = cms.int32(0),
#  metSrc    = cms.InputTag('patMETsPF'),
  metSrc    = cms.InputTag('patPFMETsTypeIcorrected'), 
  corrmetSrc = cms.InputTag('patPFMETsTypeIcorrected'),
  mvametSrc = cms.InputTag('patPFMETsTypeIcorrected')  
)
