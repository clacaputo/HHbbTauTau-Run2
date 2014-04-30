import FWCore.ParameterSet.Config as cms

tauBlock = cms.EDAnalyzer("TauBlock",
    patTauSrc = cms.InputTag('patTausTriggerMatch'),
)
