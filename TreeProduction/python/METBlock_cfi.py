import FWCore.ParameterSet.Config as cms

metBlock = cms.EDAnalyzer("METBlock",
    metSrc = cms.PSet(
        METs = cms.InputTag('patMETs'),
        METsPF = cms.InputTag('patMETsPF'),
        METsTC = cms.InputTag('patMETsTC'),
        METsMVAmuTau = cms.InputTag('patMETsMVAmuTau'),
        METsMVAeTau = cms.InputTag('patMETsMVAeTau'),
        METsMVAtauTau = cms.InputTag('patMETsMVAtauTau')
    )
)
