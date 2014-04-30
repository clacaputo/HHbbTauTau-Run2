import FWCore.ParameterSet.Config as cms

jetBlock = cms.EDAnalyzer("JetBlock",
    jetSrc = cms.InputTag('patJets'),
)
