import FWCore.ParameterSet.Config as cms

genJetBlock = cms.EDAnalyzer("GenJetBlock",
  verbosity = cms.int32(0),
  genJetSrc = cms.InputTag('ak5GenJets')
)
