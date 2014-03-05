import FWCore.ParameterSet.Config as cms

jetBlock = cms.EDAnalyzer("JetBlock",
  verbosity        = cms.int32(0),
  jetSrc           = cms.InputTag('patJets'),
  jecUncertainty   = cms.string('CondFormats/JetMETObjects/data/Spring10_Uncertainty_AK5PF.txt'),
  applyResidualJEC = cms.bool(False),
  residualJEC      = cms.string('CondFormats/JetMETObjects/data/Spring10_L2L3Residual_AK5PF.txt')
)
