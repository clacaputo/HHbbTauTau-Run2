import FWCore.ParameterSet.Config as cms

genParticleBlock = cms.EDAnalyzer("GenParticleBlock",
  verbosity      = cms.int32(0),
  genParticleSrc = cms.InputTag('genParticles')
)
