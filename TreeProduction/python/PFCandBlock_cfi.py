import FWCore.ParameterSet.Config as cms

pfCandBlock = cms.EDAnalyzer("PFCandBlock",
    srcPFCandidates = cms.InputTag('particleFlow'),
)
