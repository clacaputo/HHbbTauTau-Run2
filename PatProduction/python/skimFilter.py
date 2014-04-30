## patTaus.py - configuration file that defines parameters related to PAT Tau objects.

import FWCore.ParameterSet.Config as cms

def applySkim(process):
    process.bbttSkim   = cms.EDFilter("SkimFilter",
                                      vertexTag  = cms.untracked.InputTag('patVertices'),
                                      muonTag  = cms.untracked.InputTag('patMuonsWithEmbeddedVariables'),
                                      electronTag=cms.untracked.InputTag("patElectronsWithEmbeddedVariables"),
                                      tauTag  = cms.untracked.InputTag("patTausTriggerMatch"),
                                      filter = cms.bool(True)
                                      )
    return
