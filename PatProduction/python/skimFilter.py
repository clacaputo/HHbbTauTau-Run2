## patTaus.py - configuration file that defines parameters related to PAT Tau objects.

import FWCore.ParameterSet.Config as cms

def applySkim(process):
    process.bbttSkim   = cms.EDFilter("SkimFilter",
                                      vertexSrc  = cms.untracked.InputTag('patVertices'),
                                      muonSrc  = cms.untracked.InputTag('patMuonsWithEmbeddedVariables'),
                                      electronSrc=cms.untracked.InputTag("patElectronsWithEmbeddedVariables"),
                                      tauSrc  = cms.untracked.InputTag("patTausTriggerMatch")
                                      )
    return
