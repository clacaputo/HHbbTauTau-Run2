## patTaus.py - configuration file that defines parameters related to PAT Tau objects.

import FWCore.ParameterSet.Config as cms

def applyVertexParameters(process):
    process.patVertices = cms.EDProducer('PatVertexProducer',
        inputTag = cms.InputTag('offlinePrimaryVerticesWithBS')
    )

    return
