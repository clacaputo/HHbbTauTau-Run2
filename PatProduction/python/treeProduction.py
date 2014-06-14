## treeProduction.py - configuration file that defines tree sequence.

import FWCore.ParameterSet.Config as cms

def addTreeSequence(process, includeSim, treeOutput):
    process.TFileService = cms.Service("TFileService", fileName = cms.string(treeOutput) )
    process.load("HHbbTauTau.TreeProduction.TreeContentConfig_cff")

    process.mainTreeContentSequence = cms.Sequence(
        process.eventBlock
      + process.vertexBlock
      + process.electronBlock
      + process.jetBlock
      + process.metBlock
      + process.muonBlock
      + process.tauBlock
      + process.pfCandBlock
      + process.triggerBlock
      + process.triggerObjectBlock
    )

    process.simTreeContentSequence = cms.Sequence()
    if includeSim:
        process.simTreeContentSequence = cms.Sequence(process.genParticleBlock + process.genMETBlock)

    return
