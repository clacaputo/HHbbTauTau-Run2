## patTaus.py - configuration file that defines parameters related to PAT Tau objects.

def applyTauParameters(proccess):
    # Include high Pt tau reconstruction sequence.
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    from PhysicsTools.PatAlgos.tools.tauTools import *
    switchToPFTauHPS(process)

    process.patTaus.embedLeadPFCand = cms.bool(True)
    process.patTaus.embedSignalPFCands = cms.bool(True)
    process.patTaus.embedIsolationPFCands = cms.bool(True)
    #process.patTaus.embedLeadTrack = cms.bool(True)
    #process.patTaus.embedSignalTracks = cms.bool(True)
    #process.patTaus.embedIsolationTracks = cms.bool(True)
    process.patTaus.embedIsolationPFChargedHadrCands = cms.bool(True)
    process.patTaus.embedIsolationPFNeutralHadrCands = cms.bool(True)
    process.patTaus.embedIsolationPFGammaCands = cms.bool(True)
    process.patTaus.embedSignalPFChargedHadrCands = cms.bool(True)
    process.patTaus.embedSignalPFNeutralHadrCands = cms.bool(True)
    process.patTaus.embedSignalPFGammaCands = cms.bool(True)
    process.patTaus.embedLeadPFChargedHadrCand = cms.bool(True)
    process.patTaus.embedLeadPFNeutralCand = cms.bool(True)

    return
