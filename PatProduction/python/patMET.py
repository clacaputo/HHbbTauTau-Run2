## patMET.py - configuration file that defines parameters related to PAT MET objects.

import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.metTools import *

def applyMETParameters(process, isMC):
    # MET corrections.
    addPfMET(process, 'PF')
    process.patMETsPF.metSource = cms.InputTag("pfMet")

    process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
    process.pfJetMETcorr.offsetCorrLabel = cms.string("ak5PFL1Fastjet")
    if isMC:
        process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
    else:
        process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

    process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
    process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
    process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfMETcorrType0'),
        cms.InputTag('pfJetMETcorr', 'type1')
    )

    process.patPFMETsTypeIcorrected = process.patMETs.clone(
        metSource = cms.InputTag('pfType1CorrectedMet'),
        addMuonCorrections = cms.bool(False),
        genMETSource = cms.InputTag('genMetTrue'),
        addGenMET = cms.bool(False)
    )

    return

def applyMVAMETParamteres(process, isMC):
    #process.load('RecoMET.METPUSubtraction.mvaPFMET_leptons_cff')
    process.load('HHbbTauTau.PatProduction.mvaPFMET_leptons_cff')
    if isMC:
        process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3")
    else:
        process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual")

    process.patPFMetMVA = process.patMETs.clone(
        metSource = cms.InputTag('pfMEtMVA'),
        addMuonCorrections = cms.bool(False),
        addGenMET = cms.bool(False),
        genMETSource = cms.InputTag('genMetTrue')
    )

    return
