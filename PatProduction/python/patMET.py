## patMET.py - configuration file that defines parameters related to PAT MET objects.

import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.metTools import *

def applyMETParameters(process, isMC):
    process.patMETs.addGenMET = cms.bool(False)
    addTcMET(process, 'TC')
    addPfMET(process, 'PF')
    process.patMETsPF.metSource = cms.InputTag("pfMet")

    process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
    process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")

    process.pfJetMETcorr.offsetCorrLabel = cms.string("ak5PFL1Fastjet")
    if isMC:
        process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
    else:
        process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

    process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
    process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfMETcorrType0'),
        cms.InputTag('pfJetMETcorr', 'type1')
    )

    return

def applyMVAMETParamteres(process, isMC):
    #process.load('RecoMET.METPUSubtraction.mvaPFMET_leptons_cff')
    process.load('HHbbTauTau.PatProduction.mvaPFMET_leptons_cff')
    if isMC:
        process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3")
    else:
        process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual")

    process.patMETsMVAmuTau = process.patMETs.clone( metSource = cms.InputTag('pfMEtMVAmuTau') )
    process.patMETsMVAeTau = process.patMETs.clone( metSource = cms.InputTag('pfMEtMVAeTau') )
    process.patMETsMVAtauTau = process.patMETs.clone( metSource = cms.InputTag('pfMEtMVAtauTau') )

    process.patMETsMVAsequence = cms.Sequence(
        process.patMETsMVAmuTau *
        process.patMETsMVAeTau *
        process.patMETsMVAtauTau
    )

    return
