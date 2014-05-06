## patElectrons.py - configuration file that defines parameters related to PAT Jet objects.

import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.jetTools import *

def applyJetParameters(process, isMC):
    # Jet corrections
    jec = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
    if not isMC:
        jec.extend([ 'L2L3Residual' ])
    switchJetCollection(process, cms.InputTag('ak5PFJets'),
         doJTA        = True,
         doBTagging   = True,
         jetCorrLabel = ('AK5PF', cms.vstring(jec)),
         doType1MET   = True,
         genJetCollection=cms.InputTag("ak5GenJets"),
         doJetID      = True,
         jetIdLabel   = 'ak5'
    )

    process.load("RecoJets.JetProducers.PileupJetID_cfi")
    process.pileupJetIdProducer.jets = cms.InputTag('ak5PFJets')
    process.pileupJetIdProducerChs.jets = cms.InputTag('ak5PFJets')
    process.pileupJetIdProducer.applyJec = cms.bool(True)
    process.pileupJetIdProducer.inputIsCorrected = cms.bool(False)

    # Embed into PAT jets as userdata
    process.patJets.userData.userFloats.src = cms.VInputTag(
        cms.InputTag('pileupJetIdProducer', 'fullDiscriminant'),
        cms.InputTag('pileupJetIdProducer', 'cutbasedDiscriminant'),
        cms.InputTag('pileupJetIdProducer', 'philv1Discriminant'))
    process.patJets.userData.userInts.src = cms.VInputTag(
        cms.InputTag('pileupJetIdProducer', 'fullId'),
        cms.InputTag('pileupJetIdProducer', 'cutbasedId'),
        cms.InputTag('pileupJetIdProducer', 'philv1Id'))

    return
