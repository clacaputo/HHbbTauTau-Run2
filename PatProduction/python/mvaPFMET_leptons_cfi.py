## mvaPFMET_leptons_cfi.py - adapted copy of mvaPFMET_leptons_cfi.py from RecoMET/METPUSubtraction package.
## Original file:https://github.com/cms-met/cmssw/blob/53X-MVaNoPuMET-20131217-01/RecoMET/METPUSubtraction/python/mvaPFMET_leptons_cfi.py

import FWCore.ParameterSet.Config as cms

# Single muon for Wjets
isomuons = cms.EDFilter(
    "MuonSelector",
    src = cms.InputTag('muons'),
    cut = cms.string(
        "abs(eta) < 2.1 && pt > 17."
        ),
    filter = cms.bool(False)
)

isoelectrons = cms.EDFilter(
    "GsfElectronSelector",
    src = cms.InputTag('gsfElectrons'),
    cut = cms.string(
        'abs(eta) < 2.1 && pt > 20.0 ' +
        '&& gsfTrack.trackerExpectedHitsInner.numberOfHits == 0 '# +
#        '&& (abs(superCluster.eta)<0.8 && mvaNonTrigV0 > 0.925) || (abs(leg2().sourcePtr().superCluster().eta())>0.8 && abs(leg2().sourcePtr().superCluster().eta())<1.479 && leg2().mvaNonTrigV0() > 0.975) || (abs(leg2().sourcePtr().superCluster().eta())>1.479 && leg2().mvaNonTrigV0() > 0.985)'
    ),
    filter = cms.bool(False)
)

from RecoTauTag.Configuration.HPSPFTaus_cff import hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits

hpsPFTauDiscriminationByCombinedIsolationDBSumPtCorr3HitsMVAmet = hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits.clone(
    maximumSumPtCut = 10.0
)

isotausTT = cms.EDFilter(
    "PFTauSelector",
    src = cms.InputTag('hpsPFTauProducer'),
    BooleanOperator = cms.string("and"),
    discriminators = cms.VPSet(
        cms.PSet(
            discriminator = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
            selectionCut = cms.double(0.5)
        ),
        #cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByMVAIsolation"), selectionCut=cms.double(0.5)),
        cms.PSet(
            discriminator = cms.InputTag("hpsPFTauDiscriminationByCombinedIsolationDBSumPtCorr3HitsMVAmet"),
            selectionCut = cms.double(0.5)
        )
#        cms.PSet(
#            discriminator = cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"),
#            selectionCut = cms.double(0.5)
#        ),
#        cms.PSet(
#            discriminator = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection2"),
#            selectionCut=cms.double(0.5)
#        )
    ),
    cut = cms.string('abs(eta) < 2.1 && pt > 40.0'),
    filter = cms.bool(False)
)

isotausET = isotausTT.clone(cut = cms.string("abs(eta) < 2.3 && pt > 15.0 "))
isotausMT = isotausTT.clone(cut = cms.string("abs(eta) < 2.3 && pt > 15.0 "),
                            discriminators = cms.VPSet(
                                cms.PSet(
                                    discriminator = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
                                    selectionCut = cms.double(0.5)
                                ),
		                        cms.PSet(
                                    discriminator = cms.InputTag("hpsPFTauDiscriminationByCombinedIsolationDBSumPtCorr3HitsMVAmet"),
                                    selectionCut = cms.double(0.5)
                                ),
                                cms.PSet(
                                    discriminator = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection"),
                                    selectionCut=cms.double(0.5)
                                )
                            )
                           )

isomuonseq     = cms.Sequence(isomuons)
isoelectronseq = cms.Sequence(isoelectrons)
isotauseq      = cms.Sequence(
    hpsPFTauDiscriminationByCombinedIsolationDBSumPtCorr3HitsMVAmet *
    isotausTT *
    isotausET *
    isotausMT
)
