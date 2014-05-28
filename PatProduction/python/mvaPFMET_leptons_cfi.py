## mvaPFMET_leptons_cfi.py - adapted copy of mvaPFMET_leptons_cfi.py from RecoMET/METPUSubtraction package.
## Original file:https://github.com/cms-met/cmssw/blob/53X-MVaNoPuMET-20131217-01/RecoMET/METPUSubtraction/python/mvaPFMET_leptons_cfi.py

import FWCore.ParameterSet.Config as cms

# Single muon for Wjets
isomuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag('patMuonsWithEmbeddedVariables'),
    cut = cms.string(
        #"abs(eta) < 2.1 && pt > 17."
        "abs(eta) < 2.6 && pt > 7."
        ),
    filter = cms.bool(False)
)

isoelectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag('patElectronsWithEmbeddedVariables'),
    cut = cms.string(
        'abs(eta) < 2.1 && pt > 20.0 ' +
        '&& gsfTrack.trackerExpectedHitsInner.numberOfHits == 0 ' +
        '&& (abs(superCluster.eta)<0.8 && userFloat("mvaPOGNonTrig") > 0.925) ' +
        '|| (abs(superCluster.eta)>0.8 && abs(superCluster.eta)<1.479 && userFloat("mvaPOGNonTrig") > 0.975) ' +
        '|| (abs(superCluster.eta)>1.479 && userFloat("mvaPOGNonTrig") > 0.985) '
    ),
    filter = cms.bool(False)
)

isotausTT = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag('patTaus'),
    cut = cms.string(
        'abs(eta) < 2.1 && pt > 40.0 ' +
        '&& tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 10.0 ' +
        '&& tauID("decayModeFinding") > 0.5 '
    ),
    filter = cms.bool(False)
)

isotausET = isotausTT.clone(cut = cms.string(
                                'abs(eta) < 2.3 && pt > 15.0 ' +
                                '&& tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 10.0 ' +
                                '&& tauID("decayModeFinding") > 0.5 '
                                ))
isotausMT = isotausTT.clone(cut = cms.string(
                                #'abs(eta) < 2.3 && pt > 15.0 ' +
                                'abs(eta) < 2.6 && pt > 18.0 ' +
                                #'&& tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 10.0 ' +
                                '&& tauID("decayModeFinding") > 0.5 ' +
                                #'&& tauID("againstMuonTight") > 0.5 '
                                ))


isomuonseq     = cms.Sequence(isomuons)
isoelectronseq = cms.Sequence(isoelectrons)
isotauseq      = cms.Sequence(
    isotausTT *
    isotausET *
    isotausMT
)
