## mvaPFMET_leptons_cff.py - adapted copy of mvaPFMET_leptons_cff.py from RecoMET/METPUSubtraction package.
## Original file: https://github.com/cms-met/cmssw/blob/53X-MVaNoPuMET-20131217-01/RecoMET/METPUSubtraction/python/mvaPFMET_leptons_cff.py

import FWCore.ParameterSet.Config as cms

from HHbbTauTau.PatProduction.mvaPFMET_leptons_cfi import *

calibratedAK5PFJetsForPFMEtMVA = cms.EDProducer('PFJetCorrectionProducer',
    src = cms.InputTag('ak5PFJets'),
    correctors = cms.vstring("ak5PFL1FastL2L3") # NOTE: use "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
)

pfMEtMVAmuTau = cms.EDProducer("PFMETProducerMVA",
    srcCorrJets = cms.InputTag('calibratedAK5PFJetsForPFMEtMVA'),
    srcUncorrJets = cms.InputTag('ak5PFJets'),
    srcPFCandidates = cms.InputTag('particleFlow'),
    srcVertices = cms.InputTag('offlinePrimaryVertices'),
    srcLeptons = cms.VInputTag("isomuons","isotausMT"),
    minNumLeptons = cms.int32(0),
    srcRho = cms.InputTag('kt6PFJets','rho'),
    globalThreshold = cms.double(-1.),#pfMet.globalThreshold,
    minCorrJetPt = cms.double(-1.),
    corrector = cms.string("ak5PFL1Fastjet"),
    useType1  = cms.bool(True),
    useOld42  = cms.bool(False),
    verbosity = cms.int32(0),
    # For PFMETAlgorithmMVA
    loadMVAfromDB = cms.bool(False),
    inputFileNames = cms.PSet(
        DPhi = cms.FileInPath('HHbbTauTau/PatProduction/data/gbrmetphi_53_Dec2012.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_53_Dec2012.root'),
        U = cms.FileInPath('HHbbTauTau/PatProduction/data/gbrmet_53_Dec2012.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_53_Dec2012.root')
    ),
    is42 = cms.bool(False), # CV: set this flag to true if you are running mvaPFMET in CMSSW_4_2_x
    dZcut     = cms.double(0.1),
    # For PileupJetIdAlgo
    tmvaVariables = cms.vstring('nvtx',
        'dZ',
        'beta',
        'betaStar',
        'nCharged',
        'nNeutrals',
        'dR2Mean',
        'ptD',
        'frac01',
        'frac02',
        'frac03',
        'frac04',
        'frac05'),
    tmvaMethod = cms.string('JetIDMVAHighPt'),
    cutBased = cms.bool(False),
    tmvaWeights = cms.string('RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_53X_Dec2012.weights.xml'),
    tmvaSpectators = cms.vstring('jetPt',
        'jetEta',
        'jetPhi'),
    version = cms.int32(-1),
    JetIdParams = cms.PSet(
        Pt2030_Tight = cms.vdouble(0.73, 0.05, -0.26, -0.42),
        Pt2030_Loose = cms.vdouble(-0.63, -0.6, -0.55, -0.45),
        Pt3050_Medium = cms.vdouble(0.1, -0.36, -0.54, -0.54),
        Pt1020_MET = cms.vdouble(0.3, -0.2, -0.4, -0.4),
        Pt2030_Medium = cms.vdouble(0.1, -0.36, -0.54, -0.54),
        Pt010_Tight = cms.vdouble(-0.83, -0.81, -0.74, -0.81),
        Pt1020_Tight = cms.vdouble(-0.83, -0.81, -0.74, -0.81),
        Pt3050_MET = cms.vdouble(0.0, 0.0, -0.1, -0.2),
        Pt010_MET = cms.vdouble(0.0, -0.6, -0.4, -0.4),
        Pt1020_Loose = cms.vdouble(-0.95, -0.96, -0.94, -0.95),
        Pt010_Medium = cms.vdouble(-0.83, -0.92, -0.9, -0.92),
        Pt1020_Medium = cms.vdouble(-0.83, -0.92, -0.9, -0.92),
        Pt2030_MET = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        Pt010_Loose = cms.vdouble(-0.95, -0.96, -0.94, -0.95),
        Pt3050_Loose = cms.vdouble(-0.63, -0.6, -0.55, -0.45),
        Pt3050_Tight = cms.vdouble(0.73, 0.05, -0.26, -0.42)
    ),
    impactParTkThreshold = cms.double(1.0)
)

pfMEtMVAeTau = pfMEtMVAmuTau.clone( srcLeptons = cms.VInputTag("patElectronsWithEmbeddedVariables","isotausET") )
pfMEtMVAtauTau = pfMEtMVAmuTau.clone( srcLeptons = cms.VInputTag("isotausTT") )

pfMEtMVAsequence  = cms.Sequence(
    (isomuonseq + isotauseq + isoelectronseq) *
    calibratedAK5PFJetsForPFMEtMVA *
    pfMEtMVAmuTau *
    pfMEtMVAeTau *
    pfMEtMVAtauTau
)

