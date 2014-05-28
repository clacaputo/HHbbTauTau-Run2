## mvaPFMET_leptons_cff.py - adapted copy of mvaPFMET_leptons_cff.py from RecoMET/METPUSubtraction package.
## Original file: https://github.com/cms-met/cmssw/blob/53X-MVaNoPuMET-20131217-01/RecoMET/METPUSubtraction/python/mvaPFMET_leptons_cff.py

import FWCore.ParameterSet.Config as cms

from HHbbTauTau.PatProduction.mvaPFMET_leptons_cfi import *
from RecoJets.JetProducers.PileupJetIDParams_cfi import JetIdParams

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
        'jetPt',
        'jetEta',
        'jetPhi',
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
    tmvaMethod = cms.string('JetIDMVAMET'),
    cutBased = cms.bool(False),
    tmvaWeights = cms.string('RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml'),
    tmvaSpectators = cms.vstring(),
    label = cms.string('met53x'),
    version = cms.int32(-1),
    JetIdParams = cms.PSet(
        Pt2030_Tight = cms.vdouble(-2, -2, -2, -2, -2),
        Pt2030_Loose = cms.vdouble(-2, -2, -2, -2, -2),
        Pt3050_Medium = cms.vdouble(-2, -2, -2, -2, -2),
        Pt1020_MET = cms.vdouble(-0.2, -0.2, -0.5, -0.3),
        Pt2030_Medium = cms.vdouble(-2, -2, -2, -2, -2),
        Pt010_Tight = cms.vdouble(-2, -2, -2, -2, -2),
        Pt1020_Tight = cms.vdouble(-2, -2, -2, -2, -2),
        Pt3050_MET = cms.vdouble(-0.2, -0.2, 0.0, 0.2),
        Pt010_MET = cms.vdouble(-0.2, -0.3, -0.5, -0.5),
        Pt1020_Loose = cms.vdouble(-2, -2, -2, -2, -2),
        Pt010_Medium = cms.vdouble(-2, -2, -2, -2, -2),
        Pt1020_Medium = cms.vdouble(-2, -2, -2, -2, -2),
        Pt2030_MET = cms.vdouble(-0.2, -0.2, -0.2, 0.1),
        Pt010_Loose = cms.vdouble(-2, -2, -2, -2, -2),
        Pt3050_Loose = cms.vdouble(-2, -2, -2, -2, -2),
        Pt3050_Tight = cms.vdouble(-2, -2, -2, -2, -2)
    ),
    impactParTkThreshold = cms.double(1.0)

)

pfMEtMVA = cms.EDProducer("PFMETProducerMVA",
    srcCorrJets = cms.InputTag('calibratedAK5PFJetsForPFMEtMVA'),
    srcUncorrJets = cms.InputTag('ak5PFJets'),
    srcPFCandidates = cms.InputTag('particleFlow'),
    srcVertices = cms.InputTag('offlinePrimaryVertices'),
    srcLeptons = cms.VInputTag(),#"isomuons","isoelectrons","isotaus") # NOTE: you need to set this to collections of electrons, muons and tau-jets
                                 #                                             passing the lepton reconstruction & identification criteria applied in your analysis
    minNumLeptons = cms.int32(0),
    srcRho = cms.InputTag('kt6PFJets','rho'),
    globalThreshold = cms.double(-1.),#pfMet.globalThreshold,
    minCorrJetPt = cms.double(-1.),
    inputFileNames = cms.PSet(
        DPhi = cms.FileInPath('HHbbTauTau/PatProduction/data/gbrmetphi_53_Dec2012.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_53_Dec2012.root'),
        U = cms.FileInPath('HHbbTauTau/PatProduction/data/gbrmet_53_Dec2012.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_53_Dec2012.root')
    ),
    loadMVAfromDB = cms.bool(False),
    is42 = cms.bool(False), # CV: set this flag to true if you are running mvaPFMET in CMSSW_4_2_x
    corrector = cms.string("ak5PFL1Fastjet"),
    useType1  = cms.bool(True),
    useOld42  = cms.bool(False),
    dZcut     = cms.double(0.1),
    impactParTkThreshold = cms.double(0.),
    tmvaWeights = cms.string("RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml"),
    tmvaMethod = cms.string("JetID"),
    version = cms.int32(-1),
    cutBased = cms.bool(False),
    tmvaVariables = cms.vstring(
        "nvtx",
        "jetPt",
        "jetEta",
        "jetPhi",
        "dZ",
        "beta",
        "betaStar",
        "nCharged",
        "nNeutrals",
        "dR2Mean",
        "ptD",
        "frac01",
        "frac02",
        "frac03",
        "frac04",
        "frac05",
    ),
    tmvaSpectators = cms.vstring(),
    JetIdParams = JetIdParams,
    verbosity = cms.int32(0)
)

pfMEtMVAeTau = pfMEtMVAmuTau.clone( srcLeptons = cms.VInputTag("isoelectrons","isotausET") )
pfMEtMVAtauTau = pfMEtMVA.clone( srcLeptons = cms.VInputTag("isotausTT") )

pfMEtMVAsequence  = cms.Sequence(
    (isomuonseq + isotauseq + isoelectronseq) *
    calibratedAK5PFJetsForPFMEtMVA *
    pfMEtMVAmuTau *
    pfMEtMVAeTau *
    pfMEtMVAtauTau
)

