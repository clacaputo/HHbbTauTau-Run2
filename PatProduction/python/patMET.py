## patMET.py - configuration file that defines parameters related to PAT MET objects.

from PhysicsTools.PatAlgos.tools.metTools import *
#from RecoMET.METPUSubtraction.mvaPFMET_leptons_cfi import *
from RecoJets.JetProducers.PileupJetIDParams_cfi import JetIdParams

def applyMETParameters(cms, process, isMC):
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

from RecoTauTag.RecoTau.PFRecoTauDiscriminationByHPSSelection_cfi import hpsSelectionDiscriminator

#MVA MET recommended
'''
def applyMVAMETParamteres(cms, process, isMC):
	process.calibratedAK5PFJetsForPFMEtMVA = cms.EDProducer('PFJetCorrectionProducer',
		src = cms.InputTag('ak5PFJets'),
	    correctors = cms.vstring("ak5PFL1FastL2L3") # NOTE: use "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
	)
	
	process.pfMEtMVA = cms.EDProducer("PFMETProducerMVA",
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
			U     = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_53_June2013_type1.root'),
			DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'),
			CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_53_Dec2012.root'),
			CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_53_Dec2012.root')
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
	process.isomuons = cms.EDFilter(
		"MuonSelector",
		src = cms.InputTag('muons'),
		cut = cms.string(    "(isTrackerMuon) && abs(eta) < 2.5 && pt > 9.5"+#17. "+
			                 "&& isPFMuon"+
				             "&& globalTrack.isNonnull"+
					         "&& innerTrack.hitPattern.numberOfValidPixelHits > 0"+
						     "&& innerTrack.normalizedChi2 < 10"+
							 "&& numberOfMatches > 0"+
							 "&& innerTrack.hitPattern.numberOfValidTrackerHits>5"+
	                         "&& globalTrack.hitPattern.numberOfValidHits>0"+
		                     "&& (pfIsolationR03.sumChargedHadronPt+pfIsolationR03.sumNeutralHadronEt+pfIsolationR03.sumPhotonEt)/pt < 0.3"+
			                 "&& abs(innerTrack().dxy)<2.0"
				             ),
	    filter = cms.bool(False)
    )

	process.isoelectrons = cms.EDFilter(
	    "GsfElectronSelector",
		src = cms.InputTag('gsfElectrons'),
		cut = cms.string(
			"abs(eta) < 2.5 && pt > 9.5"                               +
	        "&& gsfTrack.trackerExpectedHitsInner.numberOfHits == 0"   +
#	         "&& (pfIsolationVariables.chargedHadronIso+pfIsolationVariables.neutralHadronIso)/et     < 0.3"  +
			"&& (isolationVariables03.tkSumPt)/et              < 0.2"  +
			"&& ((abs(eta) < 1.4442  "                                 +
	        "&& abs(deltaEtaSuperClusterTrackAtVtx)            < 0.007"+
		    "&& abs(deltaPhiSuperClusterTrackAtVtx)            < 0.8"  +
			"&& sigmaIetaIeta                                  < 0.01" +
			"&& hcalOverEcal                                   < 0.15" +
	        "&& abs(1./superCluster.energy - 1./p)             < 0.05)"+
		    "|| (abs(eta)  > 1.566 "+
			"&& abs(deltaEtaSuperClusterTrackAtVtx)            < 0.009"+
	        "&& abs(deltaPhiSuperClusterTrackAtVtx)            < 0.10" +
		    "&& sigmaIetaIeta                                  < 0.03" +
			"&& hcalOverEcal                                   < 0.10" +
			"&& abs(1./superCluster.energy - 1./p)             < 0.05))" 
	    ),
		filter = cms.bool(False)
    )

	process.hpsPFTauDiscriminationByDecayModeFinding = hpsSelectionDiscriminator.clone(
        PFTauProducer = cms.InputTag('hpsPFTauProducer')
    )

	process.requireDecayMode = cms.PSet(
        BooleanOperator = cms.string("and"),
        decayMode = cms.PSet(
            Producer = cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding'),
            cut = cms.double(0.5)
        )
    )

	process.hpsPFTauDiscriminationAgainstMuon2 = cms.EDProducer("PFRecoTauDiscriminationAgainstMuon2",
		PFTauProducer = cms.InputTag('hpsPFTauProducer'),
        Prediscriminants = process.requireDecayMode.clone(),
		discriminatorOption = cms.string('loose'), # available options are: 'loose', 'medium', 'tight'
		HoPMin = cms.double(0.2)
    )

	process.isotaus = cms.EDFilter(
		"PFTauSelector",
		src = cms.InputTag('hpsPFTauProducer'),
		BooleanOperator = cms.string("and"),
		discriminators = cms.VPSet(
			cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),       selectionCut=cms.double(0.5)),
			#cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByMVAIsolation"),           selectionCut=cms.double(0.5)),
			cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"),           selectionCut=cms.double(0.5)),
			cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"), selectionCut=cms.double(0.5)),
			cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationAgainstMuon2"),             selectionCut=cms.double(0.5)) 
		),
		cut = cms.string("abs(eta) < 2.3 && pt > 19.0 "),
		filter = cms.bool(False)
    )

	process.isomuonseq = cms.Sequence(process.isomuons)
	process.isoelectronseq = cms.Sequence(process.isoelectrons)
	process.isotauseq      = cms.Sequence(
#		hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits*
		#kt6PFJetsForRhoComputationVoronoiMet*
		#hpsPFTauDiscriminationByMVAIsolation*
		process.hpsPFTauDiscriminationAgainstMuon2*
		process.isotaus
    )

	process.leptonSelection = cms.PSet(
		SelectEvents = cms.untracked.PSet(
			SelectEvents = cms.vstring(
			    'isomuonseq',
			    'isoelectronseq',
			    'isotauseq'
			)
		)
    )

#	process.load('RecoMET.METPUSubtraction.mvaPFMET_leptons_cfi')
	process.pfMEtMVAsequence  = cms.Sequence(
#	    process.isoelectronseq*
	    (process.isomuonseq+process.isotauseq+process.isoelectronseq)*
		process.calibratedAK5PFJetsForPFMEtMVA*
	    process.pfMEtMVA
    )
	return
'''

def applyMVAMETParamteres(cms, process, isMC):
#	process.load('RecoMET.METPUSubtraction.mvaPFMET_leptons_cff')
	process.load('HHbbTauTau.PatProduction.mvaPFMET_leptons_cff')
#    if isMC:
#            process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3")
#    else:
#            process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual")

#    process.patPFMetMVA = process.patMETs.clone(
#            metSource = cms.InputTag('pfMEtMVA'),
#            addMuonCorrections = cms.bool(False),
#            addGenMET = cms.bool(False),
#            genMETSource = cms.InputTag('genMetTrue')
#    )
	return
