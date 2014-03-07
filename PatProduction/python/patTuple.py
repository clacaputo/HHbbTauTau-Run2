## Import skeleton process.
from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.options.wantSummary = False

process.source.dropDescendantsOfDroppedBranches = cms.untracked.bool(False)
process.source.inputCommands = cms.untracked.vstring(
        'keep *',
        'drop recoPFTaus_*_*_*'
    )

## Parse and apply options.
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('globalTag',
                  'START53_V7A::All',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Global Tag to use. Default: START53_V7A::All")
options.register ('isMC',
                  True,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Sample Type: MC or data")
options.register ('runOnCrab',
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Indicates if script will be executed on CRAB.")
options.register ('fileList',
                  'fileList.py',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "List of root files to process.")                  
options.register ('includeSim',
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Include Sim. Default: False")

options.parseArguments()

process.GlobalTag.globaltag = options.globalTag

if not options.runOnCrab:
    fileList = cms.untracked.vstring()
    process.source.fileNames = fileList
    execfile(options.fileList)
    process.out.fileName = options.outputFile
    process.maxEvents.input = options.maxEvents


## Muons
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFMuonIso

process.muIsoSequence       = setupPFMuonIso(process,'muons')

from CommonTools.ParticleFlow.pfParticleSelection_cff import pfParticleSelectionSequence
process.pfParticleSelectionSequence = pfParticleSelectionSequence


process.patMuons.isoDeposits = cms.PSet(
    pfAllParticles   = cms.InputTag("muPFIsoDepositPUPFIso"),      # all PU   CH+MU+E
    pfChargedHadrons = cms.InputTag("muPFIsoDepositChargedPFIso"), # all noPU CH

    pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutralPFIso"), # all NH
    pfPhotons        = cms.InputTag("muPFIsoDepositGammaPFIso"),   # all PH

    user = cms.VInputTag(
    cms.InputTag("muPFIsoDepositChargedAllPFIso"),                 # all noPU CH+MU+E
    )
    )

process.patMuons.isolationValues = cms.PSet(
    pfAllParticles   = cms.InputTag("muPFIsoValuePU04PFIso"),
    pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04PFIso"),

    pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04PFIso"),
    pfPhotons        = cms.InputTag("muPFIsoValueGamma04PFIso"),

    user = cms.VInputTag(
    cms.InputTag("muPFIsoValueChargedAll04PFIso"),
    )
    )

process.patMuons.embedTrack = cms.bool(True)

process.patMuonsWithEmbeddedVariables = cms.EDProducer('MuonsUserEmbedded',

        muonTag = cms.InputTag("patMuons"),
        vertexTag = cms.InputTag("offlinePrimaryVerticesWithBS"),

)

## Taus
## Include high Pt tau reconstruction sequence.
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

## Electrons
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso

process.electronIsoSequence = setupPFElectronIso(process,'gsfElectrons')

process.patElectrons.isoDeposits = cms.PSet(
    pfAllParticles   = cms.InputTag("elPFIsoDepositPUPFIso"),      # all PU   CH+MU+E

    pfChargedHadrons = cms.InputTag("elPFIsoDepositChargedPFIso"), # all noPU CH
    pfNeutralHadrons = cms.InputTag("elPFIsoDepositNeutralPFIso"), # all NH

    pfPhotons        = cms.InputTag("elPFIsoDepositGammaPFIso"),   # all PH
    user = cms.VInputTag(
    cms.InputTag("elPFIsoDepositChargedAllPFIso"),                 # all noPU CH+MU+E

    )
    )
process.patElectrons.isolationValues = cms.PSet(
    pfAllParticles   = cms.InputTag("elPFIsoValuePU04PFIdPFIso"),
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged04PFIdPFIso"),

    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral04PFIdPFIso"),
    pfPhotons        = cms.InputTag("elPFIsoValueGamma04PFIdPFIso"),

    user = cms.VInputTag(
    cms.InputTag("elPFIsoValueChargedAll04PFIdPFIso"),
    cms.InputTag("elPFIsoValueChargedAll04NoPFIdPFIso"),

    cms.InputTag("elPFIsoValuePU04NoPFIdPFIso"),
    cms.InputTag("elPFIsoValueCharged04NoPFIdPFIso"),
    cms.InputTag("elPFIsoValueGamma04NoPFIdPFIso"),

    cms.InputTag("elPFIsoValueNeutral04NoPFIdPFIso")
    )

    )

process.patElectronsWithEmbeddedVariables = cms.EDProducer('ElectronsUserEmbedder',
        electronTag = cms.InputTag("patElectrons"),
        vertexTag = cms.InputTag("offlinePrimaryVerticesWithBS"),
        isMC = cms.bool(options.isMC),
        doMVAPOG = cms.bool(True),

        inputFileName0v2 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_TrigV0_Cat1.weights.xml'),
        inputFileName1v2 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_TrigV0_Cat2.weights.xml'),
        inputFileName2v2 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_TrigV0_Cat3.weights.xml'),
        inputFileName3v2 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_TrigV0_Cat4.weights.xml'),
        inputFileName4v2 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_TrigV0_Cat5.weights.xml'),
        inputFileName5v2 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_TrigV0_Cat6.weights.xml'),

        inputFileName0v3 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat1.weights.xml'),
        inputFileName1v3 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat2.weights.xml'),
        inputFileName2v3 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat3.weights.xml'),
        inputFileName3v3 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat4.weights.xml'),
        inputFileName4v3 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat5.weights.xml'),
        inputFileName5v3 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat6.weights.xml'),

)

process.patElectrons.embedTrack = cms.bool(True)
process.patElectrons.embedGsfTrack = cms.bool(True)

process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.patElectrons.electronIDSources = cms.PSet(mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"))


## MET corrections.
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')
process.patMETsPF.metSource = cms.InputTag("pfMet")

process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.pfJetMETcorr.offsetCorrLabel = cms.string("ak5PFL1Fastjet")
if options.isMC:
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

##Jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *

jec = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
if not options.isMC:
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

## remove MC matching from the default sequence
if not options.isMC:
        removeMCMatching(process, ['All'])
        #removeMCMatching(process, ['METs'], "TC")
        removeMCMatching(process, ['METs'], "PF")
        process.patDefaultSequence.remove(process.patJetPartonMatch)
        #process.patDefaultSequence.remove(process.patJetPartonMatchAK5PF)
        #process.patDefaultSequence.remove(process.patJetGenJetMatchAK5PF)
        process.patDefaultSequence.remove(process.patJetFlavourId)
        process.patDefaultSequence.remove(process.patJetPartons)
        process.patDefaultSequence.remove(process.patJetPartonAssociation)
        #process.patDefaultSequence.remove(process.patJetPartonAssociationAK5PF)
        process.patDefaultSequence.remove(process.patJetFlavourAssociation)
        #process.patDefaultSequence.remove(process.patJetFlavourAssociationAK5PF)
        runOnData(process)


## Define process path.
process.p = cms.Path(
    process.PFTau
  * process.pfParticleSelectionSequence
  * process.muIsoSequence
  * process.electronIsoSequence
  * process.mvaTrigV0
  * process.mvaNonTrigV0
  * process.type0PFMEtCorrection
  * process.producePFMETCorrections
  * process.patPFMETsTypeIcorrected
  * process.patDefaultSequence
  * process.patMuonsWithEmbeddedVariables
  * process.patElectronsWithEmbeddedVariables
  )


## Output selection.
process.out.outputCommands = [
    'drop *',
    'keep patElectrons_patElectronsWithEmbeddedVariables_*_*',
    'keep patJets_patJets_*_*',
    'keep patMETs_patPFMETsTypeIcorrected_*_*',
    'keep patMuons_patMuonsWithEmbeddedVariables_*_*',
    'keep patTaus_patTaus_*_*',
    'keep recoVertexs_offlinePrimaryVerticesWithBS_*_*',
    'keep PileupSummaryInfos_*_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_TriggerResults_*_HLT',
    'keep *_genMetTrue_*_*',
    'keep recoTracks_generalTracks_*_*',
    'keep L1GlobalTriggerReadoutRecord_*_*_*',
                             ]

if options.includeSim:
        process.out.outputCommands.extend(['keep recoGenParticles_genParticles_*_*'])
