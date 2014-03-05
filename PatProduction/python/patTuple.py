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

options.register ('fileList',
                  'fileList.py',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "List of root files to process.")

options.parseArguments()

process.GlobalTag.globaltag = options.globalTag
process.maxEvents.input = options.maxEvents

fileList = cms.untracked.vstring()
process.source.fileNames = fileList
execfile(options.fileList)
process.out.fileName = options.outputFile

## Muons
process.patMuons.embedTrack = cms.bool(True)

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


## Define process path.
process.p = cms.Path(
    process.PFTau
  * process.type0PFMEtCorrection
  * process.producePFMETCorrections
  * process.patPFMETsTypeIcorrected
  * process.patDefaultSequence
  )


## Output selection.
process.out.outputCommands = [
    'drop *',
    'keep patElectrons_patElectrons_*_*',
    'keep patJets_patJets_*_*',
    'keep patMETs_patPFMETsTypeIcorrected_*_*',
    'keep patMuons_patMuons_*_*',
    'keep patTaus_patTaus_*_*',
    'keep recoGenParticles_genParticles_*_*',
    'keep recoVertexs_offlinePrimaryVerticesWithBS_*_*',
    'keep PileupSummaryInfos_*_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_TriggerResults_*_HLT',
    'keep *_genMetTrue_*_*',
    'keep recoTracks_generalTracks_*_*',
    'keep L1GlobalTriggerReadoutRecord_*_*_*',
                             ]
