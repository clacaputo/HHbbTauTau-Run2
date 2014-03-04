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

## Include high Pt tau reconstruction sequence.
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

## Define process path.
process.p = cms.Path( process.PFTau * process.patDefaultSequence )

## MET corrections.
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')
process.patMETsPF.metSource = cms.InputTag("pfMet")

## Output selection.
process.out.outputCommands = [
    'drop *',
    'keep patElectrons_patElectrons_*_*',
    'keep patJets_patJets_*_*',
    'keep patMETs_patMETsPF_*_*',
    'keep patMuons_patMuons_*_*',
    'keep patTaus_patTaus_*_*',
    'keep recoGenParticles_genParticles_*_*',
    'keep recoVertexs_offlinePrimaryVertices_*_*',
    'keep PileupSummaryInfos_*_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_TriggerResults_*_HLT',
    'keep *_genMetTrue_*_*'
                             ]
