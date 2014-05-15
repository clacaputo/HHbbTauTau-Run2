## patTuple.py - configuration file to produce PATtooples for HHbbtautau analysis.

## Import skeleton process.
from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options.wantSummary = False

## Drop input reco taus to ensure that the new tau collection will be used.
process.source.dropDescendantsOfDroppedBranches = cms.untracked.bool(False)
process.source.inputCommands = cms.untracked.vstring(
        'keep *',
        'drop recoPFTaus_*_*_*'
    )

from HHbbTauTau.PatProduction.patOptions import *
parseAndApplyOptions(process)

from HHbbTauTau.PatProduction.patMuons import *
applyMuonParameters(process)

from HHbbTauTau.PatProduction.patTaus import *
applyTauParameters(process)

from HHbbTauTau.PatProduction.patElectrons import *
applyElectronParameters(process, options.isMC)

from HHbbTauTau.PatProduction.patJets import *
applyJetParameters(process, options.isMC)

from HHbbTauTau.PatProduction.patMET import *
applyMETParameters(process, options.isMC)
applyMVAMETParamteres(process, options.isMC)

from HHbbTauTau.PatProduction.patVertices import *
applyVertexParameters(process)

from HHbbTauTau.PatProduction.patTrigger import *
applyTriggerParameters(process)

from HHbbTauTau.PatProduction.skimFilter import *
applySkim(process)

from HHbbTauTau.PatProduction.treeProduction import *
addTreeSequence(process, options.includeSim, options.treeOutput)

## Remove MC matching from the default sequence
if not options.isMC:
    removeMCMatching(process, ['All'])
#    removeMCMatching(process, ['METs'], "TC")
#    removeMCMatching(process, ['METs'], "PF")
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
    process.PFTau *
    process.pfParticleSelectionSequence *
    process.muIsoSequence *
    process.electronIsoSequence *
    process.mvaTrigV0 *
    process.mvaNonTrigV0 *
    process.type0PFMEtCorrection *
    process.producePFMETCorrections *
    process.pfMEtMVAsequence *
    process.puJetId *
    process.pileupJetIdProducer *
    process.patDefaultSequence *
    process.patMuonsWithEmbeddedVariables *
    process.patElectronsWithEmbeddedVariables *
    process.patMETsMVAsequence *
    process.patVertices *
    process.bbttSkim *
    process.mainTreeContentSequence *
    process.simTreeContentSequence
)

## Output selection.
process.out.outputCommands = [
    'drop *',
    'keep patElectrons_patElectronsWithEmbeddedVariables_*_*',
    'keep patJets_patJetsTriggerMatch_*_*',
    'keep patMETs_*_*_*',
    'keep patMuons_patMuonsWithEmbeddedVariables_*_*',
    'keep patTaus_patTausTriggerMatch_*_*',
    'keep patVertexs_patVertices_*_*',
    'keep recoVertexs_offlinePrimaryVerticesWithBS_*_*',
    'keep PileupSummaryInfos_*_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_TriggerResults_*_HLT',
    'keep L1GlobalTriggerReadoutRecord_*_*_*',
    'keep GenFilterInfo_*_*_*',
    'keep *_patTrigger_*_*',
    'keep patTriggerEvent_*_*_*'
                             ]

if options.includeSim:
    process.out.outputCommands.extend(['keep recoGenParticles_genParticles_*_*'])
    process.out.outputCommands.extend(['keep *_genMetTrue_*_*'])

#print process.dumpPython()
