## patTuple.py - configuration file to produce PATtooples for HHbbtautau analysis.

## Import skeleton process.
from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.options.wantSummary = False

## Drop input reco taus
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

from HHbbTauTau.PatProduction.patMET import *
applyMETParameters(process, options.isMC)

from HHbbTauTau.PatProduction.patJets import *
applyJetParameters(process, options.isMC)

## Remove MC matching from the default sequence
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
    'keep GenFilterInfo_*_*_*',
                             ]

if options.includeSim:
        process.out.outputCommands.extend(['keep recoGenParticles_genParticles_*_*'])
