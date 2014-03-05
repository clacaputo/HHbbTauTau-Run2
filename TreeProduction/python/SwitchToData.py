import FWCore.ParameterSet.Config as cms
import copy

import PhysicsTools.PatAlgos.tools.coreTools as coreTools

def _setattr_ifexists(obj, attrName, attrValue):
	if hasattr(obj, attrName):
		setattr(obj, attrName, attrValue)

def switchToData(process):
  #remove MC matching from standard PAT sequences
  #from HHbbTauTau.TreeProduction.SwitchToData import switchToData
  # remove MC matching

  coreTools.runOnData(process, outputInProcess = False)      
  #coreTools.removeMCMatching(process, ["All"], outputInProcess = False)
  coreTools.removeMCMatching(process, ['METs'], "TC", outputInProcess = False)
  coreTools.removeMCMatching(process, ['METs'], "PF", outputInProcess = False)
  coreTools.removeMCMatching(process, ['METs'], "AK5Calo", outputInProcess = False)

  process.patDefaultSequence.remove(process.patJetPartonMatch)
  process.patDefaultSequence.remove(process.patJetPartonMatchAK5PF)
  process.patDefaultSequence.remove(process.patJetPartonMatchAK5Calo)
  process.patDefaultSequence.remove(process.patJetGenJetMatchAK5PF)
  process.patDefaultSequence.remove(process.patJetGenJetMatchAK5Calo)
  process.patDefaultSequence.remove(process.patJetFlavourId)
  process.patDefaultSequence.remove(process.patJetPartons)
  process.patDefaultSequence.remove(process.patJetPartonAssociation)
  process.patDefaultSequence.remove(process.patJetPartonAssociationAK5PF)
  process.patDefaultSequence.remove(process.patJetPartonAssociationAK5Calo)  
  process.patDefaultSequence.remove(process.patJetFlavourAssociation)
  process.patDefaultSequence.remove(process.patJetFlavourAssociationAK5PF)
  process.patDefaultSequence.remove(process.patJetFlavourAssociationAK5Calo)  
