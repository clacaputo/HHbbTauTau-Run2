import FWCore.ParameterSet.Config as cms

from HHbbTauTau.TreeProduction.EventBlock_cfi import eventBlock
from HHbbTauTau.TreeProduction.VertexBlock_cfi import vertexBlock
from HHbbTauTau.TreeProduction.ElectronBlock_cfi import electronBlock
from HHbbTauTau.TreeProduction.MuonBlock_cfi import muonBlock
from HHbbTauTau.TreeProduction.TauBlock_cfi import tauBlock
from HHbbTauTau.TreeProduction.JetBlock_cfi import jetBlock
from HHbbTauTau.TreeProduction.METBlock_cfi import metBlock
from HHbbTauTau.TreeProduction.TriggerBlock_cfi import triggerBlock
from HHbbTauTau.TreeProduction.TriggerObjectBlock_cfi import triggerObjectBlock

treeContentSequence = cms.Sequence(
   eventBlock
 + vertexBlock
 + electronBlock
 + muonBlock
 + tauBlock
 + jetBlock  
 + metBlock
 + triggerBlock
 + triggerObjectBlock
)
