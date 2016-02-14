## @package patTuple
#  Configuration file to produce PAT-tuples and ROOT-tuples for X->HH->bbTauTau analysis.
#
#  \author Claudio Caputo
#
#  Copyright 2015
#
#  This file is part of X->HH->bbTauTau.
#
#  X->HH->bbTauTau is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  X->HH->bbTauTau is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.


import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('isData',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Include Sim. Default: False")
options.register ('sampleType',
                  'Fall15MC',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Indicates the sample type: Spring15MC, Run2015B, Run2015C, Run2015D")
options.register ('computeHT',
                  'False',
                   VarParsing.multiplicity.singleton,
                   VarParsing.varType.bool,
                  "Compute HT variable and HT binning")

options.parseArguments()

process = cms.Process("USER")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v12')

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

inputSignal_v2 = cms.untracked.vstring("file:768F5AFB-D771-E511-9ABD-B499BAABD280.root")

DYSample = cms.untracked.vstring("/store/user/ccaputo/HHbbtautau/Run2/DYSample_forHT.root")

Diboson = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/1A826380-CB6D-E511-BCFA-0025901D4D6E.root')

##inputSignal_v2 = cms.untracked.vstring(
##                '/store/mc/RunIISpring15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/10563B6E-D871-E511-9513-B499BAABD280.root')
## Input files
process.source = cms.Source("PoolSource",
    fileNames = Diboson
)

## Output file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

from HHbbTauTau.microAODProduction.skimmedBranches_cff import *

if options.computeHT and not options.isData:
    skimmedBranches = cms.untracked.vstring(BaseMCBranches+
                                            ['keep LHEEventProduct_externalLHEProducer__LHE'])
if not options.computeHT and not options.isData:
    skimmedBranches = cms.untracked.vstring(BaseMCBranches)

if options.isData:
    skimmedBranches = cms.untracked.vstring(BaseDATABranches)

allBranches = cms.untracked.vstring(['keep *'])

process.load("RecoMET.METProducers.METSignificance_cfi")
process.load("RecoMET.METProducers.METSignificanceParams_cfi")

process.bbttSkim   = cms.EDFilter("SkimFilterMiniAOD",
                                  vertexSrc  = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
                                  muonSrc  = cms.untracked.InputTag('slimmedMuons'),
                                  electronSrc=cms.untracked.InputTag("slimmedElectrons"),
                                  tauSrc  = cms.untracked.InputTag("slimmedTaus")
                                  )
## Load module for Electron MVA ID
## It will append a value maps the miniAOD, that it's accesible throught a well Handle
## Example code here:
##  https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_7.4.12/ElectronNtupler/plugins/ElectronNtuplerVIDwithMVADemo.cc#L99
## process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
##-------------
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat) ##also compute a maps with the electrons that pass an MVA cut

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
##------------
#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService", fileName = cms.string("syncTree.root") )

#-------------
# SyncTree Producer
#-------------
process.syncNtupler = cms.EDAnalyzer('SyncTreeProducer',

                                 genParticles = cms.InputTag("genParticles"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #

                                 tauSrc           = cms.InputTag("slimmedTaus"),
                                 muonSrc          = cms.InputTag("slimmedMuons"),
                                 vtxSrc           = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 jetSrc           = cms.InputTag("slimmedJets"),
                                 ##pfMETSrc       = cms.InputTag("slimmedMETsNoHF"),
                                 pfMETSrc         = cms.InputTag("slimmedMETs"),
                                 bits             = cms.InputTag("TriggerResults","","HLT"),
                                 prescales        = cms.InputTag("patTrigger"),
                                 objects          = cms.InputTag("selectedPatTrigger"),
                                 lheEventProducts = cms.InputTag("externalLHEProducer"),
                                 genEventInfoProduct = cms.InputTag("generator"),
                                 HTBinning        = cms.bool(options.computeHT),
                                 sampleType = cms.string(options.sampleType),
                                )

process.p = cms.Path(
             #process.electronMVAValueMapProducer*
             process.METSignificance*
             process.egmGsfElectronIDSequence*
             process.bbttSkim*
             process.syncNtupler
	   	    )


process.out = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string('microAOD.root'),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = skimmedBranches,
    #outputCommands = HTBinBranches,
    #outputCommands = allBranches,
    dropMetaData = cms.untracked.string('ALL'),
    fastCloning = cms.untracked.bool(False),
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)

process.endpath= cms.EndPath(process.out)


