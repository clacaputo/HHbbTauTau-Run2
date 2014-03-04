## Import skeleton process.
from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.options.wantSummary = False

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


## Customize output.
#from PhysicsTools.PatAlgos.tools.coreTools import *
## remove MC matching from the default sequence
# removeMCMatching(process, ['Muons'])
# runOnData(process)

## remove certain objects from the default sequence
# removeAllPATObjectsBut(process, ['Muons'])
# removeSpecificPATObjects(process, ['Electrons', 'Muons', 'Taus'])
