## patOptions.py - configuration file that defines command-line options for PATtoople production for HHbbtautau analysis.

from FWCore.ParameterSet.VarParsing import VarParsing
from HHbbTauTau.RunTools.readFileList import *

import sys

options = VarParsing('analysis')

def parseAndApplyOptions(process) :
    options.register ('globalTag', 'START53_V7A::All', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "Global Tag to use.")
    options.register ('isMC', True, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Sample Type: MC or data.")
    options.register ('runOnCrab', False, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Indicates if script will be executed on CRAB.")
    options.register ('fileList', 'fileList.txt', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "List of root files to process.")
    options.register ('fileNamePrefix', '', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "Prefix to add to input file names.")
    options.register ('includeSim', False, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Include Sim information.")
    options.register ('treeOutput', 'Tree.root', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "Tree root file.")

    if len(sys.argv) > 0:
        last = sys.argv.pop()
        sys.argv.extend(last.split(","))

    options.parseArguments()

    process.GlobalTag.globaltag = options.globalTag

    if not options.runOnCrab:
        readFileList(process.source.fileNames, options.fileList, options.fileNamePrefix)
        process.out.fileName = options.outputFile
        process.maxEvents.input = options.maxEvents

    return
