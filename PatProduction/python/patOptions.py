## @package patOptions
#  Configuration file that defines command-line options for PATtoople production for X->HH->bbTauTau analysis.
#
#  \author Konstantin Androsov (Siena University, INFN Pisa)
#  \author Maria Teresa Grippo (Siena University, INFN Pisa)
#
#  Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
#                 Maria Teresa Grippo <grippomariateresa@gmail.com>
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

from FWCore.ParameterSet.VarParsing import VarParsing
from HHbbTauTau.RunTools.readFileList import *

import sys

options = VarParsing('analysis')

def parseAndApplyOptions(process) :
    options.register ('globalTag', 'START53_V7A::All', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "Global Tag to use.")
    options.register ('isMC', True, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Sample Type: MC or data.")
    options.register ('isEmbedded', False, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Indicates if sample is made using embedding technique.")
    options.register ('runOnCrab', False, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Indicates if script will be executed on CRAB.")
    options.register ('fileList', 'fileList.txt', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "List of root files to process.")
    options.register ('fileNamePrefix', '', VarParsing.multiplicity.singleton,
                      VarParsing.varType.string, "Prefix to add to input file names.")
    options.register ('includeSim', False, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Include Sim information.")
    options.register ('keepPat', False, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Keep PAT information in the output file.")
    options.register ('runTree', True, VarParsing.multiplicity.singleton,
                      VarParsing.varType.bool, "Run TREE production sequence.")
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
