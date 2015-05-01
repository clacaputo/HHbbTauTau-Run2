## @package treeProducer
#  Configuration file to produce ROOT-tuples for X->HH->bbTauTau analysis.
#
#  \author Subir Sarkar
#  \author Rosamaria Venditti (INFN Bari, Bari University)
#  \author Konstantin Androsov (University of Siena, INFN Pisa)
#  \author Maria Teresa Grippo (University of Siena, INFN Pisa)
#
#  Copyright 2011-2013 Subir Sarkar, Rosamaria Venditti (INFN Bari, Bari University)
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

import FWCore.ParameterSet.Config as cms
process = cms.Process("HTauTauTree")
#------------------------
# Message Logger Settings
#------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500


from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('globalTag',
                  'START53_V21::All',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Global Tag to use. Default: START53_V21::All")
options.register ('fileList',
                  'fileList.txt',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "List of root files to process.")
options.register ('fileNamePrefix',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Prefix to add to input file names.")
options.register ('includeSim',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Include Sim. Default: False")

options.parseArguments()

#--------------------------------------
# Event Source & # of Events to process
#---------------------------------------

from HHbbTauTau.RunTools.readFileList import *

process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring())
readFileList(process.source.fileNames, options.fileList, options.fileNamePrefix)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

#-----------------------------
# Geometry
#-----------------------------
process.load("Configuration.Geometry.GeometryIdeal_cff")
#-----------------------------
# Magnetic Field
#-----------------------------
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#-------------
# Global Tag
#-------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = options.globalTag

#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile) )
#--------------------------------------------------
# VHTauTau Tree Specific
#--------------------------------------------------
process.load("HHbbTauTau.TreeProduction.TreeContentConfig_cff")

#-------------------------------------------------------
# PAT 
#------------------------------------------------------
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")

import PhysicsTools.PatAlgos.tools.tauTools as tauTools
#tauTools.switchToPFTauHPS(process) # For HPS Taus

## --
## Switch on PAT trigger
## --
import PhysicsTools.PatAlgos.tools.trigTools as trigTools
trigTools.switchOnTrigger( process, outputModule='' ) # This is optional and can be omitted.

process.mainTreeContentSequence = cms.Sequence(
    process.eventBlock
  + process.vertexBlock
  + process.electronBlock
  + process.jetBlock
  + process.metBlock
  + process.muonBlock
  + process.tauBlock
  + process.triggerBlock
  + process.triggerObjectBlock
)

process.simTreeContentSequence = cms.Sequence()
if options.includeSim:
    process.simTreeContentSequence = cms.Sequence(process.genParticleBlock + process.genMETBlock)


process.p = cms.Path(
  process.mainTreeContentSequence +
  process.simTreeContentSequence
)
