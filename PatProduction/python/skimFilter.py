## @package skimFilter
#  Configuration file that defines a filter to skim events for X->HH->bbTauTau analysis.
#
#  \author Rosamaria Venditti (INFN Bari, Bari University)
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

import FWCore.ParameterSet.Config as cms

def applySkim(process):
    process.bbttSkim   = cms.EDFilter("SkimFilter",
                                      vertexSrc  = cms.untracked.InputTag('patVertices'),
                                      muonSrc  = cms.untracked.InputTag('patMuonsWithEmbeddedVariables'),
                                      electronSrc=cms.untracked.InputTag("patElectronsWithEmbeddedVariables"),
                                      tauSrc  = cms.untracked.InputTag("patTausTriggerMatch")
                                      )
    return
