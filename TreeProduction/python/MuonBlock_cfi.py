## @package MuonBlock_cfi
#  Configuration file that defines the producer of ROOT-tuple for muons.
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

muonBlock = cms.EDAnalyzer("MuonBlock",
  muonSrc         = cms.InputTag('patMuonsWithEmbeddedVariables'),
  vertexSrc       = cms.InputTag('patVertices'),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
  beamSpotCorr    = cms.bool(True),
  muonID          = cms.string('GlobalMuonPromptTight')
)
