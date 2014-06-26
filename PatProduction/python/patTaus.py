## @package patTaus
#  Configuration file that defines parameters related to PAT Tau objects.
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

import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.tauTools import *

## Should work for both standard and high pt taus
def applyTauParameters(process):
    # Include tau reconstruction sequence.
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    switchToPFTauHPS(process)

    process.patTaus.embedLeadPFCand = cms.bool(True)
    process.patTaus.embedSignalPFCands = cms.bool(True)
    process.patTaus.embedIsolationPFCands = cms.bool(True)
    process.patTaus.embedLeadTrack = cms.bool(True)
    process.patTaus.embedSignalTracks = cms.bool(True)
    process.patTaus.embedIsolationTracks = cms.bool(True)
    process.patTaus.embedIsolationPFChargedHadrCands = cms.bool(True)
    process.patTaus.embedIsolationPFNeutralHadrCands = cms.bool(True)
    process.patTaus.embedIsolationPFGammaCands = cms.bool(True)
    process.patTaus.embedSignalPFChargedHadrCands = cms.bool(True)
    process.patTaus.embedSignalPFNeutralHadrCands = cms.bool(True)
    process.patTaus.embedSignalPFGammaCands = cms.bool(True)
    process.patTaus.embedLeadPFChargedHadrCand = cms.bool(True)
    process.patTaus.embedLeadPFNeutralCand = cms.bool(True)

    return
