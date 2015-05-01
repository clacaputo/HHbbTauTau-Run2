## @package patMuons
#  Configuration file that defines parameters related to PAT Muon objects.
#
#  \author Konstantin Androsov (University of Siena, INFN Pisa)
#  \author Maria Teresa Grippo (University of Siena, INFN Pisa)
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

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFMuonIso
from CommonTools.ParticleFlow.pfParticleSelection_cff import pfParticleSelectionSequence

def applyMuonParameters(process):
    process.muIsoSequence       = setupPFMuonIso(process,'muons')
    process.pfParticleSelectionSequence = pfParticleSelectionSequence

    process.patMuons.isoDeposits = cms.PSet(
        pfAllParticles   = cms.InputTag("muPFIsoDepositPUPFIso"),              # all PU   CH+MU+E
        pfChargedHadrons = cms.InputTag("muPFIsoDepositChargedPFIso"),         # all noPU CH
        pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutralPFIso"),         # all NH
        pfPhotons        = cms.InputTag("muPFIsoDepositGammaPFIso"),           # all PH
        user = cms.VInputTag( cms.InputTag("muPFIsoDepositChargedAllPFIso") )  # all noPU CH+MU+E
        )

    process.patMuons.isolationValues = cms.PSet(
        pfAllParticles   = cms.InputTag("muPFIsoValuePU04PFIso"),
        pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04PFIso"),
        pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04PFIso"),
        pfPhotons        = cms.InputTag("muPFIsoValueGamma04PFIso"),
        user = cms.VInputTag( cms.InputTag("muPFIsoValueChargedAll04PFIso") )
        )

    process.patMuons.embedTrack = cms.bool(True)

    process.patMuonsWithEmbeddedVariables = cms.EDProducer('MuonsUserEmbedded',
        muonTag = cms.InputTag("patMuonsTriggerMatch"),
        vertexTag = cms.InputTag("offlinePrimaryVerticesWithBS")
        )

    return
