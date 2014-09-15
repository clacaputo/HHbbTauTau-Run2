/*!
 * \file MCfinalState.h
 * \brief Definition of MCfinalState class for analysis.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-03-19 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "GenParticle.h"

namespace analysis {

namespace finalState {

struct bbTauTau{

    const GenParticle* resonance;
    const GenParticle* Higgs_TauTau;
    const GenParticle* Higgs_BB;
    VisibleGenObjectVector b_jets;
    VisibleGenObjectVector taus;
    VisibleGenObjectVector hadronic_taus;

    bbTauTau() : resonance(nullptr), Higgs_TauTau(nullptr), Higgs_BB(nullptr){}

    virtual ~bbTauTau() {}
    virtual void Reset()
    {
        resonance = Higgs_TauTau = Higgs_BB = nullptr;
        b_jets.clear();
        taus.clear();
        hadronic_taus.clear();
    }
};

struct bbETaujet : public bbTauTau {
    const GenParticle* electron;
    const VisibleGenObject* tau_jet;

    bbETaujet() : electron(nullptr), tau_jet(nullptr) {}
    virtual void Reset() override
    {
        bbTauTau::Reset();
        electron = nullptr;
        tau_jet = nullptr;
    }
};

struct bbMuTaujet : public bbTauTau{
    const GenParticle* muon;
    const VisibleGenObject* tau_jet;

    bbMuTaujet() : muon(nullptr), tau_jet(nullptr) {}
    virtual void Reset() override
    {
        bbTauTau::Reset();
        muon = nullptr;
        tau_jet = nullptr;
    }
};

struct bbTaujetTaujet : public bbTauTau {
    const VisibleGenObject* leading_tau_jet;
    const VisibleGenObject* subleading_tau_jet;

    bbTaujetTaujet() : leading_tau_jet(nullptr), subleading_tau_jet(nullptr) {}
    virtual void Reset() override
    {
        bbTauTau::Reset();
        leading_tau_jet = subleading_tau_jet = nullptr;
    }
};

struct TauTau {
    const GenParticle* resonance;
    GenParticlePtrVector taus;
    GenParticlePtrVector hadronic_taus;
    TauTau() : resonance(nullptr) {}
    virtual ~TauTau() {}
};

struct ETaujet : public TauTau {
    const GenParticle* electron;
    const GenParticle* tau_jet;

    ETaujet() : electron(nullptr), tau_jet(nullptr) {}
    ETaujet(const TauTau& other) : TauTau(other), electron(nullptr), tau_jet(nullptr) {}
};

struct MuTaujet : public TauTau{
    const GenParticle* muon;
    const GenParticle* tau_jet;

    MuTaujet() : muon(nullptr), tau_jet(nullptr) {}
    MuTaujet(const TauTau& other) : TauTau(other), muon(nullptr), tau_jet(nullptr) {}
};

struct TaujetTaujet : public TauTau {
    const GenParticle* leading_tau_jet;
    const GenParticle* subleading_tau_jet;

    TaujetTaujet() : leading_tau_jet(nullptr), subleading_tau_jet(nullptr) {}
    TaujetTaujet(const TauTau& other) : TauTau(other), leading_tau_jet(nullptr), subleading_tau_jet(nullptr) {}
};

} // finalState
} // analysis
