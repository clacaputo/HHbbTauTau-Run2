/*!
 * \file MCfinalState.h
 * \brief Definition of MCfinalState class for analysis.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-19 created
 */

#pragma once
#include "GenParticle.h"


namespace analysis {

namespace finalState {

struct bbTauTau{

    const GenParticle* resonance;
    const GenParticle* Higgs_TauTau;
    const GenParticle* Higgs_BB;
    GenParticleVector b_jets;
    GenParticleVector taus;

    bbTauTau() : resonance(nullptr), Higgs_TauTau(nullptr), Higgs_BB(nullptr){}

    virtual ~bbTauTau() {}
};

struct bbETaujet : public bbTauTau {
    const GenParticle* electron;
    const GenParticle* tau_jet;

    bbETaujet() : electron(nullptr), tau_jet(nullptr) {}
};

struct bbMuTaujet : public bbTauTau{
    const GenParticle* muon;
    const GenParticle* tau_jet;

    bbMuTaujet() : muon(nullptr), tau_jet(nullptr) {}
};

struct bbTaujetTaujet : public bbTauTau {};

} // finalState
} // analysis
