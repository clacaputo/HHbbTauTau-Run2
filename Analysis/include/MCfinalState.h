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
    const GenParticle* muon;
    const GenParticle* tau_jet;
    GenParticleVector b_jets;

    bbTauTau() : Higgs_TauTau(nullptr), Higgs_BB(nullptr), muon(nullptr), tau_jet(nullptr) {}


};

}


}
