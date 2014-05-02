/*!
 * \file SVfit.h
 * \brief Definition of wrappers for SVfit_standalone code.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-02 created
 */

#pragma once

#include "HHbbTauTau/SVfit_standalone/src/LikelihoodFunctions.cc"
#include "HHbbTauTau/SVfit_standalone/src/SVfitStandaloneAlgorithm.cc"
#include "HHbbTauTau/SVfit_standalone/src/svFitStandaloneAuxFunctions.cc"
#include "HHbbTauTau/SVfit_standalone/src/SVfitStandaloneLikelihood.cc"
#include "HHbbTauTau/SVfit_standalone/src/SVfitStandaloneMarkovChainIntegrator.cc"

#include "Candidate.h"
#include "HHbbTauTau/TreeProduction/interface/MET.h"

namespace analysis {
    Candidate CorrectMassBySVfit(const Candidate& higgsCandidate, const ntuple::MET& met)
    {
        if(higgsCandidate.type != Candidate::Higgs)
            throw std::runtime_error("Invalid candidate type for SVfit");
        Candidate correctedCandidate(higgsCandidate);
//        const TLorentzVector& momentum = higgsCandidate.momentum;
//        met.
//        svFitStandalone::Vector met_direction()
//        correctedCandidate.momentum.SetPtEtaPhiM();
        return correctedCandidate;
    }
}
