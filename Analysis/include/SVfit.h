/*!
 * \file SVfit.h
 * \brief Definition of wrappers for SVfit_standalone code.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-05-02 created
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

#include "SVfit/source/generalAuxFunctions.cc"
#include "SVfit/source/LikelihoodFunctions.cc"
#include "SVfit/source/MarkovChainIntegrator.cc"
#include "SVfit/source/NSVfitStandaloneAlgorithm.cc"
#include "SVfit/source/NSVfitStandaloneLikelihood.cc"
#include "SVfit/source/svFitAuxFunctions.cc"

#include "TreeProduction/interface/MET.h"

#include "AnalysisBase/include/Candidate.h"

namespace analysis {

namespace sv_fit {

enum class FitAlgorithm { Vegas, MarkovChain };
std::ostream& operator<< (std::ostream& s, FitAlgorithm algo)
{
    static const std::map<FitAlgorithm, std::string> algo_names = { { FitAlgorithm::Vegas, "VEGAS" },
                                                                    { FitAlgorithm::MarkovChain, "Markov chain" } };
    s << algo_names.at(algo);
    return s;
}

struct FitResults {
    static constexpr double default_value = std::numeric_limits<double>::lowest();

    bool has_valid_mass;
    double mass;

    bool has_valid_pt;
    double pt;

    FitResults() : has_valid_mass(false), mass(default_value), has_valid_pt(false), pt(default_value) {}
};

struct FitResultsWithUncertainties {
    FitResults fit_vegas;
    FitResults fit_vegas_down;
    FitResults fit_vegas_up;

    FitResults fit_mc;
    FitResults fit_mc_down;
    FitResults fit_mc_up;
};

inline FitResults Fit(FitAlgorithm fitAlgorithm, const Candidate& higgsCandidate, const ntuple::MET& met,
                      double tauESfactor)
{
    static const bool debug = false;
    FitResults result;

    if(higgsCandidate.type != Candidate::Higgs)
        throw std::runtime_error("Invalid candidate type for SVfit");
    if (higgsCandidate.finalStateDaughters.size() != 2)
        throw std::runtime_error("Invalid candidate type for SVfit - it doesn't have 2 daughters");

    // setup the MET coordinates and significance
    const NSVfitStandalone::Vector measuredMET(met.pt * std::cos(met.phi), met.pt * std::sin(met.phi), 0.0);
    const TMatrixD covMET = ntuple::VectorToSignificanceMatrix(met.significanceMatrix);

    // setup measure tau lepton vectors
    std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
    for (const Candidate& daughter : higgsCandidate.finalStateDaughters){
        NSVfitStandalone::LorentzVector lepton(daughter.momentum.Px(), daughter.momentum.Py(),
                                              daughter.momentum.Pz(), daughter.momentum.E());
        NSVfitStandalone::kDecayType decayType;
        if (daughter.type == Candidate::Electron || daughter.type == Candidate::Mu)
            decayType = NSVfitStandalone::kLepDecay;
        else if (daughter.type == Candidate::Tau){
            decayType = NSVfitStandalone::kHadDecay;
            lepton *= tauESfactor;
        }
        else
            throw std::runtime_error("final state daughters are not compatible with a leptonic or hadronic tau decay");
        measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(decayType, lepton));
    }

    NSVfitStandalone::NSVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMET, covMET, 0);
    algo.addLogM(false);

    if(fitAlgorithm == FitAlgorithm::Vegas)
        algo.integrateVEGAS();
    else if(fitAlgorithm == FitAlgorithm::MarkovChain)
        algo.integrateMarkovChain();
    else
        throw std::runtime_error("Unsupported algorithm.");

    if(algo.isValidSolution()) {
        result.mass = algo.mass();
        result.has_valid_mass = true;
        if(fitAlgorithm == FitAlgorithm::MarkovChain) {
            result.pt = algo.pt();
            result.has_valid_pt = true;
        }
    } else
        std::cerr << "Can't fit with " << fitAlgorithm << std::endl;

    if(debug) {
        std::cout << "Original mass = " << higgsCandidate.momentum.M()
                  << "\nOriginal momentum = " << higgsCandidate.momentum
                  << "\nFirst daughter momentum = " << higgsCandidate.daughters.at(0).momentum
                  << "\nSecond daughter momentum = " << higgsCandidate.daughters.at(1).momentum
                  << "\nSVfit algorithm = " << fitAlgorithm;
        if(result.has_valid_mass)
            std::cout << "\nSVfit mass = " << result.mass;
        if(result.has_valid_pt)
            std::cout << "\nSVfit pt = " << result.pt;
        std::cout << std::endl;
    }

    return result;
}

inline FitResultsWithUncertainties FitWithUncertainties(const Candidate& higgsCandidate, const ntuple::MET& met,
                                                        double tau_energy_uncertainty, bool fitWithVegas,
                                                        bool fitWithMarkovChain)
{
    FitResultsWithUncertainties result;
    if(fitWithVegas) {
        result.fit_vegas = Fit(FitAlgorithm::Vegas, higgsCandidate, met, 1.);
        result.fit_vegas_up = Fit(FitAlgorithm::Vegas, higgsCandidate, met, 1 + tau_energy_uncertainty);
        result.fit_vegas_down = Fit(FitAlgorithm::Vegas, higgsCandidate, met, 1 - tau_energy_uncertainty);
    }
    if(fitWithMarkovChain) {
        result.fit_mc = Fit(FitAlgorithm::MarkovChain, higgsCandidate, met, 1.);
        result.fit_mc_up = Fit(FitAlgorithm::MarkovChain, higgsCandidate, met, 1 + tau_energy_uncertainty);
        result.fit_mc_down = Fit(FitAlgorithm::MarkovChain, higgsCandidate, met, 1 - tau_energy_uncertainty);
    }
    return result;
}

} // namespace sv_fit

} // namespace analysis

