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

struct SVFitResults
{
    double mass_Vegas;
    TLorentzVector momentum_MarkovChain;
};

struct SVFitResultsUncertainties
{
    SVFitResults svfitResults;
    SVFitResults svfitResults_Up;
    SVFitResults svfitResults_Down;
};

inline SVFitResults CorrectMassBySVfit(const Candidate& higgsCandidate, const ntuple::MET& met, double tauESfactor)
{
    static const bool debug = false;
    SVFitResults svfitresults;

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

    // construct the class object from the minimal necesarry information
    NSVfitStandalone::NSVfitStandaloneAlgorithm algoVegas(measuredTauLeptons, measuredMET, covMET, 0);
    algoVegas.addLogM(false);
    algoVegas.integrateVEGAS();

    if (!algoVegas.isValidSolution()) {
        std::cerr << "Can't fit\n";
        svfitresults.mass_Vegas = higgsCandidate.momentum.M();
        //throw std::runtime_error("SVFit algo is not valid");
    }

    svfitresults.mass_Vegas = algoVegas.mass();

    NSVfitStandalone::NSVfitStandaloneAlgorithm algoMC(measuredTauLeptons, measuredMET, covMET, 0);
    algoMC.integrateMarkovChain();
    if (!algoMC.isValidSolution()) {
        std::cerr << "Can't fit\n";
        svfitresults.momentum_MarkovChain.SetPtEtaPhiM(higgsCandidate.momentum.Pt(), higgsCandidate.momentum.Eta(),
                                                       higgsCandidate.momentum.Phi(), higgsCandidate.momentum.M());
        //throw std::runtime_error("SVFit algo is not valid");
    }


    if(debug) {
        std::cout << "original mass = " << higgsCandidate.momentum.M()
                  << "\nOriginal momentum = " << higgsCandidate.momentum
                  << "\nFirst daughter momentum = " << higgsCandidate.daughters.at(0).momentum
                  << "\nSecond daughter momentum = " << higgsCandidate.daughters.at(1).momentum
                     << "\nSVfit mass Vegas = " << algoVegas.mass()
                  << "\nSVfit mass MC = " << algoMC.mass() << ", SVfit pt MC = " << algoMC.pt() << std::endl;
    }
    svfitresults.momentum_MarkovChain.SetPtEtaPhiM(algoMC.pt(), higgsCandidate.momentum.Eta(),
                                                   higgsCandidate.momentum.Phi(), algoMC.mass());

    return svfitresults;
}

inline SVFitResultsUncertainties CollectSVFitResults(const Candidate& higgsCandidate, const ntuple::MET& met)
{
    SVFitResultsUncertainties svfitResultsUncertainties;
    svfitResultsUncertainties.svfitResults = CorrectMassBySVfit(higgsCandidate,met,1.);
    svfitResultsUncertainties.svfitResults_Up = CorrectMassBySVfit(higgsCandidate,met,1.03);
    svfitResultsUncertainties.svfitResults_Down = CorrectMassBySVfit(higgsCandidate,met,0.97);
    return svfitResultsUncertainties;
}

} // analysis

