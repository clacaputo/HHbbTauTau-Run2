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
    if (higgsCandidate.finalStateDaughters.size() != 2)
        throw std::runtime_error("Invalid candidate type for SVfit - it doesn't have 2 daughters");

    // setup the MET coordinates and significance
    const svFitStandalone::Vector measuredMET(met.pt * std::sin(met.phi), met.pt * std::cos(met.phi), 0.0);
    const TMatrixD covMET = ntuple::VectorToSignificanceMatrix(met.significanceMatrix);

    // setup measure tau lepton vectors
    std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
    for (const Candidate* daughter : higgsCandidate.finalStateDaughters){
        svFitStandalone::LorentzVector lepton(daughter->momentum.Px(), daughter->momentum.Py(),
                                              daughter->momentum.Pz(), daughter->momentum.E());
        svFitStandalone::kDecayType decayType;
        if (daughter->type == Candidate::Electron || daughter->type == Candidate::Mu)
            decayType = svFitStandalone::kLepDecay;
        else if (daughter->type == Candidate::Tau)
            decayType = svFitStandalone::kHadDecay;
        else
            throw std::runtime_error("final state daughters are not compatible with a leptonic or hadronic tau decay");
        measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(decayType, lepton));
    }

    // construct the class object from the minimal necesarry information
    SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMET, covMET, 0);
    algo.addLogM(false);
    algo.integrateMarkovChain();
    if (!algo.isValidSolution()) {
        std::cerr << "Can't fit\n";
        return Candidate(higgsCandidate);
        //throw std::runtime_error("SVFit algo is not valid");
    }

    Candidate correctedCandidate(higgsCandidate);
    correctedCandidate.momentum.SetPtEtaPhiM(algo.pt(),
                                             higgsCandidate.momentum.Eta(),
                                             higgsCandidate.momentum.Phi(),
                                             algo.getMass());
    return correctedCandidate;
}

template<typename Histogram>
CandidateVector CorrectMassBySVfit(const CandidateVector& higgsCandidates, const ntuple::MET& met, Histogram& hist,
                                   double weight = 1.0)
{
    CandidateVector result;
    for(const auto& candidate : higgsCandidates) {
        const Candidate& correctedCandidate = CorrectMassBySVfit(candidate, met);
        hist.Fill(correctedCandidate.momentum.M(), weight);
        result.push_back(correctedCandidate);
    }
    return result;
}

} // analysis
