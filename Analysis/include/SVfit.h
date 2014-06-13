/*!
 * \file SVfit.h
 * \brief Definition of wrappers for SVfit_standalone code.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-02 created
 */

#pragma once

//#include "HHbbTauTau/SVfit_standalone/src/LikelihoodFunctions.cc"
//#include "HHbbTauTau/SVfit_standalone/src/SVfitStandaloneAlgorithm.cc"
//#include "HHbbTauTau/SVfit_standalone/src/svFitStandaloneAuxFunctions.cc"
//#include "HHbbTauTau/SVfit_standalone/src/SVfitStandaloneLikelihood.cc"
//#include "HHbbTauTau/SVfit_standalone/src/SVfitStandaloneMarkovChainIntegrator.cc"

#include "HHbbTauTau/SVfit/source/generalAuxFunctions.cc"
#include "HHbbTauTau/SVfit/source/LikelihoodFunctions.cc"
#include "HHbbTauTau/SVfit/source/NSVfitStandaloneAlgorithm.cc"
#include "HHbbTauTau/SVfit/source/NSVfitStandaloneLikelihood.cc"
#include "HHbbTauTau/SVfit/source/svFitAuxFunctions.cc"

#include "Candidate.h"
#include "HHbbTauTau/TreeProduction/interface/MET.h"

namespace analysis {
double CorrectMassBySVfit(const Candidate& higgsCandidate, const ntuple::MET& met, double tauESfactor)
{
    static const bool debug = false;

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
    NSVfitStandalone::NSVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMET, covMET, 0);
    algo.addLogM(false);
    algo.integrateVEGAS();
    if (!algo.isValidSolution()) {
        std::cerr << "Can't fit\n";
        return higgsCandidate.momentum.M();
        //throw std::runtime_error("SVFit algo is not valid");
    }

//    Candidate correctedCandidate(higgsCandidate);
//    correctedCandidate.momentum.SetPtEtaPhiM(algo.pt(),
//                                             higgsCandidate.momentum.Eta(),
//                                             higgsCandidate.momentum.Phi(),
//                                             algo.mass());

    if(debug) {
        std::cout << "SVfit mass = " << algo.mass() << ", SVfit pt = " << algo.pt() << std::endl;
    }
    return algo.mass();
//    return correctedCandidate;
}

} // analysis
