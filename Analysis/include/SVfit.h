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
        const TLorentzVector& momentum = higgsCandidate.momentum;
        const CandidatePtrVector& higgsDaughters = higgsCandidate.finalStateDaughters;
        if (higgsDaughters.size() != 2)
            throw std::runtime_error("Invalid candidate type for SVfit - it doesn't have 2 daughters");
        svFitStandalone::Vector measuredMET(met.pt * TMath::Sin(met.phi), met.pt * TMath::Cos(met.phi), 0);
        if (met.significanceMatrix.size() != 4)
            throw std::runtime_error("not all elements inside covariance matrix");
        // setup the MET significance
        TMatrixD covMET(2,2);
        covMET[0][0] = met.significanceMatrix.at(0);
        covMET[0][1] = met.significanceMatrix.at(1);
        covMET[1][0] = met.significanceMatrix.at(2);
        covMET[1][1] = met.significanceMatrix.at(3);
        // setup measure tau lepton vectors

        std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
        for (const Candidate* daughter : higgsDaughters){
            svFitStandalone::LorentzVector lepton(daughter->momentum.Px(), daughter->momentum.Py(),
                                              daughter->momentum.Pz(), daughter->momentum.E());
            if (firstDaughter->type == (Candidate::Electron || Candidate::Mu))
                measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kLepDecay, lepton));
            else if (firstDaughter->type == Candidate::Tau)
                measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kHadDecay, lepton));
            else
                throw std::runtime_error("final state daughters are not compatible with a leptonic or hadronic tau decay");
        }

        // construct the class object from the minimal necesarry information
        SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMET, covMET, 0);
        algo.addLogM(false);
        algo.integrateMarkovChain();
        if (!algo.isValidSolution())
            throw std::runtime_error("SVFit algo is not valid");
        double diTauMass = algo.getMass();
        //double diTauMassErr = algo.massUncert(); // mass uncertainty and Pt of Z/Higgs are new features of the Markov Chain integration
        double diTauPt = algo.pt();
        //double diTauPtErr = algo.ptUncert();

        correctedCandidate.momentum.SetPtEtaPhiM(diTauPt,higgsCandidate.momentum.Eta(), higgsCandidate.momentum.Phi(),diTauMass);
        return correctedCandidate;
    }
}
