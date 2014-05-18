/*!
 * \file AnalysisTools.h
 * \brief Common tools and definitions suitable for analysis purposes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-07 created
 */

#pragma once

#include <vector>
#include <string>
#include <memory>

#include <TFile.h>
#include <TH1D.h>

#include "Particles.h"
#include "H_BaseAnalyzer.h"

namespace analysis {

static const particles::ParticleCodes TauMuonicDecay = { particles::mu, particles::nu_mu, particles::nu_tau };
static const particles::ParticleCodes TauElectronDecay = { particles::e, particles::nu_e, particles::nu_tau };

inline bool HaveTriggerMatched(const std::vector<std::string>& objectMatchedPaths,
                               const std::vector<std::string>& interestinghltPaths, size_t& n)
{
    for (; n < objectMatchedPaths.size(); ++n){
        for (const std::string& interestingPath : interestinghltPaths){
            const std::string& objectMatchedPath = objectMatchedPaths.at(n);
            const size_t found = objectMatchedPath.find(interestingPath);
            if (found != std::string::npos) return true;
        }
    }
    return false;
}

inline bool HaveTriggerMatched(const std::vector<std::string>& objectMatchedPaths,
                               const std::vector<std::string>& interestinghltPaths)
{
    size_t n = 0;
    return HaveTriggerMatched(objectMatchedPaths, interestinghltPaths, n);
}

inline bool HaveTriggerMatched(const ntuple::TriggerObjectVector& triggerObjects,
                               const std::string& interestingPath,
                               const analysis::Candidate& candidate)
{
    if(candidate.finalStateDaughters.size()) {
        for(const Candidate& daughter : candidate.finalStateDaughters) {
            if(!HaveTriggerMatched(triggerObjects, interestingPath, daughter))
                return false;
        }
        return true;
    }

    for (const ntuple::TriggerObject& triggerObject : triggerObjects){
        TLorentzVector triggerObjectMomentum;
        triggerObjectMomentum.SetPtEtaPhiE(triggerObject.pt, triggerObject.eta, triggerObject.phi, triggerObject.energy);
        for (unsigned n = 0; n < triggerObject.pathNames.size(); ++n){
            const std::string& objectMatchedPath = triggerObject.pathNames.at(n);
            const size_t found = objectMatchedPath.find(interestingPath);
            if (found != std::string::npos && triggerObject.pathValues.at(n) == 1 &&
                    triggerObjectMomentum.DeltaR(candidate.momentum) < 0.5)
                return true;
        }
    }
    return false;
}

inline bool HaveTriggerMatched(const ntuple::TriggerObjectVector& triggerObjects,
                               const std::vector<std::string>& interestingPaths,
                               const analysis::Candidate& candidate)
{
    for (const std::string& interestinPath : interestingPaths){
        if (HaveTriggerMatched(triggerObjects,interestinPath,candidate)) return true;
    }
    return false;
}

inline std::shared_ptr<TH1D> LoadPUWeights(const std::string& reweightFileName, std::shared_ptr<TFile> outputFile )
{
    std::shared_ptr<TFile> reweightFile(new TFile(reweightFileName.c_str(),"READ"));
    if(reweightFile->IsZombie()){
        std::ostringstream ss;
        ss << "reweight file " << reweightFileName << " not found." ;
        throw std::runtime_error(ss.str());
    }
    TObject* originalWeights = reweightFile->Get("weights");
    if (!originalWeights)
        throw std::runtime_error("histograms with weights not found");
    if(outputFile)
        outputFile->cd();
    return std::shared_ptr<TH1D>(static_cast<TH1D*>(originalWeights->Clone("PUweights")));
}

//see AN-13-178
inline double Calculate_MT(const TLorentzVector& lepton_momentum, double met_pt, double met_phi)
{
    const double delta_phi = TVector2::Phi_mpi_pi(lepton_momentum.Phi() - met_phi);
    return std::sqrt( 2.0 * lepton_momentum.Pt() * met_pt * ( 1.0 - std::cos(delta_phi) ) );
}

} // analysis
