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

namespace analysis {

static const particles::ParticleCodes TauMuonicDecay = { particles::mu, particles::nu_mu, particles::nu_tau };
static const particles::ParticleCodes TauElectronDecay = { particles::e, particles::nu_e, particles::nu_tau };

inline bool HaveTriggerMatched(const std::vector<std::string>& objectMatchedPaths,
                               const std::vector<std::string>& interestingPaths, size_t& n)
{
    for (; n < objectMatchedPaths.size(); ++n){
        for (size_t k = 0; k < interestingPaths.size(); ++k){
            const std::string& objectMatchedPath = objectMatchedPaths.at(n);
            const std::string& interestingPath = interestingPaths.at(k);
            const size_t found = objectMatchedPath.find(interestingPath);
            if (found != std::string::npos) return true;
        }
    }
    return false;
}

inline bool HaveTriggerMatched(const std::vector<std::string>& objectMatchedPaths,
                               const std::vector<std::string>& interestingPaths)
{
    size_t n = 0;
    return HaveTriggerMatched(objectMatchedPaths, interestingPaths, n);
}

inline bool HaveTriggerMatched(const )
{

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

} // analysis
