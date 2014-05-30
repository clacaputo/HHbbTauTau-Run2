/*!
 * \file MvaMet.h
 * \brief Definition of wrapper for MVA MET code to apply MVA corrections for pfMet.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-29 created
 */

#pragma once

#include "HHbbTauTau/METPUSubtraction/source/GBRForest.cxx"
#include "HHbbTauTau/METPUSubtraction/source/GBRTree.cxx"
#include "HHbbTauTau/METPUSubtraction/source/mvaMEtUtilities.cc"
#include "HHbbTauTau/METPUSubtraction/source/PFMETAlgorithmMVA.cc"

#include "Candidate.h"

namespace analysis {

class MvaMetProducer {
public:
    MvaMetProducer(double dZcut, const std::string& inputFileNameU, const std::string& inputFileNameDPhi,
                   const std::string& inputFileNameCovU1, const std::string& inputFileNameCovU2)
        : metAlgo(dZcut)
    {
        metAlgo.initialize(inputFileNameU, inputFileNameDPhi, inputFileNameCovU1, inputFileNameCovU2);
    }

    ntuple::MET ComputeMvaMet(const Candidate& signalCandidate, const ntuple::PFCandidateVector& pfCandidates,
                              const ntuple::JetVector& jets, const Vertex& selectedVertex,
                              const ntuple::VertexVector& vertices, const ntuple::TauVector& taus)
    {
        const auto leptonInfo = ComputeLeptonInfo(signalCandidate);
        auto pfCandidateInfo = ComputePFCandidateInfo(pfCandidates, selectedVertex.position);
        const auto vertexInfo = ComputeVertexInfo(vertices);
        const auto jetInfo = ComputeJetInfo(jets, vertexInfo, selectedVertex, leptonInfo, pfCandidateInfo);
        metAlgo.setInput(leptonInfo, jetInfo, pfCandidateInfo, vertexInfo);
        metAlgo.setHasPhotons(false);
        metAlgo.evaluateMVA();
        const TLorentzVector& metMomentum = metAlgo.getMEt();
        const TMatrix& metCov = metAlgo.getMEtCov();
        ntuple::MET mvaMET;
        mvaMET.eta = metMomentum.Eta();
        mvaMET.phi = metMomentum.Phi();
        mvaMET.pt = metMomentum.Pt();
        mvaMET.significanceMatrix = ntuple::SignificanceMatrixToVector(metCov);
        return mvaMET;
    }

private:
    std::vector<mvaMEtUtilities::leptonInfo> ComputeLeptonInfo(const Candidate& signalCandidate)
    {
        std::vector<mvaMEtUtilities::leptonInfo> leptonInfos;
        for(const Candidate& daughter : signalCandidate.finalStateDaughters) {
            mvaMEtUtilities::leptonInfo info;
            info.p4_ = daughter.momentum;
            info.
            leptonInfo.push_back();
            daughter.
        }
        return leptonInfos;
    }

    std::vector<mvaMEtUtilities::pfCandInfo> ComputePFCandidateInfo(const ntuple::PFCandidateVector& pfCandidates,
                                                                    const TVector3& selectedVertex)
    {

    }

    std::vector<TVector3> ComputeVertexInfo(const ntuple::VertexVector& vertices)
    {

    }

    std::vector<mvaMEtUtilities::JetInfo> ComputeJetInfo(const ntuple::JetVector& jets,
                                                         const std::vector<TVector3>& vertices,
                                                         const TVector3& selectedVertex,
                                                         const std::vector<mvaMEtUtilities::leptonInfo>& signalLeptons,
                                                         std::vector<mvaMEtUtilities::pfCandInfo>& pfCandidates)
    {

    }

//    double ComputeChargedFraction(const Candidate& candidate)
//    {
//        if(candidate.type == Candidate::Mu || candidate.type == Candidate::Electron)
//            return 1.0;
//        double ptTotal = 0.0, ptCharged = 0.0;
//        for() {

//        }
//    }

private:
    PFMETAlgorithmMVA metAlgo;
};

} // analysis
