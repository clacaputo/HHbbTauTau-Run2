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
#include "AnalysisTools.h"

namespace analysis {

class MvaMetProducer {
public:
    MvaMetProducer(double dZcut, const std::string& inputFileNameU, const std::string& inputFileNameDPhi,
                   const std::string& inputFileNameCovU1, const std::string& inputFileNameCovU2)
        : metAlgo(dZcut), useType1Correction(true), minCorrJetPt(-1)
    {
        metAlgo.initialize(inputFileNameU, inputFileNameDPhi, inputFileNameCovU1, inputFileNameCovU2);
    }

    ntuple::MET ComputeMvaMet(const Candidate& signalCandidate, const ntuple::PFCandVector& pfCandidates,
                              const ntuple::JetVector& jets, const Vertex& selectedVertex,
                              const std::vector<Vertex>& goodVertices, const ntuple::TauVector& taus)
    {
        const auto leptonInfo = ComputeLeptonInfo(signalCandidate, taus);
        auto pfCandidateInfo = ComputePFCandidateInfo(pfCandidates, selectedVertex.position);
        const auto vertexInfo = ComputeVertexInfo(goodVertices);
        const auto jetInfo = ComputeJetInfo(jets, leptonInfo, pfCandidateInfo);
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
    static double DefaultDeltaZ() { return -999.; }

    std::vector<mvaMEtUtilities::leptonInfo> ComputeLeptonInfo(const Candidate& signalCandidate,
                                                               const ntuple::TauVector& taus)
    {
        std::vector<mvaMEtUtilities::leptonInfo> leptonInfos;
        for(const Candidate& daughter : signalCandidate.finalStateDaughters) {
            mvaMEtUtilities::leptonInfo info;
            info.p4_ = daughter.momentum;
            info.chargedFrac_ = ComputeChargedFraction(daughter, taus);
            leptonInfos.push_back(info);
        }
        return leptonInfos;
    }

    std::vector<mvaMEtUtilities::pfCandInfo> ComputePFCandidateInfo(const ntuple::PFCandVector& pfCandidates,
                                                                    const TVector3& selectedVertex)
    {
        std::vector<mvaMEtUtilities::pfCandInfo> candInfos;
        for(const ntuple::PFCand& candidate : pfCandidates) {
            mvaMEtUtilities::pfCandInfo info;
            info.p4_.SetPtEtaPhiM(candidate.pt, candidate.eta, candidate.phi, candidate.mass);
            if(candidate.haveTrackInfo) {
                const TVector3 trkV(candidate.trk_vx, candidate.trk_vy, candidate.trk_vz);
                TVector3 trkP;
                trkP.SetPtEtaPhi(candidate.trk_pt, candidate.trk_eta, candidate.trk_phi);
                info.dZ_ = Calculate_dz(trkV, selectedVertex, trkP);
            }
            else
                info.dZ_ = DefaultDeltaZ();
            candInfos.push_back(info);
        }
        return candInfos;
    }

    std::vector<TVector3> ComputeVertexInfo(const std::vector<Vertex>& goodVertices)
    {
        std::vector<TVector3> vertexInfos;
        for(const Vertex& vertex : goodVertices)
            vertexInfos.push_back(vertex.position);
        return vertexInfos;
    }

    std::vector<mvaMEtUtilities::JetInfo> ComputeJetInfo(const ntuple::JetVector& jets,
                                                         const std::vector<mvaMEtUtilities::leptonInfo>& signalLeptons,
                                                         std::vector<mvaMEtUtilities::pfCandInfo>& pfCandidates)
    {
        static const double MinDeltaRtoSignalObjects = 0.5;
        static const double MinJetPtForPFcandCreation = 10.0;
        static const double MaxJetEtaForNeutralEnFrac = 2.5;

        std::vector<mvaMEtUtilities::JetInfo> jetInfos;
        for(const ntuple::Jet& jet : jets) {
            if(!jet.passLooseID) continue;
            mvaMEtUtilities::JetInfo jetInfo;
            jetInfo.p4_.SetPtEtaPhiE(jet.pt, jet.eta, jet.phi, jet.energy);
            double lType1Corr = 0;
            if(useType1Correction) {
                double pCorr = jet.correction;
                lType1Corr = jet.pt - pCorr * jet.pt_raw;
                TLorentzVector pType1Corr;
                pType1Corr.SetPtEtaPhiM(lType1Corr, 0, jet.phi, 0);
                bool pOnLepton = false;
                for(const mvaMEtUtilities::leptonInfo& signal : signalLeptons) {
                    if(signal.p4_.DeltaR(jetInfo.p4_) < MinDeltaRtoSignalObjects) pOnLepton = true;
                }
                if(jet.pt > MinJetPtForPFcandCreation && !pOnLepton) {
                    mvaMEtUtilities::pfCandInfo candInfo;
                    candInfo.p4_ = pType1Corr;
                    candInfo.dZ_ = DefaultDeltaZ();
                    pfCandidates.push_back(candInfo);
                }
                lType1Corr = pCorr*jet.pt_raw - jet.pt_raw;
                lType1Corr /= jet.pt;
            }
            if(jetInfo.p4_.Pt() <= minCorrJetPt) continue;
            jetInfo.mva_ = jet.puIdMVA_met;
            jetInfo.neutralEnFrac_ = jet.neutralEmEnergyFraction + jet.neutralHadronEnergyFraction;
            if(jet.eta > MaxJetEtaForNeutralEnFrac) jetInfo.neutralEnFrac_ = 1.0;
            if(useType1Correction) jetInfo.neutralEnFrac_ -= lType1Corr*jetInfo.neutralEnFrac_;
            jetInfos.push_back(jetInfo);
        }
        return jetInfos;
    }

    double ComputeChargedFraction(const Candidate& candidate, const ntuple::TauVector& taus)
    {
        if(candidate.type == Candidate::Mu || candidate.type == Candidate::Electron)
            return 1.0;
        if(candidate.type != Candidate::Tau)
            throw std::runtime_error("Unsupported candidate type to compute charged fraction.");
        const ntuple::Tau& tau = taus.at(candidate.index);
        double ptTotal = 0.0, ptCharged = 0.0;
        for(double pt : tau.signalChHadCand_Pt) {
            ptCharged += pt;
            ptTotal += pt;
        }
        for(double pt : tau.signalNeutrHadCand_Pt)
            ptTotal += pt;
        for(double pt : tau.signalGammaCand_Pt)
            ptTotal += pt;
        if(ptTotal) ptTotal = 1.0;
        return ptCharged / ptTotal;
    }

private:
    PFMETAlgorithmMVA metAlgo;
    bool useType1Correction;
    double minCorrJetPt;
};

} // analysis
