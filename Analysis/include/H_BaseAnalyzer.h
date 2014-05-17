/*!
 * \file H_BaseAnalyzer.h
 * \brief Definition of H_BaseAnalyzer class which is the base class for all H->tautau analyzers.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-07 created
 */

#pragma once

#include "BaseAnalyzer.h"
#include "SyncTree.h"

namespace analysis {

class H_BaseAnalyzer : public BaseAnalyzer {
public:
    H_BaseAnalyzer(const std::string& inputFileName, const std::string& outputFileName,
                 const std::string& _prefix = "none", size_t _maxNumberOfEvents = 0, bool _useMCtruth = false,
                 const std::string& reweightFileName = "none")
        : BaseAnalyzer(inputFileName, outputFileName, _prefix, _maxNumberOfEvents, _useMCtruth, reweightFileName),
          syncTree("syncTree")
    {}


protected:
    virtual analysis::Candidate SelectBackgroundElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                                         root_ext::AnalyzerData& _anaData,
                                                         const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::electronVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event.electrons().at(id);

        cut(true, ">0 ele cand");
        cut(X(pt) > pt, "pt");
        const double eta = std::abs( X(eta) );
        cut(eta < eta_high, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = (ele_vertex - primaryVertex.position).Perp(); // same as dB
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");
        const size_t pt_index = object.pt < ref_pt ? 0 : 1;
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[pt_index][eta_index], "mva");

        return analysis::Candidate(analysis::Candidate::Electron, id, object,object.charge);
    }

    virtual analysis::Candidate SelectBackgroundMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                                     root_ext::AnalyzerData& _anaData,
                                                     const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::muonVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Muon& object = event.muons().at(id);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = (mu_vertex - primaryVertex.position).Perp();
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(isGlobalMuonPromptTight) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon) == isPFMuon, "PF");
        cut(X(nMatchedStations) > nMatched_Stations, "stations");
        cut(X(pixHits) > pixHits, "pix_hits");
        cut(X(trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");
        cut(X(pfRelIso) < pfRelIso, "pFRelIso");

        return analysis::Candidate(analysis::Candidate::Mu, id, object,object.charge);
    }

    analysis::CandidateVector CollectBJets(const Candidate& higgs)
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> analysis::Candidate
            { return SelectBjet(id, _objectSelector, _anaData, _selection_label, higgs); };
        return CollectObjects<analysis::Candidate>("bjets", base_selector, event.jets().size());
    }

    virtual analysis::Candidate SelectBjet(size_t id, cuts::ObjectSelector* objectSelector,
                                           root_ext::AnalyzerData& _anaData, const std::string& selection_label,
                                           const Candidate& higgs)
    {
        using namespace cuts::Htautau_Summer13::btag;
        cuts::Cutter cut(objectSelector);
        const ntuple::Jet& object = event.jets().at(id);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(combinedSecondaryVertexBJetTags) > CSV, "CSV");
        cut(X(passLooseID) == pfLooseID, "pfLooseID");
        const bool passPUlooseID = (object.puIdBits & (1 << ntuple::JetID_MVA::kLoose)) != 0;
        cut(Y(passPUlooseID) == puLooseID, "puLooseId");
        const Candidate bjet(Candidate::Bjet, id, object);
        for (const Candidate& daughter : higgs.finalStateDaughters){
            const double DeltaR = bjet.momentum.DeltaR(daughter.momentum);
            cut(Y(DeltaR) > deltaR_signalObjects, "DR_signalLeptons");
        }

        return bjet;
    }

    analysis::CandidateVector CollectJets()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> analysis::Candidate
            { return SelectJet(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<analysis::Candidate>("jets", base_selector, event.jets().size());
    }

    virtual analysis::Candidate SelectJet(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::jetID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Jet& object = event.jets().at(id);

        cut(true, ">0 jet cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(passLooseID) == pfLooseID, "pfLooseID");
        const bool passPUlooseID = (object.puIdBits & (1 << ntuple::JetID_MVA::kLoose)) != 0;
        cut(Y(passPUlooseID) == puLooseID, "puLooseID");

        return Candidate(analysis::Candidate::Jet, id, object);
    }

    void ApplyTauCorrections(const finalState::TauTau& mcFinalState, const ntuple::MET& metMVA)
    {
        using namespace cuts::Htautau_Summer13::tauCorrections;

        if(!useMCtruth) {
            correctedTaus = event.taus();
            correctedMET = metMVA;
            return;
        }
        correctedTaus.clear();

        TLorentzVector sumCorrectedTaus, sumTaus;
        for(const ntuple::Tau& tau : event.taus()) {
            TLorentzVector momentum;
            momentum.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
            sumTaus += momentum;

            const bool hasMCmatch = FindMatchedParticles(momentum, mcFinalState.taus, deltaR).size() != 0;
            const double scaleFactor = MomentumScaleFactor(hasMCmatch, momentum.Pt(),
                                                         static_cast<ntuple::tau_id::hadronicDecayMode>(tau.decayMode));
            const TLorentzVector correctedMomentum = momentum * scaleFactor;
            ntuple::Tau correctedTau(tau);
            correctedTau.pt = correctedMomentum.Pt();
            correctedTau.eta = correctedMomentum.Eta();
            correctedTau.phi = correctedMomentum.Phi();
            correctedTau.energy = correctedMomentum.E();
            correctedTaus.push_back(correctedTau);
            sumCorrectedTaus += correctedMomentum;
        }

        TLorentzVector met, metCorrected;
        met.SetPtEtaPhiM(metMVA.pt, metMVA.eta, metMVA.phi, 0.);
        metCorrected = met + sumTaus - sumCorrectedTaus;
        correctedMET = metMVA;
        correctedMET.pt = metCorrected.Pt();
        correctedMET.eta = metCorrected.Eta();
        correctedMET.phi = metCorrected.Phi();
    }

    Candidate ApplyCorrections(const Candidate& higgs, const GenParticle* resonance, const size_t njets)
    {
        if (useMCtruth){
            postRecoilMET = ApplyPostRecoilCorrection(correctedMET, higgs.momentum, resonance->momentum, njets);
        }
        else postRecoilMET = correctedMET;
        return CorrectMassBySVfit(higgs, postRecoilMET);
    }

    const analysis::Candidate& SelectSemiLeptonicHiggs(const analysis::CandidateVector& higgses)
    {
        if(!higgses.size())
            throw std::runtime_error("no available higgs candidate to select");
        const auto higgsSelector = [&] (const analysis::Candidate& first, const analysis::Candidate& second) -> bool
        {
            const double first_Pt1 = first.daughters.at(0).momentum.Pt();
            const double first_Pt2 = first.daughters.at(1).momentum.Pt();
            const double first_sumPt = first_Pt1 + first_Pt2;
            const double second_Pt1 = second.daughters.at(0).momentum.Pt();
            const double second_Pt2 = second.daughters.at(1).momentum.Pt();
            const double second_sumPt = second_Pt1 + second_Pt2;

            return first_sumPt < second_sumPt;
        };
        return *std::max_element(higgses.begin(), higgses.end(), higgsSelector);
    }


    bool FindAnalysisFinalState(analysis::finalState::TauTau& final_state)
    {
        static const particles::ParticleCodes resonanceCodes = { particles::Higgs, particles::Z };
        static const particles::ParticleCodes resonanceDecay = { particles::tau, particles::tau };

        genEvent.Initialize(event.genParticles());

        const analysis::GenParticleSet resonances = genEvent.GetParticles(resonanceCodes);
        if (resonances.size() != 1)
            throw std::runtime_error("not one resonance per event");

        final_state.resonance = *resonances.begin();

        analysis::GenParticlePtrVector resonanceDecayProducts;
        if(!analysis::FindDecayProducts(*final_state.resonance, resonanceDecay,resonanceDecayProducts)) {

//            genEvent.PrintChain(final_state.resonance);
//            throw std::runtime_error("Resonance does not decay into 2 taus");
            std::cerr << "event id = " << event.eventId().eventId
                      << " Resonance does not decayed into 2 taus" << std::endl;
            return false;
        }

        final_state.taus = resonanceDecayProducts;

        return true;
    }

    void FillSyncTree(const Candidate& higgs, const Candidate& higgs_corr, const CandidateVector& jets,
                      const CandidateVector& bjets, const VertexVector& vertices, const Candidate& tau)
    {
        static const double default_value = -10000;
        syncTree.run() = event.eventInfo().run;
        syncTree.lumi() = event.eventInfo().lumis;
        syncTree.evt() = event.eventInfo().EventId;
        syncTree.npv() = vertices.size();
        for (unsigned n = 0; n < event.eventInfo().bunchCrossing.size(); ++n){
            if (event.eventInfo().bunchCrossing.at(n) == 0){
                syncTree.npu() = event.eventInfo().trueNInt.at(n); //only in-time PU
            }
        }
        syncTree.puweight() = weight;

        syncTree.mvis() = higgs.momentum.M();
        syncTree.m_sv() = higgs_corr.momentum.M();
        syncTree.pt_sv() = higgs_corr.momentum.Pt();
        syncTree.eta_sv() = higgs_corr.momentum.Eta();
        syncTree.phi_sv() = higgs_corr.momentum.Phi();

        const ntuple::Tau& ntuple_tau = correctedTaus.at(tau.index);
        //const ntuple::Tau& ntuple_tau = event.taus().at(tau.index);

        //tau - not corrected
        syncTree.pt_2() = ntuple_tau.pt;
        syncTree.eta_2() = ntuple_tau.eta;
        syncTree.phi_2() = ntuple_tau.phi;
        TLorentzVector tau_momentum;
        tau_momentum.SetPtEtaPhiE(ntuple_tau.pt,ntuple_tau.eta,ntuple_tau.phi,ntuple_tau.energy);
        syncTree.m_2() = tau_momentum.M();
        syncTree.q_2() = ntuple_tau.charge;
        syncTree.iso_2() = ntuple_tau.byIsolationMVAraw;
        const TVector3 tau_vertex(ntuple_tau.vx, ntuple_tau.vy, ntuple_tau.vz);
        syncTree.d0_2() = (tau_vertex - primaryVertex.position).Perp();
        syncTree.dZ_2() = std::abs(ntuple_tau.vz - primaryVertex.position.Z());

        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = ntuple_tau.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        syncTree.againstElectronMVA3raw_2() = ntuple_tau.againstElectronMVA3raw;
        syncTree.againstElectronMVA3category_2() = ntuple_tau.againstElectronMVA3category;
        syncTree.byIsolationMVA2raw_2() = ntuple_tau.byIsolationMVA2raw;
        syncTree.againstMuonLoose_2() = ntuple_tau.againstMuonLoose;
        syncTree.againstMuonLoose2_2() = ntuple_tau.againstMuonLoose2;
        syncTree.againstMuonMedium2_2() = ntuple_tau.againstMuonMedium2;
        syncTree.againstMuonTight2_2() = ntuple_tau.againstMuonTight2;
        syncTree.againstElectronLooseMVA3_2() = ntuple_tau.againstElectronLooseMVA3;
        syncTree.againstElectronLoose_2() = ntuple_tau.againstElectronLoose;

        syncTree.met() = event.metPF().pt; //raw
        syncTree.metphi() = event.metPF().phi; //raw
        syncTree.mvamet() = postRecoilMET.pt;
        syncTree.mvametphi() = postRecoilMET.phi;
        syncTree.metcov00() = event.metPF().significanceMatrix.at(0);
        syncTree.metcov01() = event.metPF().significanceMatrix.at(1);
        syncTree.metcov10() = event.metPF().significanceMatrix.at(2);
        syncTree.metcov11() = event.metPF().significanceMatrix.at(3);
        syncTree.mvacov00() = postRecoilMET.significanceMatrix.at(0);
        syncTree.mvacov01() = postRecoilMET.significanceMatrix.at(1);
        syncTree.mvacov10() = postRecoilMET.significanceMatrix.at(2);
        syncTree.mvacov11() = postRecoilMET.significanceMatrix.at(3);

        syncTree.njets() = jets.size();
        syncTree.nbtag() = bjets.size();
        if (jets.size() >= 1){
            const ntuple::Jet& ntuple_jet = event.jets().at(jets.at(0).index);
            syncTree.jpt_1() = jets.at(0).momentum.Pt();
            syncTree.jeta_1() = jets.at(0).momentum.Eta();
            syncTree.jphi_1() = jets.at(0).momentum.Phi();
            syncTree.jmva_1() = ntuple_jet.puIdMVA;
        }
        else {
            syncTree.jpt_1() = default_value;
            syncTree.jeta_1() = default_value;
            syncTree.jphi_1() = default_value;
            syncTree.jmva_1() = default_value;
        }
        if (jets.size() >= 2){
            const ntuple::Jet& ntuple_jet = event.jets().at(jets.at(1).index);
            syncTree.jpt_2() = jets.at(1).momentum.Pt();
            syncTree.jeta_2() = jets.at(1).momentum.Eta();
            syncTree.jphi_2() = jets.at(1).momentum.Phi();
            syncTree.jmva_2() = ntuple_jet.puIdMVA;
        }
        else {
            syncTree.jpt_2() = default_value;
            syncTree.jeta_2() = default_value;
            syncTree.jphi_2() = default_value;
            syncTree.jmva_2() = default_value;
        }
        if (bjets.size() >= 1){
            syncTree.bpt() = bjets.at(0).momentum.Pt();
            syncTree.beta() = bjets.at(0).momentum.Eta();
            syncTree.bphi() = bjets.at(0).momentum.Phi();
        }
        else {
            syncTree.bpt() = default_value;
            syncTree.beta() = default_value;
            syncTree.bphi() = default_value;
        }
    }

protected:
    ntuple::SyncTree syncTree;
    ntuple::TauVector correctedTaus;
    ntuple::MET correctedMET;
    ntuple::MET postRecoilMET;
};
} // analysis
