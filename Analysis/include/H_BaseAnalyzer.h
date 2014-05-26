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

    virtual void CalculateTriggerWeights(const Candidate& candidate) override
    {
        triggerWeights.clear();
        triggerWeights.push_back(1);
        triggerWeights.push_back(1);
    }

    virtual void CalculateIsoWeights(const Candidate& candidate) override
    {
        IsoWeights.clear();
        IsoWeights.push_back(1);
        IsoWeights.push_back(1);
    }

    virtual void CalculateIdWeights(const Candidate& candidate) override
    {
        IDweights.clear();
        IDweights.push_back(1);
        IDweights.push_back(1);
    }

    virtual void CalculateDMWeights(const Candidate& candidate) override
    {
        DMweights.clear();
        DMweights.push_back(1);
        DMweights.push_back(1);
    }

    virtual analysis::Candidate SelectBackgroundElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                                         root_ext::AnalyzerData& _anaData,
                                                         const std::string& selection_label) override
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
        cut(X(missingHits) < missingHits, "missingHits");
        cut(X(hasMatchedConversion) == hasMatchedConversion, "conversion");

        return analysis::Candidate(analysis::Candidate::Electron, id, object);
    }

    virtual analysis::Candidate SelectBackgroundMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                                     root_ext::AnalyzerData& _anaData,
                                                     const std::string& selection_label) override
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

        return analysis::Candidate(analysis::Candidate::Mu, id, object);
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

    void ApplyTauCorrections(const finalState::TauTau& mcFinalState, const ntuple::MET& metMVA,
                             bool useLegacyCorrections)
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
                                   ntuple::tau_id::ConvertToHadronicDecayMode(tau.decayMode), useLegacyCorrections);
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

    analysis::CandidateVector ApplyTriggerMatch(const analysis::CandidateVector& higgses,
                                                        const std::vector<std::string>& hltPaths,
                                                        bool useStandardTriggerMatch)
    {
        analysis::CandidateVector triggeredHiggses;
        for (const auto& higgs : higgses){
            if(!useStandardTriggerMatch && analysis::HaveTriggerMatched(event.triggerObjects(), hltPaths, higgs))
                triggeredHiggses.push_back(higgs);
            if (useStandardTriggerMatch && analysis::HaveTriggerMatched(event,hltPaths,higgs))
                triggeredHiggses.push_back(higgs);
        }
        return triggeredHiggses;

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
                      const CandidateVector& bjets, const VertexVector& vertices, const Candidate& leg1,
                      const Candidate& leg2)
    {
        static const double default_value = ntuple::DefaultFillValueForSyncTree();
        syncTree.run() = event.eventInfo().run;
        syncTree.lumi() = event.eventInfo().lumis;
        syncTree.evt() = event.eventInfo().EventId;

        syncTree.npv() = vertices.size();
        const size_t bxIndex = tools::find_index(event.eventInfo().bunchCrossing, 0);
        if(bxIndex >= event.eventInfo().bunchCrossing.size())
            throw std::runtime_error("in-time BX not found");
        syncTree.npu() = event.eventInfo().trueNInt.at(bxIndex);
        //syncTree.rho();

        //syncTree.mcweight();
        syncTree.puweight() = PUweight;
        syncTree.trigweight_1() = triggerWeights.at(0);
        syncTree.trigweight_2() = triggerWeights.at(1);
        syncTree.idweight_1() = IDweights.at(0);
        syncTree.idweight_2() = IDweights.at(1);
        syncTree.isoweight_1() = IsoWeights.at(0);
        syncTree.isoweight_2() = IsoWeights.at(1);
        //syncTree.fakeweight();
        //syncTree.effweight();
        syncTree.weight() = eventWeight;
        //syncTree.embeddedWeight();
        //syncTree.signalWeight();

        syncTree.mvis() = higgs.momentum.M();
        syncTree.m_sv() = higgs_corr.momentum.M();
        syncTree.pt_sv() = higgs_corr.momentum.Pt();
        syncTree.eta_sv() = higgs_corr.momentum.Eta();
        syncTree.phi_sv() = higgs_corr.momentum.Phi();
        //syncTree.m_sv_Up();
        //syncTree.m_sv_Down();

        syncTree.pt_1() = leg1.momentum.Pt();
        syncTree.phi_1() = leg1.momentum.Phi();
        syncTree.eta_1() = leg1.momentum.Eta();
        syncTree.m_1() = leg1.momentum.M();
        syncTree.q_1() = leg1.charge;
        syncTree.mt_1() = analysis::Calculate_MT(leg1.momentum, correctedMET.pt, correctedMET.phi);
        syncTree.d0_1() = (leg1.vertexPosition - primaryVertex.position).Perp();
        syncTree.dZ_1() = std::abs(leg1.vertexPosition.Z() - primaryVertex.position.Z());

        // leg1 lepton specific variable should be filled outside. Here all them set to the default value.
        syncTree.iso_1() = default_value;
        syncTree.mva_1() = default_value;
        //syncTree.passid_1() = default_value;
        //syncTree.passiso_1() = default_value;
        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = default_value;
        syncTree.againstElectronMVA3raw_1() = default_value;
        syncTree.byIsolationMVA2raw_1() = default_value;
        syncTree.againstMuonLoose2_1() = default_value;
        syncTree.againstMuonMedium2_1() = default_value;
        syncTree.againstMuonTight2_1() = default_value;

        syncTree.pt_2() = leg2.momentum.Pt();
        syncTree.phi_2() = leg2.momentum.Phi();
        syncTree.eta_2() = leg2.momentum.Eta();
        syncTree.m_2() = leg2.momentum.M();
        syncTree.q_2() = leg2.charge;
        syncTree.mt_2() = analysis::Calculate_MT(leg2.momentum, correctedMET.pt, correctedMET.phi);
        syncTree.d0_2() = (leg2.vertexPosition - primaryVertex.position).Perp();
        syncTree.dZ_2() = std::abs(leg2.vertexPosition.Z() - primaryVertex.position.Z());

        const ntuple::Tau& ntuple_tau_leg2 = correctedTaus.at(leg2.index);
        syncTree.iso_2() = ntuple_tau_leg2.byIsolationMVAraw;
        syncTree.mva_2() = ntuple_tau_leg2.againstElectronMVA3raw;
        //syncTree.passid_2();
        //syncTree.passiso_2();
        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = ntuple_tau_leg2.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        syncTree.againstElectronMVA3raw_2() = ntuple_tau_leg2.againstElectronMVA3raw;
        syncTree.byIsolationMVA2raw_2() = ntuple_tau_leg2.byIsolationMVA2raw;
        syncTree.againstMuonLoose2_2() = ntuple_tau_leg2.againstMuonLoose2;
        syncTree.againstMuonMedium2_2() = ntuple_tau_leg2.againstMuonMedium2;
        syncTree.againstMuonTight2_2() = ntuple_tau_leg2.againstMuonTight2;

        syncTree.pt_tt() = (leg1.momentum + leg2.momentum).Pt();

        syncTree.met() = event.metPF().pt;
        syncTree.metphi() = event.metPF().phi;
        syncTree.mvamet() = postRecoilMET.pt;
        syncTree.mvametphi() = postRecoilMET.phi;
        //syncTree.pzetavis();
        //syncTree.pzetamiss();
        const TMatrixD metPFcov = ntuple::VectorToSignificanceMatrix(event.metPF().significanceMatrix);
        syncTree.metcov00() = metPFcov[0][0];
        syncTree.metcov01() = metPFcov[0][1];
        syncTree.metcov10() = metPFcov[1][0];
        syncTree.metcov11() = metPFcov[1][1];
        const TMatrixD metMVAcov = ntuple::VectorToSignificanceMatrix(postRecoilMET.significanceMatrix);
        syncTree.mvacov00() = metMVAcov[0][0];
        syncTree.mvacov01() = metMVAcov[0][1];
        syncTree.mvacov10() = metMVAcov[1][0];
        syncTree.mvacov11() = metMVAcov[1][1];

        syncTree.njets() = jets.size();
        //syncTree.njetspt20();

        if (jets.size() >= 1) {
            const Candidate& jet = jets.at(0);
            const ntuple::Jet& ntuple_jet = event.jets().at(jet.index);
            syncTree.jpt_1() = jet.momentum.Pt();
            syncTree.jeta_1() = jet.momentum.Eta();
            syncTree.jphi_1() = jet.momentum.Phi();
            syncTree.jptraw_1() = ntuple_jet.pt_raw;
            //syncTree.jptunc_1();
            syncTree.jmva_1() = ntuple_jet.puIdMVA;
            //syncTree.jlrm_1();
            //syncTree.jctm_1();
            syncTree.jpass_1() = ntuple::JetID_MVA::PassLooseId(ntuple_jet.puIdBits);
        } else {
            syncTree.jpt_1() = default_value;
            syncTree.jeta_1() = default_value;
            syncTree.jphi_1() = default_value;
            syncTree.jptraw_1() = default_value;
            syncTree.jptunc_1() = default_value;
            syncTree.jmva_1() = default_value;
            syncTree.jlrm_1() = default_value;
            syncTree.jctm_1() = default_value;
            syncTree.jpass_1() = default_value;
        }

        if (jets.size() >= 2) {
            const Candidate& jet = jets.at(1);
            const ntuple::Jet& ntuple_jet = event.jets().at(jet.index);
            syncTree.jpt_2() = jet.momentum.Pt();
            syncTree.jeta_2() = jet.momentum.Eta();
            syncTree.jphi_2() = jet.momentum.Phi();
            syncTree.jptraw_2() = ntuple_jet.pt_raw;
            //syncTree.jptunc_2();
            syncTree.jmva_2() = ntuple_jet.puIdMVA;
            //syncTree.jlrm_2();
            //syncTree.jctm_2();
            syncTree.jpass_2() = ntuple::JetID_MVA::PassLooseId(ntuple_jet.puIdBits);
        } else {
            syncTree.jpt_2() = default_value;
            syncTree.jeta_2() = default_value;
            syncTree.jphi_2() = default_value;
            syncTree.jptraw_2() = default_value;
            syncTree.jptunc_2() = default_value;
            syncTree.jmva_2() = default_value;
            syncTree.jlrm_2() = default_value;
            syncTree.jctm_2() = default_value;
            syncTree.jpass_2() = default_value;
        }

        syncTree.nbtag() = bjets.size();

        if (bjets.size() >= 1) {
            const Candidate& bjet = bjets.at(0);
            const ntuple::Jet& ntuple_bjet = event.jets().at(bjet.index);
            syncTree.bpt_1() = bjet.momentum.Pt();
            syncTree.beta_1() = bjet.momentum.Eta();
            syncTree.bphi_1() = bjet.momentum.Phi();
            syncTree.bcsv_1() = ntuple_bjet.combinedSecondaryVertexBJetTags;
        } else {
            syncTree.bpt_1() = default_value;
            syncTree.beta_1() = default_value;
            syncTree.bphi_1() = default_value;
            syncTree.bcsv_1() = default_value;
        }

        if (bjets.size() >= 2) {
            const Candidate& bjet = bjets.at(1);
            const ntuple::Jet& ntuple_bjet = event.jets().at(bjet.index);
            syncTree.bpt_2() = bjet.momentum.Pt();
            syncTree.beta_2() = bjet.momentum.Eta();
            syncTree.bphi_2() = bjet.momentum.Phi();
            syncTree.bcsv_2() = ntuple_bjet.combinedSecondaryVertexBJetTags;
        } else {
            syncTree.bpt_2() = default_value;
            syncTree.beta_2() = default_value;
            syncTree.bphi_2() = default_value;
            syncTree.bcsv_2() = default_value;
        }

        if (bjets.size() >= 3){
            const Candidate& bjet = bjets.at(2);
            const ntuple::Jet& ntuple_bjet = event.jets().at(bjet.index);
            syncTree.bpt_3() = bjet.momentum.Pt();
            syncTree.beta_3() = bjet.momentum.Eta();
            syncTree.bphi_3() = bjet.momentum.Phi();
            syncTree.bcsv_3() = ntuple_bjet.combinedSecondaryVertexBJetTags;
        } else {
            syncTree.bpt_3() = default_value;
            syncTree.beta_3() = default_value;
            syncTree.bphi_3() = default_value;
            syncTree.bcsv_3() = default_value;
        }
    }

protected:
    ntuple::SyncTree syncTree;
    ntuple::TauVector correctedTaus;
    ntuple::MET correctedMET;
    ntuple::MET postRecoilMET;
};
} // analysis
