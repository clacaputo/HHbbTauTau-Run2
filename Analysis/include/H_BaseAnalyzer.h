/*!
 * \file H_BaseAnalyzer.h
 * \brief Definition of H_BaseAnalyzer class which is the base class for all H->tautau analyzers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-05-07 created
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

#include "AnalysisBase/include/SyncTree.h"

#include "BaseAnalyzer.h"

namespace analysis {

class H_BaseAnalyzer : public BaseAnalyzer {
public:
    H_BaseAnalyzer(const std::string& inputFileName, const std::string& outputFileName,
                   const std::string& configFileName, const std::string& _prefix = "none",
                   size_t _maxNumberOfEvents = 0)
        : BaseAnalyzer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
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

    virtual void CalculateFakeWeights(const Candidate& candidate) override
    {
        fakeWeights.clear();
        fakeWeights.push_back(1);
        fakeWeights.push_back(1);
    }

    virtual analysis::Candidate SelectBackgroundElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                                         root_ext::AnalyzerData& _anaData,
                                                         const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::electronVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event->electrons().at(id);
        const analysis::Candidate electron(analysis::Candidate::Electron, id, object);

        cut(true, ">0 ele cand");
        cut(X(pt) > pt, "pt");
        const double eta = std::abs( X(eta) );
        cut(eta < eta_high, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(ele_vertex,primaryVertex.position,electron.momentum); // same as dB
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");
        const size_t pt_index = object.pt < ref_pt ? 0 : 1;
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[pt_index][eta_index], "mva");
        cut(X(missingHits) < missingHits, "missingHits");
        cut(X(hasMatchedConversion) == hasMatchedConversion, "conversion");

        return electron;
    }

    virtual analysis::Candidate SelectBackgroundMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                                     root_ext::AnalyzerData& _anaData,
                                                     const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::muonVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Muon& object = event->muons().at(id);
        const analysis::Candidate muon(analysis::Candidate::Mu, id, object);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(mu_vertex,primaryVertex.position,muon.momentum);
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(isGlobalMuonPromptTight) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon) == isPFMuon, "PF");
        cut(X(nMatchedStations) > nMatched_Stations, "stations");
        cut(X(pixHits) > pixHits, "pix_hits");
        cut(X(trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");
        cut(X(pfRelIso) < pfRelIso, "pFRelIso");

        return muon;
    }

    CandidateVector CollectBJets(const CandidateVector& looseJets, bool doReTag)
    {
        using namespace cuts::Htautau_Summer13::btag;
        analysis::CandidateVector bjets;
        for(const Candidate& looseJetCandidate : looseJets) {
            const ntuple::Jet& looseJet = event->jets().at(looseJetCandidate.index);
            if(looseJet.pt <= pt || std::abs(looseJet.eta) >= eta) continue;
            if(doReTag && !analysis::btag::ReTag(looseJet, btag::payload::EPS13, btag::tagger::CSVM, 0, 0, CSV))
                continue;
            else if(!doReTag && looseJet.combinedSecondaryVertexBJetTags <= CSV)
                continue;

            bjets.push_back(looseJetCandidate);
        }

        const auto bjetsSelector = [&] (const analysis::Candidate& first, const analysis::Candidate& second) -> bool
        {
            const ntuple::Jet& first_bjet = event->jets().at(first.index);
            const ntuple::Jet& second_bjet = event->jets().at(second.index);

            return first_bjet.combinedSecondaryVertexBJetTags > second_bjet.combinedSecondaryVertexBJetTags;
        };

        std::sort(bjets.begin(), bjets.end(), bjetsSelector);
        return bjets;
    }

    CandidateVector CollectLooseJets()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> Candidate
            { return SelectLooseJet(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<Candidate>("jets", base_selector, event->jets().size());
    }

    virtual Candidate SelectLooseJet(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::jetID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Jet& object = event->jets().at(id);

        cut(true, ">0 jet cand");
        cut(X(pt) > pt_loose, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const bool passLooseID = analysis::passPFLooseId(object);
        cut(Y(passLooseID) == pfLooseID, "pfLooseID");
        const bool passPUlooseID = ntuple::JetID_MVA::PassLooseId(object.puIdBits);
        cut(Y(passPUlooseID) == puLooseID, "puLooseID");

        return Candidate(analysis::Candidate::Jet, id, object);
    }

    CandidateVector CollectJets(const CandidateVector& looseJets)
    {
        using namespace cuts::Htautau_Summer13::jetID;
        analysis::CandidateVector jets;
        for(const Candidate& looseJet : looseJets) {
            if(looseJet.momentum.Pt() > pt)
                jets.push_back(looseJet);
        }
        return jets;
    }

    void ApplyTauCorrections(const GenParticlePtrVector& hadronic_taus, bool useLegacyCorrections)
    {
        using namespace cuts::Htautau_Summer13::tauCorrections;

        correctedTaus.clear();

        for(const ntuple::Tau& tau : event->taus()) {
            TLorentzVector momentum;
            momentum.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass);

            const bool hasMCmatch = FindMatchedParticles(momentum, hadronic_taus, deltaR).size() != 0;
            const double scaleFactor = MomentumScaleFactor(hasMCmatch, momentum.Pt(),
                                   ntuple::tau_id::ConvertToHadronicDecayMode(tau.decayMode), useLegacyCorrections);
            const TLorentzVector correctedMomentum = momentum * scaleFactor;
            ntuple::Tau correctedTau(tau);
            correctedTau.pt = correctedMomentum.Pt();
            correctedTau.eta = correctedMomentum.Eta();
            correctedTau.phi = correctedMomentum.Phi();
            correctedTau.mass = correctedMomentum.M();
            correctedTaus.push_back(correctedTau);
        }
    }

    void ApplyTauCorrectionsToMVAMET(const ntuple::MET& metMVA)
    {

        TLorentzVector sumCorrectedTaus, sumTaus;
        for(const ntuple::Tau& tau : event->taus()) {
            TLorentzVector momentum;
            momentum.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass);
            sumTaus += momentum;
        }

        for (const ntuple::Tau& correctedTau : correctedTaus) {
            TLorentzVector correctedMomentum;
            correctedMomentum.SetPtEtaPhiM(correctedTau.pt, correctedTau.eta, correctedTau.phi,correctedTau.mass);
            sumCorrectedTaus += correctedMomentum;
        }

        TLorentzVector met, metCorrected;
        met.SetPtEtaPhiM(metMVA.pt, 0, metMVA.phi, 0.);
        metCorrected = met + sumTaus - sumCorrectedTaus;
        correctedMET = metMVA;
        correctedMET.pt = metCorrected.Pt();
        correctedMET.phi = metCorrected.Phi();

    }

    analysis::CandidateVector ApplyTriggerMatch(const analysis::CandidateVector& higgses,
                                                        const std::vector<std::string>& hltPaths,
                                                        bool useStandardTriggerMatch)
    {
        analysis::CandidateVector triggeredHiggses;
        for (const auto& higgs : higgses){
            if(!useStandardTriggerMatch && analysis::HaveTriggerMatched(event->triggerObjects(), hltPaths, higgs,
                                                                        cuts::Htautau_Summer13::DeltaR_triggerMatch))
                triggeredHiggses.push_back(higgs);
            if (useStandardTriggerMatch && analysis::HaveTriggerMatched(*event,hltPaths,higgs))
                triggeredHiggses.push_back(higgs);
        }
        return triggeredHiggses;

    }

    void ApplyRecoilCorrections(const Candidate& higgs, const GenParticle* resonance, const size_t njets)
    {
        //for W - not implemented
        if (config.ApplyRecoilCorrection()){
            if(!resonance)
                resonance = FindWboson();
            if(resonance)
                postRecoilMET =
                        recoilCorrectionProducer->ApplyCorrection(correctedMET, higgs.momentum, resonance->momentum, njets);
            else
                postRecoilMET = correctedMET;
        }
        else
            postRecoilMET = correctedMET;

    }

    const GenParticle* FindWboson()
    {
        static const particles::ParticleCodes resonanceCodes = { particles::W_plus };
        static const particles::ParticleCodes resonanceDecay = { particles::tau, particles::nu_tau };

//        genEvent.Print();
        const analysis::GenParticleSet all_resonances = genEvent.GetParticles(resonanceCodes);

        analysis::GenParticleSet resonances;
        for(const GenParticle* w : all_resonances) {
            if(w->mothers.size() == 1 && w->mothers.at(0)->pdg.Code != particles::tau)
                resonances.insert(w);
        }

        if (resonances.size() == 0) return nullptr;

        if (resonances.size() > 2)
            throw exception("more than 2 W in the event");

        if (resonances.size() == 1){
            const GenParticle* Wboson = *resonances.begin();

            analysis::GenParticlePtrVector resonanceDecayProducts;
            if(analysis::FindDecayProducts(*Wboson, resonanceDecay,resonanceDecayProducts)){
                return Wboson;
            }
            return nullptr;
        }


        const GenParticle* first_resonance = *resonances.begin();
        const GenParticle* second_resonance = *std::next(resonances.begin());
        const GenParticle* W_mother;
        const GenParticle* W_daughter;


        if (first_resonance->daughters.size() == 1 && first_resonance->daughters.at(0) == second_resonance){
            W_mother = first_resonance;
            W_daughter = second_resonance;
        }
        else if (second_resonance->daughters.size() == 1 && second_resonance->daughters.at(0) == first_resonance){
            W_mother = second_resonance;
            W_daughter = first_resonance;
        }
        else
            throw exception("two resonances are not in relationship");

        analysis::GenParticlePtrVector resonanceDecayProducts;
        if(analysis::FindDecayProducts(*W_mother, resonanceDecay,resonanceDecayProducts)){
            return W_daughter;
        }

        return nullptr;
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

        genEvent.Initialize(event->genParticles());
//        genEvent.Print();

        const analysis::GenParticleSet resonances = genEvent.GetParticles(resonanceCodes);

        analysis::GenParticlePtrVector resonancesToTauTau;

        for (const GenParticle* resonance : resonances){
            analysis::GenParticlePtrVector resonanceDecayProducts;
            if(analysis::FindDecayProducts(*resonance, resonanceDecay,resonanceDecayProducts)){
                final_state.taus = resonanceDecayProducts;
                final_state.resonance = resonance;
                resonancesToTauTau.push_back(resonance);
            }
        }

        if (resonancesToTauTau.size() > 1)
            throw exception("more than one resonance to tautau per event");

        if(resonancesToTauTau.size() == 0){
            if(config.ExpectedAtLeastOneSMResonanceToTauTauOrToBB()) /* ExpectedOneResonanceToTauTau */
                throw exception("resonance to tautau not found");
            else
                return false;
        }

        for (const analysis::GenParticle* tau_MC : final_state.taus) {
            analysis::GenParticlePtrVector TauProducts;
            if (!analysis::FindDecayProducts(*tau_MC, analysis::TauMuonicDecay, TauProducts)
                    && !analysis::FindDecayProducts(*tau_MC, analysis::TauElectronDecay, TauProducts))
                final_state.hadronic_taus.push_back(tau_MC);
        }
        return true;
    }

    void FillSyncTree(const Candidate& higgs, double m_sv,
                      const CandidateVector& jets, const CandidateVector& jetsPt20,
                      const CandidateVector& bjets, const analysis::CandidateVector& retagged_bjets,
                      const VertexVector& vertices, const Candidate& leg1, const Candidate& leg2,
                      const ntuple::MET& pfMET)
    {
        static const double default_value = ntuple::DefaultFillValueForSyncTree();
        syncTree.run() = event->eventInfo().run;
        syncTree.lumi() = event->eventInfo().lumis;
        syncTree.evt() = event->eventInfo().EventId;

        syncTree.npv() = vertices.size();
        if (config.ApplyPUreweight()){
            const size_t bxIndex = tools::find_index(event->eventInfo().bunchCrossing, 0);
            if(bxIndex >= event->eventInfo().bunchCrossing.size())
                throw std::runtime_error("in-time BX not found");
            syncTree.npu() = event->eventInfo().trueNInt.at(bxIndex);
        }
        //syncTree.rho();

        //syncTree.mcweight();
        syncTree.puweight() = PUweight;
        syncTree.trigweight_1() = triggerWeights.at(0);
        syncTree.trigweight_2() = triggerWeights.at(1);
        syncTree.idweight_1() = IDweights.at(0);
        syncTree.idweight_2() = IDweights.at(1);
        syncTree.isoweight_1() = IsoWeights.at(0);
        syncTree.isoweight_2() = IsoWeights.at(1);
        syncTree.fakeweight() = fakeWeights.at(1);
        //syncTree.effweight();
        syncTree.weight() = eventWeight;
        //syncTree.embeddedWeight();
        //syncTree.signalWeight();

        syncTree.mvis() = higgs.momentum.M();
        syncTree.m_sv() = m_sv;
        syncTree.pt_sv() = default_value;
        syncTree.eta_sv() = default_value;
        syncTree.phi_sv() = default_value;
        syncTree.m_sv_Up() = default_value;
        syncTree.m_sv_Down() = default_value;

        syncTree.pt_1() = leg1.momentum.Pt();
        syncTree.phi_1() = leg1.momentum.Phi();
        syncTree.eta_1() = leg1.momentum.Eta();
        syncTree.m_1() = leg1.momentum.M();
        syncTree.q_1() = leg1.charge;
        syncTree.mt_1() = analysis::Calculate_MT(leg1.momentum, postRecoilMET.pt, postRecoilMET.phi);
        syncTree.d0_1() = analysis::Calculate_dxy(leg1.vertexPosition, primaryVertex.position,leg1.momentum);
        syncTree.dZ_1() = leg1.vertexPosition.Z() - primaryVertex.position.Z();

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
        syncTree.mt_2() = analysis::Calculate_MT(leg2.momentum, postRecoilMET.pt, postRecoilMET.phi);
        syncTree.d0_2() = analysis::Calculate_dxy(leg2.vertexPosition, primaryVertex.position,leg2.momentum);
        syncTree.dZ_2() = leg2.vertexPosition.Z() - primaryVertex.position.Z();

        const ntuple::Tau& ntuple_tau_leg2 = correctedTaus.at(leg2.index);
        syncTree.iso_2() = ntuple_tau_leg2.byIsolationMVAraw;
        //syncTree.passid_2();
        //syncTree.passiso_2();
        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = ntuple_tau_leg2.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        syncTree.againstElectronMVA3raw_2() = ntuple_tau_leg2.againstElectronMVA3raw;
        syncTree.byIsolationMVA2raw_2() = ntuple_tau_leg2.byIsolationMVA2raw;
        syncTree.againstMuonLoose2_2() = ntuple_tau_leg2.againstMuonLoose2;
        syncTree.againstMuonMedium2_2() = ntuple_tau_leg2.againstMuonMedium2;
        syncTree.againstMuonTight2_2() = ntuple_tau_leg2.againstMuonTight2;

        TLorentzVector postRecoilMetMomentum;
        postRecoilMetMomentum.SetPtEtaPhiM(postRecoilMET.pt, 0, postRecoilMET.phi, 0.);
        syncTree.pt_tt() = (leg1.momentum + leg2.momentum + postRecoilMetMomentum).Pt();

        syncTree.met() = pfMET.pt;
        syncTree.metphi() = pfMET.phi;
        syncTree.mvamet() = postRecoilMET.pt;
        syncTree.mvametphi() = postRecoilMET.phi;
        //syncTree.pzetavis();
        //syncTree.pzetamiss();
        if(pfMET.significanceMatrix.size()) {
            const TMatrixD metPFcov = ntuple::VectorToSignificanceMatrix(pfMET.significanceMatrix);
            syncTree.metcov00() = metPFcov[0][0];
            syncTree.metcov01() = metPFcov[0][1];
            syncTree.metcov10() = metPFcov[1][0];
            syncTree.metcov11() = metPFcov[1][1];
        }
        const TMatrixD metMVAcov = ntuple::VectorToSignificanceMatrix(postRecoilMET.significanceMatrix);
        syncTree.mvacov00() = metMVAcov[0][0];
        syncTree.mvacov01() = metMVAcov[0][1];
        syncTree.mvacov10() = metMVAcov[1][0];
        syncTree.mvacov11() = metMVAcov[1][1];

        syncTree.njets() = jets.size();
        syncTree.njetspt20() = jetsPt20.size();

        if (jets.size() >= 1) {
            const Candidate& jet = jets.at(0);
            const ntuple::Jet& ntuple_jet = event->jets().at(jet.index);
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
            //syncTree.jptunc_1() = default_value;
            syncTree.jmva_1() = default_value;
//            syncTree.jlrm_1() = default_value;
//            syncTree.jctm_1() = default_value;
            syncTree.jpass_1() = default_value;
        }

        if (jets.size() >= 2) {
            const Candidate& jet = jets.at(1);
            const ntuple::Jet& ntuple_jet = event->jets().at(jet.index);
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
//            syncTree.jptunc_2() = default_value;
            syncTree.jmva_2() = default_value;
//            syncTree.jlrm_2() = default_value;
//            syncTree.jctm_2() = default_value;
            syncTree.jpass_2() = default_value;
        }

        syncTree.nbtag() = retagged_bjets.size();

        if (bjets.size() >= 1) {
            const Candidate& bjet = bjets.at(0);
            const ntuple::Jet& ntuple_bjet = event->jets().at(bjet.index);
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
            const Candidate& lead_bjet = bjets.at(0);
            const Candidate& subLead_bjet = bjets.at(1);
            const ntuple::Jet& ntuple_bjet = event->jets().at(subLead_bjet.index);
            syncTree.bpt_2() = subLead_bjet.momentum.Pt();
            syncTree.beta_2() = subLead_bjet.momentum.Eta();
            syncTree.bphi_2() = subLead_bjet.momentum.Phi();
            syncTree.bcsv_2() = ntuple_bjet.combinedSecondaryVertexBJetTags;
            const Candidate higgs_bb(Candidate::Higgs,lead_bjet,subLead_bjet);
            const Candidate resonance(Candidate::Resonance,higgs,higgs_bb);
            syncTree.m_bb() = higgs_bb.momentum.M();
            syncTree.m_ttbb() = (resonance.momentum + postRecoilMetMomentum).M();
        } else {
            syncTree.bpt_2() = default_value;
            syncTree.beta_2() = default_value;
            syncTree.bphi_2() = default_value;
            syncTree.bcsv_2() = default_value;
            syncTree.m_bb() = default_value;
            syncTree.m_ttbb() = default_value;
        }

        if (bjets.size() >= 3){
            const Candidate& bjet = bjets.at(2);
            const ntuple::Jet& ntuple_bjet = event->jets().at(bjet.index);
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
