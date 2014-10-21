/*!
 * \file BaseFlatTreeProducer.h
 * \brief Definition of BaseFlatTreeProducer class, the base class for flat tree producers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-07-11 created
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

#include "AnalysisBase/include/FlatTree.h"

#include "BaseAnalyzer.h"

namespace analysis {

class BaseFlatTreeProducer : public BaseAnalyzer {
public:
    BaseFlatTreeProducer(const std::string& inputFileName, const std::string& outputFileName,
                         const std::string& configFileName, const std::string& _prefix = "none",
                         size_t _maxNumberOfEvents = 0,
                         std::shared_ptr<ntuple::FlatTree> _flatTree = std::shared_ptr<ntuple::FlatTree>())
        : BaseAnalyzer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          flatTree(_flatTree), writeFlatTree(!flatTree)
    {
        if(!flatTree)
            flatTree = std::shared_ptr<ntuple::FlatTree>(new ntuple::FlatTree("flatTree"));
    }

    virtual ~BaseFlatTreeProducer() override
    {
        if(writeFlatTree)
            flatTree->Write();
    }


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

    CandidateVector CollectBJets(const CandidateVector& looseJets, bool doReTag, bool applyCsvCut)
    {
        using namespace cuts::Htautau_Summer13::btag;
        analysis::CandidateVector bjets;
        for(const Candidate& looseJetCandidate : looseJets) {
            const ntuple::Jet& looseJet = event->jets().at(looseJetCandidate.index);
            if(looseJet.pt <= pt || std::abs(looseJet.eta) >= eta) continue;
            if(doReTag && !analysis::btag::ReTag(looseJet, btag::payload::EPS13, btag::tagger::CSVM, 0, 0, CSV))
                continue;
            else if(!doReTag && applyCsvCut && looseJet.combinedSecondaryVertexBJetTags <= CSV)
                continue;

            bjets.push_back(looseJetCandidate);
        }

        const auto bjetsSelector = [&] (const Candidate& first, const Candidate& second) -> bool
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
        const bool passLooseID = passPFLooseId(object);
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

    bool FindAnalysisFinalState(finalState::bbTauTau& final_state)
    {
        static const particles::ParticleCodes resonanceCodes = { particles::MSSM_H, particles::MSSM_A };
        static const particles::ParticleCodes resonanceDecay_1 = { particles::Higgs, particles::Higgs };
        static const particles::ParticleCodes resonanceDecay_2 = { particles::Z, particles::Higgs };
        static const particles::ParticleCodes SM_ResonanceCodes = { particles::Higgs, particles::Z,
                                                                    particles::MSSM_H, particles::MSSM_A };
        static const particles::ParticleCodes SM_ResonanceDecay_1 = { particles::tau, particles::tau };
        static const particles::ParticleCodes SM_ResonanceDecay_2 = { particles::b, particles::b };
        static const particles::ParticleCodes2D HiggsDecays = { SM_ResonanceDecay_1, SM_ResonanceDecay_2 };

        genEvent.Initialize(event->genParticles());
        final_state.Reset();

        const GenParticleSet resonances = genEvent.GetParticles(resonanceCodes);

        if (resonances.size() > 1)
            throw exception("more than 1 resonance per event");

        if (resonances.size() == 1) {
            final_state.resonance = *resonances.begin();

            bool doubleHiggsSignal = true;
            GenParticlePtrVector HiggsBosons;
            if(!FindDecayProducts(*final_state.resonance, resonanceDecay_1, HiggsBosons) ||
                    !FindDecayProducts(*final_state.resonance, resonanceDecay_2, HiggsBosons))
                doubleHiggsSignal = false;

            if(doubleHiggsSignal) {
                GenParticleVector2D HiggsDecayProducts;
                GenParticleIndexVector HiggsIndexes;
                if(!FindDecayProducts2D(HiggsBosons, HiggsDecays, HiggsDecayProducts, HiggsIndexes))
                    throw exception("NOT HH -> bb tautau");

                for(const GenParticle* tau : HiggsDecayProducts.at(0)) {
                    const VisibleGenObject tau_products(tau);
                    final_state.taus.push_back(tau_products);
                    if(tau_products.finalStateChargedHadrons.size() != 0)
                        final_state.hadronic_taus.push_back(tau_products);
                }
                for(const GenParticle* b : HiggsDecayProducts.at(1))
                    final_state.b_jets.push_back(VisibleGenObject(b));

                final_state.Higgs_TauTau = HiggsBosons.at(HiggsIndexes.at(0));
                final_state.Higgs_BB = HiggsBosons.at(HiggsIndexes.at(1));
                return true;
            }
        }

        if(config.ExpectedOneNonSMResonance())
            throw exception("Non-SM resonance not found.");

        //search H->bb, H->tautau
        const GenParticleSet SM_particles = genEvent.GetParticles(SM_ResonanceCodes);

        GenParticlePtrVector SM_ResonanceToTauTau;
        GenParticlePtrVector SM_ResonanceToBB;

        for (const GenParticle* SM_particle : SM_particles){
            GenParticlePtrVector resonanceDecayProducts;
            if(FindDecayProducts(*SM_particle, SM_ResonanceDecay_1,resonanceDecayProducts)){
                for(const GenParticle* tau : resonanceDecayProducts) {
                    const VisibleGenObject tau_products(tau);
                    final_state.taus.push_back(tau_products);
                    if(tau_products.finalStateChargedHadrons.size() != 0)
                        final_state.hadronic_taus.push_back(tau_products);
                }
                final_state.Higgs_TauTau = SM_particle;
                SM_ResonanceToTauTau.push_back(SM_particle);
            }
            else if (FindDecayProducts(*SM_particle, SM_ResonanceDecay_2,resonanceDecayProducts)){
                for(const GenParticle* b : resonanceDecayProducts)
                    final_state.b_jets.push_back(VisibleGenObject(b));
                final_state.Higgs_BB = SM_particle;
                SM_ResonanceToBB.push_back(SM_particle);
            }
        }

        if (SM_ResonanceToTauTau.size() > 1)
            throw exception("more than one SM resonance to tautau per event");
        if (SM_ResonanceToBB.size() > 1)
            throw exception("more than one SM resonance to bb per event");

        if (SM_ResonanceToTauTau.size() + SM_ResonanceToBB.size() == 0) {
            if(config.ExpectedAtLeastOneSMResonanceToTauTauOrToBB())
                throw exception("SM resonance to tautau or to bb not found.");
            return false;
        }

        return true;
    }


    ntuple::EventType DoEventCategorization(const analysis::Candidate& higgs)
    {
        using namespace cuts::Htautau_Summer13::DrellYannCategorization;
        if(!config.DoZEventCategorization())
            return ntuple::EventType::Unknown;

        static const particles::ParticleCodes Zcode = { particles::Z };
        static const particles::ParticleCodes ZDecay_electrons = { particles::e, particles::e };
        static const particles::ParticleCodes ZDecay_muons = { particles::mu, particles::mu };
        static const particles::ParticleCodes ZDecay_taus = { particles::tau, particles::tau };

        const GenParticleSet Zparticles_all = genEvent.GetParticles(Zcode);

        GenParticleSet Zparticles;
        for(const GenParticle* z : Zparticles_all) {
            if(z->mothers.size() == 1) {
                const GenParticle* mother = z->mothers.front();
                if(mother->pdg.Code == particles::Z && mother->status == particles::HardInteractionProduct)
                    Zparticles.insert(z);
            }
         }

        if (Zparticles.size() > 1 || Zparticles.size() == 0)
            throw exception("not 1 Z per event");

        const GenParticle* Z_mc = *Zparticles.begin();
        while(Z_mc->daughters.size() == 1 && Z_mc->daughters.front()->pdg.Code == particles::Z)
            Z_mc = Z_mc->daughters.front();


        GenParticlePtrVector ZProducts;
        const bool z_tt = FindDecayProducts(*Z_mc, ZDecay_taus, ZProducts, true);
        if (!z_tt && !FindDecayProducts(*Z_mc, ZDecay_electrons, ZProducts, true)
                 && !FindDecayProducts(*Z_mc, ZDecay_muons, ZProducts, true))
            throw exception("not leptonic Z decay");

        const auto matchedParticles_first =
                FindMatchedParticles(higgs.daughters.at(0).momentum, ZProducts, deltaR_matchGenParticle);
        const auto matchedParticles_second =
                FindMatchedParticles(higgs.daughters.at(1).momentum, ZProducts, deltaR_matchGenParticle);

        if (matchedParticles_first.size() >= 1 && matchedParticles_second.size() >= 1)
            return z_tt ? ntuple::EventType::ZTT : ntuple::EventType::ZL;
        return z_tt ? ntuple::EventType::ZTT_no_match : ntuple::EventType::ZJ;
    }

    bool GenFilterForZevents(const finalState::bbTauTau& final_state)
    {
        using namespace cuts::Htautau_Summer13::DYEmbedded;
        if (final_state.taus.size() != 2)
            throw exception("not 2 taus in the event at Gen Level");
        const GenParticle* firstTau = final_state.taus.at(0).GetOriginGenParticle();
        const GenParticle* secondTau = final_state.taus.at(1).GetOriginGenParticle();
        if ((firstTau->momentum + secondTau->momentum).M() > invariantMassCut) return true;
        return false;
    }

    bool MatchTausFromHiggsWithGenTaus(const analysis::Candidate& higgs, const finalState::bbTauTau& final_state)
    {
        using namespace cuts::Htautau_Summer13::DYEmbedded;
        std::set<const GenParticle *> all_matched_particles;
        for(const auto& daughter : higgs.daughters) {
            const VisibleGenObjectVector matched_particles =
                    FindMatchedObjects(daughter.momentum, final_state.taus, deltaR_betweenTaus);
            if(!matched_particles.size())
                return false;
            for(const auto& gen_object : matched_particles)
                all_matched_particles.insert(gen_object.origin);
        }
        return all_matched_particles.size() >= higgs.daughters.size();
    }

    analysis::kinematic_fit::FitResultsMap RunKinematicFit(const CandidateVector& bjets, const Candidate& higgs_to_taus,
                                                           const ntuple::MET& met, bool fit_two_bjets,
                                                           bool fit_four_body)
    {
        using namespace analysis::kinematic_fit;
        using namespace cuts::Htautau_Summer13;

        FitResultsMap result;
        TLorentzVector met_momentum;
        met_momentum.SetPtEtaPhiM(met.pt, 0, met.phi, 0);
        const TMatrix met_cov = ntuple::VectorToSignificanceMatrix(met.significanceMatrix);

        for(size_t n = 0; n < bjets.size(); ++n) {
            for(size_t k = n + 1; k < bjets.size(); ++k) {
                const FitId id(n, k);
                const four_body::FitInput input(bjets.at(n).momentum, bjets.at(k).momentum,
                                                higgs_to_taus.daughters.at(0).momentum,
                                                higgs_to_taus.daughters.at(1).momentum,
                                                met_momentum, met_cov);
                result[id] = FitWithUncertainties(input, cuts::jetCorrections::energyUncertainty,
                                                  tauCorrections::energyUncertainty, fit_two_bjets, fit_four_body);

            }
        }
        return result;
    }

    void FillHistogramsForMCstudies(const finalState::bbTauTau& final_state, const CandidateVector& bjets_all)
    {
        if(final_state.b_jets.size() >= 2 && final_state.taus.size() >= 2) {

            const VisibleGenObject& bjet1_visible = final_state.b_jets.at(0);
            const GenParticle* bjet1_MC = bjet1_visible.origin;

            const VisibleGenObject& bjet2_visible = final_state.b_jets.at(1);
            const GenParticle* bjet2_MC = bjet2_visible.origin;

            const VisibleGenObject tau_1 = final_state.taus.at(0);
            const VisibleGenObject tau_2 = final_state.taus.at(1);
            const GenParticle* tau_MC_1 = tau_1.origin;
            const GenParticle* tau_MC_2 = tau_2.origin;

            if (bjets_all.size() > 1 && bjet1_MC->momentum.Pt() > 20 && bjet2_MC->momentum.Pt() > 20 &&
                    std::abs(bjet1_MC->momentum.Eta()) < 2.4 && std::abs(bjet2_MC->momentum.Eta()) < 2.4 ){

                double deltaRmin_firstCouple =
                        std::min(bjet1_MC->momentum.DeltaR(tau_MC_1->momentum),bjet1_MC->momentum.DeltaR(tau_MC_2->momentum));
                double deltaRmin_secondCouple =
                        std::min(bjet2_MC->momentum.DeltaR(tau_MC_1->momentum),bjet2_MC->momentum.DeltaR(tau_MC_2->momentum));
                double deltaRmin_MC = std::min(deltaRmin_firstCouple,deltaRmin_secondCouple);

                //std::cout << "deltaRmin_MC=" << deltaRmin_MC << std::endl;
                GetAnaData().deltaRmin_MC().Fill(deltaRmin_MC);
                GetAnaData().DeltaRbjets_MC().Fill(bjet1_MC->momentum.DeltaR(bjet2_MC->momentum));
                GetAnaData().MinPtBjetsMC().Fill(std::min(bjet1_MC->momentum.Pt(),bjet2_MC->momentum.Pt()));
                TLorentzVector bb_MC = bjet1_MC->momentum + bjet2_MC->momentum;
                TLorentzVector bb_MC_visible = bjet1_visible.visibleMomentum + bjet2_visible.visibleMomentum;

                GetAnaData().MassBB_MC().Fill(bb_MC.M());
                GetAnaData().MassBB_MCvis().Fill(bb_MC_visible.M());

                double deltaRmin1_original = std::numeric_limits<double>::max();
                double deltaRmin2_original = std::numeric_limits<double>::max();
                double deltaRmin1_visible = std::numeric_limits<double>::max();
                double deltaRmin2_visible = std::numeric_limits<double>::max();
                unsigned index_bjet1 = 0;
                unsigned index_bjet2 = 0;
                unsigned index_bjet1_vis = 0;
                unsigned index_bjet2_vis = 0;
                for (unsigned i = 0; i < bjets_all.size(); ++i){
                    const Candidate& bjet = bjets_all.at(i);
    //                double deltaPt1 = std::abs(bjet.momentum.Pt() - bjet1_MC->momentum.Pt())/bjet1_MC->momentum.Pt();
    //                double deltaPt2 = std::abs(bjet.momentum.Pt() - bjet2_MC->momentum.Pt())/bjet2_MC->momentum.Pt();
                    if (bjet.momentum.DeltaR(bjet1_MC->momentum) < deltaRmin1_original /*&& deltaPt1 < 0.4*/){
                        deltaRmin1_original = bjet.momentum.DeltaR(bjet1_MC->momentum);
                        index_bjet1 = i;
                    }

                    if (bjet.momentum.DeltaR(bjet2_MC->momentum) < deltaRmin2_original /*&& deltaPt2 < 0.4*/){
                        deltaRmin2_original = bjet.momentum.DeltaR(bjet2_MC->momentum);
                        index_bjet2 = i;
                    }

                    if (bjet.momentum.DeltaR(bjet1_visible.visibleMomentum) < deltaRmin1_visible){
                        deltaRmin1_visible = bjet.momentum.DeltaR(bjet1_visible.visibleMomentum);
                        index_bjet1_vis = i;
                    }

                    if (bjet.momentum.DeltaR(bjet2_visible.visibleMomentum) < deltaRmin2_visible){
                        deltaRmin2_visible = bjet.momentum.DeltaR(bjet2_visible.visibleMomentum);
                        index_bjet2_vis = i;
                    }

                }

                GetAnaData().DeltaRmin1_original().Fill(deltaRmin1_original);
                GetAnaData().DeltaRmin2_original().Fill(deltaRmin2_original);

                double deltaRmax_original = std::max(deltaRmin1_original,deltaRmin2_original);
                GetAnaData().deltaRmax_original().Fill(deltaRmax_original);

                GetAnaData().DeltaRmin1_visible().Fill(deltaRmin1_visible);
                GetAnaData().DeltaRmin2_visible().Fill(deltaRmin2_visible);

                double deltaRmax_visible = std::max(deltaRmin1_visible,deltaRmin2_visible);
                GetAnaData().deltaRmax_visible().Fill(deltaRmax_visible);

                double deltaPtMax;
                if (deltaRmax_original == std::numeric_limits<double>::max()){
                    deltaPtMax = std::numeric_limits<double>::max();
                }
                else {
                    const Candidate& selectedBjets1 = bjets_all.at(index_bjet1);
                    const Candidate& selectedBjets2 = bjets_all.at(index_bjet2);
                    double deltaPt1 = std::abs(selectedBjets1.momentum.Pt() - bjet1_MC->momentum.Pt());
                    double deltaPt2 = std::abs(selectedBjets2.momentum.Pt() - bjet2_MC->momentum.Pt());
                    deltaPtMax = std::max(deltaPt1/bjet1_MC->momentum.Pt(),deltaPt2/bjet2_MC->momentum.Pt());
                }
                GetAnaData().deltaPtMax().Fill(deltaPtMax);

                double deltaPtMax_vis;
                if (deltaRmax_visible == std::numeric_limits<double>::max()){
                    deltaPtMax_vis = std::numeric_limits<double>::max();
                }
                else {
                    const Candidate& selectedBjets1 = bjets_all.at(index_bjet1_vis);
                    const Candidate& selectedBjets2 = bjets_all.at(index_bjet2_vis);
                    double deltaPt1 = std::abs(selectedBjets1.momentum.Pt() - bjet1_visible.visibleMomentum.Pt());
                    double deltaPt2 = std::abs(selectedBjets2.momentum.Pt() - bjet2_visible.visibleMomentum.Pt());
                    deltaPtMax_vis =
                            std::max(deltaPt1/bjet1_visible.visibleMomentum.Pt(),deltaPt2/bjet2_visible.visibleMomentum.Pt());
                }
                GetAnaData().deltaPtMax_vis().Fill(deltaPtMax_vis);
            }
        }
    }

    void FillFlatTree(const Candidate& higgs, const analysis::sv_fit::FitResultsWithUncertainties& svfitResults,
                      const analysis::kinematic_fit::FitResultsMap& kinfitResults,
                      const CandidateVector& jets , const CandidateVector& jetsPt20,
                      const CandidateVector& bjets_all, const CandidateVector& retagged_bjets,
                      const VertexVector& vertices, const Candidate& leg1, const Candidate& leg2,
                      const ntuple::MET& pfMET, const finalState::bbTauTau& final_state_MC)
    {
        static const float default_value = ntuple::DefaultFloatFillValueForFlatTree();

        // Event
        flatTree->run() = event->eventInfo().run;
        flatTree->lumi() = event->eventInfo().lumis;
        flatTree->evt() = event->eventInfo().EventId;
        flatTree->eventType() = static_cast<int>(DoEventCategorization(higgs));

        flatTree->npv() = vertices.size();
        if (config.ApplyPUreweight()){
            const size_t bxIndex = tools::find_index(event->eventInfo().bunchCrossing, 0);
            if(bxIndex >= event->eventInfo().bunchCrossing.size())
                throw std::runtime_error("in-time BX not found");
            flatTree->npu() = event->eventInfo().trueNInt.at(bxIndex);
        }

        // Weights
        flatTree->puweight()        = PUweight;
        flatTree->trigweight_1()    = triggerWeights.at(0);
        flatTree->trigweight_2()    = triggerWeights.at(1);
        flatTree->idweight_1()      = IDweights.at(0);
        flatTree->idweight_2()      = IDweights.at(1);
        flatTree->isoweight_1()     = IsoWeights.at(0);
        flatTree->isoweight_2()     = IsoWeights.at(1);
        // flatTree->fakeweight()      = fakeWeights.at(1); // what's this?
        flatTree->weight()          = eventWeight;
        flatTree->embeddedWeight()  =
                config.isDYEmbeddedSample() ? event->genEvent().embeddedWeight : 1.; // embedded weight

        // HTT candidate
        flatTree->mvis() = higgs.momentum.M();
        flatTree->m_sv_vegas() = svfitResults.fit_vegas.has_valid_mass ? svfitResults.fit_vegas.mass : default_value;
        flatTree->m_sv_up_vegas() = svfitResults.fit_vegas_up.has_valid_mass ? svfitResults.fit_vegas_up.mass : default_value;
        flatTree->m_sv_down_vegas() = svfitResults.fit_vegas_down.has_valid_mass ? svfitResults.fit_vegas_down.mass : default_value;
        flatTree->m_sv_MC() = svfitResults.fit_mc.has_valid_mass ? svfitResults.fit_mc.mass : default_value;
        flatTree->pt_sv_MC() = svfitResults.fit_mc.has_valid_pt ? svfitResults.fit_mc.pt : default_value;
        flatTree->m_sv_up_MC() = svfitResults.fit_mc_up.has_valid_mass ? svfitResults.fit_mc_up.mass : default_value;
        flatTree->pt_sv_up_MC() = svfitResults.fit_mc_up.has_valid_pt ? svfitResults.fit_mc_up.pt : default_value;
        flatTree->m_sv_down_MC() = svfitResults.fit_mc_down.has_valid_mass ? svfitResults.fit_mc_down.mass : default_value;
        flatTree->pt_sv_down_MC() = svfitResults.fit_mc_down.has_valid_pt ? svfitResults.fit_mc_down.pt : default_value;
        flatTree->DeltaR_leptons() = leg1.momentum.DeltaR(leg2.momentum) ;
        flatTree->pt_tt()          = (leg1.momentum + leg2.momentum).Pt();

        // Kinematic fit
        for(const auto& fit_result_entry : kinfitResults) {
            const analysis::kinematic_fit::four_body::FitResults& result_4body = fit_result_entry.second.fit_bb_tt;
            flatTree->kinfit_bb_tt_mass().push_back(result_4body.mass);
            flatTree->kinfit_bb_tt_convergence().push_back(result_4body.convergence);
            flatTree->kinfit_bb_tt_chi2().push_back(result_4body.chi2);
            flatTree->kinfit_bb_tt_pull_balance().push_back(result_4body.pull_balance);

//            const analysis::kinematic_fit::two_body::FitResults& result_2body = fit_result_entry.second.fit_bb;
//            TLorentzVector ditau_momentum;
//            ditau_momentum.SetPtEtaPhiM(svfitResults.fit_mc.pt, higgs.momentum.Eta(), higgs.momentum.Phi(),
//                                        svfitResults.fit_mc.mass);
//            const TLorentzVector combined_momentum = ditau_momentum + result_2body.bjet_momentums.at(0)
//                                                                    + result_2body.bjet_momentums.at(1);
//            flatTree->kinfit_bb_sv_mass().push_back(combined_momentum.M());
//            flatTree->kinfit_bb_convergence().push_back(result_2body.convergence);
//            flatTree->kinfit_bb_chi2().push_back(result_2body.chi2);
        }

        // Hhh generator info candidate
        if(final_state_MC.resonance) {
            const TLorentzVector& momentum = final_state_MC.resonance->momentum;
            flatTree->pt_resonance_MC()    = momentum.Pt()  ;
            flatTree->eta_resonance_MC()   = momentum.Eta() ;
            flatTree->phi_resonance_MC()   = momentum.Phi() ;
            flatTree->mass_resonance_MC()  = momentum.M()   ;
            flatTree->pdgId_resonance_MC() = final_state_MC.resonance->pdg.ToInteger();
        } else {
            flatTree->pt_resonance_MC()    = default_value  ;
            flatTree->eta_resonance_MC()   = default_value  ;
            flatTree->phi_resonance_MC()   = default_value  ;
            flatTree->mass_resonance_MC()  = default_value  ;
            flatTree->pdgId_resonance_MC() = particles::NONEXISTENT.RawCode() ;
        }

        if(final_state_MC.Higgs_TauTau) {
            const TLorentzVector& momentum = final_state_MC.Higgs_TauTau->momentum;
            flatTree->pt_Htt_MC()    = momentum.Pt()  ;
            flatTree->eta_Htt_MC()   = momentum.Eta() ;
            flatTree->phi_Htt_MC()   = momentum.Phi() ;
            flatTree->mass_Htt_MC()  = momentum.M()   ;
            flatTree->pdgId_Htt_MC() = final_state_MC.Higgs_TauTau->pdg.ToInteger();
        } else {
            flatTree->pt_Htt_MC()    = default_value  ;
            flatTree->eta_Htt_MC()   = default_value  ;
            flatTree->phi_Htt_MC()   = default_value  ;
            flatTree->mass_Htt_MC()  = default_value  ;
            flatTree->pdgId_Htt_MC() = particles::NONEXISTENT.RawCode() ;
        }

        if(final_state_MC.Higgs_BB) {
            const TLorentzVector& momentum = final_state_MC.Higgs_BB->momentum;
            flatTree->pt_Hbb_MC()    = momentum.Pt()  ;
            flatTree->eta_Hbb_MC()   = momentum.Eta() ;
            flatTree->phi_Hbb_MC()   = momentum.Phi() ;
            flatTree->mass_Hbb_MC()  = momentum.M()   ;
            flatTree->pdgId_Hbb_MC() = final_state_MC.Higgs_BB->pdg.ToInteger();
        } else {
            flatTree->pt_Hbb_MC()    = default_value  ;
            flatTree->eta_Hbb_MC()   = default_value  ;
            flatTree->phi_Hbb_MC()   = default_value  ;
            flatTree->mass_Hbb_MC()  = default_value  ;
            flatTree->pdgId_Hbb_MC() = particles::NONEXISTENT.RawCode()  ;
        }

        flatTree->n_extraJets_MC()     = default_value ; // needs to be filles with NUP! https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauTauAnalyzer.py#L51

        // Leg 1, lepton
        flatTree->pt_1()     = leg1.momentum.Pt()  ;
        flatTree->phi_1()    = leg1.momentum.Phi() ;
        flatTree->eta_1()    = leg1.momentum.Eta() ;
        flatTree->m_1()      = leg1.momentum.M()   ;
        flatTree->energy_1() = leg1.momentum.E()   ;
        flatTree->q_1()      = leg1.charge         ;
        flatTree->mt_1()     = analysis::Calculate_MT(leg1.momentum, postRecoilMET.pt, postRecoilMET.phi);
        flatTree->d0_1()     = analysis::Calculate_dxy(leg1.vertexPosition, primaryVertex.position,leg1.momentum);
        flatTree->dZ_1()     = leg1.vertexPosition.Z() - primaryVertex.position.Z();

        // Leg 2, tau
        flatTree->pt_2()     = leg2.momentum.Pt();
        flatTree->phi_2()    = leg2.momentum.Phi();
        flatTree->eta_2()    = leg2.momentum.Eta();
        flatTree->m_2()      = leg2.momentum.M();
        flatTree->energy_2() = leg2.momentum.E();
        flatTree->q_2()      = leg2.charge;
        flatTree->mt_2()     = analysis::Calculate_MT(leg2.momentum, postRecoilMET.pt, postRecoilMET.phi);
        flatTree->d0_2()     = analysis::Calculate_dxy(leg2.vertexPosition, primaryVertex.position,leg2.momentum);
        flatTree->dZ_2()     = leg2.vertexPosition.Z() - primaryVertex.position.Z();

        // RM: for the three channels, mt, et, tt this leg is always a tau
        const ntuple::Tau& ntuple_tau_leg2 = correctedTaus.at(leg2.index);
        flatTree->decayMode_2()                                = ntuple_tau_leg2.decayMode;
        flatTree->againstElectronLooseMVA_2() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau_leg2.againstElectronMVA3category, ntuple_tau_leg2.againstElectronMVA3raw, 0);
        flatTree->againstElectronMediumMVA_2() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau_leg2.againstElectronMVA3category, ntuple_tau_leg2.againstElectronMVA3raw, 1);
        flatTree->againstElectronTightMVA_2() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau_leg2.againstElectronMVA3category, ntuple_tau_leg2.againstElectronMVA3raw, 2);
        flatTree->againstElectronVTightMVA_2() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau_leg2.againstElectronMVA3category, ntuple_tau_leg2.againstElectronMVA3raw, 3);
        flatTree->againstElectronLoose_2()                     = ntuple_tau_leg2.againstElectronLoose  ;
        flatTree->againstElectronMedium_2()                    = ntuple_tau_leg2.againstElectronMedium ;
        flatTree->againstElectronTight_2()                     = ntuple_tau_leg2.againstElectronTight  ;
        flatTree->againstMuonLoose_2()                         = ntuple_tau_leg2.againstMuonLoose      ;
        flatTree->againstMuonMedium_2()                        = ntuple_tau_leg2.againstMuonMedium     ;
        flatTree->againstMuonTight_2()                         = ntuple_tau_leg2.againstMuonTight      ;
        flatTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = ntuple_tau_leg2.byCombinedIsolationDeltaBetaCorrRaw3Hits ;

        // MET
        TLorentzVector postRecoilMetMomentum;
        postRecoilMetMomentum.SetPtEtaPhiM(postRecoilMET.pt, 0, postRecoilMET.phi, 0.);
        flatTree->pt_tt_MET() = (leg1.momentum + leg2.momentum + postRecoilMetMomentum).Pt();

        flatTree->met()       = pfMET.pt          ;
        flatTree->metphi()    = pfMET.phi         ;
        flatTree->mvamet()    = postRecoilMET.pt  ;
        flatTree->mvametphi() = postRecoilMET.phi ;
        //flatTree->pzetavis();
        //flatTree->pzetamiss();
        if(pfMET.significanceMatrix.size()) {
            const TMatrixD metPFcov = ntuple::VectorToSignificanceMatrix(pfMET.significanceMatrix);
            flatTree->metcov00() = metPFcov[0][0];
            flatTree->metcov01() = metPFcov[0][1];
            flatTree->metcov10() = metPFcov[1][0];
            flatTree->metcov11() = metPFcov[1][1];
        }
        const TMatrixD metMVAcov = ntuple::VectorToSignificanceMatrix(postRecoilMET.significanceMatrix);
        flatTree->mvacov00() = metMVAcov[0][0];
        flatTree->mvacov01() = metMVAcov[0][1];
        flatTree->mvacov10() = metMVAcov[1][0];
        flatTree->mvacov11() = metMVAcov[1][1];

        // Jets
        flatTree->njets()     = jets.size()           ;
        flatTree->njetspt20() = jetsPt20.size()       ;
        flatTree->nBjets()    = retagged_bjets.size() ;

        for (const Candidate& jet : bjets_all) {
            const ntuple::Jet& ntuple_jet = event->jets().at(jet.index);

            flatTree->pt_Bjets()      .push_back( jet.momentum.Pt() );
            flatTree->eta_Bjets()     .push_back( jet.momentum.Eta() );
            flatTree->phi_Bjets()     .push_back( jet.momentum.Phi() );
            flatTree->energy_Bjets()  .push_back( jet.momentum.E() );
            flatTree->chargedHadronEF_Bjets().push_back( ntuple_jet.chargedHadronEnergyFraction );
            flatTree->neutralHadronEF_Bjets()  .push_back( ntuple_jet.neutralHadronEnergyFraction );
            flatTree->photonEF_Bjets()         .push_back( ntuple_jet.photonEnergyFraction );
            flatTree->muonEF_Bjets()  .push_back( ntuple_jet.muonEnergyFraction );
            flatTree->electronEF_Bjets()  .push_back( ntuple_jet.electronEnergyFraction );
            flatTree->csv_Bjets()     .push_back( ntuple_jet.combinedSecondaryVertexBJetTags );
            // inspect the flavour of the gen jet
            const VisibleGenObjectVector matched_bjets_MC = FindMatchedObjects(jet.momentum, final_state_MC.b_jets,
                                                                               cuts::DeltaR_MC_Match);
            const bool isJet_MC_Bjet = matched_bjets_MC.size() != 0;
            const bool isJet_MC_Bjet_withLeptonicDecay = isJet_MC_Bjet
                    && matched_bjets_MC.at(0).finalStateChargedLeptons.size() != 0;
            flatTree->isBjet_MC_Bjet()                  .push_back( isJet_MC_Bjet );
            flatTree->isBjet_MC_Bjet_withLeptonicDecay().push_back( isJet_MC_Bjet_withLeptonicDecay );
        }

        flatTree->x_PV() = primaryVertex.position.x();
        flatTree->y_PV() = primaryVertex.position.y();
        flatTree->z_PV() = primaryVertex.position.z();
    }

protected:
    std::shared_ptr<ntuple::FlatTree> flatTree;
    bool writeFlatTree;
    ntuple::TauVector correctedTaus;
    ntuple::MET correctedMET;
    ntuple::MET postRecoilMET;
};
} // analysis
