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
        static const particles::ParticleCodes ZDecay_1 = { particles::tau, particles::tau };
        static const particles::ParticleCodes ZDecay_2 = { particles::e, particles::e };
        static const particles::ParticleCodes ZDecay_3 = { particles::mu, particles::mu };

        static const particles::ParticleCodes particlesCodes = { particles::e, particles::mu };

        const GenParticleSet Zparticles = genEvent.GetParticles(Zcode);

        if (Zparticles.size() > 1 || Zparticles.size() == 0)
            throw exception("not 1 Z per event");

        const GenParticle* Z_mc = *Zparticles.begin();
        GenParticlePtrVector ZProducts;
        if(FindDecayProducts(*Z_mc, ZDecay_1,ZProducts))
            return ntuple::EventType::ZTT;

        if (!FindDecayProducts(*Z_mc, ZDecay_2,ZProducts) &&
                !FindDecayProducts(*Z_mc, ZDecay_3,ZProducts))
            throw exception("not leptonic decay");


        const GenParticleSet selectedParticles = genEvent.GetParticles(particlesCodes,minimal_genParticle_pt);

        const analysis::Candidate& firstDaughter = higgs.daughters.at(0);
        const analysis::Candidate& secondDaughter = higgs.daughters.at(1);
        const GenParticleSet matchedParticles_first =
                FindMatchedParticles(firstDaughter.momentum,selectedParticles,deltaR_matchGenParticle);
        const GenParticleSet matchedParticles_second =
                FindMatchedParticles(secondDaughter.momentum,selectedParticles,deltaR_matchGenParticle);

        if (matchedParticles_first.size() == 1 && matchedParticles_second.size() == 1)
            return ntuple::EventType::ZL;
        else if (matchedParticles_first.size() == 0 && matchedParticles_second.size() == 0)
            return ntuple::EventType::ZJ;
        else
            return ntuple::EventType::OtherZ;
    }

    bool GenFilterForZevents(const finalState::bbTauTau& final_state)
    {
        if (final_state.taus.size() != 2)
            throw exception("not 2 taus in the event at Gen Level");
        const GenParticle* firstTau = final_state.taus.at(0).GetOriginGenParticle();
        const GenParticle* secondTau = final_state.taus.at(1).GetOriginGenParticle();
        if ((firstTau->momentum + secondTau->momentum).M() > 50) return true;
        return false;
    }

    bool MatchTausFromHiggsWithGenTaus(const analysis::Candidate& higgs, const finalState::bbTauTau& final_state)
    {
        using namespace cuts::Htautau_Summer13::Embedded;
        const analysis::Candidate& firstDaughter = higgs.daughters.at(0);
        const analysis::Candidate& secondDaughter = higgs.daughters.at(1);
        const VisibleGenObjectVector matchedParticles_first =
                FindMatchedObjects(firstDaughter.momentum,final_state.taus,deltaR_betweenTaus);
        const VisibleGenObjectVector matchedParticles_second =
                FindMatchedObjects(secondDaughter.momentum,final_state.taus,deltaR_betweenTaus);

        if (matchedParticles_first.size() == 1 && matchedParticles_second.size() == 1)
            return true;
        return false;
    }


    analysis::kinematic_fit::FitResultsWithUncertainties RunKinematicFit(const CandidateVector& bjets,
                                                                         const Candidate& higgs_to_taus,
                                                                         const ntuple::MET& met, bool fit_two_bjets,
                                                                         bool fit_four_body)
    {
        using namespace analysis::kinematic_fit;
        using namespace cuts::Htautau_Summer13;
        if(bjets.size() < 2)
            return FitResultsWithUncertainties();
        FitResultsWithUncertainties result;
        if(bjets.size() >= 2) {
            TLorentzVector met_momentum;
            met_momentum.SetPtEtaPhiM(met.pt, 0, met.phi, 0);
            const TMatrix met_cov = ntuple::VectorToSignificanceMatrix(met.significanceMatrix);
            const FourBodyFitInput input(bjets.at(0).momentum, bjets.at(1).momentum,
                                         higgs_to_taus.daughters.at(0).momentum, higgs_to_taus.daughters.at(1).momentum,
                                         met_momentum, met_cov);
            result = FitWithUncertainties(input, cuts::jetCorrections::energyUncertainty,
                                          tauCorrections::energyUncertainty, fit_two_bjets, fit_four_body);
        }
        return result;
    }

    void FillFlatTree(const Candidate& higgs, const analysis::sv_fit::FitResultsWithUncertainties& svfitResults,
                      const analysis::kinematic_fit::FitResultsWithUncertainties& kinfitResults,
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
                config.isEmbeddedSample() ? event->genEvent().embeddedWeight : 1.; // embedded weight

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
        //flatTree->m_kinfit_bb_tt() = kinfitResults.fit_bb_tt.has_valid_mass ? kinfitResults.fit_bb_tt.mass : default_value;
        flatTree->m_kinfit_bb_down_tt_down() = kinfitResults.fit_bb_down_tt_down.has_valid_mass ? kinfitResults.fit_bb_down_tt_down.mass : default_value;
        flatTree->m_kinfit_bb_down_tt_up() = kinfitResults.fit_bb_down_tt_up.has_valid_mass ? kinfitResults.fit_bb_down_tt_up.mass : default_value;
        flatTree->m_kinfit_bb_up_tt_down() = kinfitResults.fit_bb_up_tt_down.has_valid_mass ? kinfitResults.fit_bb_up_tt_down.mass : default_value;
        flatTree->m_kinfit_bb_up_tt_up() = kinfitResults.fit_bb_up_tt_up.has_valid_mass ? kinfitResults.fit_bb_up_tt_up.mass : default_value;
        flatTree->m_kinfit_bb() = kinfitResults.fit_bb.has_valid_mass ? kinfitResults.fit_bb.mass : default_value;
        flatTree->m_kinfit_bb_down() = kinfitResults.fit_bb_down.has_valid_mass ? kinfitResults.fit_bb_down.mass : default_value;
        flatTree->m_kinfit_bb_up() = kinfitResults.fit_bb_up.has_valid_mass ? kinfitResults.fit_bb_up.mass : default_value;

        flatTree->m_kinfit_bb_tt() = kinfitResults.fit_bb_tt.mass;
        flatTree->m_kinfit_bb_tt_convergence() = kinfitResults.fit_bb_tt.convergence;
        flatTree->m_kinfit_bb_tt_chi2() = kinfitResults.fit_bb_tt.chi2;
        flatTree->m_kinfit_bb_tt_pull_balance() = kinfitResults.fit_bb_tt.pull_balance;

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

        // RM: SAVE ALL THE JETS WITH PT > 20 (AND SAVE THE CSV DISCRIMINATOR)
        CandidateVector pt_sorted_jetsPt20 = jetsPt20;
        std::sort( pt_sorted_jetsPt20.begin(), pt_sorted_jetsPt20.end(),
                   []( const Candidate& a, const Candidate& b ) { return a.momentum.Pt() > b.momentum.Pt(); } ) ;

        // fill the jet collections
        for (const Candidate& jet : pt_sorted_jetsPt20) {
            const ntuple::Jet& ntuple_jet = event->jets().at(jet.index);
            flatTree->pt_jets()      .push_back( jet.momentum.Pt() );
            flatTree->eta_jets()     .push_back( jet.momentum.Eta() );
            flatTree->phi_jets()     .push_back( jet.momentum.Phi() );
            flatTree->mass_jets()    .push_back( jet.momentum.M() );
            flatTree->energy_jets()  .push_back( jet.momentum.E() );
            flatTree->csv_jets()     .push_back( ntuple_jet.combinedSecondaryVertexBJetTags );
            // inspect the flavour of the gen jet
            const VisibleGenObjectVector matched_bjets_MC = FindMatchedObjects(jet.momentum, final_state_MC.b_jets,
                                                                               cuts::DeltaR_MC_Match);
            const bool isJet_MC_Bjet = matched_bjets_MC.size() != 0;
            const bool isJet_MC_Bjet_withLeptonicDecay = isJet_MC_Bjet
                    && matched_bjets_MC.at(0).finalStateChargedLeptons.size() != 0;
            flatTree->isJet_MC_Bjet()                  .push_back( isJet_MC_Bjet );
            flatTree->isJet_MC_Bjet_withLeptonicDecay().push_back( isJet_MC_Bjet_withLeptonicDecay );
        }

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
