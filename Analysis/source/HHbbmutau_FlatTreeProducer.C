/*!
 * \file HmutauFlatTreeProducer.C
 * \brief Generate flat-tree for Htautau analysis using looser selection.
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

#include "Analysis/include/BaseFlatTreeProducer.h"

class HHbbmutau_FlatTreeProducer : public analysis::BaseFlatTreeProducer {
public:
    HHbbmutau_FlatTreeProducer(const std::string& inputFileName, const std::string& outputFileName,
                               const std::string& configFileName, const std::string& _prefix = "none",
                               size_t _maxNumberOfEvents = 0,
                               std::shared_ptr<ntuple::FlatTree> _flatTree = std::shared_ptr<ntuple::FlatTree>())
        : BaseFlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents, _flatTree),
          anaData(*outputFile)
    {
        anaData.getOutputFile().cd();
    }

    virtual analysis::BaseAnalyzerData& GetAnaData() override { return anaData; }

    virtual void ProcessEvent(std::shared_ptr<const analysis::EventDescriptor> _event) override
    {
        BaseFlatTreeProducer::ProcessEvent(_event);
        using namespace analysis;
        using namespace cuts::Htautau_Summer13;
        using namespace cuts::Htautau_Summer13::MuTau;


        if (!FindAnalysisFinalState(muTau_MC) && config.RequireSpecificFinalState()) return;



        cuts::Cutter cut(&anaData.Selection("event"));
        cut(true, "total");

        cut(HaveTriggerFired(trigger::hltPaths), "trigger");

        const VertexVector vertices = CollectVertices();
        cut(vertices.size(), "vertex");
        primaryVertex = vertices.front();

        const auto z_muons = CollectZmuons();
        const auto z_muon_candidates = FindCompatibleObjects(z_muons, ZmumuVeto::deltaR, Candidate::Z, "Z_mu_mu", 0);
        cut(!z_muon_candidates.size(), "z_mumu_veto");

        const auto electrons_bkg = CollectBackgroundElectrons();
        cut(!electrons_bkg.size(), "no_electron");

        const auto muons = CollectSignalMuons();
        cut(muons.size(), "muon_cand");
        cut(muons.size() == 1, "one_muon_cand");

        const auto muons_extra = CollectBackgroundMuons();
        const bool have_bkg_muon = muons_extra.size() > 1 ||
                ( muons_extra.size() == 1 && muons_extra.front() != muons.front() );
        cut(!have_bkg_muon, "no_bkg_muon");

        const auto allmuons = CollectMuons();

        correctedTaus = config.ApplyTauESCorrection() ? ApplyTauCorrections(muTau_MC.hadronic_taus,false) : event->taus();

        const auto alltaus = CollectTaus();
        cut(alltaus.size(), "tau_cand");

        const auto signaltaus = CollectSignalTaus() ;

        analysis::CandidateVector higgses ;

        // First check OS, isolated higgs candidates
        const auto higgses_sig = FindCompatibleObjects(muons, signaltaus, DeltaR_betweenSignalObjects,
                                                   Candidate::Higgs, "H_mu_tau",0);

        // If no OS candidate, keep any higgs-ish candidates for bkg estimation (don't cut on sign nor isolation)
        if (!higgses_sig.size())
        {
          const auto higgses_bkg = FindCompatibleObjects(allmuons, alltaus, DeltaR_betweenSignalObjects,
                                                     Candidate::Higgs, "H_mu_tau");
          higgses = higgses_bkg ;
        }
        else
        {
          higgses = higgses_sig ;
        }

        cut(higgses.size(), "mu_tau");

        const auto higgsTriggered = ApplyTriggerMatch(higgses,trigger::hltPaths,false);
        cut(higgsTriggered.size(), "trigger obj match");

        const Candidate higgs = SelectSemiLeptonicHiggs(higgsTriggered);

        const ntuple::MET mvaMet = mvaMetProducer.ComputeMvaMet(higgs,event->pfCandidates(),
                                                                event->jets(),primaryVertex,
                                                                vertices,event->taus());

        correctedMET = config.ApplyTauESCorrection() ? ApplyTauCorrectionsToMVAMET(mvaMet, correctedTaus) : mvaMet;

        const auto looseJets = CollectLooseJets();

        const auto filteredLooseJets = FilterCompatibleObjects(looseJets, higgs,
                                                               cuts::Htautau_Summer13::jetID::deltaR_signalObjects);


        const auto jets = CollectJets(filteredLooseJets);
        const auto bjets_all = CollectBJets(filteredLooseJets, false, false);
        const auto retagged_bjets = CollectBJets(filteredLooseJets, config.isMC(), true);


        postRecoilMET = ApplyRecoilCorrections(higgs, muTau_MC.resonance, jets.size(), correctedMET);

        const auto svfitResults = analysis::sv_fit::FitWithUncertainties(higgs, postRecoilMET,
                                                                         tauCorrections::energyUncertainty,
                                                                         true, true);
//        std::cout << "Event ID: " << event->eventInfo().EventId << std::endl;
        const auto kinfitResults = RunKinematicFit(bjets_all, higgs, postRecoilMET, true, true);

        CalculateFullEventWeight(higgs);

        const ntuple::MET pfMET = config.isMC() ? event->metPF() : mvaMetProducer.ComputePFMet(event->pfCandidates(), primaryVertex);

        FillFlatTree(higgs, svfitResults, kinfitResults, jets, filteredLooseJets, bjets_all, retagged_bjets, vertices, pfMET);
    }

protected:


    virtual analysis::Candidate SelectMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                           root_ext::AnalyzerData& _anaData,
                                           const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::muonID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Muon& object = event->muons().at(id);
        const analysis::Candidate muon(analysis::Candidate::Mu, id, object);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(isGlobalMuonPromptTight) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon) == isPFMuon, "PF");
        cut(X(nMatchedStations) > nMatched_Stations, "stations");
        cut(X(pixHits) > pixHits, "pix_hits");
        cut(X(trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double dB_PV = analysis::Calculate_dxy(mu_vertex,primaryVertex.position,muon.momentum);
        cut(std::abs( Y(dB_PV) ) < dB, "dB");
        //cut(X(pfRelIso) < cuts::skim::MuTau::pFRelIso, "pFRelIso");

        return muon;
    }

    virtual analysis::Candidate SelectSignalMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                                 root_ext::AnalyzerData& _anaData,
                                                 const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::muonID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Muon& object = event->muons().at(id);
        const analysis::Candidate muon(analysis::Candidate::Mu, id, object);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(isGlobalMuonPromptTight) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon) == isPFMuon, "PF");
        cut(X(nMatchedStations) > nMatched_Stations, "stations");
        cut(X(pixHits) > pixHits, "pix_hits");
        cut(X(trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double dB_PV = analysis::Calculate_dxy(mu_vertex,primaryVertex.position,muon.momentum);
        cut(std::abs( Y(dB_PV) ) < dB, "dB");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");

        return muon;
    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData,
                                          const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::tauID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Tau& object = correctedTaus.at(id);

        cut(true, ">0 tau cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(againstMuonLoose) > cuts::skim::MuTau::tauID::againstMuonLoose, "vs_mu_loose");
        cut(X(againstElectronLoose) > againstElectronLoose, "vs_e_loose");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits) < cuts::skim::MuTau::tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits, "looseIso3Hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");


        const analysis::Candidate tau(analysis::Candidate::Tau, id, object);

        return tau;
    }

    virtual analysis::Candidate SelectSignalTau(size_t id, cuts::ObjectSelector* objectSelector,
                                                root_ext::AnalyzerData& _anaData,
                                                const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::tauID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Tau& object = correctedTaus.at(id);

        cut(true, ">0 tau cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(againstMuonTight) > againstMuonTight, "vs_mu_tight");
        cut(X(againstElectronLoose) > againstElectronLoose, "vs_e_loose");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits) < byCombinedIsolationDeltaBetaCorrRaw3Hits, "looseIso3Hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");


        const analysis::Candidate tau(analysis::Candidate::Tau, id, object);

        return tau;
    }

    analysis::CandidateVector CollectZmuons()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> analysis::Candidate
            { return SelectZmuon(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<analysis::Candidate>("z_muons", base_selector, event->muons().size());
    }

    virtual analysis::Candidate SelectZmuon(size_t id, cuts::ObjectSelector* objectSelector,
                                            root_ext::AnalyzerData& _anaData,
                                            const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::MuTau::ZmumuVeto;
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
        cut(X(isTrackerMuon) == isTrackerMuon, "trackerMuon");
        cut(X(isPFMuon) == isPFMuon, "PFMuon");
        cut(X(pfRelIso) < pfRelIso, "pFRelIso");

        return muon;
    }

    bool FindAnalysisFinalState(analysis::finalState::bbMuTaujet& final_state)
    {
        const bool base_result = BaseFlatTreeProducer::FindAnalysisFinalState(final_state);
        if(!base_result)
            return base_result;

        for(const analysis::VisibleGenObject& tau_MC : final_state.taus) {
            if(tau_MC.finalStateChargedLeptons.size() == 1
                    && tau_MC.finalStateChargedLeptons.at(0)->pdg.Code == particles::mu)
                final_state.muon = tau_MC.finalStateChargedLeptons.at(0);
            else if(tau_MC.finalStateChargedHadrons.size() >= 1)
                final_state.tau_jet = &tau_MC;
        }

        if (!final_state.muon || !final_state.tau_jet) return false;
        return true;
    }

    virtual void CalculateTriggerWeights(const analysis::Candidate& higgs) override
    {
        using namespace analysis::Htautau_Summer13::trigger::Run2012ABCD::MuTau;
        const analysis::Candidate& mu = higgs.GetDaughter(analysis::Candidate::Mu);
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        triggerWeights = CalculateWeights(mu.momentum, tau.momentum);
    }

    virtual void CalculateIsoWeights(const analysis::Candidate& higgs) override
    {
        using namespace cuts::Htautau_Summer13::MuTau::muonISOscaleFactor;
        IsoWeights.clear();
        const analysis::Candidate& mu = higgs.GetDaughter(analysis::Candidate::Mu);
        const double mu_pt = mu.momentum.Pt(), mu_eta = std::abs(mu.momentum.Eta());
        if(mu_pt < pt.at(0))
            throw std::runtime_error("No information about ISO. Muon pt is too small");
        const size_t pt_bin = mu_pt < pt.at(1) ? 0 : 1;
        if(mu_eta >= eta.at(2))
            throw std::runtime_error("No information about ISO. Muon eta is too big");
        const size_t eta_bin = mu_eta < eta.at(0) ? 0 : ( mu_eta < eta.at(1) ? 1 : 2 );
        const double scale = scaleFactors.at(pt_bin).at(eta_bin);
        // first mu, second tau
        IsoWeights.push_back(scale);
        IsoWeights.push_back(1);
    }

    virtual void CalculateIdWeights(const analysis::Candidate& higgs) override
    {
        using namespace cuts::Htautau_Summer13::MuTau::muonIDscaleFactor;
        IDweights.clear();
        const analysis::Candidate& mu = higgs.GetDaughter(analysis::Candidate::Mu);
        const double mu_pt = mu.momentum.Pt(), mu_eta = std::abs(mu.momentum.Eta());
        if(mu_pt < pt.at(0))
            throw std::runtime_error("No information about ID. Muon pt is too small");
        const size_t pt_bin = mu_pt < pt.at(1) ? 0 : 1;
        if(mu_eta >= eta.at(2))
            throw std::runtime_error("No information about ID. Muon eta is too big");
        const size_t eta_bin = mu_eta < eta.at(0) ? 0 : ( mu_eta < eta.at(1) ? 1 : 2 );
        const double scale = scaleFactors.at(pt_bin).at(eta_bin);
        // first mu, second tau
        IDweights.push_back(scale);
        IDweights.push_back(1);
    }

    virtual void CalculateFakeWeights(const analysis::Candidate& higgs) override
    {
        fakeWeights.clear();
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        const double tau_pt = tau.momentum.Pt();
        const double fake_weight =
                (1.15743)-(0.00736136*tau_pt)+(4.3699e-05*tau_pt*tau_pt)-(1.188e-07*tau_pt*tau_pt*tau_pt);
        // first mu, second tau
        fakeWeights.push_back(1);
        fakeWeights.push_back(fake_weight);
    }

    virtual void CalculateDMWeights(const analysis::Candidate& higgs) override
    {
        DMweights.clear();
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        const ntuple::Tau& tau_leg = correctedTaus.at(tau.index);
        const double weight = tau_leg.decayMode == ntuple::tau_id::kOneProng0PiZero ?
                                  cuts::Htautau_Summer13::tauCorrections::DecayModeWeight : 1;
        // first mu, second tau
        DMweights.push_back(1);
        DMweights.push_back(weight);
    }

    void FillFlatTree(const analysis::Candidate& higgs,
                      const analysis::sv_fit::FitResultsWithUncertainties& svfitResults,
                      const analysis::kinematic_fit::FitResultsWithUncertainties& kinfitResults,
                      const analysis::CandidateVector& jets, const analysis::CandidateVector& jetsPt20,
                      const analysis::CandidateVector& bjets_all, const analysis::CandidateVector& retagged_bjets,
                      const analysis::VertexVector& vertices, const ntuple::MET& pfMET)
    {
        static const float default_value = ntuple::DefaultFloatFillValueForFlatTree();
        static const float default_int_value = ntuple::DefaultIntegerFillValueForFlatTree();

        const analysis::Candidate& muon = higgs.GetDaughter(analysis::Candidate::Mu);
        const ntuple::Muon& ntuple_muon = event->muons().at(muon.index);
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);

        BaseFlatTreeProducer::FillFlatTree(higgs, svfitResults, kinfitResults, jets, jetsPt20, bjets_all,
                                           retagged_bjets, vertices, muon, tau, pfMET, muTau_MC);

        flatTree->channel() = static_cast<int>(ntuple::Channel::MuTau);
        flatTree->pfRelIso_1() = ntuple_muon.pfRelIso;
        flatTree->mva_1() = 0;
        flatTree->passid_1() = true;
        flatTree->passiso_1() = ntuple_muon.pfRelIso < cuts::Htautau_Summer13::MuTau::muonID::pFRelIso;
        flatTree->decayMode_1() = default_int_value;
        flatTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = default_value;
        flatTree->againstElectronLooseMVA_1() = false;
        flatTree->againstElectronMediumMVA_1() = false;
        flatTree->againstElectronTightMVA_1() = false;
        flatTree->againstElectronVTightMVA_1() = false;
        flatTree->againstElectronLoose_1() = false;
        flatTree->againstElectronMedium_1() = false;
        flatTree->againstElectronTight_1() = false;
        flatTree->againstMuonLoose_1() = false;
        flatTree->againstMuonMedium_1() = false;
        flatTree->againstMuonTight_1() = false;

        const bool muon_matched = analysis::HasMatchWithMCParticle(muon.momentum, muTau_MC.muon, cuts::DeltaR_MC_Match);
        if(muon_matched) {
            const TLorentzVector& momentum = muTau_MC.muon->momentum;
            flatTree->pt_1_MC   () = momentum.Pt()  ;
            flatTree->phi_1_MC  () = momentum.Phi() ;
            flatTree->eta_1_MC  () = momentum.Eta() ;
            flatTree->m_1_MC    () = momentum.M()   ;
            flatTree->pt_1_visible_MC   () = momentum.Pt()  ;
            flatTree->phi_1_visible_MC  () = momentum.Phi() ;
            flatTree->eta_1_visible_MC  () = momentum.Eta() ;
            flatTree->m_1_visible_MC    () = momentum.M()   ;
            flatTree->pdgId_1_MC() = muTau_MC.muon->pdg.ToInteger();
        } else {
            flatTree->pt_1_MC   () = default_value ;
            flatTree->phi_1_MC  () = default_value ;
            flatTree->eta_1_MC  () = default_value ;
            flatTree->m_1_MC    () = default_value ;
            flatTree->pt_1_visible_MC   () = default_value ;
            flatTree->phi_1_visible_MC  () = default_value ;
            flatTree->eta_1_visible_MC  () = default_value ;
            flatTree->m_1_visible_MC    () = default_value ;
            flatTree->pdgId_1_MC() = particles::NONEXISTENT.RawCode();
        }

        const bool tau_matched = analysis::HasMatchWithMCObject(tau.momentum, muTau_MC.tau_jet, cuts::DeltaR_MC_Match);
        if(tau_matched) {
            const TLorentzVector& momentum = muTau_MC.tau_jet->origin->momentum;
            flatTree->pt_2_MC   () = momentum.Pt()  ;
            flatTree->phi_2_MC  () = momentum.Phi() ;
            flatTree->eta_2_MC  () = momentum.Eta() ;
            flatTree->m_2_MC    () = momentum.M()   ;
            const TLorentzVector& visible_momentum = muTau_MC.tau_jet->visibleMomentum;
            flatTree->pt_2_visible_MC   () = visible_momentum.Pt()  ;
            flatTree->phi_2_visible_MC  () = visible_momentum.Phi() ;
            flatTree->eta_2_visible_MC  () = visible_momentum.Eta() ;
            flatTree->m_2_visible_MC    () = visible_momentum.M()   ;
            flatTree->pdgId_2_MC() = muTau_MC.tau_jet->origin->pdg.ToInteger();
        } else {
            flatTree->pt_2_MC   () = default_value ;
            flatTree->phi_2_MC  () = default_value ;
            flatTree->eta_2_MC  () = default_value ;
            flatTree->m_2_MC    () = default_value ;
            flatTree->pt_2_visible_MC   () = default_value ;
            flatTree->phi_2_visible_MC  () = default_value ;
            flatTree->eta_2_visible_MC  () = default_value ;
            flatTree->m_2_visible_MC    () = default_value ;
            flatTree->pdgId_2_MC() = particles::NONEXISTENT.RawCode();
        }

        flatTree->Fill();
    }

private:
    analysis::BaseAnalyzerData anaData;
    analysis::finalState::bbMuTaujet muTau_MC;
};

#include "METPUSubtraction/interface/GBRProjectDict.cxx"
