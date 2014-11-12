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

namespace analysis {
struct SelectionResults_mutau : public SelectionResults {
    finalState::bbMuTaujet muTau_MC;
    const Candidate& GetMuon() const { return higgs.GetDaughter(Candidate::Mu); }
    const Candidate& GetTau() const { return higgs.GetDaughter(Candidate::Tau); }

    virtual const Candidate& GetLeg1() const override { return GetMuon(); }
    virtual const Candidate& GetLeg2() const override { return GetTau(); }
    virtual const finalState::bbTauTau& GetFinalStateMC() const override { return muTau_MC; }
};
} // namespace analysis

class HHbbmutau_FlatTreeProducer : public virtual analysis::BaseFlatTreeProducer {
public:
    HHbbmutau_FlatTreeProducer(const std::string& inputFileName, const std::string& outputFileName,
                               const std::string& configFileName, const std::string& _prefix = "none",
                               size_t _maxNumberOfEvents = 0,
                               std::shared_ptr<ntuple::FlatTree> _flatTree = std::shared_ptr<ntuple::FlatTree>())
        : BaseFlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents, _flatTree),
          baseAnaData(*outputFile)
    {
        baseAnaData.getOutputFile().cd();
        if(config.ApplyRecoilCorrection())
            recoilCorrectionProducer_mutau = std::shared_ptr<analysis::RecoilCorrectionProducer>(
                        new analysis::RecoilCorrectionProducer(config.RecoilCorrection_fileCorrectTo_MuTau(),
                                                               config.RecoilCorrection_fileZmmData_MuTau(),
                                                               config.RecoilCorrection_fileZmmMC_MuTau()));
    }

    virtual analysis::BaseAnalyzerData& GetAnaData() override { return baseAnaData; }
    virtual analysis::RecoilCorrectionProducer& GetRecoilCorrectionProducer() override
    {
        return *recoilCorrectionProducer_mutau;
    }

protected:
    virtual analysis::SelectionResults& ApplyBaselineSelection() override
    {
        using namespace analysis;
        using namespace cuts::Htautau_Summer13;
        using namespace cuts::Htautau_Summer13::MuTau;

        selection = SelectionResults_mutau();

        cuts::Cutter cut(&GetAnaData().Selection("event"));
        cut(true, "total");

        cut(FindAnalysisFinalState(selection.muTau_MC) || !config.RequireSpecificFinalState(), "spec_final_state");
        cut(!config.isDYEmbeddedSample() || GenFilterForZevents(selection.muTau_MC), "genFilter");

        const auto& selectedTriggerPath = config.isDYEmbeddedSample()
                ? DYEmbedded::trigger::hltPaths : trigger::hltPaths;
        cut(HaveTriggerFired(selectedTriggerPath), "trigger");

        selection.vertices = CollectVertices();
        cut(selection.vertices.size(), "vertex");
        primaryVertex = selection.vertices.front();

        const auto z_muons = CollectZmuons();
        const auto z_muon_candidates = FindCompatibleObjects(z_muons, ZmumuVeto::deltaR, Candidate::Z, "Z_mu_mu", 0);
        cut(!z_muon_candidates.size(), "z_mumu_veto");

        const auto electrons_bkg = CollectBackgroundElectrons();
        cut(!electrons_bkg.size(), "no_electron");

        const auto muons = CollectSignalMuons();

        const auto muons_extra = CollectBackgroundMuons();
        const bool have_bkg_muon = muons_extra.size() > 1 || muons.size() > 1 ||
                ( muons_extra.size() == 1 && muons.size() == 1 && muons_extra.front() != muons.front() );
        cut(!have_bkg_muon, "no_bkg_muon");

        const auto allmuons = CollectMuons();
        cut(allmuons.size(), "muon_cand");

        correctedTaus = config.ApplyTauESCorrection()
                ? ApplyTauCorrections(selection.muTau_MC.hadronic_taus,false) : event->taus();

        const auto alltaus = CollectTaus();
        cut(alltaus.size(), "tau_cand");

        const auto signaltaus = CollectSignalTaus() ;

        // First check OS, isolated higgs candidates
        // If no OS candidate, keep any higgs-ish candidates for bkg estimation (don't cut on sign nor isolation)
        auto higgses = FindCompatibleObjects(muons, signaltaus, DeltaR_betweenSignalObjects, Candidate::Higgs,
                                             "H_mu_tau", 0);
        if(!higgses.size())
            higgses = FindCompatibleObjects(allmuons, alltaus, DeltaR_betweenSignalObjects,
                                            Candidate::Higgs, "H_mu_tau");

        cut(higgses.size(), "mu_tau");

        const auto higgsTriggered = config.isDYEmbeddedSample() ? higgses :
                                                                ApplyTriggerMatch(higgses,trigger::hltPaths,false);

        cut(higgsTriggered.size(), "trigger obj match");

        selection.higgs = SelectSemiLeptonicHiggs(higgsTriggered);
        selection.eventType = DoEventCategorization(selection.higgs);

        cut(!config.isDYEmbeddedSample() || selection.eventType == ntuple::EventType::ZTT, "tau match with MC truth");

        CalculateFullEventWeight(selection.higgs);

        if (!config.isMC() || config.isDYEmbeddedSample()){
            selection.pfMET = mvaMetProducer.ComputePFMet(event->pfCandidates(), primaryVertex);
        }
        else
            selection.pfMET = event->metPF();


        const ntuple::MET mvaMet = mvaMetProducer.ComputeMvaMet(selection.higgs, event->pfCandidates(),
                                                                event->jets(), primaryVertex,
                                                                selection.vertices, event->taus());

        const ntuple::MET correctedMET = config.ApplyTauESCorrection()
                ? ApplyTauCorrectionsToMVAMET(mvaMet, correctedTaus) : mvaMet;

        const auto looseJets = CollectLooseJets();
        selection.jetsPt20 = FilterCompatibleObjects(looseJets, selection.higgs, jetID::deltaR_signalObjects);


        selection.jets = CollectJets(selection.jetsPt20);
        selection.bjets_all = CollectBJets(selection.jetsPt20, false, false);
        selection.retagged_bjets = CollectBJets(selection.jetsPt20, config.isMC(), true);


        selection.MET_with_recoil_corrections = ApplyRecoilCorrections(selection.higgs, selection.muTau_MC.resonance,
                                                                       selection.jets.size(), correctedMET);
        return selection;
    }

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
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits)
            < cuts::skim::MuTau::tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits, "looseIso3Hits");
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
//            if(tau_MC.finalStateChargedLeptons.size() == 1
//                    && (*tau_MC.finalStateChargedLeptons.begin())->pdg.Code == particles::mu)
//                final_state.muon = *tau_MC.finalStateChargedLeptons.begin();
            analysis::GenParticlePtrVector tauProducts;
            if (analysis::FindDecayProducts(*tau_MC.origin,analysis::TauMuonicDecay,tauProducts,false)){
                const analysis::GenParticle* muon = tauProducts.at(0);
                final_state.muon = muon;
            }
//            else if(tau_MC.finalStateChargedHadrons.size() >= 1)
            else if(!analysis:: IsLeptonicTau(tau_MC.origin)){
                final_state.tau_jet = &tau_MC;
            }
        }

        if (!final_state.muon || !final_state.tau_jet) return false;
        return true;
    }

    virtual void CalculateTriggerWeights(const analysis::Candidate& higgs) override
    {
        using namespace analysis::Htautau_Summer13::trigger::Run2012ABCD::MuTau;
        const analysis::Candidate& mu = higgs.GetDaughter(analysis::Candidate::Mu);
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        triggerWeights = config.isDYEmbeddedSample() ? CalculateTurnOnCurveData(mu.momentum,tau.momentum) :
                                                     CalculateWeights(mu.momentum, tau.momentum);
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
        double fakeJetToTauWeight = 1;
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        if (config.ApplyJetToTauFakeRate())
            fakeJetToTauWeight = cuts::Htautau_Summer13::jetToTauFakeRateWeight::CalculateJetToTauFakeWeight(tau);
        // first mu, second tau
        fakeWeights.push_back(1);
        fakeWeights.push_back(fakeJetToTauWeight);

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

    virtual void FillFlatTree(const analysis::SelectionResults& /*selection*/) override
    {
        static const float default_value = ntuple::DefaultFloatFillValueForFlatTree();
        static const float default_int_value = ntuple::DefaultIntegerFillValueForFlatTree();

        BaseFlatTreeProducer::FillFlatTree(selection);

        const analysis::Candidate& muon = selection.GetMuon();
        const ntuple::Muon& ntuple_muon = event->muons().at(muon.index);
        const analysis::Candidate& tau = selection.GetTau();

        flatTree->channel() = static_cast<int>(analysis::Channel::MuTau);
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

        const bool muon_matched = analysis::HasMatchWithMCParticle(muon.momentum, selection.muTau_MC.muon,
                                                                   cuts::DeltaR_MC_Match);
        if(muon_matched) {
            const TLorentzVector& momentum = selection.muTau_MC.muon->momentum;
            flatTree->pt_1_MC   () = momentum.Pt()  ;
            flatTree->phi_1_MC  () = momentum.Phi() ;
            flatTree->eta_1_MC  () = momentum.Eta() ;
            flatTree->m_1_MC    () = momentum.M()   ;
            flatTree->pt_1_visible_MC   () = momentum.Pt()  ;
            flatTree->phi_1_visible_MC  () = momentum.Phi() ;
            flatTree->eta_1_visible_MC  () = momentum.Eta() ;
            flatTree->m_1_visible_MC    () = momentum.M()   ;
            flatTree->pdgId_1_MC() = selection.muTau_MC.muon->pdg.ToInteger();
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

        const bool tau_matched = analysis::HasMatchWithMCObject(tau.momentum, selection.muTau_MC.tau_jet,
                                                                cuts::DeltaR_MC_Match);
        if(tau_matched) {
            const TLorentzVector& momentum = selection.muTau_MC.tau_jet->origin->momentum;
            flatTree->pt_2_MC   () = momentum.Pt()  ;
            flatTree->phi_2_MC  () = momentum.Phi() ;
            flatTree->eta_2_MC  () = momentum.Eta() ;
            flatTree->m_2_MC    () = momentum.M()   ;
            const TLorentzVector& visible_momentum = selection.muTau_MC.tau_jet->visibleMomentum;
            flatTree->pt_2_visible_MC   () = visible_momentum.Pt()  ;
            flatTree->phi_2_visible_MC  () = visible_momentum.Phi() ;
            flatTree->eta_2_visible_MC  () = visible_momentum.Eta() ;
            flatTree->m_2_visible_MC    () = visible_momentum.M()   ;
            flatTree->pdgId_2_MC() = selection.muTau_MC.tau_jet->origin->pdg.ToInteger();
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

protected:
    analysis::BaseAnalyzerData baseAnaData;
    analysis::SelectionResults_mutau selection;
    std::shared_ptr<analysis::RecoilCorrectionProducer> recoilCorrectionProducer_mutau;
};

#include "METPUSubtraction/interface/GBRProjectDict.cxx"
