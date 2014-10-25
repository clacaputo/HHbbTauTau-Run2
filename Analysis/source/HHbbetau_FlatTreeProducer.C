/*!
 * \file HetauFlatTreeProducer.C
 * \brief Generate flat-tree for H->tautau->e_taujet analysis using looser selection.
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
struct SelectionResults_etau : public SelectionResults {
    finalState::bbETaujet eTau_MC;
    const Candidate& GetElectron() const { return higgs.GetDaughter(Candidate::Electron); }
    const Candidate& GetTau() const { return higgs.GetDaughter(Candidate::Tau); }

    virtual const Candidate& GetLeg1() const override { return GetElectron(); }
    virtual const Candidate& GetLeg2() const override { return GetTau(); }
    virtual const finalState::bbTauTau& GetFinalStateMC() const override { return eTau_MC; }
};
} // namespace analysis

class HHbbetau_FlatTreeProducer : public virtual analysis::BaseFlatTreeProducer {
public:
    HHbbetau_FlatTreeProducer(const std::string& inputFileName, const std::string& outputFileName,
                              const std::string& configFileName, const std::string& _prefix = "none",
                              size_t _maxNumberOfEvents = 0,
                              std::shared_ptr<ntuple::FlatTree> _flatTree = std::shared_ptr<ntuple::FlatTree>())
        : BaseFlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents, _flatTree),
          baseAnaData(*outputFile)
    {
        baseAnaData.getOutputFile().cd();
    }

    virtual analysis::BaseAnalyzerData& GetAnaData() override { return baseAnaData; }

protected:
    virtual analysis::SelectionResults& ApplyBaselineSelection() override
    {
        using namespace analysis;
        using namespace cuts::Htautau_Summer13;
        using namespace cuts::Htautau_Summer13::ETau;

        selection = SelectionResults_etau();

        cuts::Cutter cut(&GetAnaData().Selection("event"));
        cut(true, "total");

        cut(FindAnalysisFinalState(selection.eTau_MC) || !config.RequireSpecificFinalState(), "spec_final_state");
        cut(!config.isDYEmbeddedSample() || GenFilterForZevents(selection.eTau_MC), "genFilter");

        const auto& selectedTriggerPath = config.isDYEmbeddedSample()
                ? DYEmbedded::trigger::hltPaths : trigger::hltPaths;
        cut(HaveTriggerFired(selectedTriggerPath), "trigger");

        selection.vertices = CollectVertices();
        cut(selection.vertices.size(), "vertex");
        primaryVertex = selection.vertices.front();

        const auto z_electrons = CollectZelectrons();
        const auto z_electron_candidates = FindCompatibleObjects(z_electrons, ZeeVeto::deltaR,Candidate::Z, "Z_e_e", 0);
        cut(!z_electron_candidates.size(), "z_ee_veto");

        const auto muons_bkg = CollectBackgroundMuons();
        cut(!muons_bkg.size(), "no_muon");

        const auto electrons = CollectSignalElectrons();

        const auto electrons_bkg = CollectBackgroundElectrons();
        const bool have_bkg_electron = electrons_bkg.size() > 1 || electrons.size() > 1 ||
                ( electrons_bkg.size() == 1 && electrons.size() == 1 && electrons_bkg.front() != electrons.front() );
        cut(!have_bkg_electron, "no_bkg_electron");

        const auto allelectrons = CollectElectrons();
        cut(allelectrons.size(), "electron_cand");

        correctedTaus = config.ApplyTauESCorrection()
                ? ApplyTauCorrections(selection.eTau_MC.hadronic_taus,false) : event->taus();

        const auto alltaus = CollectTaus();
        cut(alltaus.size(), "tau_cand");

        const auto signaltaus = CollectSignalTaus() ;

        // First check OS, isolated higgs candidates
        // If no OS candidate, keep any higgs-ish candidates for bkg estimation (don't cut on sign nor isolation)
        auto higgses = FindCompatibleObjects(electrons, signaltaus, DeltaR_betweenSignalObjects,
                                             Candidate::Higgs, "H_e_tau", 0);
        if (!higgses.size())
          higgses = FindCompatibleObjects(allelectrons, alltaus, DeltaR_betweenSignalObjects,
                                          Candidate::Higgs, "H_e_tau");

        cut(higgses.size(), "e_tau");

        const auto higgsTriggered = config.isDYEmbeddedSample() ? higgses :
                                                                ApplyTriggerMatch(higgses,trigger::hltPaths,false);

        cut(higgsTriggered.size(), "trigger obj match");


        selection.higgs = SelectSemiLeptonicHiggs(higgsTriggered);

        cut(!config.isDYEmbeddedSample() || MatchTausFromHiggsWithGenTaus(selection.higgs, selection.eTau_MC),
            "tau match with MC truth");

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


        selection.MET_with_recoil_corrections = ApplyRecoilCorrections(selection.higgs, selection.eTau_MC.resonance,
                                                                       selection.jets.size(), correctedMET);
        return selection;
    }

    virtual analysis::Candidate SelectElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                               root_ext::AnalyzerData& _anaData,
                                               const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::electronID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event->electrons().at(id);
        const analysis::Candidate electron(analysis::Candidate::Electron, id, object);

        cut(true, ">0 ele cand");
        cut(X(pt) > pt, "pt");
        const double eta = std::abs( X(eta) );
        cut(eta < eta_high, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        cut(X(missingHits) < missingHits, "missingHits");
        cut(X(hasMatchedConversion) == hasMatchedConversion, "conversion");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(ele_vertex,primaryVertex.position,electron.momentum); // same as dB
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[eta_index], "mva");
        //cut(X(pfRelIso) < cuts::skim::MuTau::pFRelIso , "pFRelIso");

//        const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
//        cut(Y(haveTriggerMatch, 2, -0.5, 1.5), "triggerMatch");
        return analysis::Candidate(analysis::Candidate::Electron, id, object);
    }

    virtual analysis::Candidate SelectSignalElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                               root_ext::AnalyzerData& _anaData,
                                               const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::electronID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event->electrons().at(id);
        const analysis::Candidate electron(analysis::Candidate::Electron, id, object);

        cut(true, ">0 ele cand");
        cut(X(pt) > pt, "pt");
        const double eta = std::abs( X(eta) );
        cut(eta < eta_high, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        cut(X(missingHits) < missingHits, "missingHits");
        cut(X(hasMatchedConversion) == hasMatchedConversion, "conversion");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(ele_vertex,primaryVertex.position,electron.momentum); // same as dB
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[eta_index], "mva");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");
//        const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
//        cut(Y(haveTriggerMatch, 2, -0.5, 1.5), "triggerMatch");
        return analysis::Candidate(analysis::Candidate::Electron, id, object);
    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData,
                                          const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::tauID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Tau& object = correctedTaus.at(id);

        cut(true, ">0 tau cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(againstMuonLoose) > againstMuonLoose, "vs_mu_loose");
        const bool againstElectron =  againstElectronMediumMVA3_Custom(object);
        cut(Y(againstElectron), "vs_e_mediumMVA");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits) <
            cuts::skim::ETau::tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits, "looseIso3Hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");

        const analysis::Candidate tau(analysis::Candidate::Tau, id, object);

        return tau;
    }

    virtual analysis::Candidate SelectSignalTau(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData,
                                          const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::tauID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Tau& object = correctedTaus.at(id);

        cut(true, ">0 tau cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(againstMuonLoose) > againstMuonLoose, "vs_mu_loose");
        const bool againstElectron =  againstElectronMediumMVA3_Custom(object);
        cut(Y(againstElectron), "vs_e_mediumMVA");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits) < byCombinedIsolationDeltaBetaCorrRaw3Hits, "looseIso3Hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");

        return analysis::Candidate(analysis::Candidate::Tau, id, object);
    }

    analysis::CandidateVector CollectZelectrons()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> analysis::Candidate
            { return SelectZelectron(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<analysis::Candidate>("z_electrons", base_selector, event->electrons().size());
    }

    virtual analysis::Candidate SelectZelectron(size_t id, cuts::ObjectSelector* objectSelector,
                                                root_ext::AnalyzerData& _anaData,
                                                const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::ETau::ZeeVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event->electrons().at(id);
        const analysis::Candidate electron(analysis::Candidate::Electron, id, object);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(ele_vertex,primaryVertex.position,electron.momentum); // same as dB
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(pfRelIso) < pfRelIso, "pFRelIso");
        const size_t eta_index = std::abs(object.eta) <= barrel_eta_high ? barrel_index : endcap_index;
        cut(X(sigmaIEtaIEta) < sigma_ieta_ieta[eta_index], "sigmaIetaIeta");
        cut(X(deltaEtaTrkSC) < delta_eta[eta_index], "deltaEtaSC");
        cut(X(deltaPhiTrkSC) < delta_phi[eta_index], "deltaPhiSC");
        cut(X(eSuperClusterOverP) < HoverE[eta_index], "HoverE");

        return electron;
    }

    bool againstElectronMediumMVA3_Custom(const ntuple::Tau& tau)
    {
        using namespace cuts::Htautau_Summer13::ETau::tauID;
        const int icut = std::round(tau.againstElectronMVA3category);
        if(icut < 0) return false;
		const size_t ucut = (size_t)icut;
        if(ucut >= againstElectronMediumMVA3_customValues.size()) return true;
        return tau.againstElectronMVA3raw > againstElectronMediumMVA3_customValues.at(ucut);
    }

    bool FindAnalysisFinalState(analysis::finalState::bbETaujet& final_state)
    {
        const bool base_result = BaseFlatTreeProducer::FindAnalysisFinalState(final_state);
        if(!base_result)
            return base_result;

        for (const analysis::VisibleGenObject& tau_MC : final_state.taus) {
            if(tau_MC.finalStateChargedLeptons.size() == 1
                    && tau_MC.finalStateChargedLeptons.at(0)->pdg.Code == particles::e)
                final_state.electron = tau_MC.finalStateChargedLeptons.at(0);
            else if(tau_MC.finalStateChargedHadrons.size() >= 1)
                final_state.tau_jet = &tau_MC;
        }

        if (!final_state.electron || !final_state.tau_jet) return false;
        return true;
    }

    virtual void CalculateTriggerWeights(const analysis::Candidate& higgs) override
    {
        using namespace analysis::Htautau_Summer13::trigger::Run2012ABCD::ETau;
        const analysis::Candidate& ele = higgs.GetDaughter(analysis::Candidate::Electron);
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        triggerWeights = config.isDYEmbeddedSample() ? CalculateTurnOnCurveData(ele.momentum, tau.momentum) :
                                                     CalculateWeights(ele.momentum, tau.momentum);
    }

    virtual void CalculateIsoWeights(const analysis::Candidate& higgs) override
    {
        using namespace cuts::Htautau_Summer13::ETau::electronISOscaleFactor;
        IsoWeights.clear();
        const analysis::Candidate& ele = higgs.GetDaughter(analysis::Candidate::Electron);
        const double ele_pt = ele.momentum.Pt(), ele_eta = std::abs(ele.momentum.Eta());
        if(ele_pt < pt.at(0))
            throw std::runtime_error("No information about ISO. Electron pt is too small");
        const size_t pt_bin = ele_pt < pt.at(1) ? 0 : 1;
        if(ele_eta >= eta.at(1))
            throw std::runtime_error("No information about ISO. Electron eta is too big");
        const size_t eta_bin = ele_eta < eta.at(0) ? 0 : 1;
        const double scale = scaleFactors.at(pt_bin).at(eta_bin);
        // first ele, second tau
        IsoWeights.push_back(scale);
        IsoWeights.push_back(1);
    }

    virtual void CalculateIdWeights(const analysis::Candidate& higgs) override
    {
        using namespace cuts::Htautau_Summer13::ETau::electronIDscaleFactor;
        IDweights.clear();
        const analysis::Candidate& ele = higgs.GetDaughter(analysis::Candidate::Electron);
        const double ele_pt = ele.momentum.Pt(), ele_eta = std::abs(ele.momentum.Eta());
        if(ele_pt < pt.at(0))
            throw std::runtime_error("No information about ID. Electron pt is too small");
        const size_t pt_bin = ele_pt < pt.at(1) ? 0 : 1;
        if(ele_eta >= eta.at(1))
            throw std::runtime_error("No information about ID. Electron eta is too big");
        const size_t eta_bin = ele_eta < eta.at(0) ? 0 : 1;
        const double scale = scaleFactors.at(pt_bin).at(eta_bin);
        // first ele, second tau
        IDweights.push_back(scale);
        IDweights.push_back(1);
    }

    virtual void CalculateDMWeights(const analysis::Candidate& higgs) override
    {
        DMweights.clear();
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        const ntuple::Tau& tau_leg = correctedTaus.at(tau.index);
        const double weight = tau_leg.decayMode == ntuple::tau_id::kOneProng0PiZero ?
                                  cuts::Htautau_Summer13::tauCorrections::DecayModeWeight : 1;
        // first electron, second tau
        DMweights.push_back(1);
        DMweights.push_back(weight);
    }

    virtual void FillFlatTree(const analysis::SelectionResults& /*selection*/) override
    {
        static const float default_value = ntuple::DefaultFloatFillValueForFlatTree();
        static const float default_int_value = ntuple::DefaultIntegerFillValueForFlatTree();

        const analysis::Candidate& electron = selection.GetElectron();
        const ntuple::Electron& ntuple_electron = event->electrons().at(electron.index);
        const analysis::Candidate& tau = selection.GetTau();

        BaseFlatTreeProducer::FillFlatTree(selection);

        flatTree->channel() = static_cast<int>(analysis::Channel::ETau);
        flatTree->pfRelIso_1() = ntuple_electron.pfRelIso;
        flatTree->mva_1() = ntuple_electron.mvaPOGNonTrig;
        flatTree->passid_1() = true;
        flatTree->passiso_1() = ntuple_electron.pfRelIso < cuts::Htautau_Summer13::ETau::electronID::pFRelIso;
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

        const bool electron_matched = analysis::HasMatchWithMCParticle(electron.momentum, selection.eTau_MC.electron,
                                                                       cuts::DeltaR_MC_Match);
        if(electron_matched) {
            const TLorentzVector& momentum = selection.eTau_MC.electron->momentum;
            flatTree->pt_1_MC   () = momentum.Pt()  ;
            flatTree->phi_1_MC  () = momentum.Phi() ;
            flatTree->eta_1_MC  () = momentum.Eta() ;
            flatTree->m_1_MC    () = momentum.M()   ;
            flatTree->pt_1_visible_MC   () = momentum.Pt()  ;
            flatTree->phi_1_visible_MC  () = momentum.Phi() ;
            flatTree->eta_1_visible_MC  () = momentum.Eta() ;
            flatTree->m_1_visible_MC    () = momentum.M()   ;
            flatTree->pdgId_1_MC() = selection.eTau_MC.electron->pdg.ToInteger();
        } else {
            flatTree->pt_1_MC   () = default_value ;
            flatTree->phi_1_MC  () = default_value ;
            flatTree->eta_1_MC  () = default_value ;
            flatTree->m_1_MC    () = default_value ;
            flatTree->pt_1_visible_MC   () = default_value ;
            flatTree->phi_1_visible_MC  () = default_value ;
            flatTree->eta_1_visible_MC  () = default_value ;
            flatTree->m_1_visible_MC    () = default_value ;
            flatTree->pdgId_1_MC() = default_int_value;
        }

        const bool tau_matched = analysis::HasMatchWithMCObject(tau.momentum, selection.eTau_MC.tau_jet,
                                                                cuts::DeltaR_MC_Match);
        if(tau_matched) {
            const TLorentzVector& momentum = selection.eTau_MC.tau_jet->origin->momentum;
            flatTree->pt_2_MC   () = momentum.Pt()  ;
            flatTree->phi_2_MC  () = momentum.Phi() ;
            flatTree->eta_2_MC  () = momentum.Eta() ;
            flatTree->m_2_MC    () = momentum.M()   ;
            const TLorentzVector& visible_momentum = selection.eTau_MC.tau_jet->visibleMomentum;
            flatTree->pt_2_visible_MC   () = visible_momentum.Pt()  ;
            flatTree->phi_2_visible_MC  () = visible_momentum.Phi() ;
            flatTree->eta_2_visible_MC  () = visible_momentum.Eta() ;
            flatTree->m_2_visible_MC    () = visible_momentum.M()   ;
            flatTree->pdgId_2_MC() = selection.eTau_MC.tau_jet->origin->pdg.ToInteger();
        } else {
            flatTree->pt_2_MC   () = default_value ;
            flatTree->phi_2_MC  () = default_value ;
            flatTree->eta_2_MC  () = default_value ;
            flatTree->m_2_MC    () = default_value ;
            flatTree->pt_2_visible_MC   () = default_value ;
            flatTree->phi_2_visible_MC  () = default_value ;
            flatTree->eta_2_visible_MC  () = default_value ;
            flatTree->m_2_visible_MC    () = default_value ;
            flatTree->pdgId_2_MC() = default_int_value;
        }

        flatTree->Fill();
    }

protected:
    analysis::BaseAnalyzerData baseAnaData;
    analysis::SelectionResults_etau selection;
};

#include "METPUSubtraction/interface/GBRProjectDict.cxx"
