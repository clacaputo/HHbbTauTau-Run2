/*!
 * \file FlatTreeProducer_tautau.C
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
struct SelectionResults_tautau : public SelectionResults {
    finalState::bbTaujetTaujet tauTau_MC;
    CandidatePtr GetLeadingTau() const { return higgs->GetLeadingDaughter(Candidate::Type::Tau); }
    CandidatePtr GetSubleadingTau() const { return higgs->GetSubleadingDaughter(Candidate::Type::Tau); }

    virtual CandidatePtr GetLeg1() const override { return GetLeadingTau(); }
    virtual CandidatePtr GetLeg2() const override { return GetSubleadingTau(); }
    virtual const finalState::bbTauTau& GetFinalStateMC() const override { return tauTau_MC; }
};
} // namespace analysis

class FlatTreeProducer_tautau : public virtual analysis::BaseFlatTreeProducer {
public:
    typedef std::map<analysis::CandidatePtr, analysis::CandidatePtrVector> Higgs_JetsMap;
    typedef std::vector< std::pair<std::string, bool> > TriggerPathVector;
    typedef std::pair<analysis::CandidatePtr, TriggerPathVector> HiggsWithTriggerPath;
    typedef std::map<analysis::CandidatePtr, TriggerPathVector> Higgs_TriggerPathMap;

    FlatTreeProducer_tautau(const std::string& inputFileName, const std::string& outputFileName,
                            const std::string& configFileName, const std::string& _prefix = "none",
                            size_t _maxNumberOfEvents = 0,
                            std::shared_ptr<ntuple::FlatTree> _flatTree = std::shared_ptr<ntuple::FlatTree>())
        : BaseFlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents, _flatTree),
          baseAnaData(*outputFile)
    {
        baseAnaData.getOutputFile().cd();
        if(config.ApplyRecoilCorrection())
            recoilCorrectionProducer_tautau = std::shared_ptr<analysis::RecoilCorrectionProducer>(
                        new analysis::RecoilCorrectionProducer(config.RecoilCorrection_fileCorrectTo_TauTau(),
                                                               config.RecoilCorrection_fileZmmData_TauTau(),
                                                               config.RecoilCorrection_fileZmmMC_TauTau()));
    }

    virtual analysis::BaseAnalyzerData& GetAnaData() override { return baseAnaData; }
    virtual analysis::RecoilCorrectionProducer& GetRecoilCorrectionProducer() override
    {
        return *recoilCorrectionProducer_tautau;
    }

protected:
    virtual analysis::SelectionResults& ApplyBaselineSelection() override
    {
        using namespace analysis;
        using namespace cuts::Htautau_Summer13;
        using namespace cuts::Htautau_Summer13::TauTau;

        selection = SelectionResults_tautau();

        cuts::Cutter cut(&GetAnaData().Selection("event"));
        cut(true, "total");

        cut(FindAnalysisFinalState(selection.tauTau_MC) || !config.RequireSpecificFinalState(), "spec_final_state");
        cut(!config.isDYEmbeddedSample() || GenFilterForZevents(selection.tauTau_MC), "genFilter");

        const auto& selectedTriggerPath = config.isDYEmbeddedSample()
                ? DYEmbedded::trigger::hltPaths : trigger::hltPaths;
        cut(HaveTriggerFired(selectedTriggerPath), "trigger");

        selection.vertices = CollectVertices();
        cut(selection.vertices.size(), "vertex");
        primaryVertex = selection.GetPrimaryVertex();

        const auto electrons_bkg = CollectBackgroundElectrons();
        cut(!electrons_bkg.size(), "no_electron");

        const auto muons_bkg = CollectBackgroundMuons();
        cut(!muons_bkg.size(), "no_muon");

        correctedTaus = config.ApplyTauESCorrection()
                ? ApplyTauCorrections(selection.tauTau_MC.hadronic_taus,false) : event->taus();

        const auto taus = CollectTaus();
        cut(taus.size(), "tau_cand");
        cut(taus.size() >= 2, "at least 2 taus");

        const auto higgses = FindCompatibleObjects(taus, DeltaR_betweenSignalObjects, Candidate::Type::Higgs, "H_2tau");

        cut(higgses.size(), "DeltaR taus");

        const auto looseJets = CollectLooseJets();
        const auto jets = CollectJets(looseJets);

        const Higgs_JetsMap higgs_JetsMap = MatchedHiggsAndJets(higgses, jets);
        const Higgs_JetsMap higgs_looseJetsMap = MatchedHiggsAndJets(higgses, looseJets);

        const auto higgsTriggered = config.isDYEmbeddedSample()
                ? ApplyTriggerMatchForEmbedded(higgs_JetsMap) : ApplyTriggerMatch(higgs_JetsMap,false);
        cut(higgsTriggered.size(), "trigger obj match");

        selectedHiggsWithTriggerPath = SelectFullyHadronicHiggs(higgsTriggered);
        selection.higgs = selectedHiggsWithTriggerPath.first;
        selection.eventType = DoEventCategorization(selection.higgs);

        cut(!config.isDYEmbeddedSample() || selection.eventType == ntuple::EventType::ZTT, "tau match with MC truth");

        CalculateFullEventWeight(selection.higgs);

        if (!config.isMC() || config.isDYEmbeddedSample()){
            selection.pfMET = mvaMetProducer.ComputePFMet(event->pfCandidates(), primaryVertex);
        }
        else
            selection.pfMET = event->metPF();

        const ntuple::MET mvaMet = ComputeMvaMet(selection.higgs, selection.vertices);

        const ntuple::MET correctedMET = config.ApplyTauESCorrection()
                ? ApplyTauCorrectionsToMVAMET(mvaMet, correctedTaus) : mvaMet;

        selection.jetsPt20 = higgs_looseJetsMap.at(selection.higgs);
        selection.jets = higgs_JetsMap.at(selection.higgs);
        selection.bjets_all = CollectBJets(selection.jetsPt20, false, false);
        selection.retagged_bjets = CollectBJets(selection.jetsPt20, config.isMC(), true);


        selection.MET_with_recoil_corrections = ApplyRecoilCorrections(selection.higgs, selection.tauTau_MC.Higgs_TauTau,
                                                                       selection.jets.size(), correctedMET);

        return selection;
    }

    virtual void SelectTau(const analysis::CandidatePtr& tau, analysis::SelectionManager& selectionManager,
                           cuts::Cutter& cut) override
    {
        using namespace cuts::Htautau_Summer13::TauTau;
        using namespace cuts::Htautau_Summer13::TauTau::tauID;
        const ntuple::Tau& object = tau->GetNtupleObject<ntuple::Tau>();

        cut(true, ">0 tau cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex->GetPosition().Z());
        cut(Y(DeltaZ)  < dz, "dz");
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(againstMuonLoose) > againstMuonLoose, "vs_mu_loose");
        cut(X(againstElectronLoose) > againstElectronLoose, "vs_e_loose");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits) <
            cuts::skim::TauTau::tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits, "relaxed_Iso3Hits");
    }

    Higgs_JetsMap MatchedHiggsAndJets(const analysis::CandidatePtrVector& higgses,
                                            const analysis::CandidatePtrVector& jets)
    {
        Higgs_JetsMap higgs_JetsMap;
        for (const analysis::CandidatePtr& higgs : higgses){
            analysis::CandidatePtrVector goodJets;
            for (const analysis::CandidatePtr& jet : jets){
                const double deltaR_1 = jet->GetMomentum().DeltaR(higgs->GetFinalStateDaughters().at(0)->GetMomentum());
                const double deltaR_2 = jet->GetMomentum().DeltaR(higgs->GetFinalStateDaughters().at(1)->GetMomentum());
                if (deltaR_1 > cuts::Htautau_Summer13::jetID::deltaR_signalObjects &&
                        deltaR_2 > cuts::Htautau_Summer13::jetID::deltaR_signalObjects)
                    goodJets.push_back(jet);
            }
            higgs_JetsMap[higgs] = goodJets;
        }
        return higgs_JetsMap;
    }

    Higgs_TriggerPathMap ApplyTriggerMatchForEmbedded(const Higgs_JetsMap& higgs_JetsMap)
    {
        using namespace cuts::Htautau_Summer13;
        const auto firedPaths = CollectPathsForTriggerFired(DYEmbedded::trigger::hltPaths);
        Higgs_TriggerPathMap triggeredHiggses;
        for (const auto& firedPath : firedPaths){
            for (const auto& higgs_iter : higgs_JetsMap){
                const TriggerPathVector::value_type path(firedPath, false);
                triggeredHiggses[higgs_iter.first].push_back(path);
            }
        }
        return triggeredHiggses;
    }

    Higgs_TriggerPathMap ApplyTriggerMatch(const Higgs_JetsMap& higgs_JetsMap, bool useStandardTriggerMatch)
    {
        Higgs_TriggerPathMap triggeredHiggses;
        for (const auto& higgs_iter : higgs_JetsMap) {
            const analysis::CandidatePtr& higgs = higgs_iter.first;
            const analysis::CandidatePtrVector& jets = higgs_iter.second;
            for (const auto& interestingPathIter : cuts::Htautau_Summer13::TauTau::trigger::hltPathsMap) {
                const std::string& interestingPath = interestingPathIter.first;
                const bool jetTriggerRequest = interestingPathIter.second;

                if(!useStandardTriggerMatch && !analysis::HaveTriggerMatched(event->triggerObjects(), interestingPath,
                                                                    *higgs, cuts::Htautau_Summer13::DeltaR_triggerMatch))
                    continue;

                if (useStandardTriggerMatch && !analysis::HaveTriggerMatched(interestingPath,*higgs))
                    continue;

                bool jetMatched = false;
                if(jetTriggerRequest) {
                    for (const auto& jet : jets){
                        if (!useStandardTriggerMatch && analysis::HaveTriggerMatched(event->triggerObjects(),
                                                   interestingPath, *jet, cuts::Htautau_Summer13::DeltaR_triggerMatch)) {
                            jetMatched = true;
                            break;
                        }
                        if (useStandardTriggerMatch && analysis::HaveTriggerMatched(interestingPath, *jet)){
                            jetMatched = true;
                            break;
                        }
                    }
                }

                if(!jetTriggerRequest || jetMatched) {
                    const TriggerPathVector::value_type path(interestingPath, jetTriggerRequest);
                    triggeredHiggses[higgs].push_back(path);
                }
            }
        }
        return triggeredHiggses;
    }

    const Higgs_TriggerPathMap::value_type& SelectFullyHadronicHiggs(const Higgs_TriggerPathMap& higgses)
    {
        if(!higgses.size())
            throw std::runtime_error("no available higgs candidate to select");
        const auto higgsSelector = [&] (const Higgs_TriggerPathMap::value_type& first,
                const Higgs_TriggerPathMap::value_type& second) -> bool
        {
            const ntuple::Tau& first_tau1 = first.first->GetDaughters().at(0)->GetNtupleObject<ntuple::Tau>();
            const ntuple::Tau& first_tau2 = first.first->GetDaughters().at(1)->GetNtupleObject<ntuple::Tau>();
            const ntuple::Tau& second_tau1 = second.first->GetDaughters().at(0)->GetNtupleObject<ntuple::Tau>();
            const ntuple::Tau& second_tau2 = second.first->GetDaughters().at(1)->GetNtupleObject<ntuple::Tau>();
            const double first_iso = std::max(first_tau1.byCombinedIsolationDeltaBetaCorrRaw3Hits,
                                              first_tau2.byCombinedIsolationDeltaBetaCorrRaw3Hits);
            const double second_iso = std::max(second_tau1.byCombinedIsolationDeltaBetaCorrRaw3Hits,
                                               second_tau2.byCombinedIsolationDeltaBetaCorrRaw3Hits);
            return first_iso < second_iso;
        };
        return *std::min_element(higgses.begin(), higgses.end(), higgsSelector);
    }

    bool FindAnalysisFinalState(analysis::finalState::bbTaujetTaujet& final_state)
    {
        const bool base_result = BaseFlatTreeProducer::FindAnalysisFinalState(final_state);
        if(!base_result)
            return base_result;

        if(final_state.hadronic_taus.size() != 2) return false;

        if(final_state.hadronic_taus.at(0).visibleMomentum.Pt() > final_state.hadronic_taus.at(1).visibleMomentum.Pt()) {
            final_state.leading_tau_jet = &final_state.hadronic_taus.at(0);
            final_state.subleading_tau_jet = &final_state.hadronic_taus.at(1);
        } else {
            final_state.leading_tau_jet = &final_state.hadronic_taus.at(1);
            final_state.subleading_tau_jet = &final_state.hadronic_taus.at(0);
        }
        return true;
    }

    virtual void CalculateTriggerWeights(const analysis::CandidatePtr& higgs) override
    {
        using namespace analysis::Htautau_Summer13::trigger::Run2012ABCD::TauTau;
        if(higgs != selectedHiggsWithTriggerPath.first)
            throw analysis::exception("Inconsistet higgs selection");
        const analysis::CandidatePtr& leadTau = higgs->GetLeadingDaughter(analysis::Candidate::Type::Tau);
        const analysis::CandidatePtr& subLeadTau = higgs->GetSubleadingDaughter(analysis::Candidate::Type::Tau);

        bool useDiTauJetWeight = true;
        for(const auto& path : selectedHiggsWithTriggerPath.second) {
            if(!path.second) {
                useDiTauJetWeight = false;
                break;
            }
        }

        if(config.isDYEmbeddedSample())
            triggerWeights = DiTau::CalculateTurnOnCurveData(leadTau->GetMomentum(), subLeadTau->GetMomentum());
        else if(useDiTauJetWeight)
            triggerWeights = DiTauJet::CalculateWeights(leadTau->GetMomentum(), subLeadTau->GetMomentum());
        else
            triggerWeights = DiTau::CalculateWeights(leadTau->GetMomentum(), subLeadTau->GetMomentum());
    }

    virtual void CalculateDMWeights(const analysis::CandidatePtr& higgs) override
    {
        DMweights.clear();
        const analysis::CandidatePtr& leadTau = higgs->GetLeadingDaughter(analysis::Candidate::Type::Tau);
        const analysis::CandidatePtr& subLeadTau = higgs->GetSubleadingDaughter(analysis::Candidate::Type::Tau);
        const ntuple::Tau& leg1 = leadTau->GetNtupleObject<ntuple::Tau>();
        const ntuple::Tau& leg2 = subLeadTau->GetNtupleObject<ntuple::Tau>();
        const double leadWeight = leg1.decayMode == ntuple::tau_id::kOneProng0PiZero ?
                                  cuts::Htautau_Summer13::tauCorrections::DecayModeWeight : 1;
        const double subLeadWeight = leg2.decayMode == ntuple::tau_id::kOneProng0PiZero ?
                                  cuts::Htautau_Summer13::tauCorrections::DecayModeWeight : 1;
        // first leadTau, second subLeadTau
        DMweights.push_back(leadWeight);
        DMweights.push_back(subLeadWeight);
    }

    virtual void CalculateFakeWeights(const analysis::CandidatePtr& /*higgs*/) override
    {
        fakeWeights.clear();
        // first leadTau, second subLeadTau
        fakeWeights.push_back(1);
        fakeWeights.push_back(1);
    }

    virtual void FillFlatTree(const analysis::SelectionResults& /*selection*/) override
    {
        static const float default_value = ntuple::DefaultFloatFillValueForFlatTree();
        static const float default_int_value = ntuple::DefaultIntegerFillValueForFlatTree();

        const analysis::CandidatePtr& leadTau = selection.GetLeadingTau();
        const analysis::CandidatePtr& subLeadTau = selection.GetSubleadingTau();
        const ntuple::Tau& ntuple_tau_leg1 = leadTau->GetNtupleObject<ntuple::Tau>();

        BaseFlatTreeProducer::FillFlatTree(selection);

        flatTree->channel() = static_cast<int>(analysis::Channel::TauTau);
        flatTree->pfRelIso_1() = default_value;
        flatTree->mva_1() = default_int_value;
        flatTree->passid_1() = false;
        flatTree->passiso_1() = false;
        flatTree->decayMode_1() = ntuple_tau_leg1.decayMode;
        flatTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = ntuple_tau_leg1.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        flatTree->againstElectronLooseMVA_1() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau_leg1, 0);
        flatTree->againstElectronMediumMVA_1() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau_leg1, 1);
        flatTree->againstElectronTightMVA_1() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau_leg1, 2);
        flatTree->againstElectronVTightMVA_1() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau_leg1, 3);
        flatTree->againstElectronLoose_1() = ntuple_tau_leg1.againstElectronLoose;
        flatTree->againstElectronMedium_1() = ntuple_tau_leg1.againstElectronMedium;
        flatTree->againstElectronTight_1() = ntuple_tau_leg1.againstElectronTight;
        flatTree->againstMuonLoose_1() = ntuple_tau_leg1.againstMuonLoose;
        flatTree->againstMuonMedium_1() = ntuple_tau_leg1.againstMuonMedium;
        flatTree->againstMuonTight_1() = ntuple_tau_leg1.againstMuonTight;

        const auto leadTau_matches = analysis::FindMatchedObjects(leadTau->GetMomentum(),
                                                                  selection.tauTau_MC.hadronic_taus,
                                                                  cuts::DeltaR_MC_Match);
        if(leadTau_matches.size() != 0) {
            const TLorentzVector& momentum = leadTau_matches.at(0).origin->momentum;
            flatTree->pt_1_MC   () = momentum.Pt()  ;
            flatTree->phi_1_MC  () = momentum.Phi() ;
            flatTree->eta_1_MC  () = momentum.Eta() ;
            flatTree->m_1_MC    () = momentum.M()   ;
            const TLorentzVector& visible_momentum = leadTau_matches.at(0).visibleMomentum;
            flatTree->pt_1_visible_MC   () = visible_momentum.Pt()  ;
            flatTree->phi_1_visible_MC  () = visible_momentum.Phi() ;
            flatTree->eta_1_visible_MC  () = visible_momentum.Eta() ;
            flatTree->m_1_visible_MC    () = visible_momentum.M()   ;
            flatTree->pdgId_1_MC() = leadTau_matches.at(0).origin->pdg.ToInteger();
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

        const auto subLeadTau_matches = analysis::FindMatchedObjects(subLeadTau->GetMomentum(),
                                                                     selection.tauTau_MC.hadronic_taus,
                                                                     cuts::DeltaR_MC_Match);
        if(subLeadTau_matches.size() != 0) {
            const TLorentzVector& momentum = subLeadTau_matches.at(0).origin->momentum;
            flatTree->pt_2_MC   () = momentum.Pt()  ;
            flatTree->phi_2_MC  () = momentum.Phi() ;
            flatTree->eta_2_MC  () = momentum.Eta() ;
            flatTree->m_2_MC    () = momentum.M()   ;
            const TLorentzVector& visible_momentum = subLeadTau_matches.at(0).visibleMomentum;
            flatTree->pt_2_visible_MC   () = visible_momentum.Pt()  ;
            flatTree->phi_2_visible_MC  () = visible_momentum.Phi() ;
            flatTree->eta_2_visible_MC  () = visible_momentum.Eta() ;
            flatTree->m_2_visible_MC    () = visible_momentum.M()   ;
            flatTree->pdgId_2_MC() = subLeadTau_matches.at(0).origin->pdg.ToInteger();
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
    analysis::SelectionResults_tautau selection;
    HiggsWithTriggerPath selectedHiggsWithTriggerPath;
    std::shared_ptr<analysis::RecoilCorrectionProducer> recoilCorrectionProducer_tautau;
};

#include "METPUSubtraction/interface/GBRProjectDict.cxx"
