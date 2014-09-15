/*!
 * \file HtautauFlatTreeProducer.C
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

class HHbbtautau_FlatTreeProducer : public analysis::BaseFlatTreeProducer {
public:
    typedef std::map<analysis::Candidate, analysis::CandidateVector> Higgs_JetsMap;

    HHbbtautau_FlatTreeProducer(const std::string& inputFileName, const std::string& outputFileName,
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
        using namespace cuts::Htautau_Summer13::TauTau;

        if (!FindAnalysisFinalState(tauTau_MC) && config.RequireSpecificFinalState()) return;

        cuts::Cutter cut(&anaData.Selection("event"));
        cut(true, "total");

        cut(HaveTriggerFired(trigger::hltPaths), "trigger");

        const VertexVector vertices = CollectVertices();
        cut(vertices.size(), "vertex");
        primaryVertex = vertices.front();

        const auto electrons_bkg = CollectBackgroundElectrons();
        cut(!electrons_bkg.size(), "no_electron");

        const auto muons_bkg = CollectBackgroundMuons();
        cut(!muons_bkg.size(), "no_muon");

        correctedTaus = config.ApplyTauESCorrection() ? ApplyTauCorrections(tauTau_MC.hadronic_taus,false) : event->taus();

        const auto taus = CollectTaus();
        cut(taus.size(), "tau_cand");
        cut(taus.size() >= 2, "at least 2 taus");

        const auto higgses = FindCompatibleObjects(taus,
                   cuts::Htautau_Summer13::DeltaR_betweenSignalObjects, analysis::Candidate::Higgs, "H_2tau");

        cut(higgses.size(), "DeltaR taus");

        const auto looseJets = CollectLooseJets();
        const auto jets = CollectJets(looseJets);

        const Higgs_JetsMap higgs_JetsMap = MatchedHiggsAndJets(higgses,jets);
        const Higgs_JetsMap higgs_looseJetsMap = MatchedHiggsAndJets(higgses, looseJets);

        const auto higgsTriggered = ApplyTriggerMatch(higgs_JetsMap,false);
        cut(higgsTriggered.size(), "trigger obj match");

        const Candidate higgs = SelectFullyHadronicHiggs(higgsTriggered);

        const ntuple::MET mvaMet = mvaMetProducer.ComputeMvaMet(higgs,event->pfCandidates(),
                                                                event->jets(),primaryVertex,
                                                                vertices,event->taus());

        correctedMET = config.ApplyTauESCorrection() ? ApplyTauCorrectionsToMVAMET(mvaMet, correctedTaus) : mvaMet;

        const auto bjets = CollectBJets(higgs_looseJetsMap.at(higgs), false);
        const auto retagged_bjets = CollectBJets(higgs_looseJetsMap.at(higgs), config.isMC());


        postRecoilMET = ApplyRecoilCorrections(higgs, tauTau_MC.resonance, jets.size(), correctedMET);

        const double m_sv      = CorrectMassBySVfit(higgs, postRecoilMET,1   );
        const double m_sv_Up   = CorrectMassBySVfit(higgs, postRecoilMET,1.03);
        const double m_sv_Down = CorrectMassBySVfit(higgs, postRecoilMET,0.97);

        CalculateFullEventWeight(higgs);

        const ntuple::MET pfMET = config.isMC() ? event->metPF() : mvaMetProducer.ComputePFMet(event->pfCandidates(), primaryVertex);

        FillFlatTree(higgs, m_sv, m_sv_Up, m_sv_Down, higgs_JetsMap.at(higgs), higgs_looseJetsMap.at(higgs),
                     bjets, retagged_bjets, vertices, pfMET);
    }

protected:

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData,
                                          const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::TauTau;
        using namespace cuts::Htautau_Summer13::TauTau::tauID;

        cuts::Cutter cut(objectSelector);
        const ntuple::Tau& object = correctedTaus.at(id);

        cut(true, ">0 tau cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(againstMuonLoose) > againstMuonLoose, "vs_mu_loose");
        cut(X(againstElectronLoose) > againstElectronLoose, "vs_e_loose");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits) <
            cuts::skim::TauTau::tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits, "looseIso3Hits");

        return analysis::Candidate(analysis::Candidate::Tau, id, object);
    }

    analysis::CandidateVector ApplyTauFullSelection(const analysis::CandidateVector& higgses)
    {
        using namespace analysis;
        using namespace cuts::Htautau_Summer13::TauTau::tauID;
        CandidateVector result;
        for(const Candidate& higgs : higgses) {
            const Candidate& subleading_tau = higgs.GetSubleadingDaughter(Candidate::Tau);
            const ntuple::Tau& ntuple_subleadingTau = correctedTaus.at(subleading_tau.index);
            if(ntuple_subleadingTau.againstElectronLooseMVA3 > againstElectronLooseMVA3)
                result.push_back(higgs);
        }
        return result;
    }

    Higgs_JetsMap MatchedHiggsAndJets(const analysis::CandidateVector& higgses,
                                            const analysis::CandidateVector& jets)
    {
        Higgs_JetsMap higgs_JetsMap;
        for (const analysis::Candidate& higgs : higgses){
            analysis::CandidateVector goodJets;
            for (const analysis::Candidate& jet : jets){
                const double deltaR_1 = jet.momentum.DeltaR(higgs.finalStateDaughters.at(0).momentum);
                const double deltaR_2 = jet.momentum.DeltaR(higgs.finalStateDaughters.at(1).momentum);
                if (deltaR_1 > cuts::Htautau_Summer13::jetID::deltaR_signalObjects &&
                        deltaR_2 > cuts::Htautau_Summer13::jetID::deltaR_signalObjects)
                    goodJets.push_back(jet);
            }
            higgs_JetsMap[higgs] = goodJets;
        }
        return higgs_JetsMap;
    }

    analysis::CandidateVector ApplyTriggerMatch(const Higgs_JetsMap& higgs_JetsMap, bool useStandardTriggerMatch)
    {
        analysis::CandidateVector triggeredHiggses;
        for (const auto& higgs_iter : higgs_JetsMap) {
            const analysis::Candidate& higgs = higgs_iter.first;
            const analysis::CandidateVector& jets = higgs_iter.second;
            for (const auto& interestingPathIter : cuts::Htautau_Summer13::TauTau::trigger::hltPathsMap) {
                const std::string& interestingPath = interestingPathIter.first;
                const bool jetTriggerRequest = interestingPathIter.second;

                if(!useStandardTriggerMatch && !analysis::HaveTriggerMatched(event->triggerObjects(), interestingPath,
                                                                    higgs, cuts::Htautau_Summer13::DeltaR_triggerMatch))
                    continue;

                if (useStandardTriggerMatch && !analysis::HaveTriggerMatched(*event,interestingPath,higgs))
                    continue;

                bool jetMatched = false;
                if(jetTriggerRequest) {
                    for (const auto& jet : jets){
                        if (!useStandardTriggerMatch && analysis::HaveTriggerMatched(event->triggerObjects(),
                                                   interestingPath, jet, cuts::Htautau_Summer13::DeltaR_triggerMatch)) {
                            jetMatched = true;
                            break;
                        }
                        if (useStandardTriggerMatch && analysis::HaveTriggerMatched(*event, interestingPath, jet)){
                            jetMatched = true;
                            break;
                        }
                    }
                }

                if(!jetTriggerRequest || jetMatched) {
                    triggeredHiggses.push_back(higgs);
                    break;
                }
            }
        }
        return triggeredHiggses;
    }

    const analysis::Candidate& SelectFullyHadronicHiggs(const analysis::CandidateVector& higgses)
    {
        if(!higgses.size())
            throw std::runtime_error("no available higgs candidate to select");
        const auto higgsSelector = [&] (const analysis::Candidate& first, const analysis::Candidate& second) -> bool
        {
            const ntuple::Tau& first_tau1 = correctedTaus.at(first.daughters.at(0).index);
            const ntuple::Tau& first_tau2 = correctedTaus.at(first.daughters.at(1).index);
            const ntuple::Tau& second_tau1 = correctedTaus.at(second.daughters.at(0).index);
            const ntuple::Tau& second_tau2 = correctedTaus.at(second.daughters.at(1).index);
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

    virtual void CalculateTriggerWeights(const analysis::Candidate& higgs) override
    {
        using namespace analysis::Htautau_Summer13::trigger::Run2012ABCD::TauTau;
        const analysis::Candidate& leadTau = higgs.GetLeadingDaughter(analysis::Candidate::Tau);
        const analysis::Candidate& subLeadTau = higgs.GetSubleadingDaughter(analysis::Candidate::Tau);
        triggerWeights = CalculateWeights(leadTau.momentum, subLeadTau.momentum);
    }

    virtual void CalculateDMWeights(const analysis::Candidate& higgs) override
    {
        DMweights.clear();
        const analysis::Candidate& leadTau = higgs.GetLeadingDaughter(analysis::Candidate::Tau);
        const analysis::Candidate& subLeadTau = higgs.GetSubleadingDaughter(analysis::Candidate::Tau);
        const ntuple::Tau& leg1 = correctedTaus.at(leadTau.index);
        const ntuple::Tau& leg2 = correctedTaus.at(subLeadTau.index);
        const double leadWeight = leg1.decayMode == ntuple::tau_id::kOneProng0PiZero ?
                                  cuts::Htautau_Summer13::tauCorrections::DecayModeWeight : 1;
        const double subLeadWeight = leg2.decayMode == ntuple::tau_id::kOneProng0PiZero ?
                                  cuts::Htautau_Summer13::tauCorrections::DecayModeWeight : 1;
        // first leadTau, second subLeadTau
        DMweights.push_back(leadWeight);
        DMweights.push_back(subLeadWeight);
    }

    void FillFlatTree(const analysis::Candidate& higgs, double m_sv,
                      double m_sv_Up, double m_sv_Down,
                      const analysis::CandidateVector& jets, const analysis::CandidateVector& jetsPt20,
                      const analysis::CandidateVector& bjets,  const analysis::CandidateVector& retagged_bjets,
                      const analysis::VertexVector& vertices, const ntuple::MET& pfMET)
    {
        static const float default_value = ntuple::DefaultFloatFillValueForFlatTree();
        static const float default_int_value = ntuple::DefaultIntegerFillValueForFlatTree();

        const analysis::Candidate& leadTau = higgs.GetLeadingDaughter(analysis::Candidate::Tau);
        const analysis::Candidate& subLeadTau = higgs.GetSubleadingDaughter(analysis::Candidate::Tau);
        const ntuple::Tau& ntuple_tau_leg1 = correctedTaus.at(leadTau.index);

        BaseFlatTreeProducer::FillFlatTree(higgs, m_sv, m_sv_Up, m_sv_Down, jets, jetsPt20,
                                           bjets, retagged_bjets, vertices, leadTau, subLeadTau,pfMET, tauTau_MC);

        flatTree->channel() = static_cast<int>(ntuple::Channel::TauTau);
        flatTree->pfRelIso_1() = default_value;
        flatTree->mva_1() = default_int_value;
        flatTree->passid_1() = false;
        flatTree->passiso_1() = false;
        flatTree->decayMode_1() = ntuple_tau_leg1.decayMode;
        flatTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = ntuple_tau_leg1.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        flatTree->againstElectronLooseMVA_1() = ComputeAntiElectronMVA3New(ntuple_tau_leg1.againstElectronMVA3category,
                                                                           ntuple_tau_leg1.againstElectronMVA3raw, 0);
        flatTree->againstElectronMediumMVA_1() = ComputeAntiElectronMVA3New(ntuple_tau_leg1.againstElectronMVA3category,
                                                                            ntuple_tau_leg1.againstElectronMVA3raw, 1);
        flatTree->againstElectronTightMVA_1() = ComputeAntiElectronMVA3New(ntuple_tau_leg1.againstElectronMVA3category,
                                                                           ntuple_tau_leg1.againstElectronMVA3raw, 2);
        flatTree->againstElectronVTightMVA_1() = ComputeAntiElectronMVA3New(ntuple_tau_leg1.againstElectronMVA3category,
                                                                            ntuple_tau_leg1.againstElectronMVA3raw, 3);
        flatTree->againstElectronLoose_1() = ntuple_tau_leg1.againstElectronLoose;
        flatTree->againstElectronMedium_1() = ntuple_tau_leg1.againstElectronMedium;
        flatTree->againstElectronTight_1() = ntuple_tau_leg1.againstElectronTight;
        flatTree->againstMuonLoose_1() = ntuple_tau_leg1.againstMuonLoose;
        flatTree->againstMuonMedium_1() = ntuple_tau_leg1.againstMuonMedium;
        flatTree->againstMuonTight_1() = ntuple_tau_leg1.againstMuonTight;

        const analysis::VisibleGenObjectVector leadTau_matches =
                analysis::FindMatchedObjects(leadTau.momentum, tauTau_MC.hadronic_taus, cuts::DeltaR_MC_Match);
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

        const analysis::VisibleGenObjectVector subLeadTau_matches =
                analysis::FindMatchedObjects(subLeadTau.momentum, tauTau_MC.hadronic_taus, cuts::DeltaR_MC_Match);

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

private:
    analysis::BaseAnalyzerData anaData;
    analysis::finalState::bbTaujetTaujet tauTau_MC;

};

#include "METPUSubtraction/interface/GBRProjectDict.cxx"
