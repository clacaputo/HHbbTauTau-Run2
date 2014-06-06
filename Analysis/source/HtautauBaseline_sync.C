/*!
 * \file HtautauBaseline_sync.C
 * \brief Generate sync-tree for Htautau analysis using baseline selection.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-05 created
 */

#include "../include/H_BaseAnalyzer.h"
#include "../include/SyncTree.h"

class BaselineAnalyzerData : public analysis::BaseAnalyzerData {
public:
    BaselineAnalyzerData(TFile& outputFile) : BaseAnalyzerData(outputFile) {}

};

class HtautauBaseline_sync : public analysis::H_BaseAnalyzer {
public:
    HtautauBaseline_sync(const std::string& inputFileName, const std::string& outputFileName,
                          const std::string& _prefix = "none", Long64_t _maxNumberOfEvents = 0,
                          bool _useMCtruth = false, const std::string& reweightFileName = "none")
        : H_BaseAnalyzer(inputFileName,outputFileName,_prefix,_maxNumberOfEvents,_useMCtruth, reweightFileName),
          anaData(*outputFile)
    {
        anaData.getOutputFile().cd();
    }

    virtual ~HtautauBaseline_sync() override
    {
        anaData.getOutputFile().cd();
        syncTree.Write();
    }

protected:
    typedef std::map<analysis::Candidate, analysis::CandidateVector> Higgs_JetsMap;
    virtual analysis::BaseAnalyzerData& GetAnaData() override { return anaData; }

    virtual void ProcessEvent() override
    {
        H_BaseAnalyzer::ProcessEvent();
        using namespace analysis;
        using namespace cuts::Htautau_Summer13;
        using namespace cuts::Htautau_Summer13::TauTau;
        finalState::TaujetTaujet tauTau;
        if (useMCtruth && !FindAnalysisFinalState(tauTau, false)) return;

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

        //ApplyTauCorrections(tauTau, event.metMVAtauTau(), false);

        //ApplyTauCorrections(tauTau,false);
        correctedTaus = event.taus(); //no tau corrections

        const auto taus = CollectTaus();
        cut(taus.size(), "tau_cand");
        cut(taus.size() >= 2, "at least 2 taus");

        const auto higgses = FindCompatibleObjects(taus,
                   cuts::Htautau_Summer13::DeltaR_betweenSignalObjects, analysis::Candidate::Higgs, "H_2tau");

        cut(higgses.size(), "DeltaR taus");

        const auto higgses_tautau = ApplyTauFullSelection(higgses);

        cut(higgses_tautau.size(), "tau_tau selection");

        const auto jetsPt20 = CollectJetsPt20();
        const auto jets = CollectJets();

        const Higgs_JetsMap higgs_JetsMap = MatchedHiggsAndJets(higgses_tautau,jets);

        const auto higgsTriggered = ApplyTriggerMatch(higgs_JetsMap,false);
        cut(higgsTriggered.size(), "trigger obj match");

        const Candidate higgs_withoutTauCorrections = SelectFullyHadronicHiggs(higgsTriggered);

        const ntuple::MET mvaMet = mvaMetProducer.ComputeMvaMet(higgs_withoutTauCorrections,event.pfCandidates(),event.jets(),primaryVertex,
                                                                vertices,event.taus());

        ApplyTauCorrections(tauTau,false);
        const Candidate higgs = ApplyTauCorrectionsToMVAMETandHiggs(mvaMet,higgs_withoutTauCorrections);

        const auto bjets = CollectBJets(higgs);

        //ApplyPostRecoilCorrections(higgs, tauTau.resonance, higgs_JetsMap.at(higgs_withoutTauCorrections).size());
        postRecoilMET = correctedMET; //with tau corrections

        const Candidate higgs_sv = CorrectMassBySVfit(higgs, postRecoilMET,1);
        const Candidate higgs_sv_up = CorrectMassBySVfit(higgs, postRecoilMET,1.03);
        const Candidate higgs_sv_down = CorrectMassBySVfit(higgs, postRecoilMET,0.97);

        CalculateFullEventWeight(higgs_sv);

        FillSyncTree(higgs, higgs_sv, higgs_sv_up, higgs_sv_down, higgs_JetsMap.at(higgs_withoutTauCorrections), jetsPt20, bjets, vertices);

        //postRecoilMET = mvaMet;
//        postRecoilMET = correctedMET; //with tau corrections

//        CalculateFullEventWeight(higgs);
//        FillSyncTree(higgs, higgs, higgs, higgs, higgs_JetsMap.at(higgs_withoutTauCorrections), jetsPt20, bjets, vertices);
    }

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
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(againstMuonLoose) > againstMuonLoose, "vs_mu_loose");
        cut(X(againstElectronLoose) > againstElectronLoose, "vs_e_loose");
        cut(X(byMediumCombinedIsolationDeltaBetaCorr3Hits) > byMediumCombinedIsolationDeltaBetaCorr3Hits, "mediumIso3Hits");

//        const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
        const analysis::Candidate tau(analysis::Candidate::Tau, id, object);
//        const bool haveTriggerMatch = analysis::HaveTriggerMatched(event.triggerObjects(), trigger::hltPaths,tau);
//        cut(Y(haveTriggerMatch), "triggerMatch");

        return tau;
    }

    analysis::CandidateVector ApplyTauFullSelection(const analysis::CandidateVector& higgses)
    {
        using namespace analysis;
        using namespace cuts::Htautau_Summer13::TauTau::tauID;
        CandidateVector result;
        for(const Candidate& higgs : higgses) {
//            const Candidate& leading_tau = higgs.GetLeadingDaughter(Candidate::Tau);
//            std::cout << "leading tau mom " << leading_tau.momentum << std::endl;
            const Candidate& subleading_tau = higgs.GetSubleadingDaughter(Candidate::Tau);
//            std::cout << "subleading tau mom " << subleading_tau.momentum << std::endl;
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

                if(!useStandardTriggerMatch && !analysis::HaveTriggerMatched(event.triggerObjects(), interestingPath, higgs))
                    continue;

                if (useStandardTriggerMatch && !analysis::HaveTriggerMatched(event,interestingPath,higgs))
                    continue;

                bool jetMatched = false;
                if(jetTriggerRequest) {
                    for (const auto& jet : jets){
                        if (!useStandardTriggerMatch && analysis::HaveTriggerMatched(event.triggerObjects(), interestingPath, jet)){
                            jetMatched = true;
                            break;
                        }
                        if (useStandardTriggerMatch && analysis::HaveTriggerMatched(event, interestingPath, jet)){
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

//    bool FindAnalysisFinalState(analysis::finalState::TaujetTaujet& final_state)
//    {
//        if(!H_BaseAnalyzer::FindAnalysisFinalState(final_state))
//            return false;

//        unsigned n_hadronic_taus = 0;

//        if(final_state.taus.size() != 2)
//            throw std::runtime_error("bad tautau MC final state");

//        for (const analysis::GenParticle* tau_MC : final_state.taus) {
//            analysis::GenParticlePtrVector TauProducts;
//            if (!analysis::FindDecayProducts(*tau_MC, analysis::TauMuonicDecay, TauProducts)
//                    && !analysis::FindDecayProducts(*tau_MC, analysis::TauElectronDecay, TauProducts))
//                ++n_hadronic_taus;
//        }

//        if (n_hadronic_taus != 2) return false;
//        if(final_state.taus.at(0)->momentum.Pt() > final_state.taus.at(1)->momentum.Pt()) {
//            final_state.leading_tau_jet = final_state.taus.at(0);
//            final_state.subleading_tau_jet = final_state.taus.at(1);
//        } else {
//            final_state.leading_tau_jet = final_state.taus.at(1);
//            final_state.subleading_tau_jet = final_state.taus.at(0);
//        }
//        return true;
//    }

    virtual void CalculateTriggerWeights(const analysis::Candidate& higgs) override
    {
        triggerWeights.clear();
        const analysis::Candidate& leadTau = higgs.GetLeadingDaughter(analysis::Candidate::Tau);
        const analysis::Candidate& subLeadTau = higgs.GetSubleadingDaughter(analysis::Candidate::Tau);
        analysis::Htautau_Summer13::TriggerEfficiency efficiency;
        const double eff_data_leadTau = efficiency.eff2012IsoParkedTau19fb_Simone(leadTau.momentum.Pt(), leadTau.momentum.Eta());
        const double eff_data_subLeadTau = efficiency.eff2012IsoParkedTau19fb_Simone(subLeadTau.momentum.Pt(), subLeadTau.momentum.Eta());
        const double eff_mc_leadTau = efficiency.eff2012IsoParkedTau19fbMC_Simone(leadTau.momentum.Pt(), leadTau.momentum.Eta());
        const double eff_mc_subLeadTau = efficiency.eff2012IsoParkedTau19fbMC_Simone(subLeadTau.momentum.Pt(), subLeadTau.momentum.Eta());
        // first leadTau, second subLeadTau
        triggerWeights.push_back(eff_data_leadTau/eff_mc_leadTau);
        triggerWeights.push_back(eff_data_subLeadTau/eff_mc_subLeadTau);
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

    void FillSyncTree(const analysis::Candidate& higgs, const analysis::Candidate& higgs_sv,
                      const analysis::Candidate& higgs_sv_up, const analysis::Candidate& higgs_sv_down,
                      const analysis::CandidateVector& jets, const analysis::CandidateVector& jetsPt20,
                      const analysis::CandidateVector& bjets, const analysis::VertexVector& vertices)
    {
        const analysis::Candidate& leadTau = higgs.GetLeadingDaughter(analysis::Candidate::Tau);
        const analysis::Candidate& subLeadTau = higgs.GetSubleadingDaughter(analysis::Candidate::Tau);
        const ntuple::Tau& ntuple_tau_leg1 = correctedTaus.at(leadTau.index);
//        const ntuple::Tau& ntuple_tau_leg2 = correctedTaus.at(subLeadTau.index);

        H_BaseAnalyzer::FillSyncTree(higgs, higgs_sv, higgs_sv_up, higgs_sv_down, jets, jetsPt20, bjets, vertices, leadTau, subLeadTau);

        syncTree.iso_1() = ntuple_tau_leg1.byIsolationMVAraw;
        syncTree.mva_1() = 0;
        //syncTree.passid_1();
        //syncTree.passiso_1();
        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = ntuple_tau_leg1.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        syncTree.againstElectronMVA3raw_1() = ntuple_tau_leg1.againstElectronMVA3raw;
        syncTree.byIsolationMVA2raw_1() = ntuple_tau_leg1.byIsolationMVA2raw;
        syncTree.againstMuonLoose2_1() = ntuple_tau_leg1.againstMuonLoose2;
        syncTree.againstMuonMedium2_1() = ntuple_tau_leg1.againstMuonMedium2;
        syncTree.againstMuonTight2_1() = ntuple_tau_leg1.againstMuonTight2;
        //syncTree.mva_2();
        syncTree.Fill();
    }

private:
    BaselineAnalyzerData anaData;
};

#include "METPUSubtraction/interface/GBRProjectDict.cxx"
