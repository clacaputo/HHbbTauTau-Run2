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

    virtual ~HtautauBaseline_sync()
    {
        anaData.getOutputFile().cd();
        syncTree.Write();
    }

protected:
    virtual analysis::BaseAnalyzerData& GetAnaData() { return anaData; }

    virtual void ProcessEvent()
    {
        H_BaseAnalyzer::ProcessEvent();
        using namespace analysis;
        using namespace cuts::Htautau_Summer13;
        using namespace cuts::Htautau_Summer13::TauTau;
        finalState::TaujetTaujet tauTau;
        if (useMCtruth && !FindAnalysisFinalState(tauTau)) return;



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

        ApplyTauCorrections(tauTau);
        const auto taus = CollectTaus();
        cut(taus.size(), "tau_cand");
        cut(taus.size() >= 2, "at least 2 taus");



        const auto higgses = FindCompatibleObjects(taus,
                   cuts::Htautau_Summer13::DeltaR_betweenSignalObjects, analysis::Candidate::Higgs, "H_2tau");

        cut(higgses.size(), "DeltaR higgses");

        const auto higgses_tautau = ApplyTauFullSelection(higgses);

        cut(higgses_tautau.size(), "tau_tau selection");

        const Candidate higgs = SelectFullyHadronicHiggs(higgses_tautau);


        const auto jets = CollectJets(higgs);
        const auto bjets = CollectBJets(higgs);

        //const Candidate higgs_corr = ApplyCorrections(higgs, tauTau.resonance, jets.size());
        //FillSyncTree(higgs, higgs_corr, jets, bjets, vertices);
        postRecoilMET = correctedMET;
        FillSyncTree(higgs, higgs, jets, bjets, vertices);
        //std::cout << "I filled the tree" << std::endl;

    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::TauTau;
        using namespace cuts::Htautau_Summer13::TauTau::tauID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Tau& object = correctedTaus.at(id);

        cut(true, ">0 tau cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(decayModeFinding, 2, -0.5, 1.5) > decayModeFinding, "decay_mode");
        cut(X(againstMuonLoose, 2, -0.5, 1.5) > againstMuonLoose, "vs_mu_loose");
        cut(X(againstElectronLoose, 2, -0.5, 1.5) > againstElectronLoose, "vs_e_loose");
        cut(X(byMediumCombinedIsolationDeltaBetaCorr3Hits, 2, -0.5, 1.5)
            > byMediumCombinedIsolationDeltaBetaCorr3Hits, "mediumIso3Hits");

        //const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
        const analysis::Candidate tau(analysis::Candidate::Tau, id, object,object.charge);
        const bool haveTriggerMatch = analysis::HaveTriggerMatched(event.triggerObjects(), trigger::hltPaths,tau);
        cut(Y(haveTriggerMatch, 2, -0.5, 1.5), "triggerMatch");

        return tau;
    }

    analysis::CandidateVector ApplyTauFullSelection(const analysis::CandidateVector& higgses)
    {
        using namespace analysis;
        using namespace cuts::Htautau_Summer13::TauTau::tauID;
        CandidateVector result;
        for(const Candidate& higgs : higgses) {
            //const Candidate& leading_tau = higgs.GetLeadingDaughter(Candidate::Tau);
            const Candidate& subleading_tau = higgs.GetSubleadingDaughter(Candidate::Tau);
            const ntuple::Tau& ntuple_subleadingTau = correctedTaus.at(subleading_tau.index);
            if(ntuple_subleadingTau.againstElectronLooseMVA3 > againstElectronLooseMVA3)
                result.push_back(higgs);
        }
        return result;
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

    void FillSyncTree(const analysis::Candidate& higgs, const analysis::Candidate& higgs_corr,
                      const analysis::CandidateVector& jets, const analysis::CandidateVector& bjets,
                      const analysis::VertexVector& vertices)
    {
        Double_t DMweight = 1;
        analysis::Candidate leadTau, subLeadTau;
        if (higgs.finalStateDaughters.at(0).momentum.Pt() < higgs.finalStateDaughters.at(1).momentum.Pt()){
            leadTau = higgs.finalStateDaughters.at(1);
            subLeadTau = higgs.finalStateDaughters.at(0);
        }
        else {
            leadTau = higgs.finalStateDaughters.at(0);
            subLeadTau = higgs.finalStateDaughters.at(1);
        }
        const ntuple::Tau& leg1 = correctedTaus.at(leadTau.index);
        const ntuple::Tau& leg2 = correctedTaus.at(subLeadTau.index);
//        const ntuple::Tau& leg1 = event.taus().at(leadTau.index);
//        const ntuple::Tau& leg2 = event.taus().at(subLeadTau.index);
        if (leg1.decayMode == 0)
            DMweight *= 0.88;
        if (leg2.decayMode == 0)
            DMweight *= 0.88;
        syncTree.decaymodeweight() = DMweight;

        H_BaseAnalyzer::FillSyncTree(higgs, higgs_corr, jets, bjets, vertices, subLeadTau);

        //leadTau - not corrected!
        syncTree.pt_1() = leg1.pt;
        syncTree.eta_1() = leg1.eta;
        syncTree.phi_1() = leg1.phi;
        TLorentzVector leg1_momentum;
        leg1_momentum.SetPtEtaPhiE(leg1.pt,leg1.eta,leg1.phi,leg1.energy);
        syncTree.m_1() = leg1_momentum.M();
        syncTree.q_1() = leg1.charge;
        syncTree.iso_1() = leg1.byIsolationMVAraw;
        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = leg1.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        syncTree.againstElectronMVA3raw_1() = leg1.againstElectronMVA3raw;
        syncTree.againstElectronMVA3category_1() = leg1.againstElectronMVA3category;
        syncTree.byIsolationMVA2raw_1() = leg1.byIsolationMVA2raw;
        syncTree.againstMuonLoose_1() = leg1.againstMuonLoose;
        syncTree.againstMuonLoose2_1() = leg1.againstMuonLoose2;
        syncTree.againstMuonMedium2_1() = leg1.againstMuonMedium2;
        syncTree.againstMuonTight2_1() = leg1.againstMuonTight2;
        syncTree.againstElectronLooseMVA3_1() = leg1.againstElectronLooseMVA3;
        syncTree.againstElectronLoose_1() = leg1.againstElectronLoose;

        syncTree.Fill();
    }

private:
    BaselineAnalyzerData anaData;
};

