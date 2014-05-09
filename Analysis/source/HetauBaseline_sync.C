/*!
 * \file HtautauBaseline_sync.C
 * \brief Generate sync-tree for H->tautau->e_taujet analysis using baseline selection.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-05 created
 */

#include "../include/H_BaseAnalyzer.h"

class BaselineAnalyzerData : public analysis::BaseAnalyzerData {
public:
    BaselineAnalyzerData(TFile& outputFile) : BaseAnalyzerData(outputFile) {}

};

class HetauBaseline_sync : public analysis::H_BaseAnalyzer {
public:
    HetauBaseline_sync(const std::string& inputFileName, const std::string& outputFileName,
                       const std::string& _prefix = "none", Long64_t _maxNumberOfEvents = 0,
                       bool _useMCtruth = false, const std::string& reweightFileName = "none")
        : H_BaseAnalyzer(inputFileName,outputFileName,_prefix,_maxNumberOfEvents,_useMCtruth, reweightFileName),
          anaData(*outputFile)
    {
        anaData.getOutputFile().cd();
    }

    ~HetauBaseline_sync()
    {
        anaData.getOutputFile().cd();
        syncTree.Write();
    }

protected:
    virtual analysis::SignalAnalyzerData& GetAnaData() { return anaData; }

    virtual void ProcessEvent()
    {
        H_BaseAnalyzer::ProcessEvent();
        using namespace analysis;
        using namespace cuts::Htautau_Summer13::ETau;
        finalState::ETaujet eTau;
        if (useMCtruth && !FindAnalysisFinalState(eTau)) return;

        ApplyTauCorrections(eTau);

        cuts::Cutter cut(&anaData.Selection("event"));
        cut(true, "total");

        cut(HaveTriggerFired(trigger::hltPaths), "trigger");

        const VertexVector vertices = CollectVertices();
        cut(vertices.size(), "vertex");
        primaryVertex = vertices.back();

        const auto z_electrons = CollectZelectrons();
        const auto z_electron_candidates = FindCompatibleObjects(z_electrons, ZeeVeto::deltaR,
                                                                 Candidate::Z, "Z_e_e", 0);
        cut(!z_electron_candidates.size(), "z_ee_veto");

        const auto muons_bkg = CollectBackgroundMuons();
        cut(!muons_bkg.size(), "no_muon");

        const auto electrons = CollectElectrons();
        cut(electrons.size(), "electron_cand");
        cut(electrons.size() == 1, "one_electron_cand");

        const auto electrons_bkg = CollectBackgroundElectrons();
        const bool have_bkg_electron = electrons_bkg.size() > 1 ||
                ( electrons_bkg.size() == 1 && electrons_bkg.front() != electrons.front() );
        cut(!have_bkg_electron, "no_bkg_electron");

        const auto taus = CollectTaus();
        cut(taus.size(), "tau_cand");

        const auto higgses = FindCompatibleObjects(electrons, taus, DeltaR_betweenSignalObjects,
                                                   Candidate::Higgs, "H_e_tau", 0);
        cut(higgses.size(), "e_tau");

        const Candidate higgs = SelectSemiLeptonicHiggs(higgses);

        const auto jets = CollectJets(higgs);
        const auto bjets = CollectBJets(higgs);
        const Candidate higgs_corr = ApplyCorrections(higgs, eTau.resonance, jets.size());

        FillSyncTree(higgs, higgs_corr, jets, bjets, vertices);
    }

    virtual analysis::Candidate SelectElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                               root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::electronID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event.electrons().at(id);

        cut(true, ">0 ele cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        const double eta = std::abs( X(eta, 120, -6.0, 6.0) );
        cut(eta < eta_high, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");
        cut(X(missingHits, 20, 0.0, 20.0) < missingHits, "missingHits");
        cut(X(hasMatchedConversion, 2, -0.5, 1.5) == hasMatchedConversion, "conversion");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        const double dB_PV = (ele_vertex - primaryVertex.position).Perp();
        cut(std::abs( Y(dB_PV, 50, 0.0, 0.5) ) < dB, "dB");
        cut(X(pfRelIso, 1000, 0.0, 100.0) < pFRelIso, "pFRelIso");
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig, 300, -1.5, 1.5) > MVApogNonTrig[eta_index], "mva");

        const bool haveTriggerMatch = HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
        cut(Y(haveTriggerMatch, 2, -0.5, 1.5), "triggerMatch");
        return analysis::Candidate(analysis::Candidate::Electron, id, object,object.charge);
    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::tauID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Tau& object = correctedTaus.at(id);

        cut(true, ">0 tau cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(decayModeFinding, 2, -0.5, 1.5) > decayModeFinding, "decay_mode");
        cut(X(againstMuonTight, 2, -0.5, 1.5) > againstMuonTight, "vs_mu_tight");
        cut(X(againstElectronLoose, 2, -0.5, 1.5) > againstElectronLoose, "vs_e_loose");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits, 100, 0, 10)
            < byCombinedIsolationDeltaBetaCorrRaw3Hits, "looseIso3Hits");

        const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
        cut(Y(haveTriggerMatch, 2, -0.5, 1.5), "triggerMatch");

        return analysis::Candidate(analysis::Candidate::Tau, id, object,object.charge);
    }

    analysis::CandidateVector CollectZelectrons()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> analysis::Candidate
            { return SelectZelectron(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<analysis::Candidate>("z_electrons", base_selector, event.electrons().size());
    }

    virtual analysis::Candidate SelectZelectron(size_t id, cuts::ObjectSelector* objectSelector,
                                                root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::MuTau::ZmumuVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Muon& object = event.muons().at(id);

        cut(true, ">0 mu cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = (mu_vertex - primaryVertex.position).Perp();
        cut(std::abs( Y(d0_PV, 50, 0.0, 0.5) ) < d0, "d0");
        cut(X(isTrackerMuon, 2, -0.5, 1.5) == isTrackerMuon, "trackerMuon");
        cut(X(isPFMuon, 2, -0.5, 1.5) == isPFMuon, "PFMuon");
        cut(X(pfRelIso, 1000, 0.0, 100.0) < pfRelIso, "pFRelIso");

        return analysis::Candidate(analysis::Candidate::Mu, id, object,object.charge);
    }

    bool FindAnalysisFinalState(analysis::finalState::ETaujet& final_state)
    {
        if(!H_BaseAnalyzer::FindAnalysisFinalState(final_state))
            return false;

        final_state.electron = final_state.tau_jet = nullptr;
        for (const analysis::GenParticle* tau_MC : final_state.taus) {
            analysis::GenParticlePtrVector tauProducts;
            if (analysis::FindDecayProducts(*tau_MC, analysis::TauElectronDecay, tauProducts))
                final_state.electron = tauProducts.at(0);
            else if (!analysis::FindDecayProducts(*tau_MC, analysis::TauMuonicDecay, tauProducts))
                final_state.tau_jet = tau_MC;
        }

        if (!final_state.electron || !final_state.tau_jet) return false;
        return true;
    }

    void FillSyncTree(ntuple::SyncTree& syncTree, const analysis::Candidate& higgs)
    {
        syncTree.run() = event.eventInfo().run;
        syncTree.lumi() = event.eventInfo().lumis;
        syncTree.evt() = event.eventInfo().EventId;
        for (unsigned n = 0; n < event.eventInfo().bunchCrossing.size(); ++n){
            if (event.eventInfo().bunchCrossing.at(n) == 0){
                syncTree.npu() = event.eventInfo().nPU.at(n); //only in-time PU
            }
        }
        syncTree.puweight() = weight;
        Double_t DMweight = 1;
        analysis::Candidate cand1 = higgs.finalStateDaughters.at(0);
        analysis::Candidate cand2 = higgs.finalStateDaughters.at(1);
        if (cand1.momentum.Pt() < cand2.momentum.Pt()){
            cand1 = higgs.finalStateDaughters.at(1);
            cand2 = higgs.finalStateDaughters.at(0);
        }
        const ntuple::Tau& leg1 = correctedTaus.at(cand1.index);
        const ntuple::Tau& leg2 = correctedTaus.at(cand2.index);
        if (leg1.decayMode == 0)
            DMweight *= 0.88;
        if (leg2.decayMode == 0)
            DMweight *= 0.88;
        syncTree.decaymodeweight() = DMweight;
        syncTree.mvis() = higgs.momentum.M();

        syncTree.pt_1() = cand1.momentum.Pt();
        syncTree.eta_1() = cand1.momentum.Eta();
        syncTree.phi_1() = cand1.momentum.Phi();
        syncTree.m_1() = cand1.momentum.M();
        syncTree.q_1() = leg1.charge;
        syncTree.iso_1() = leg1.byIsolationMVAraw;
        syncTree.passid_1() = 1;
        syncTree.passiso_1() = 1;
        syncTree.mt_1() = cand1.momentum.Mt();
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

        syncTree.pt_2() = cand2.momentum.Pt();
        syncTree.eta_2() = cand2.momentum.Eta();
        syncTree.phi_2() = cand2.momentum.Phi();
        syncTree.m_2() = cand2.momentum.M();
        syncTree.q_2() = leg2.charge;
        syncTree.iso_2() = leg2.byIsolationMVAraw;
        syncTree.passid_2() = 1;
        syncTree.passiso_2() = 1;
        syncTree.mt_2() = cand2.momentum.Mt();
        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = leg2.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        syncTree.againstElectronMVA3raw_2() = leg2.againstElectronMVA3raw;
        syncTree.againstElectronMVA3category_2() = leg2.againstElectronMVA3category;
        syncTree.byIsolationMVA2raw_2() = leg2.byIsolationMVA2raw;
        syncTree.againstMuonLoose_2() = leg2.againstMuonLoose;
        syncTree.againstMuonLoose2_2() = leg2.againstMuonLoose2;
        syncTree.againstMuonMedium2_2() = leg2.againstMuonMedium2;
        syncTree.againstMuonTight2_2() = leg2.againstMuonTight2;
        syncTree.againstElectronLooseMVA3_2() = leg2.againstElectronLooseMVA3;
        syncTree.againstElectronLoose_2() = leg2.againstElectronLoose;

        syncTree.met() = event.metPF().pt_uncorrected; //raw
        syncTree.metphi() = event.metPF().phi_uncorrected; //raw
        syncTree.mvamet() = event.metMVA().pt;
        syncTree.mvametphi() = event.metMVA().phi;
        syncTree.metcov00() = event.metPF().significanceMatrix.at(0);
        syncTree.metcov01() = event.metPF().significanceMatrix.at(1);
        syncTree.metcov10() = event.metPF().significanceMatrix.at(2);
        syncTree.metcov11() = event.metPF().significanceMatrix.at(3);
        syncTree.mvacov00() = event.metMVA().significanceMatrix.at(0);
        syncTree.mvacov01() = event.metMVA().significanceMatrix.at(1);
        syncTree.mvacov10() = event.metMVA().significanceMatrix.at(2);
        syncTree.mvacov11() = event.metMVA().significanceMatrix.at(3);

        syncTree.Fill();
    }

private:
    BaselineAnalyzerData anaData;
};

