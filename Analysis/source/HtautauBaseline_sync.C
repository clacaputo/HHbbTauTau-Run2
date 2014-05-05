/*!
 * \file HtautauBaseline_sync.C
 * \brief Generate sync-tree for Htautau analysis using baseline selection.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-05 created
 */

#include "../include/BaseAnalyzer.h"
#include "../include/SyncTree.h"

class BaselineAnalyzerData : public analysis::SignalAnalyzerData {
public:
    BaselineAnalyzerData(TFile& outputFile) : SignalAnalyzerData(outputFile) {}

};

class HtautauBaseline_sync : public analysis::BaseAnalyzer {
public:
    HtautauBaseline_sync(const std::string& inputFileName, const std::string& outputFileName,
                          const std::string& _prefix = "none", Long64_t _maxNumberOfEvents = 0,
                          bool _useMCtruth = false, const std::string& reweightFileName = "none")
        : BaseAnalyzer(inputFileName,outputFileName,_prefix,_maxNumberOfEvents,_useMCtruth, reweightFileName),
          anaData(*outputFile), etauTree("etauTree"), mutauTree("mutauTree"), tautauTree("tautauTree")
    {
        anaData.getOutputFile().cd();
    }

protected:
    virtual analysis::SignalAnalyzerData& GetAnaData() { return anaData; }

    virtual void ProcessEvent()
    {
        BaseAnalyzer::ProcessEvent();
        using namespace analysis;
        finalState::TauTau tauTau;
        if (useMCtruth && !FindAnalysisFinalState(tauTau)) return;

        cuts::Cutter cut(anaData.EventSelection());

        cut(true, "total");

        const VertexVector vertices = CollectVertices();

        cut(vertices.size(), "vertex");
        primaryVertex = vertices.back();

//        cut(HaveTriggerFired(cuts::Htautau_Summer13::trigger::MuTau::hltPaths), "trigger");


        //SIGNAL OBJECTS
        const auto electrons = CollectElectrons(true);
        const auto muons = CollectMuons(true);
        const auto taus_eTau = CollectTaus_eTau(true);
        const auto taus_muTau = CollectTaus_muTau(true);
        const auto taus_TauTau = CollectTaus_TauTau(true);

        const auto Higgses_e_tau_signal = FindCompatibleObjects(electrons, taus_eTau,
                   cuts::Htautau_Summer13::DeltaR_betweenSignalObjects,analysis::Candidate::Higgs, "H_e_tau", 0);
        const auto Higgses_mu_tau_signal = FindCompatibleObjects(muons, taus_muTau,
                   cuts::Htautau_Summer13::DeltaR_betweenSignalObjects, analysis::Candidate::Higgs, "H_mu_tau", 0);
        const auto Higgses_tautau_loose = FindCompatibleObjects(taus_TauTau,
                   cuts::Htautau_Summer13::DeltaR_betweenSignalObjects, analysis::Candidate::Higgs, "H_2tau", 0);
        const auto Higgses_tautau_signal = ApplyTauFullSelection(Higgses_tautau_loose);

        const auto Higgses_e_tau = ApplyVetos(Higgses_e_tau_signal);
        const auto Higgses_mu_tau = ApplyVetos(Higgses_mu_tau_signal);
        const auto Higgses_tautau = ApplyVetos(Higgses_tautau_signal);

        finalState::ETaujet eTaujet(tauTau);
        finalState::MuTaujet muTaujet(tauTau);
        finalState::TaujetTaujet taujetTaujet(tauTau);
        if(HaveTriggerFired(cuts::Htautau_Summer13::trigger::ETau::hltPaths) &&
                ( !useMCtruth || FindAnalysisFinalState(eTaujet) ) )
            FillSyncTree(etauTree, Higgses_e_tau);
        else if(HaveTriggerFired(cuts::Htautau_Summer13::trigger::MuTau::hltPaths) &&
                ( !useMCtruth || FindAnalysisFinalState(muTaujet) ) )
            FillSyncTree(mutauTree, Higgses_mu_tau);
        else if(HaveTriggerFired(cuts::Htautau_Summer13::trigger::TauTau::hltPaths) &&
                ( !useMCtruth || FindAnalysisFinalState(taujetTaujet) ) )
            FillSyncTree(tautauTree, Higgses_tautau);
    }

    virtual analysis::Candidate SelectMuon(size_t id, cuts::ObjectSelector& objectSelector,
                                           bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::muonID::MuTau;
        const std::string selection_label = "muon";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Muon& object = event.muons().at(id);

        cut(true, ">0 mu cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(isGlobalMuonPromptTight, 2, -0.5, 1.5) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon, 2, -0.5, 1.5) == isPFMuon, "PF");
        cut(X(nMatchedStations, 10, 0.0, 10.0) > nMatched_Stations, "stations");
        cut(X(trackerLayersWithMeasurement, 20, 0.0, 20.0) > trackerLayersWithMeasurement, "layers");
        cut(X(pixHits, 10, 0.0, 10.0) > pixHits, "pix_hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double dB_PV = (mu_vertex - primaryVertex.position).Perp();
        cut(std::abs( Y(dB_PV, 50, 0.0, 0.5) ) < dB, "dB");
        cut(X(pfRelIso, 1000, 0.0, 100.0) < pFRelIso, "pFRelIso");

        return analysis::Candidate(analysis::Candidate::Mu, id, object,object.charge);
    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector& objectSelector,
                                          bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::tauID::MuTau;
        const std::string selection_label = "tau";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Tau& object = event.taus().at(id);

        cut(true, ">0 tau cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(decayModeFinding, 2, -0.5, 1.5) > decayModeFinding, "decay_mode");
        cut(X(byLooseCombinedIsolationDeltaBetaCorr3Hits, 2, -0.5, 1.5)
            > LooseCombinedIsolationDeltaBetaCorr3Hits, "looseIso3Hits");
        cut(X(againstMuonTight, 2, -0.5, 1.5) > againstMuonTight, "vs_mu_tight");
        cut(X(againstElectronLoose, 2, -0.5, 1.5) > againstElectronLoose, "vs_e_loose");

        return analysis::Candidate(analysis::Candidate::Tau, id, object,object.charge);
    }


    bool FindAnalysisFinalState(analysis::finalState::TauTau& finalState)
    {
        return true;
    }

    void FillSyncTree(ntuple::SyncTree& syncTree, const analysis::Candidate& higgs){
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
        analysis::Candidate* cand1 = higgs.finalStateDaughters.at(0);
        analysis::Candidate* cand2 = higgs.finalStateDaughters.at(1);
        if (cand1->momentum.Pt() < cand2->momentum.Pt()){
            cand1 = higgs.finalStateDaughters.at(1);
            cand2 = higgs.finalStateDaughters.at(0);
        }
        const ntuple::Tau& leg1 = event.taus().at(cand1->index);
        const ntuple::Tau& leg2 = event.taus().at(cand2->index);
        if (leg1.decayMode == 0)
            DMweight *= 0.88;
        if (leg2.decayMode == 0)
            DMweight *= 0.88;
        syncTree.decaymodeweight() = DMweight;
        syncTree.mvis() = higgs.momentum.M();

        syncTree.pt_1() = cand1->momentum.Pt();
        syncTree.eta_1() = cand1->momentum.Eta();
        syncTree.phi_1() = cand1->momentum.Phi();
        syncTree.m_1() = cand1->momentum.M();
        syncTree.q_1() = leg1.charge;
        syncTree.iso_1() = leg1.byIsolationMVAraw;
        syncTree.passid_1() = 1;
        syncTree.passiso_1() = 1;
        syncTree.mt_1() = cand1->momentum.Mt();
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

        syncTree.pt_2() = cand2->momentum.Pt();
        syncTree.eta_2() = cand2->momentum.Eta();
        syncTree.phi_2() = cand2->momentum.Phi();
        syncTree.m_2() = cand2->momentum.M();
        syncTree.q_2() = leg2.charge;
        syncTree.iso_2() = leg2.byIsolationMVAraw;
        syncTree.passid_2() = 1;
        syncTree.passiso_2() = 1;
        syncTree.mt_2() = cand2->momentum.Mt();
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
    }

private:
    BaselineAnalyzerData anaData;
    ntuple::SyncTree etauTree, mutauTree, tautauTree;
};

