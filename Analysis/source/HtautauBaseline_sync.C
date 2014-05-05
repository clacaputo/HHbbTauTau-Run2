/*!
 * \file HtautauBaseline_sync.C
 * \brief Generate sync-tree for Htautau analysis using baseline selection.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-05 created
 */

#include "../include/BaseAnalyzer.h"

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
          anaData(*outputFile)
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
            FillSyncTree(muTauTree, Higgses_mu_tau);
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

private:
    BaselineAnalyzerData anaData;
};

