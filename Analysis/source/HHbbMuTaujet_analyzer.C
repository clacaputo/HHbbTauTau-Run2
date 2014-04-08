/*!
 * \file HHbbMuTaujet_analyzer.C
 * \brief X->HH->bbTauTau->b_jet,b_jet,mu,tau_jet analysis.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-12 created
 */

#include "../include/BaseAnalyzer.h"

class MuTauSignalAnalyzerData : public analysis::SignalAnalyzerData {
public:
    MuTauSignalAnalyzerData(TFile& outputFile) : SignalAnalyzerData(outputFile) {}

    SELECTION_ENTRY(EventSelection)


    ENTRY_1D(float, Tau_Pt_MC)
    ENTRY_1D(float, Mu_Pt_MC)
    ENTRY_1D(float, PV_sumPt)
    ENTRY_1D(float, ResonanceDeltaZ_PV)
    ENTRY_1D(float, ResonanceDeltaR_PV)

};

class HHbbMuTaujet_analyzer : public analysis::BaseAnalyzer {
public:
    HHbbMuTaujet_analyzer(const std::string& inputFileName, const std::string& outputFileName,
                          const std::string& _prefix = "none", Long64_t _maxNumberOfEvents = 0,
                          bool _useMCtruth = false)
        : BaseAnalyzer(inputFileName,outputFileName,_prefix,_maxNumberOfEvents,_useMCtruth), anaData(*outputFile)
    {
        anaData.getOutputFile().cd();
    }

protected:
    virtual analysis::SignalAnalyzerData& GetAnaData() { return anaData; }

    virtual void ProcessEvent()
    {
        using namespace analysis;
        finalState::bbMuTaujet muTauJet;
        if (useMCtruth && !FindAnalysisFinalState(muTauJet)) return;

        cuts::Cutter cut(anaData.EventSelection());

        cut(true, "total");

        const VertexVector vertices = CollectVertices();

        cut(vertices.size(), "vertex");
        primaryVertex = vertices.back();
        anaData.PV_sumPt().Fill(primaryVertex.sumPtSquared);
        if(useMCtruth) {
            const double DeltaZ_PV = muTauJet.resonance->vertex.Z() - primaryVertex.position.Z();
            anaData.ResonanceDeltaZ_PV().Fill(DeltaZ_PV);
            const double DeltaR_PV = (muTauJet.resonance->vertex - primaryVertex.position).Perp();
            anaData.ResonanceDeltaR_PV().Fill(DeltaR_PV);
        }

        //SIGNAL OBJECTS
        const auto muons = CollectMuons(true);
        cut(muons.size(), "muon");

        const auto taus = CollectTaus(true);
        cut(taus.size(), "tau");

        const auto Higgses_mu_tau = FindCompatibleObjects(muons, taus,
                   cuts::Htautau_Summer13::DeltaR_betweenSignalObjects,analysis::Candidate::Higgs, "H_mu_tau");
        cut(Higgses_mu_tau.size(), "H_mu_tau");

        const auto b_jets = CollectBJets(cuts::Htautau_Summer13::btag::CSVL, "bjets_loose");
        cut(b_jets.size() >= 2, ">=2b_loose");

        const auto Higgses_bb =FindCompatibleObjects(b_jets, cuts::Htautau_Summer13::DeltaR_betweenSignalObjects,
                                      analysis::Candidate::Higgs, "H_bb");
        cut(Higgses_bb.size(), "H_bb");

        const auto Resonances =
                FindCompatibleObjects(Higgses_mu_tau, Higgses_bb, cuts::minDeltaR_betweenHiggses,
                                      analysis::Candidate::Resonance, "resonance");
        cut(Resonances.size(), "resonance");

        //OBJECT VETO
        ApplyVetos(Resonances, cut);
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

        return analysis::Candidate(analysis::Candidate::Mu, id, object);
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

        return analysis::Candidate(analysis::Candidate::Tau, id, object);
    }


    bool FindAnalysisFinalState(analysis::finalState::bbMuTaujet& finalState)
    {
        BaseAnalyzer::FindAnalysisFinalState(finalState);

        static const analysis::ParticleCodes TauMuonicDecay = { particles::mu, particles::nu_mu, particles::nu_tau };
        static const analysis::ParticleCodes TauElectronDecay = { particles::e, particles::nu_e, particles::nu_tau };

        finalState.muon = finalState.tau_jet = nullptr;
        for (const analysis::GenParticle* tau_MC : finalState.taus) {
            analysis::GenParticlePtrVector TauProducts;
            if (analysis::FindDecayProducts(*tau_MC,TauMuonicDecay,TauProducts))
                finalState.muon = TauProducts.at(0);
            else if (!analysis::FindDecayProducts(*tau_MC,TauElectronDecay,TauProducts))
                finalState.tau_jet = tau_MC;
        }

        if (!finalState.muon || !finalState.tau_jet) return false;

        anaData.Tau_Pt_MC().Fill(finalState.tau_jet->momentum.Pt());
        anaData.Mu_Pt_MC().Fill(finalState.muon->momentum.Pt());
        return true;
    }

private:
    MuTauSignalAnalyzerData anaData;
};
