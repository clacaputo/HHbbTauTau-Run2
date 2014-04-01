/*!
 * \file HHbbETaujet_analyzer.C
 * \brief X->HH->bbTauTau->b_jet,b_jet,e,tau_jet analysis.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-04-01 created
 */

#include "../include/BaseAnalyzer.h"

class ETauSignalAnalyzerData : public analysis::SignalAnalyzerData {
public:
    ETauSignalAnalyzerData(TFile& outputFile) : SignalAnalyzerData(outputFile) {}

    SELECTION_ENTRY(EventSelection, 15, 5)


    ENTRY_1D(float, Tau_Pt_MC)
    ENTRY_1D(float, E_Pt_MC)
    ENTRY_1D(float, PV_sumPt)
    ENTRY_1D(float, ResonanceDeltaZ_PV)
    ENTRY_1D(float, ResonanceDeltaR_PV)

};

class HHbbETaujet_analyzer : public analysis::BaseAnalyzer {
public:
    HHbbETaujet_analyzer(const std::string& inputFileName, const std::string& outputFileName,
                        Long64_t _maxNumberOfEvents = 0, bool _useMCtruth = false)
        : BaseAnalyzer(inputFileName,outputFileName,_maxNumberOfEvents,_useMCtruth), anaData(*outputFile)
    {
        anaData.getOutputFile().cd();
    }

protected:
    virtual analysis::SignalAnalyzerData& GetAnaData() { return anaData; }

    virtual void ProcessEvent()
    {
        using namespace analysis;
        finalState::bbETaujet eTauJet;
        if (useMCtruth && !FindAnalysisFinalState(eTauJet)) return;

        cuts::Cutter cut(anaData.EventSelection(), anaData.EventSelection());

        cut(true, "total");

        const VertexVector vertices = CollectVertices();

        cut(vertices.size(), "vertex");
        primaryVertex = vertices.back();
        anaData.PV_sumPt().Fill(primaryVertex.sumPtSquared);
        if(useMCtruth) {
            const double DeltaZ_PV = eTauJet.resonance->vertex.Z() - primaryVertex.position.Z();
            anaData.ResonanceDeltaZ_PV().Fill(DeltaZ_PV);
            const double DeltaR_PV = (eTauJet.resonance->vertex - primaryVertex.position).Perp();
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

        const auto b_jets = CollectBJets(cuts::Htautau_Summer13::btag::CSVL, "loose");
        cut(b_jets.size() >= 2, ">=2b_loose");

        const auto Higgses_bb =FindCompatibleObjects(b_jets, cuts::Htautau_Summer13::DeltaR_betweenSignalObjects,
                                      analysis::Candidate::Higgs, "H_bb");
        cut(Higgses_bb.size(), "H_bb");

        const auto Resonances =
                FindCompatibleObjects(Higgses_mu_tau, Higgses_bb, 1,analysis::Candidate::Resonance, "resonance");
        cut(Resonances.size(), "resonance");

        //OBJECT VETO
        const auto electrons_bkg = CollectElectrons(false);
        const auto resonances_noEle = FilterBackground(Resonances,electrons_bkg,
                              cuts::Htautau_Summer13::electronID::veto::deltaR_signalObjects,"resonances_noEle");
        cut(resonances_noEle.size(), "no_electrons");

        const auto muons_bkg = CollectMuons(false);
        const auto resonances_noMu = FilterBackground(resonances_noEle,muons_bkg,
                                 cuts::Htautau_Summer13::muonID::veto::deltaR_signalObjects,"resonances_noMu");
        cut(resonances_noMu.size(), "no_muons");

        const auto bjets_bkg = CollectBJets(cuts::Htautau_Summer13::btag::CSVL, "loose",false);
        const auto resonances_noBjets = FilterBackground(resonances_noMu,bjets_bkg,
                              cuts::Htautau_Summer13::btag::veto::deltaR_signalObjects,"resonances_noBjets");
        cut(resonances_noBjets.size(), "no_bjets");

        const auto taus_bkg = CollectTaus(false);
        const auto resonances_noTau = FilterBackground(resonances_noBjets,taus_bkg,
                                 cuts::Htautau_Summer13::tauID::veto::deltaR_signalObjects,"resonances_noTau");
        cut(resonances_noTau.size(), "no_taus");

    }

    virtual analysis::Candidate SelectMuon(size_t id, bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::muonID::MuTau;
        const std::string selection_label = "";
        cuts::Cutter cut(anaData.Counter(), anaData.MuonSelection(), enabled);
        const ntuple::Muon& object = event.muons().at(id);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(isGlobalMuonPromptTight) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon) == isPFMuon, "PF");
        cut(X(nMatchedStations) > nMatched_Stations, "stations");
        cut(X(trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");
        cut(X(pixHits) > pixHits, "pix_hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double dB_PV = (mu_vertex - primaryVertex.position).Perp();
        cut(std::abs( Y(dB_PV) ) < dB, "dB");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");

        return analysis::Candidate(analysis::Candidate::Mu, id, object);
    }

    virtual analysis::Candidate SelectTau(size_t id, bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::tauID::MuTau;
        const std::string selection_label = "";
        cuts::Cutter cut(anaData.Counter(), anaData.TauSelection(), enabled);
        const ntuple::Tau& object = event.taus().at(id);

        cut(true, ">0 tau cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(byLooseCombinedIsolationDeltaBetaCorr3Hits) > LooseCombinedIsolationDeltaBetaCorr3Hits, "looseIso3Hits");
        cut(X(againstMuonTight) > againstMuonTight, "vs_mu_tight");
        cut(X(againstElectronLoose) > againstElectronLoose, "vs_e_loose");

        return analysis::Candidate(analysis::Candidate::Tau, id, object);
    }


    bool FindAnalysisFinalState(analysis::finalState::bbETaujet& finalState)
    {
        BaseAnalyzer::FindAnalysisFinalState(finalState);

        static const analysis::ParticleCodes TauMuonicDecay = { particles::mu, particles::nu_mu, particles::nu_tau };
        static const analysis::ParticleCodes TauElectronDecay = { particles::e, particles::nu_e, particles::nu_tau };

        finalState.electron = finalState.tau_jet = nullptr;
        for (const analysis::GenParticle* tau_MC : finalState.taus) {
            analysis::GenParticlePtrVector TauProducts;
            if (analysis::FindDecayProducts(*tau_MC,TauElectronDecay,TauProducts))
                finalState.electron = TauProducts.at(0);
            else if (!analysis::FindDecayProducts(*tau_MC,TauMuonicDecay,TauProducts))
                finalState.tau_jet = tau_MC;
        }

        if (!finalState.electron || !finalState.tau_jet) return false;

        anaData.Tau_Pt_MC().Fill(finalState.tau_jet->momentum.Pt());
        anaData.E_Pt_MC().Fill(finalState.electron->momentum.Pt());
        return true;
    }

private:
    MuTauSignalAnalyzerData anaData;
};
