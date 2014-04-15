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


    ENTRY_1D(float, Tau_Pt_MC)
    ENTRY_1D(float, E_Pt_MC)
//    ENTRY_1D(float, PV_sumPt)
    ENTRY_1D(float, ResonanceDeltaZ_PV)
    ENTRY_1D(float, ResonanceDeltaR_PV)

};

class HHbbETaujet_analyzer : public analysis::BaseAnalyzer {
public:
    HHbbETaujet_analyzer(const std::string& inputFileName, const std::string& outputFileName,
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
        finalState::bbETaujet eTauJet;
        if (useMCtruth && !FindAnalysisFinalState(eTauJet)) return;

        cuts::Cutter cut(anaData.EventSelection());

        cut(true, "total");

        cut(HaveTriggerFired(cuts::Htautau_Summer13::trigger::ETau::hltPaths), "trigger");

        const VertexVector vertices = CollectVertices();

        cut(vertices.size(), "vertex");
        primaryVertex = vertices.back();
//        anaData.PV_sumPt().Fill(primaryVertex.sumPtSquared);
        if(useMCtruth) {
            const double DeltaZ_PV = eTauJet.resonance->vertex.Z() - primaryVertex.position.Z();
            anaData.ResonanceDeltaZ_PV().Fill(DeltaZ_PV);
            const double DeltaR_PV = (eTauJet.resonance->vertex - primaryVertex.position).Perp();
            anaData.ResonanceDeltaR_PV().Fill(DeltaR_PV);
        }

        //SIGNAL OBJECTS
        const auto electrons = CollectElectrons(true);
        cut(electrons.size(), "electrons");

        const auto taus = CollectTaus(true);
        cut(taus.size(), "tau");

        const auto Higgses_e_tau = FindCompatibleObjects(electrons, taus,
                   cuts::Htautau_Summer13::DeltaR_betweenSignalObjects,analysis::Candidate::Higgs, "H_e_tau", 0);
        cut(Higgses_e_tau.size(), "H_e_tau");

        const auto b_jets = CollectBJets(cuts::Htautau_Summer13::btag::CSVL, "loose");
        cut(b_jets.size() >= 2, ">=2b_loose");

        const auto Higgses_bb =FindCompatibleObjects(b_jets, cuts::Htautau_Summer13::DeltaR_betweenSignalObjects,
                                      analysis::Candidate::Higgs, "H_bb", Candidate::UnknownCharge());
        cut(Higgses_bb.size(), "H_bb");

        const auto Resonances =
                FindCompatibleObjects(Higgses_e_tau, Higgses_bb, cuts::minDeltaR_betweenHiggses,
                                      analysis::Candidate::Resonance, "resonance", Candidate::UnknownCharge());
        cut(Resonances.size(), "resonance");

        //OBJECT VETO
        ApplyVetos(Resonances,cut);

    }

    virtual analysis::Candidate SelectElectron(size_t id, cuts::ObjectSelector& objectSelector,
                                               bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::electronID::ETau;
        const std::string selection_label = "electron";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Electron& object = event.electrons().at(id);

        cut(true, ">0 ele cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        const double eta = std::abs( X(eta, 120, -6.0, 6.0) );
        cut(eta < eta_high && (eta < cuts::Htautau_Summer13::electronID::eta_CrackVeto_low ||
                               eta > cuts::Htautau_Summer13::electronID::eta_CrackVeto_high), "eta");
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

        return analysis::Candidate(analysis::Candidate::Electron, id, object,object.charge);
    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector& objectSelector,
                                          bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::tauID::ETau;
        const std::string selection_label = "tau";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Tau& object = event.taus().at(id);

        cut(true, ">0 tau cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(decayModeFinding, 2, -0.5, 1.5) > decayModeFinding, "decay_mode");
        cut(X(byLooseCombinedIsolationDeltaBetaCorr3Hits, 2, -0.5, 1.5)
            > LooseCombinedIsolationDeltaBetaCorr3Hits, "looseIso3Hits");
        cut(X(againstMuonLoose, 2, -0.5, 1.5) > againstMuonLoose, "vs_mu_loose");
        cut(X(againstElectronMediumMVA5, 2, -0.5, 1.5) > againstElectronMediumMVA5, "vs_e_mediumMVA");

        return analysis::Candidate(analysis::Candidate::Tau, id, object,object.charge);
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
    ETauSignalAnalyzerData anaData;
};
