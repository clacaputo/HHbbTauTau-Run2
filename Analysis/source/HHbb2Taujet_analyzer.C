/*!
 * \file HHbb2Taujet_analyzer.C
 * \brief X->HH->bbTauTau->b_jet,b_jet,tau_jet,tau_jet analysis.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-04-01 created
 */

#include "../include/BaseAnalyzer.h"

class MuTauSignalAnalyzerData : public analysis::SignalAnalyzerData {
public:
    MuTauSignalAnalyzerData(TFile& outputFile) : SignalAnalyzerData(outputFile) {}

    ENTRY_1D(float, LeadTau_Pt_MC)
    ENTRY_1D(float, SubleadingTau_Pt_MC)
    ENTRY_1D(float, DR_tauJets_MC)
//    ENTRY_1D(float, PV_sumPtSquared)
    ENTRY_1D(float, ResonanceDeltaZ_PV)
    ENTRY_1D(float, ResonanceDeltaR_PV)

};

class HHbb2Taujet_analyzer : public analysis::BaseAnalyzer {
public:
    HHbb2Taujet_analyzer(const std::string& inputFileName, const std::string& outputFileName,
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
        finalState::bbTaujetTaujet mc_truth;
        if (useMCtruth && !FindAnalysisFinalState(mc_truth)) return;

        cuts::Cutter cut(anaData.EventSelection());

        cut(true, "total");

        const VertexVector vertices = CollectVertices();
        cut(vertices.size(), "vertex");
        primaryVertex = vertices.back();
//        anaData.PV_sumPtSquared().Fill(primaryVertex.sumPtSquared);
        if(useMCtruth) {
            const double DeltaZ_PV = mc_truth.resonance->vertex.Z() - primaryVertex.position.Z();
            anaData.ResonanceDeltaZ_PV().Fill(DeltaZ_PV);
            const double DeltaR_PV = (mc_truth.resonance->vertex - primaryVertex.position).Perp();
            anaData.ResonanceDeltaR_PV().Fill(DeltaR_PV);
        }

        //SIGNAL OBJECTS
        const auto taus = CollectTaus(true);
        cut(taus.size() >= 2, ">=2_tau");

        const auto Higgses_tautau = FindCompatibleObjects(taus, cuts::Htautau_Summer13::DeltaR_betweenSignalObjects,
                                                           analysis::Candidate::Higgs, "H_2tau");
        cut(Higgses_tautau.size(), "H_2tau");

        const auto Higgses_tautau_full = ApplyTauFullSelection(Higgses_tautau);
        cut(Higgses_tautau_full.size(), "H_2tau_full_selection");

        const auto b_jets = CollectBJets(cuts::Htautau_Summer13::btag::CSVL, "loose");
        cut(b_jets.size() >= 2, ">=2b_loose");

        const auto Higgses_bb =FindCompatibleObjects(b_jets, cuts::Htautau_Summer13::DeltaR_betweenSignalObjects,
                                                     analysis::Candidate::Higgs, "H_bb");
        cut(Higgses_bb.size(), "H_bb");

        const auto Resonances = FindCompatibleObjects(Higgses_tautau_full, Higgses_bb, cuts::minDeltaR_betweenHiggses,
                                                      analysis::Candidate::Resonance, "resonance");
        cut(Resonances.size(), "resonance");

        //OBJECT VETO
        ApplyVetos(Resonances, cut);
    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                          root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::tauID::TauTau;
        const std::string selection_label = "tau";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Tau& object = event.taus().at(id);

        cut(true, ">0 tau cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt_sublead, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(decayModeFinding, 2, -0.5, 1.5) > decayModeFinding, "decay_mode");
        cut(X(againstMuonLoose, 2, -0.5, 1.5) > againstMuonLoose, "vs_mu_tight");
        cut(X(againstElectronLoose, 2, -0.5, 1.5) > againstElectronLoose, "vs_e_loose");
        cut(X(byMediumCombinedIsolationDeltaBetaCorr3Hits, 2, -0.5, 1.5)
            > MediumCombinedIsolationDeltaBetaCorr3Hits, "looseIso3Hits");

        return analysis::Candidate(analysis::Candidate::Tau, id, object);
    }

    analysis::CandidateVector ApplyTauFullSelection(const analysis::CandidateVector& higgses)
    {
        using namespace analysis;
        using namespace cuts::Htautau_Summer13::tauID::TauTau;
        CandidateVector result;
        for(const Candidate& higgs : higgses) {
            const Candidate *leading_tau = nullptr, *subleading_tau = nullptr;
            if(higgs.daughters.size() != 2 || higgs.daughters.at(0)->type != Candidate::Tau
                    || higgs.daughters.at(1)->type != Candidate::Tau)
                throw std::runtime_error("bad higgs to tautau");
            if(higgs.daughters.at(0)->momentum.Pt() > higgs.daughters.at(1)->momentum.Pt()) {
                leading_tau = higgs.daughters.at(0);
                subleading_tau = higgs.daughters.at(1);
            } else {
                leading_tau = higgs.daughters.at(1);
                subleading_tau = higgs.daughters.at(0);
            }
            const ntuple::Tau& ntuple_subleadingTau = event.taus().at(subleading_tau->index);
            if(ntuple_subleadingTau.againstElectronLooseMVA5 > againstElectronLooseMVA5
                    && leading_tau->momentum.Pt() > pt_lead)
                result.push_back(higgs);
        }
        return result;
    }

    bool FindAnalysisFinalState(analysis::finalState::bbTaujetTaujet& finalState)
    {
        BaseAnalyzer::FindAnalysisFinalState(finalState);

        static const analysis::ParticleCodes TauMuonicDecay = { particles::mu, particles::nu_mu, particles::nu_tau };
        static const analysis::ParticleCodes TauElectronDecay = { particles::e, particles::nu_e, particles::nu_tau };

        unsigned n_hadronic_taus = 0;

        if(finalState.taus.size() != 2)
            throw std::runtime_error("bad bbtautau MC final state");

        for (const analysis::GenParticle* tau_MC : finalState.taus) {
            analysis::GenParticlePtrVector TauProducts;
            if (!analysis::FindDecayProducts(*tau_MC,TauMuonicDecay,TauProducts)
                    && !analysis::FindDecayProducts(*tau_MC,TauElectronDecay,TauProducts))
                ++n_hadronic_taus;
        }

        if (n_hadronic_taus != 2) return false;
        if(finalState.taus.at(0)->momentum.Pt() > finalState.taus.at(1)->momentum.Pt()) {
            finalState.leading_tau_jet = finalState.taus.at(0);
            finalState.subleading_tau_jet = finalState.taus.at(1);
        } else {
            finalState.leading_tau_jet = finalState.taus.at(1);
            finalState.subleading_tau_jet = finalState.taus.at(0);
        }

        anaData.LeadTau_Pt_MC().Fill(finalState.leading_tau_jet->momentum.Pt());
        anaData.SubleadingTau_Pt_MC().Fill(finalState.subleading_tau_jet->momentum.Pt());
        const double deltaR = finalState.leading_tau_jet->momentum.DeltaR(finalState.subleading_tau_jet->momentum);
        anaData.DR_tauJets_MC().Fill(deltaR);
        return true;
    }

private:
    MuTauSignalAnalyzerData anaData;
};
