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

    SELECTION_ENTRY(EventSelection, 15, 5)

    ENTRY_1D(float, Mu_tau_mass)
    ENTRY_1D(float, BB_mass)
    ENTRY_1D(float, Tau_Pt_MC)
    ENTRY_1D(float, Mu_Pt_MC)
};

class HHbbMuTaujet_analyzer : public analysis::BaseAnalyzer {
public:
    HHbbMuTaujet_analyzer(const std::string& inputFileName, const std::string& outputFileName,
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
        finalState::bbMuTaujet muTauJet;
        if (useMCtruth && !FindAnalysisFinalState(muTauJet)) return;

        cuts::Cutter cut(anaData.EventSelection(), anaData.EventSelection());

        cut(true, "total");
        const CandidateVector muons = CollectMuons(true);
        cut(muons.size(), "muon");
        const CandidateVector taus = CollectTaus(true);
        cut(taus.size(), "tau");
        const CandidateVector Higgses_mu_tau =
                FindCompatibleObjects(muons, taus, cuts::Htautau_Summer13::DeltaR_betweenSignalLeptons,
                                      analysis::Candidate::Higgs, anaData.Mu_tau_mass());
        cut(Higgses_mu_tau.size(), "mu_tau");
        const CandidateVector electrons = CollectElectrons(false);
        cut(!electrons.size(), "no_electron");

        const CandidateVector b_jets_loose = CollectBJets(cuts::Htautau_Summer13::btag::CSVL, "loose");
        const CandidateVector b_jets_medium = CollectBJets(cuts::Htautau_Summer13::btag::CSVM, "medium");
        cut.test(b_jets_loose.size() == 2, "2b_loose");
        cut.test(b_jets_loose.size() == 2 && b_jets_medium.size() >= 1, "1b_loose+1b_medium");
        if (cut.test(b_jets_medium.size() == 2, "2b_medium")){
            const Candidate Higgs_bb(Candidate::Higgs, b_jets_medium.at(0), b_jets_medium.at(1));
            anaData.BB_mass().Fill(Higgs_bb.momentum.M());
        }
        cut.test(b_jets_loose.size() >= 2, ">=2b_loose");
        cut.test(b_jets_loose.size() >= 2 && b_jets_medium.size() >= 1, ">=1b_loose+>=1b_medium");
        cut.test(b_jets_medium.size() >= 2, ">=2b_medium");
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
        cut(std::abs( X(dB) ) < dB, "dB");
        cut(std::abs( X(vtxDistZ) ) < dz, "dz");
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

    virtual analysis::Candidate SelectBackgroundElectron(size_t id, bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::electronID::veto;
        const std::string selection_label = "bkg";
        cuts::Cutter cut(anaData.Counter(), anaData.ElectronSelection(), enabled);
        const ntuple::Electron& object = event.electrons().at(id);

        cut(true, ">0 ele cand");
        cut(X(pt) > pt, "pt");
        const double eta = std::abs( X(eta) );
        cut(eta < eta_high && (eta < cuts::Htautau_Summer13::electronID::eta_CrackVeto_low ||
                               eta > cuts::Htautau_Summer13::electronID::eta_CrackVeto_high), "eta");
        cut(X(vtxDistZ) < dz, "dZ");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");
        const size_t pt_index = object.pt < ref_pt ? 0 : 1;
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[pt_index][eta_index], "mva");

        return analysis::Candidate(analysis::Candidate::Electron, id, object);
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
