/*!
 * \file HHbbETaujet_analyzer.C
 * \brief X->HH->bbTauTau->b_jet,b_jet,e,tau_jet analysis.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-20 created
 */

#include "../include/BaseAnalyzer.h"

class ETauSignalAnalyzerData : public analysis::SignalAnalyzerData {

public:
    ETauSignalAnalyzerData(TFile& outputFile) : SignalAnalyzerData(outputFile) {}

    SELECTION_ENTRY(EventSelection, 15, 5)

    ENTRY_1D(double, Ele_tau_mass)
    ENTRY_1D(double, BB_mass)
    ENTRY_1D(double,Tau_Pt_MC)
    ENTRY_1D(double,E_Pt_MC)
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
        const CandidateVector electrons = CollectElectrons();
        cut(electrons.size(), "electron");
        const CandidateVector taus = CollectTaus();
        cut(taus.size(), "tau");
        const CandidateVector Higgses_ele_tau =
                FindCompatibleObjects(electrons, taus, cuts::Htautau_Summer13::DeltaR_signalLeptons,
                                      analysis::Candidate::Higgs, anaData.Ele_tau_mass());
        cut(Higgses_ele_tau.size(), "ele_tau");
        const CandidateVector muons = CollectMuons();
        cut(!muons.size(), "no_muon");

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

    virtual analysis::Candidate SelectMuon(Int_t id, bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::muonID::MuTau;
        cuts::Cutter cut(anaData.Counter(), anaData.MuonSelection(), enabled);

        cut(true, ">0 mu cand");
        cut(X(Muon_pt) > pt, "pt");
        cut(std::abs( X(Muon_eta) ) < eta, "eta");
        cut(!isTrackerMuon || X(Muon_isTrackerMuon), "tracker");
        cut(!isGlobalMuonPromptTight || X(Muon_isGlobalMuonPromptTight), "tight");
        cut(!isPFMuon || X(Muon_isPFMuon), "PF");
        cut(X(Muon_nChambers) > nChambers, "chamers");
        cut(X(Muon_nMatchedStations) > nMatched_Stations, "stations");
        cut(X(Muon_trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");
        cut(X(Muon_pixHits) > pixHits, "pix_hits");
        cut(X(Muon_globalChi2) < globalChiSquare, "chi2");
        cut(std::abs( X(Muon_dB) ) < dB, "dB");
        cut(std::abs( X(Muon_vtxDistZ) ) < dz, "dz");
        cut(X(Muon_pfRelIso) < pFRelIso, "pFRelIso");

        TLorentzVector momentum;
        momentum.SetPtEtaPhiE(event->Muon_pt[id], event->Muon_eta[id], event->Muon_phi[id],
                                             event->Muon_energy[id]);
        return analysis::Candidate(analysis::Candidate::Mu, id, momentum);
    }

    virtual analysis::Candidate SelectTau(Int_t id, bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::tauID::ETau;
        cuts::Cutter cut(anaData.Counter(), anaData.TauSelection(), enabled);

        cut(true, ">0 tau cand");
        cut(X(Tau_pt) > pt, "pt");
        cut(std::abs( X(Tau_eta) ) < eta, "eta");
        cut(X(Tau_decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits) > LooseCombinedIsolationDeltaBetaCorr3Hits, "looseIso3Hits");
        cut(X(Tau_againstMuonLoose) > againstMuonLoose, "vs_mu_loose");
        cut(X(Tau_againstElectronMediumMVA5) > againstElectronMediumMVA5, "vs_e_medium");

        TLorentzVector momentum;
        momentum.SetPtEtaPhiE(event->Tau_pt[id], event->Tau_eta[id], event->Tau_phi[id],
                              event->Tau_energy[id]);
        return analysis::Candidate(analysis::Candidate::Tau, id, momentum);
    }

    virtual analysis::Candidate SelectElectron(Int_t id, bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::electronID;
        cuts::Cutter cut(anaData.Counter(), anaData.ElectronSelection(), enabled);

        cut(true, ">0 ele cand");
        cut(X(Electron_pt) > pt, "pt");
        const double eta = std::abs( X(Electron_eta) );
        cut(eta < eta_high && (eta < eta_CrackVeto_low || eta > eta_CrackVeto_high), "eta");
        // cut dz_pv
        cut(X(Electron_missingHits) < missingHits, "mis_hits");
        cut(X(Electron_hasMatchedConv) > hasMatchedConv, "has_conv");
        cut(X(Electron_dB) < dB, "dB");
        cut(X(Electron_vtxDistZ) < dz_pv, "dZ");
        const size_t pt_index = event->Electron_pt[id] < ref_pt ? 0 : 1;
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(Electron_mvaPOGNonTrig) > MVApogNonTrig[pt_index][eta_index], "mva");

        TLorentzVector momentum;
        momentum.SetPtEtaPhiE(event->Electron_pt[id], event->Electron_eta[id], event->Electron_phi[id],
                              event->Electron_energy[id]);
        return analysis::Candidate(analysis::Candidate::Electron, id, momentum);
    }

    bool FindAnalysisFinalState(analysis::finalState::bbETaujet& finalState)
    {
        BaseAnalyzer::FindAnalysisFinalState(finalState);

        static const analysis::ParticleCodes TauMuonicDecay = { particles::mu, particles::nu_mu, particles::nu_tau };
        static const analysis::ParticleCodes TauElectronDecay = { particles::e, particles::nu_e, particles::nu_tau };

        finalState.electron = finalState.tau_jet = nullptr;
        for (const analysis::GenParticle* tau_MC : finalState.taus) {
            analysis::GenParticleVector TauProducts;
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
