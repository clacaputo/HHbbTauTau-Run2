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

    ENTRY_1D(float, Mass)
    ENTRY_1D(float, Tau_Pt_MC)
    ENTRY_1D(float, Mu_Pt_MC)
    ENTRY_1D(float, PV_sumPt)
    ENTRY_1D(float, ResonanceDeltaZ_PV)
    ENTRY_1D(float, ResonanceDeltaR_PV)
    ENTRY_1D(unsigned, N_objects)
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

        const VertexVector vertices = CollectVertices();
        anaData.N_objects("vertices").Fill(vertices.size());
        cut(vertices.size(), "vertex");
        primaryVertex = vertices.back();
        anaData.PV_sumPt().Fill(primaryVertex.sumPt);
        if(useMCtruth) {
            const double DeltaZ_PV = muTauJet.resonance->vertex.Z() - primaryVertex.position.Z();
            anaData.ResonanceDeltaZ_PV().Fill(DeltaZ_PV);
            const double DeltaR_PV = (muTauJet.resonance->vertex - primaryVertex.position).Perp();
            anaData.ResonanceDeltaR_PV().Fill(DeltaR_PV);
        }
        const CandidateVector muons = CollectMuons(true);
        anaData.N_objects("muons").Fill(muons.size());
        cut(muons.size(), "muon");
        const CandidateVector taus = CollectTaus(true);
        anaData.N_objects("taus").Fill(taus.size());
        cut(taus.size(), "tau");
        const CandidateVector Higgses_mu_tau =
                FindCompatibleObjects(muons, taus, cuts::Htautau_Summer13::DeltaR_betweenSignalObjects,
                                      analysis::Candidate::Higgs, anaData.Mass("mu_tau"));
        anaData.N_objects("Hmutau").Fill(Higgses_mu_tau.size());
        cut(Higgses_mu_tau.size(), "mu_tau");


        const CandidateVector b_jets = CollectBJets(cuts::Htautau_Summer13::btag::CSVL, "loose");
//        const CandidateVector b_jets = CollectBJets(cuts::Htautau_Summer13::btag::CSVM, "medium");
        anaData.N_objects("bjets").Fill(b_jets.size());
        cut(b_jets.size() >= 2, ">=2b_loose");

        const CandidateVector Higgses_bb =
                FindCompatibleObjects(b_jets, cuts::Htautau_Summer13::DeltaR_betweenSignalObjects,
                                      analysis::Candidate::Higgs, anaData.Mass("bb"));
        anaData.N_objects("Hbb").Fill(Higgses_bb.size());
        cut(Higgses_bb.size(), "bb");

        const CandidateVector Resonances =
                FindCompatibleObjects(Higgses_mu_tau, Higgses_bb, 1,
                                      analysis::Candidate::Resonance, anaData.Mass("resonance"));
        anaData.N_objects("resonance").Fill(Resonances.size());
        cut(Resonances.size(), "resonance");

        const CandidateVector electrons_bkg = CollectElectrons(false);
        anaData.N_objects("electrons_bkg").Fill(electrons_bkg.size());
        const CandidateVector resonances_noEle = FilterBackground(Resonances,electrons_bkg,
                                                        cuts::Htautau_Summer13::electronID::veto::deltaR_signalObjects);
        anaData.N_objects("resonances_noEle").Fill(resonances_noEle.size());
        cut(resonances_noEle.size(), "no_electrons");

        const CandidateVector muons_bkg = CollectMuons(false);
        anaData.N_objects("muons_bkg").Fill(muons_bkg.size());
        const CandidateVector resonances_noMu = FilterBackground(resonances_noEle,muons_bkg,
                                                        cuts::Htautau_Summer13::muonID::veto::deltaR_signalObjects);
        anaData.N_objects("resonances_noMu").Fill(resonances_noMu.size());
        cut(resonances_noMu.size(), "no_muons");

        const CandidateVector bjets_bkg = CollectBJets(cuts::Htautau_Summer13::btag::CSVL, "loose",false);
        anaData.N_objects("bjets_bkg").Fill(bjets_bkg.size());
        const CandidateVector resonances_noBjets = FilterBackground(resonances_noMu,bjets_bkg,
                                                        cuts::Htautau_Summer13::btag::veto::deltaR_signalObjects);
        anaData.N_objects("resonances_noBjets").Fill(resonances_noBjets.size());
        cut(resonances_noBjets.size(), "no_bjets");

        const CandidateVector taus_bkg = CollectTaus(false);
        anaData.N_objects("taus_bkg").Fill(taus_bkg.size());
        const CandidateVector resonances_noTau = FilterBackground(resonances_noBjets,taus_bkg,
                                                        cuts::Htautau_Summer13::tauID::veto::deltaR_signalObjects);
        anaData.N_objects("resonances_noTau").Fill(resonances_noTau.size());
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

    virtual analysis::Candidate SelectBackgroundElectron(size_t id, bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::electronID::veto;
        const std::string selection_label = "bkg";
        cuts::Cutter cut(anaData.Counter(), anaData.ElectronSelectionBkg(), enabled);
        const ntuple::Electron& object = event.electrons().at(id);

        cut(true, ">0 ele cand");
        cut(X(pt) > pt, "pt");
        const double eta = std::abs( X(eta) );
        cut(eta < eta_high && (eta < cuts::Htautau_Summer13::electronID::eta_CrackVeto_low ||
                               eta > cuts::Htautau_Summer13::electronID::eta_CrackVeto_high), "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");
        const size_t pt_index = object.pt < ref_pt ? 0 : 1;
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[pt_index][eta_index], "mva");

        return analysis::Candidate(analysis::Candidate::Electron, id, object);
    }

    virtual analysis::Candidate SelectBackgroundMuon(size_t id, bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::muonID::veto;
        const std::string selection_label = "bkg";
        cuts::Cutter cut(anaData.Counter(), anaData.MuonSelectionBkg(), enabled);
        const ntuple::Muon& object = event.muons().at(id);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(isTightMuon) == isTightMuon, "tight");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");

        return analysis::Candidate(analysis::Candidate::Mu, id, object);
    }

    virtual analysis::Candidate SelectBackgroundTau(size_t id, bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::tauID::veto;
        const std::string selection_label = "bkg";
        cuts::Cutter cut(anaData.Counter(), anaData.TauSelectionBkg(), enabled);
        const ntuple::Tau& object = event.taus().at(id);

        cut(true, ">0 tau cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(byLooseCombinedIsolationDeltaBetaCorr3Hits) > LooseCombinedIsolationDeltaBetaCorr3Hits, "looseIso3Hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");

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
