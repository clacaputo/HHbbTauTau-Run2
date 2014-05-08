/*!
 * \file HtautauBaseline_sync.C
 * \brief Generate sync-tree for Htautau analysis using baseline selection.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-05 created
 */

#include "../include/H_BaseAnalyzer.h"

class BaselineAnalyzerData : public analysis::BaseAnalyzerData {
public:
    BaselineAnalyzerData(TFile& outputFile) : BaseAnalyzerData(outputFile) {}

};

class HmutauBaseline_sync : public analysis::H_BaseAnalyzer {
public:
    HmutauBaseline_sync(const std::string& inputFileName, const std::string& outputFileName,
                          const std::string& _prefix = "none", Long64_t _maxNumberOfEvents = 0,
                          bool _useMCtruth = false, const std::string& reweightFileName = "none")
        : H_BaseAnalyzer(inputFileName,outputFileName,_prefix,_maxNumberOfEvents,_useMCtruth, reweightFileName),
          anaData(*outputFile)
    {
        anaData.getOutputFile().cd();
    }

    virtual ~HmutauBaseline_sync()
    {
        anaData.getOutputFile().cd();
        syncTree.Write();
    }


protected:
    virtual analysis::BaseAnalyzerData& GetAnaData() { return anaData; }

    virtual void ProcessEvent()
    {
        H_BaseAnalyzer::ProcessEvent();
        using namespace analysis;
        using namespace cuts::Htautau_Summer13::MuTau;
        finalState::MuTaujet muTau;
        if (useMCtruth && !FindAnalysisFinalState(muTau)) return;

        ApplyTauCorrections(muTau);

        cuts::Cutter cut(&anaData.Selection("event"));
        cut(true, "total");

        cut(HaveTriggerFired(trigger::hltPaths), "trigger");

        const VertexVector vertices = CollectVertices();
        cut(vertices.size(), "vertex");
        primaryVertex = vertices.back();

        const auto z_muons = CollectZmuons();
        const auto z_muon_candidates = FindCompatibleObjects(z_muons, ZmumuVeto::deltaR, Candidate::Z, "Z_mu_mu", 0);
        cut(!z_muon_candidates.size(), "z_mumu_veto");

        const auto electrons_bkg = CollectBackgroundElectrons();
        cut(!electrons_bkg.size(), "no_electron");

        const auto muons = CollectMuons();
        cut(muons.size(), "muon_cand");
        cut(muons.size() == 1, "one_muon_cand");

        const auto muons_bkg = CollectBackgroundMuons();
        const bool have_bkg_muon = muons_bkg.size() > 1 ||
                ( muons_bkg.size() == 1 && muons_bkg.front() != muons.front() );
        cut(!have_bkg_muon, "no_bkg_muon");

        const auto taus = CollectTaus();
        cut(taus.size(), "tau_cand");

        const auto higgses = FindCompatibleObjects(muons, taus, DeltaR_betweenSignalObjects,
                                                   Candidate::Higgs, "H_mu_tau", 0);
        cut(higgses.size(), "mu_tau");

        const Candidate higgs = SelectSemiLeptonicHiggs(higgses);

        const auto jets = CollectJets(higgs);
        const auto bjets = CollectBJets(higgs);
        const Candidate higgs_corr = ApplyCorrections(higgs, muTau.resonance, jets.size());

        FillSyncTree(higgs_corr);
    }

    virtual analysis::Candidate SelectMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                           root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::muonID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Muon& object = event.muons().at(id);

        cut(true, ">0 mu cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(isGlobalMuonPromptTight, 2, -0.5, 1.5) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon, 2, -0.5, 1.5) == isPFMuon, "PF");
        cut(X(nMatchedStations, 10, 0.0, 10.0) > nMatched_Stations, "stations");
        cut(X(pixHits, 10, 0.0, 10.0) > pixHits, "pix_hits");
        cut(X(trackerLayersWithMeasurement, 20, 0.0, 20.0) > trackerLayersWithMeasurement, "layers");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double dB_PV = (mu_vertex - primaryVertex.position).Perp();
        cut(std::abs( Y(dB_PV, 50, 0.0, 0.5) ) < dB, "dB");
        cut(X(pfRelIso, 1000, 0.0, 100.0) < pFRelIso, "pFRelIso");

        const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
        cut(Y(haveTriggerMatch, 2, -0.5, 1.5), "triggerMatch");

        return analysis::Candidate(analysis::Candidate::Mu, id, object,object.charge);
    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::tauID;
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

    virtual analysis::Candidate SelectBackgroundElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                                         root_ext::AnalyzerData& _anaData,
                                                         const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::MuTau::electronVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event.electrons().at(id);

        cut(true, ">0 ele cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        const double eta = std::abs( X(eta, 120, -6.0, 6.0) );
        cut(eta < eta_high, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = (ele_vertex - primaryVertex.position).Perp(); // same as dB
        cut(std::abs( Y(d0_PV, 50, 0.0, 0.5) ) < d0, "d0");
        cut(X(pfRelIso, 1000, 0.0, 100.0) < pFRelIso, "pFRelIso");
        const size_t pt_index = object.pt < ref_pt ? 0 : 1;
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig, 300, -1.5, 1.5) > MVApogNonTrig[pt_index][eta_index], "mva");

        return analysis::Candidate(analysis::Candidate::Electron, id, object,object.charge);
    }

    virtual analysis::Candidate SelectBackgroundMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                                     root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::MuTau::muonVeto;
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
        cut(X(isGlobalMuonPromptTight, 2, -0.5, 1.5) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon, 2, -0.5, 1.5) == isPFMuon, "PF");
        cut(X(nMatchedStations, 10, 0.0, 10.0) > nMatched_Stations, "stations");
        cut(X(pixHits, 10, 0.0, 10.0) > pixHits, "pix_hits");
        cut(X(trackerLayersWithMeasurement, 20, 0.0, 20.0) > trackerLayersWithMeasurement, "layers");
        cut(X(pfRelIso, 1000, 0.0, 100.0) < pfRelIso, "pFRelIso");

        return analysis::Candidate(analysis::Candidate::Mu, id, object,object.charge);
    }

    analysis::CandidateVector CollectZmuons()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> analysis::Candidate
            { return SelectZmuon(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<analysis::Candidate>("z_muons", base_selector, event.muons().size());
    }

    virtual analysis::Candidate SelectZmuon(size_t id, cuts::ObjectSelector* objectSelector, root_ext::AnalyzerData& _anaData,
                                 const std::string& selection_label)
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


    bool FindAnalysisFinalState(analysis::finalState::MuTaujet& final_state)
    {
        if(!H_BaseAnalyzer::FindAnalysisFinalState(final_state))
            return false;

        final_state.muon = final_state.tau_jet = nullptr;
        for (const analysis::GenParticle* tau_MC : final_state.taus) {
            analysis::GenParticlePtrVector tauProducts;
            if (analysis::FindDecayProducts(*tau_MC, analysis::TauMuonicDecay, tauProducts))
                final_state.muon = tauProducts.at(0);
            else if (!analysis::FindDecayProducts(*tau_MC, analysis::TauElectronDecay, tauProducts))
                final_state.tau_jet = tau_MC;
        }

        if (!final_state.muon || !final_state.tau_jet) return false;
        return true;
    }

private:
    BaselineAnalyzerData anaData;
};

