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

    virtual ~HmutauBaseline_sync() override
    {
        anaData.getOutputFile().cd();
        syncTree.Write();
    }


protected:
    virtual analysis::BaseAnalyzerData& GetAnaData() override { return anaData; }

    virtual void ProcessEvent() override
    {
        H_BaseAnalyzer::ProcessEvent();
        using namespace analysis;
        using namespace cuts::Htautau_Summer13;
        using namespace cuts::Htautau_Summer13::MuTau;
        finalState::MuTaujet muTau;
        if (useMCtruth && !FindAnalysisFinalState(muTau, false, true)) return;

        cuts::Cutter cut(&anaData.Selection("event"));
        cut(true, "total");

        cut(HaveTriggerFired(trigger::hltPaths), "trigger");

        const VertexVector vertices = CollectVertices();
        cut(vertices.size(), "vertex");
        primaryVertex = vertices.front();

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

        ApplyTauCorrections(muTau,false);
        const auto taus = CollectTaus();
        cut(taus.size(), "tau_cand");

        const auto higgses = FindCompatibleObjects(muons, taus, DeltaR_betweenSignalObjects,
                                                   Candidate::Higgs, "H_mu_tau", 0);
        cut(higgses.size(), "mu_tau");


        const auto higgsTriggered = ApplyTriggerMatch(higgses,trigger::hltPaths,false);
        cut(higgsTriggered.size(), "trigger obj match");

        const Candidate higgs = SelectSemiLeptonicHiggs(higgsTriggered);

        const ntuple::MET mvaMet = mvaMetProducer.ComputeMvaMet(higgs,event.pfCandidates(),
                                                                event.jets(),primaryVertex,
                                                                vertices,event.taus());

        ApplyTauCorrectionsToMVAMETandHiggs(mvaMet,higgs);

        const auto looseJets = CollectLooseJets();

        const auto filteredLooseJets = FilterCompatibleObjects(looseJets, higgs,
                                                               cuts::Htautau_Summer13::jetID::deltaR_signalObjects);


        const auto jets = CollectJets(filteredLooseJets);
        const auto bjets = CollectBJets(filteredLooseJets, false);
        const auto retagged_bjets = CollectBJets(filteredLooseJets, useMCtruth);

        ApplyPostRecoilCorrections(higgs, muTau.resonance, jets.size());
        const double m_sv = CorrectMassBySVfit(higgs, postRecoilMET,1);

        CalculateFullEventWeight(higgs);

        FillSyncTree(higgs, m_sv, jets, filteredLooseJets, bjets, retagged_bjets, vertices);
    }

    virtual analysis::Candidate SelectMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                           root_ext::AnalyzerData& _anaData,
                                           const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::muonID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Muon& object = event.muons().at(id);
        const analysis::Candidate muon(analysis::Candidate::Mu, id, object);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(isGlobalMuonPromptTight) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon) == isPFMuon, "PF");
        cut(X(nMatchedStations) > nMatched_Stations, "stations");
        cut(X(pixHits) > pixHits, "pix_hits");
        cut(X(trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double dB_PV = analysis::Calculate_dxy(mu_vertex,primaryVertex.position,muon.momentum);
        cut(std::abs( Y(dB_PV) ) < dB, "dB");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");

        return muon;
    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData,
                                          const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::tauID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Tau& object = correctedTaus.at(id);

        cut(true, ">0 tau cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(againstMuonTight) > againstMuonTight, "vs_mu_tight");
        cut(X(againstElectronLoose) > againstElectronLoose, "vs_e_loose");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits) < byCombinedIsolationDeltaBetaCorrRaw3Hits, "looseIso3Hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");


        const analysis::Candidate tau(analysis::Candidate::Tau, id, object);

        return tau;
    }

    analysis::CandidateVector CollectZmuons()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> analysis::Candidate
            { return SelectZmuon(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<analysis::Candidate>("z_muons", base_selector, event.muons().size());
    }

    virtual analysis::Candidate SelectZmuon(size_t id, cuts::ObjectSelector* objectSelector,
                                            root_ext::AnalyzerData& _anaData,
                                            const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::MuTau::ZmumuVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Muon& object = event.muons().at(id);
        const analysis::Candidate muon(analysis::Candidate::Mu, id, object);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(mu_vertex,primaryVertex.position,muon.momentum);
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(isTrackerMuon) == isTrackerMuon, "trackerMuon");
        cut(X(isPFMuon) == isPFMuon, "PFMuon");
        cut(X(pfRelIso) < pfRelIso, "pFRelIso");

        return muon;
    }


    bool FindAnalysisFinalState(analysis::finalState::MuTaujet& final_state, bool oneResonanceSelection, bool inclusive)
    {
        const bool base_result = H_BaseAnalyzer::FindAnalysisFinalState(final_state, oneResonanceSelection);
        if(inclusive || !base_result)
            return base_result;

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

    virtual void CalculateTriggerWeights(const analysis::Candidate& higgs) override
    {
        using namespace analysis::Htautau_Summer13::trigger::Run2012ABCD::MuTau;
        const analysis::Candidate& mu = higgs.GetDaughter(analysis::Candidate::Mu);
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        triggerWeights = CalculateWeights(mu.momentum, tau.momentum);
    }

    virtual void CalculateIsoWeights(const analysis::Candidate& higgs) override
    {
        using namespace cuts::Htautau_Summer13::MuTau::muonISOscaleFactor;
        IsoWeights.clear();
        const analysis::Candidate& mu = higgs.GetDaughter(analysis::Candidate::Mu);
        const double mu_pt = mu.momentum.Pt(), mu_eta = std::abs(mu.momentum.Eta());
        if(mu_pt < pt.at(0))
            throw std::runtime_error("No information about ISO. Muon pt is too small");
        const size_t pt_bin = mu_pt < pt.at(1) ? 0 : 1;
        if(mu_eta >= eta.at(2))
            throw std::runtime_error("No information about ISO. Muon eta is too big");
        const size_t eta_bin = mu_eta < eta.at(0) ? 0 : ( mu_eta < eta.at(1) ? 1 : 2 );
        const double scale = scaleFactors.at(pt_bin).at(eta_bin);
        // first mu, second tau
        IsoWeights.push_back(scale);
        IsoWeights.push_back(1);
    }

    virtual void CalculateIdWeights(const analysis::Candidate& higgs) override
    {
        using namespace cuts::Htautau_Summer13::MuTau::muonIDscaleFactor;
        IDweights.clear();
        const analysis::Candidate& mu = higgs.GetDaughter(analysis::Candidate::Mu);
        const double mu_pt = mu.momentum.Pt(), mu_eta = std::abs(mu.momentum.Eta());
        if(mu_pt < pt.at(0))
            throw std::runtime_error("No information about ID. Muon pt is too small");
        const size_t pt_bin = mu_pt < pt.at(1) ? 0 : 1;
        if(mu_eta >= eta.at(2))
            throw std::runtime_error("No information about ID. Muon eta is too big");
        const size_t eta_bin = mu_eta < eta.at(0) ? 0 : ( mu_eta < eta.at(1) ? 1 : 2 );
        const double scale = scaleFactors.at(pt_bin).at(eta_bin);
        // first mu, second tau
        IDweights.push_back(scale);
        IDweights.push_back(1);
    }

    virtual void CalculateFakeWeights(const analysis::Candidate& higgs) override
    {
        fakeWeights.clear();
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        const double tau_pt = tau.momentum.Pt();
        const double fake_weight =
                (1.15743)-(0.00736136*tau_pt)+(4.3699e-05*tau_pt*tau_pt)-(1.188e-07*tau_pt*tau_pt*tau_pt);
        // first mu, second tau
        fakeWeights.push_back(1);
        fakeWeights.push_back(fake_weight);
    }

    virtual void CalculateDMWeights(const analysis::Candidate& higgs) override
    {
        DMweights.clear();
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        const ntuple::Tau& tau_leg = correctedTaus.at(tau.index);
        const double weight = tau_leg.decayMode == ntuple::tau_id::kOneProng0PiZero ?
                                  cuts::Htautau_Summer13::tauCorrections::DecayModeWeight : 1;
        // first mu, second tau
        DMweights.push_back(1);
        DMweights.push_back(weight);
    }

    void FillSyncTree(const analysis::Candidate& higgs, double m_sv,
                      const analysis::CandidateVector& jets, const analysis::CandidateVector& jetsPt20,
                      const analysis::CandidateVector& bjets, const analysis::CandidateVector& retagged_bjets,
                      const analysis::VertexVector& vertices)
    {
        const analysis::Candidate& muon = higgs.GetDaughter(analysis::Candidate::Mu);
        const ntuple::Muon& ntuple_muon = event.muons().at(muon.index);
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);

        H_BaseAnalyzer::FillSyncTree(higgs, m_sv, jets, jetsPt20, bjets, retagged_bjets, vertices, muon, tau);

        syncTree.iso_1() = ntuple_muon.pfRelIso;
        syncTree.mva_1() = 0;
        //syncTree.passid_2();
        //syncTree.passiso_2();
        //syncTree.mva_2() = ntuple_tau.againstElectronLoose;
        syncTree.Fill();
    }

private:
    BaselineAnalyzerData anaData;
};

#include "METPUSubtraction/interface/GBRProjectDict.cxx"
