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
        if (useMCtruth && !FindAnalysisFinalState(muTau,false)) return;

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

        //not more needed
        //ApplyTauCorrections(muTau,event.metMVAmuTau(), false);
        //ApplyTauCorrections(muTau,event.metPF(), false);

        ApplyTauCorrections(muTau,false);
        //correctedTaus = event.taus(); //no tau corrections

        const auto taus = CollectTaus();
        cut(taus.size(), "tau_cand");


        const auto higgses = FindCompatibleObjects(muons, taus, DeltaR_betweenSignalObjects,
                                                   Candidate::Higgs, "H_mu_tau", 0);
        cut(higgses.size(), "mu_tau");


        const auto higgsTriggered = ApplyTriggerMatch(higgses,trigger::hltPaths,false);
        cut(higgsTriggered.size(), "trigger obj match");

//        const Candidate higgs_withoutTauCorrections = SelectSemiLeptonicHiggs(higgsTriggered);

//        const ntuple::MET mvaMet = mvaMetProducer.ComputeMvaMet(higgs_withoutTauCorrections,event.pfCandidates(),
//                                                                event.jets(),primaryVertex,
//                                                                vertices,event.taus());

//        ApplyTauCorrections(muTau,false);
//        const Candidate higgs = ApplyTauCorrectionsToMVAMETandHiggs(mvaMet,higgs_withoutTauCorrections);

        const Candidate higgs = SelectSemiLeptonicHiggs(higgsTriggered);

        const ntuple::MET mvaMet = mvaMetProducer.ComputeMvaMet(higgs,event.pfCandidates(),
                                                                event.jets(),primaryVertex,
                                                                vertices,event.taus());


        ApplyTauCorrectionsToMVAMETandHiggs(mvaMet,higgs);

        const auto jetsPt20 = CollectJetsPt20();
        const auto jets = CollectJets();

        const auto filteredJets = FilterCompatibleObjects(jets,higgs,cuts::Htautau_Summer13::jetID::deltaR_signalObjects);


        const auto bjets = CollectBJets(higgs);

        //Apply Post Recoil Corrections
        ApplyPostRecoilCorrections(higgs, muTau.resonance, filteredJets.size());
        //postRecoilMET = correctedMET; //tau corrections

//        const Candidate higgs_sv = CorrectMassBySVfit(higgs, postRecoilMET,1);
//        const Candidate higgs_sv_up = CorrectMassBySVfit(higgs, postRecoilMET,1.03);
//        const Candidate higgs_sv_down = CorrectMassBySVfit(higgs, postRecoilMET,0.97);

//        CalculateFullEventWeight(higgs_sv);

//        FillSyncTree(higgs, higgs_sv, higgs_sv_up, higgs_sv_down, filteredJets, jetsPt20, bjets, vertices);

        //postRecoilMET = mvaMet;
//        postRecoilMET = correctedMET; //tau corrections

        CalculateFullEventWeight(higgs);
        FillSyncTree(higgs, higgs, higgs, higgs, filteredJets, jetsPt20, bjets, vertices);
    }

    virtual analysis::Candidate SelectMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                           root_ext::AnalyzerData& _anaData,
                                           const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::muonID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Muon& object = event.muons().at(id);

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
        const double dB_PV = (mu_vertex - primaryVertex.position).Perp();
        cut(std::abs( Y(dB_PV) ) < dB, "dB");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");

        const analysis::Candidate muon(analysis::Candidate::Mu, id, object);
//        const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
//        const bool haveTriggerMatch = analysis::HaveTriggerMatched(event.triggerObjects(), trigger::hltPaths,muon);
//        cut(Y(haveTriggerMatch), "triggerMatch");

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

        const analysis::Candidate tau(analysis::Candidate::Tau, id, object);
//        const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
//        const bool haveTriggerMatch = analysis::HaveTriggerMatched(event.triggerObjects(), trigger::hltPaths,tau);
//        cut(Y(haveTriggerMatch, 2, -0.5, 1.5), "triggerMatch");

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

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = (mu_vertex - primaryVertex.position).Perp();
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(isTrackerMuon) == isTrackerMuon, "trackerMuon");
        cut(X(isPFMuon) == isPFMuon, "PFMuon");
        cut(X(pfRelIso) < pfRelIso, "pFRelIso");

        return analysis::Candidate(analysis::Candidate::Mu, id, object);
    }


//    bool FindAnalysisFinalState(analysis::finalState::MuTaujet& final_state)
//    {
//        if(!H_BaseAnalyzer::FindAnalysisFinalState(final_state))
//            return false;

//        final_state.muon = final_state.tau_jet = nullptr;
//        for (const analysis::GenParticle* tau_MC : final_state.taus) {
//            analysis::GenParticlePtrVector tauProducts;
//            if (analysis::FindDecayProducts(*tau_MC, analysis::TauMuonicDecay, tauProducts))
//                final_state.muon = tauProducts.at(0);
//            else if (!analysis::FindDecayProducts(*tau_MC, analysis::TauElectronDecay, tauProducts))
//                final_state.tau_jet = tau_MC;
//        }

//        if (!final_state.muon || !final_state.tau_jet) return false;
//        return true;
//    }

    virtual void CalculateTriggerWeights(const analysis::Candidate& higgs) override
    {
        triggerWeights.clear();
        const analysis::Candidate& mu = higgs.GetDaughter(analysis::Candidate::Mu);
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        analysis::Htautau_Summer13::TriggerEfficiency efficiency;
        const double eff_data_Mu = efficiency.effMu_muTau_Data_2012ABCD(mu.momentum.Pt(), mu.momentum.Eta());
        const double eff_data_Tau = efficiency.effTau_muTau_Data_2012ABCD(tau.momentum.Pt(), tau.momentum.Eta());
        const double eff_mc_Mu = efficiency.effMu_muTau_MC_2012ABCD(mu.momentum.Pt(), mu.momentum.Eta());
        const double eff_mc_tau = efficiency.effTau_muTau_MC_2012ABCD(tau.momentum.Pt(), tau.momentum.Eta());
        // first mu, second tau
        triggerWeights.push_back(eff_data_Mu/eff_mc_Mu);
        triggerWeights.push_back(eff_data_Tau/eff_mc_tau);
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

    void FillSyncTree(const analysis::Candidate& higgs, const analysis::Candidate& higgs_sv,
                      const analysis::Candidate& higgs_sv_up, const analysis::Candidate& higgs_sv_down,
                      const analysis::CandidateVector& jets, const analysis::CandidateVector& jetsPt20,
                      const analysis::CandidateVector& bjets, const analysis::VertexVector& vertices)
    {
        const analysis::Candidate& muon = higgs.GetDaughter(analysis::Candidate::Mu);
        const ntuple::Muon& ntuple_muon = event.muons().at(muon.index);
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
//        const ntuple::Tau& ntuple_tau = event.taus().at(tau.index);

        H_BaseAnalyzer::FillSyncTree(higgs, higgs_sv, higgs_sv_up, higgs_sv_down, jets, jetsPt20, bjets, vertices, muon, tau);

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
