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

        FillSyncTree(higgs, higgs_corr, jets, bjets, vertices);
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

    virtual void FillSyncTree(const analysis::Candidate& higgs, const analysis::Candidate& higgs_corr,
                              const analysis::CandidateVector& jets, const analysis::CandidateVector& bjets,
                              const analysis::VertexVector& vertices)
    {
        syncTree.run() = event.eventInfo().run;
        syncTree.lumi() = event.eventInfo().lumis;
        syncTree.evt() = event.eventInfo().EventId;
        syncTree.npv() = vertices.size();
        for (unsigned n = 0; n < event.eventInfo().bunchCrossing.size(); ++n){
            if (event.eventInfo().bunchCrossing.at(n) == 0){
                syncTree.npu() = event.eventInfo().nPU.at(n); //only in-time PU
            }
        }
        syncTree.puweight() = weight;



        syncTree.mvis() = higgs.momentum.M();
        syncTree.m_sv() = higgs_corr.momentum.M();
        syncTree.pt_sv() = higgs_corr.momentum.Pt();
        syncTree.eta_sv() = higgs_corr.momentum.Eta();
        syncTree.phi_sv() = higgs_corr.momentum.Phi();

        const analysis::Candidate& muon = higgs.GetDaughter(analysis::Candidate::Mu);
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        const ntuple::Muon& ntuple_muon = event.muons().at(muon.index);
        const ntuple::Tau& ntuple_tau = correctedTaus.at(tau.index);

        //muon
        syncTree.pt_1() = muon.momentum.Pt();
        syncTree.eta_1() = muon.momentum.Eta();
        syncTree.phi_1() = muon.momentum.Phi();
        syncTree.m_1() = muon.momentum.M();
        syncTree.q_1() = muon.charge;
        syncTree.iso_1() = ntuple_muon.pfRelIso;

        const TVector3 mu_vertex(ntuple_muon.vx, ntuple_muon.vy, ntuple_muon.vz);
        syncTree.d0_1() = (mu_vertex - primaryVertex.position).Perp();
        syncTree.dZ_1() = std::abs(ntuple_muon.vz - primaryVertex.position.Z());


        //tau
        syncTree.pt_2() = tau.momentum.Pt();
        syncTree.eta_2() = tau.momentum.Eta();
        syncTree.phi_2() = tau.momentum.Phi();
        syncTree.m_2() = tau.momentum.M();
        syncTree.q_2() = tau.charge;
        syncTree.iso_2() = ntuple_tau.byIsolationMVAraw;
        const TVector3 tau_vertex(ntuple_tau.vx, ntuple_tau.vy, ntuple_tau.vz);
        syncTree.d0_2() = (tau_vertex - primaryVertex.position).Perp();
        syncTree.dZ_2() = std::abs(ntuple_tau.vz - primaryVertex.position.Z());

        Double_t DMweight = 1;
        if (ntuple_tau.decayMode == 0)
            DMweight *= 0.88;
        syncTree.decaymodeweight() = DMweight;
        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = ntuple_tau.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        syncTree.againstElectronMVA3raw_2() = ntuple_tau.againstElectronMVA3raw;
        syncTree.againstElectronMVA3category_2() = ntuple_tau.againstElectronMVA3category;
        syncTree.byIsolationMVA2raw_2() = ntuple_tau.byIsolationMVA2raw;
        syncTree.againstMuonLoose_2() = ntuple_tau.againstMuonLoose;
        syncTree.againstMuonLoose2_2() = ntuple_tau.againstMuonLoose2;
        syncTree.againstMuonMedium2_2() = ntuple_tau.againstMuonMedium2;
        syncTree.againstMuonTight2_2() = ntuple_tau.againstMuonTight2;
        syncTree.againstElectronLooseMVA3_2() = ntuple_tau.againstElectronLooseMVA3;
        syncTree.againstElectronLoose_2() = ntuple_tau.againstElectronLoose;

        syncTree.met() = event.metPF().pt; //raw
        syncTree.metphi() = event.metPF().phi; //raw
        syncTree.mvamet() = correctedMET.pt;
        syncTree.mvametphi() = correctedMET.phi;
        syncTree.metcov00() = event.metPF().significanceMatrix.at(0);
        syncTree.metcov01() = event.metPF().significanceMatrix.at(1);
        syncTree.metcov10() = event.metPF().significanceMatrix.at(2);
        syncTree.metcov11() = event.metPF().significanceMatrix.at(3);
        syncTree.mvacov00() = correctedMET.significanceMatrix.at(0);
        syncTree.mvacov01() = correctedMET.significanceMatrix.at(1);
        syncTree.mvacov10() = correctedMET.significanceMatrix.at(2);
        syncTree.mvacov11() = correctedMET.significanceMatrix.at(3);

        const TLorentzVector met;
        met.SetPtEtaPhiM(correctedMET.pt, correctedMET.eta, correctedMET.phi, 0.);
        const double Mt = std::sqrt(2*muon.momentum.Pt()*met.Pt()*(1-std::cos(muon.momentum.DeltaPhi(met)))); //see AN-13-178
        syncTree.mt_1() = Mt;
        syncTree.njets() = jets.size();
        syncTree.nbtag() = bjets.size();
        if (jets.size() >= 1){
            const ntuple::Jet& ntuple_jet = event.jets().at(jets.at(0).index);
            syncTree.jpt_1() = jets.at(0).momentum.Pt();
            syncTree.jeta_1() = jets.at(0).momentum.Eta();
            syncTree.jphi_1() = jets.at(0).momentum.Phi();
            syncTree.jmva_1() = ntuple_jet.puIdMVA;
        }
        else {
            syncTree.jpt_1() = -10000;
            syncTree.jeta_1() = -10000;
            syncTree.jphi_1() = -10000;
            syncTree.jmva_1() = -10000;
        }
        if (jets.size() >= 2){
            const ntuple::Jet& ntuple_jet = event.jets().at(jets.at(1).index);
            syncTree.jpt_2() = jets.at(1).momentum.Pt();
            syncTree.jeta_2() = jets.at(1).momentum.Eta();
            syncTree.jphi_2() = jets.at(1).momentum.Phi();
            syncTree.jmva_2() = ntuple_jet.puIdMVA;
        }
        else {
            syncTree.jpt_2() = -10000;
            syncTree.jeta_2() = -10000;
            syncTree.jphi_2() = -10000;
            syncTree.jmva_2() = -10000;
        }
        if (bjets.size() >= 1){
            syncTree.bpt() = bjets.at(0).momentum.Pt();
            syncTree.beta() = bjets.at(0).momentum.Eta();
            syncTree.bphi() = bjets.at(0).momentum.Phi();
        }
        else {
            syncTree.bpt() = -10000;
            syncTree.beta() = -10000;
            syncTree.bphi() = -10000;
        }
        syncTree.Fill();
    }

private:
    BaselineAnalyzerData anaData;
};

