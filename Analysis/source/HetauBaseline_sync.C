/*!
 * \file HtautauBaseline_sync.C
 * \brief Generate sync-tree for H->tautau->e_taujet analysis using baseline selection.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-05 created
 */

#include "../include/H_BaseAnalyzer.h"

class BaselineAnalyzerData : public analysis::BaseAnalyzerData {
public:
    BaselineAnalyzerData(TFile& outputFile) : BaseAnalyzerData(outputFile) {}

};

class HetauBaseline_sync : public analysis::H_BaseAnalyzer {
public:
    HetauBaseline_sync(const std::string& inputFileName, const std::string& outputFileName,
                       const std::string& _prefix = "none", Long64_t _maxNumberOfEvents = 0,
                       bool _useMCtruth = false, const std::string& reweightFileName = "none")
        : H_BaseAnalyzer(inputFileName,outputFileName,_prefix,_maxNumberOfEvents,_useMCtruth, reweightFileName),
          anaData(*outputFile)
    {
        anaData.getOutputFile().cd();
    }

    ~HetauBaseline_sync()
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
        using namespace cuts::Htautau_Summer13;
        using namespace cuts::Htautau_Summer13::ETau;
        finalState::ETaujet eTau;
        if (useMCtruth && !FindAnalysisFinalState(eTau)) return;

        ApplyTauCorrections(eTau);

        cuts::Cutter cut(&anaData.Selection("event"));
        cut(true, "total");

        cut(HaveTriggerFired(trigger::hltPaths), "trigger");

        const VertexVector vertices = CollectVertices();
        cut(vertices.size(), "vertex");
        primaryVertex = vertices.front();

        const auto z_electrons = CollectZelectrons();
        const auto z_electron_candidates = FindCompatibleObjects(z_electrons, ZeeVeto::deltaR,
                                                                 Candidate::Z, "Z_e_e", 0);
        cut(!z_electron_candidates.size(), "z_ee_veto");

        const auto muons_bkg = CollectBackgroundMuons();
        cut(!muons_bkg.size(), "no_muon");

        const auto electrons = CollectElectrons();
        cut(electrons.size(), "electron_cand");
        cut(electrons.size() == 1, "one_electron_cand");

        const auto electrons_bkg = CollectBackgroundElectrons();
        const bool have_bkg_electron = electrons_bkg.size() > 1 ||
                ( electrons_bkg.size() == 1 && electrons_bkg.front() != electrons.front() );
        cut(!have_bkg_electron, "no_bkg_electron");

        const auto taus = CollectTaus();
        cut(taus.size(), "tau_cand");

        const auto higgses = FindCompatibleObjects(electrons, taus, DeltaR_betweenSignalObjects,
                                                   Candidate::Higgs, "H_e_tau", 0);
        cut(higgses.size(), "e_tau");

        const Candidate higgs = SelectSemiLeptonicHiggs(higgses);

        const auto jets = CollectJets(higgs);
        const auto bjets = CollectBJets(higgs);
        const Candidate higgs_corr = ApplyCorrections(higgs, eTau.resonance, jets.size());

        postRecoilMET = correctedMET;
        FillSyncTree(higgs, higgs_corr, jets, bjets, vertices);
    }

    virtual analysis::Candidate SelectElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                               root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::electronID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event.electrons().at(id);

        cut(true, ">0 ele cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        const double eta = std::abs( X(eta, 120, -6.0, 6.0) );
        cut(eta < eta_high, "eta");
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

        const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
        cut(Y(haveTriggerMatch, 2, -0.5, 1.5), "triggerMatch");
        return analysis::Candidate(analysis::Candidate::Electron, id, object,object.charge);
    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::tauID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Tau& object = correctedTaus.at(id);

        cut(true, ">0 tau cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(decayModeFinding, 2, -0.5, 1.5) > decayModeFinding, "decay_mode");
        cut(X(againstMuonLoose, 2, -0.5, 1.5) > againstMuonLoose, "vs_mu_loose");
        cut(X(againstElectronMediumMVA3, 2, -0.5, 1.5) > againstElectronMediumMVA3, "vs_e_mediumMVA");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits, 100, 0, 10)
            < byCombinedIsolationDeltaBetaCorrRaw3Hits, "looseIso3Hits");

        const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
        cut(Y(haveTriggerMatch, 2, -0.5, 1.5), "triggerMatch");

        return analysis::Candidate(analysis::Candidate::Tau, id, object,object.charge);
    }

    analysis::CandidateVector CollectZelectrons()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> analysis::Candidate
            { return SelectZelectron(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<analysis::Candidate>("z_electrons", base_selector, event.electrons().size());
    }

    virtual analysis::Candidate SelectZelectron(size_t id, cuts::ObjectSelector* objectSelector,
                                                root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::ETau::ZeeVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event.electrons().at(id);

        cut(true, ">0 mu cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");
        const TVector3 e_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = (e_vertex - primaryVertex.position).Perp();
        cut(std::abs( Y(d0_PV, 50, 0.0, 0.5) ) < d0, "d0");
        cut(X(pfRelIso, 1000, 0.0, 100.0) < pfRelIso, "pFRelIso");
        const size_t eta_index = eta <= barrel_eta_high ? barrel_index : endcap_index;
        cut(X(sigmaIEtaIEta, 100, 0.0, 0.1) < sigma_ieta_ieta[eta_index], "sigmaIetaIeta");
        cut(X(deltaEtaTrkSC, 100, 0.0, 0.1) < delta_eta[eta_index], "deltaEtaSC");
        cut(X(deltaPhiTrkSC, 200, 0.0, 2.0) < delta_phi[eta_index], "deltaPhiSC");
        cut(X(eSuperClusterOverP, 100, 0.0, 1.0) < HoverE[eta_index], "HoverE");

        return analysis::Candidate(analysis::Candidate::Electron, id, object, object.charge);
    }

//    bool FindAnalysisFinalState(analysis::finalState::ETaujet& final_state)
//    {
//        if(!H_BaseAnalyzer::FindAnalysisFinalState(final_state))
//            return false;

//        final_state.electron = final_state.tau_jet = nullptr;
//        for (const analysis::GenParticle* tau_MC : final_state.taus) {
//            analysis::GenParticlePtrVector tauProducts;
//            if (analysis::FindDecayProducts(*tau_MC, analysis::TauElectronDecay, tauProducts))
//                final_state.electron = tauProducts.at(0);
//            else if (!analysis::FindDecayProducts(*tau_MC, analysis::TauMuonicDecay, tauProducts))
//                final_state.tau_jet = tau_MC;
//        }

//        if (!final_state.electron || !final_state.tau_jet) return false;
//        return true;
//    }

    void FillSyncTree(const analysis::Candidate& higgs, const analysis::Candidate& higgs_corr,
                      const analysis::CandidateVector& jets, const analysis::CandidateVector& bjets,
                      const analysis::VertexVector& vertices)
    {
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        const ntuple::Tau& ntuple_tau = correctedTaus.at(tau.index);
        H_BaseAnalyzer::FillSyncTree(higgs, higgs_corr, jets, bjets, vertices, tau);

        const analysis::Candidate& electron = higgs.GetDaughter(analysis::Candidate::Electron);
        const ntuple::Electron& ntuple_electron = event.electrons().at(electron.index);

        // electron
        syncTree.pt_1() = electron.momentum.Pt();
        syncTree.eta_1() = electron.momentum.Eta();
        syncTree.phi_1() = electron.momentum.Phi();
        syncTree.m_1() = electron.momentum.M();
        syncTree.q_1() = electron.charge;
        syncTree.iso_1() = ntuple_electron.pfRelIso;
        const TVector3 electron_vertex(ntuple_electron.vx, ntuple_electron.vy, ntuple_electron.vz);
        syncTree.d0_1() = (electron_vertex - primaryVertex.position).Perp();
        syncTree.dZ_1() = std::abs(ntuple_electron.vz - primaryVertex.position.Z());

        Double_t DMweight = 1;
        if (ntuple_tau.decayMode == 0)
            DMweight *= 0.88;
        syncTree.decaymodeweight() = DMweight;

        TLorentzVector met;
        met.SetPtEtaPhiM(correctedMET.pt, correctedMET.eta, correctedMET.phi, 0.);
        //see AN-2013/178
        syncTree.mt_1() = std::sqrt(2*electron.momentum.Pt()*met.Pt()*(1-std::cos(electron.momentum.DeltaPhi(met))));

        syncTree.Fill();
    }

private:
    BaselineAnalyzerData anaData;
};

