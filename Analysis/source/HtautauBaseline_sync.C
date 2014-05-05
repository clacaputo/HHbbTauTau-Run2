/*!
 * \file HtautauBaseline_sync.C
 * \brief Generate sync-tree for Htautau analysis using baseline selection.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-05 created
 */

#include "../include/BaseAnalyzer.h"
#include "../include/SyncTree.h"

class BaselineAnalyzerData : public analysis::SignalAnalyzerData {
public:
    BaselineAnalyzerData(TFile& outputFile) : SignalAnalyzerData(outputFile) {}

};

class HtautauBaseline_sync : public analysis::BaseAnalyzer {
public:
    HtautauBaseline_sync(const std::string& inputFileName, const std::string& outputFileName,
                          const std::string& _prefix = "none", Long64_t _maxNumberOfEvents = 0,
                          bool _useMCtruth = false, const std::string& reweightFileName = "none")
        : BaseAnalyzer(inputFileName,outputFileName,_prefix,_maxNumberOfEvents,_useMCtruth, reweightFileName),
          anaData(*outputFile), etauTree("etauTree"), mutauTree("mutauTree"), tautauTree("tautauTree")
    {
        anaData.getOutputFile().cd();
    }

    ~HtautauBaseline_sync()
    {
        etauTree.Write();
        mutauTree.Write();
        tautauTree.Write();
    }

protected:
    virtual analysis::SignalAnalyzerData& GetAnaData() { return anaData; }

    virtual void ProcessEvent()
    {
        BaseAnalyzer::ProcessEvent();
        using namespace analysis;
        finalState::TauTau tauTau;
        if (useMCtruth && !FindAnalysisFinalState(tauTau)) return;

        cuts::Cutter cut(anaData.EventSelection());
        cut(true, "total");

        const VertexVector vertices = CollectVertices();
        cut(vertices.size(), "vertex");
        primaryVertex = vertices.back();

        //SIGNAL OBJECTS
        const auto electrons = CollectElectrons(true);
        const auto muons = CollectMuons(true);
        const auto taus_eTau = CollectTaus_eTau();
        const auto taus_muTau = CollectTaus_muTau();
        const auto taus_tauTau = CollectTaus_tauTau();

        // BACKGROUND OBJECTS
        const auto electrons_bkg = CollectElectrons(false);
        const auto muons_bkg = CollectMuons(false);
        const auto taus_bkg = CollectTaus(false);

        const auto Higgses_e_tau_signal = FindCompatibleObjects(electrons, taus_eTau,
                   cuts::Htautau_Summer13::DeltaR_betweenSignalObjects,analysis::Candidate::Higgs, "H_e_tau", 0);
        const auto Higgses_mu_tau_signal = FindCompatibleObjects(muons, taus_muTau,
                   cuts::Htautau_Summer13::DeltaR_betweenSignalObjects, analysis::Candidate::Higgs, "H_mu_tau", 0);
        const auto Higgses_tautau_loose = FindCompatibleObjects(taus_tauTau,
                   cuts::Htautau_Summer13::DeltaR_betweenSignalObjects, analysis::Candidate::Higgs, "H_2tau", 0);
        const auto Higgses_tautau_signal = ApplyTauFullSelection(Higgses_tautau_loose);

        const auto Higgses_e_tau = ApplyVetos(Higgses_e_tau_signal, "e_tau", electrons_bkg, muons_bkg, taus_bkg);
        const auto Higgses_mu_tau = ApplyVetos(Higgses_mu_tau_signal, "mu_tau", electrons_bkg, muons_bkg, taus_bkg);
        const auto Higgses_tautau = ApplyVetos(Higgses_tautau_signal, "tautau", electrons_bkg, muons_bkg, taus_bkg);

        finalState::ETaujet eTaujet(tauTau);
        finalState::MuTaujet muTaujet(tauTau);
        finalState::TaujetTaujet taujetTaujet(tauTau);
        if(Higgses_e_tau.size() && HaveTriggerFired(cuts::Htautau_Summer13::trigger::ETau::hltPaths) &&
                ( !useMCtruth || FindAnalysisFinalState_eTau(eTaujet) ) )
            FillSyncTree(etauTree, SelectSemiLeptonicHiggs(Higgses_e_tau));
        else if(Higgses_mu_tau.size() && HaveTriggerFired(cuts::Htautau_Summer13::trigger::MuTau::hltPaths) &&
                ( !useMCtruth || FindAnalysisFinalState_muTau(muTaujet) ) )
            FillSyncTree(mutauTree, SelectSemiLeptonicHiggs(Higgses_mu_tau));
        else if(Higgses_tautau.size() && HaveTriggerFired(cuts::Htautau_Summer13::trigger::TauTau::hltPaths) &&
                ( !useMCtruth || FindAnalysisFinalState_tauTau(taujetTaujet) ) )
            FillSyncTree(tautauTree, SelectFullyHadronicHiggs(Higgses_tautau));
    }

    analysis::CandidateVector CollectTaus_eTau()
    {
        const BaseSelector base_selector = [&](unsigned id, cuts::ObjectSelector& _objectSelector,
                bool enabled, root_ext::AnalyzerData& _anaData)
                -> analysis::Candidate { return SelectTau_eTau(id, _objectSelector, enabled, _anaData); };
        return CollectObjects<analysis::Candidate>(anaData.TauSelection("eTau"), event.taus().size(),
                                                   base_selector, "taus_eTau");
    }

    analysis::CandidateVector CollectTaus_muTau()
    {
        const BaseSelector base_selector = [&](unsigned id, cuts::ObjectSelector& _objectSelector,
                bool enabled, root_ext::AnalyzerData& _anaData)
                -> analysis::Candidate { return SelectTau_muTau(id, _objectSelector, enabled, _anaData); };
        return CollectObjects<analysis::Candidate>(anaData.TauSelection("muTau"), event.taus().size(),
                                                   base_selector, "taus_muTau");
    }

    analysis::CandidateVector CollectTaus_tauTau()
    {
        const BaseSelector base_selector = [&](unsigned id, cuts::ObjectSelector& _objectSelector,
                bool enabled, root_ext::AnalyzerData& _anaData)
                -> analysis::Candidate { return SelectTau_tauTau(id, _objectSelector, enabled, _anaData); };
        return CollectObjects<analysis::Candidate>(anaData.TauSelection("tauTau"), event.taus().size(),
                                                   base_selector, "taus_tauTau");
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

    virtual analysis::Candidate SelectMuon(size_t id, cuts::ObjectSelector& objectSelector,
                                           bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::muonID::MuTau;
        const std::string selection_label = "muon";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Muon& object = event.muons().at(id);

        cut(true, ">0 mu cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(isGlobalMuonPromptTight, 2, -0.5, 1.5) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon, 2, -0.5, 1.5) == isPFMuon, "PF");
        cut(X(nMatchedStations, 10, 0.0, 10.0) > nMatched_Stations, "stations");
        cut(X(trackerLayersWithMeasurement, 20, 0.0, 20.0) > trackerLayersWithMeasurement, "layers");
        cut(X(pixHits, 10, 0.0, 10.0) > pixHits, "pix_hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double dB_PV = (mu_vertex - primaryVertex.position).Perp();
        cut(std::abs( Y(dB_PV, 50, 0.0, 0.5) ) < dB, "dB");
        cut(X(pfRelIso, 1000, 0.0, 100.0) < pFRelIso, "pFRelIso");

        return analysis::Candidate(analysis::Candidate::Mu, id, object,object.charge);
    }

    analysis::Candidate SelectTau_eTau(size_t id, cuts::ObjectSelector& objectSelector,
                                       bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::tauID::ETau;
        const std::string selection_label = "tau_eTau";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Tau& object = event.taus().at(id);

        cut(true, ">0 tau cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(decayModeFinding, 2, -0.5, 1.5) > decayModeFinding, "decay_mode");
        cut(X(byLooseCombinedIsolationDeltaBetaCorr3Hits, 2, -0.5, 1.5)
            > LooseCombinedIsolationDeltaBetaCorr3Hits, "looseIso3Hits");
        cut(X(againstMuonLoose, 2, -0.5, 1.5) > againstMuonLoose, "vs_mu_loose");
        cut(X(againstElectronMediumMVA3, 2, -0.5, 1.5) > againstElectronMediumMVA3, "vs_e_mediumMVA");

        return analysis::Candidate(analysis::Candidate::Tau, id, object,object.charge);
    }

    analysis::Candidate SelectTau_muTau(size_t id, cuts::ObjectSelector& objectSelector,
                                        bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::tauID::MuTau;
        const std::string selection_label = "tau_muTau";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Tau& object = event.taus().at(id);

        cut(true, ">0 tau cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(decayModeFinding, 2, -0.5, 1.5) > decayModeFinding, "decay_mode");
        cut(X(byLooseCombinedIsolationDeltaBetaCorr3Hits, 2, -0.5, 1.5)
            > LooseCombinedIsolationDeltaBetaCorr3Hits, "looseIso3Hits");
        cut(X(againstMuonTight, 2, -0.5, 1.5) > againstMuonTight, "vs_mu_tight");
        cut(X(againstElectronLoose, 2, -0.5, 1.5) > againstElectronLoose, "vs_e_loose");

        return analysis::Candidate(analysis::Candidate::Tau, id, object,object.charge);
    }

    analysis::Candidate SelectTau_tauTau(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                         root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::tauID::TauTau;
        const std::string selection_label = "tau_tauTau";
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

        return analysis::Candidate(analysis::Candidate::Tau, id, object,object.charge);
    }

    analysis::CandidateVector ApplyTauFullSelection(const analysis::CandidateVector& higgses)
    {
        using namespace analysis;
        using namespace cuts::Htautau_Summer13::tauID::TauTau;
        CandidateVector result;
        for(const Candidate& higgs : higgses) {
            const Candidate* leading_tau = higgs.GetLeadingDaughter(Candidate::Tau);
            const Candidate* subleading_tau = higgs.GetSubleadingDaughter(Candidate::Tau);
            const ntuple::Tau& ntuple_subleadingTau = event.taus().at(subleading_tau->index);
            if(ntuple_subleadingTau.againstElectronLooseMVA3 > againstElectronLooseMVA3
                    && leading_tau->momentum.Pt() > pt_lead)
                result.push_back(higgs);
        }
        return result;
    }

    analysis::CandidateVector ApplyVetos(const analysis::CandidateVector& resonances, const std::string& suffix,
                                         const analysis::CandidateVector& electrons_bkg,
                                         const analysis::CandidateVector& muons_bkg,
                                         const analysis::CandidateVector& taus_bkg)
    {
        const auto resonances_noEle = FilterBackground(resonances,electrons_bkg,
                           cuts::Htautau_Summer13::electronID::veto::deltaR_signalObjects, "resonances_noEle" + suffix);
        const auto resonances_noMu = FilterBackground(resonances_noEle,muons_bkg,
                           cuts::Htautau_Summer13::muonID::veto::deltaR_signalObjects, "resonances_noMu" + suffix);
        const auto resonances_noTau = FilterBackground(resonances_noMu,taus_bkg,
                           cuts::Htautau_Summer13::tauID::veto::deltaR_signalObjects, "resonances_noTau" + suffix);
        return resonances_noTau;
    }

    const analysis::Candidate& SelectSemiLeptonicHiggs(const analysis::CandidateVector& higgses)
    {
        if(!higgses.size())
            throw std::runtime_error("no available higgs candidate to select");
        return *std::max_element(higgses.begin(), higgses.end());
    }

    const analysis::Candidate& SelectFullyHadronicHiggs(const analysis::CandidateVector& higgses)
    {
        if(!higgses.size())
            throw std::runtime_error("no available higgs candidate to select");
        const auto higgsSelector = [&] (const analysis::Candidate& first, const analysis::Candidate& second) -> bool
        {
            const ntuple::Tau& first_tau1 = event.taus().at(first.daughters.at(0)->index);
            const ntuple::Tau& first_tau2 = event.taus().at(first.daughters.at(1)->index);
            const ntuple::Tau& second_tau1 = event.taus().at(second.daughters.at(0)->index);
            const ntuple::Tau& second_tau2 = event.taus().at(second.daughters.at(1)->index);
            const double first_iso = std::max(first_tau1.byCombinedIsolationDeltaBetaCorrRaw3Hits,
                                              first_tau2.byCombinedIsolationDeltaBetaCorrRaw3Hits);
            const double second_iso = std::max(second_tau1.byCombinedIsolationDeltaBetaCorrRaw3Hits,
                                               second_tau2.byCombinedIsolationDeltaBetaCorrRaw3Hits);
            return first_iso < second_iso;
        };
        return *std::min_element(higgses.begin(), higgses.end(), higgsSelector);
    }

    bool FindAnalysisFinalState(analysis::finalState::TauTau& final_state)
    {
        static const analysis::ParticleCodes resonanceCodes = { particles::Higgs, particles::Z };
        static const analysis::ParticleCodes resonanceDecay = { particles::tau, particles::tau };

        genEvent.Initialize(event.genParticles());

        const analysis::GenParticleSet resonances = genEvent.GetParticles(resonanceCodes);
        if (resonances.size() != 1)
            throw std::runtime_error("not one resonance per event");

        final_state.resonance = *resonances.begin();

        analysis::GenParticlePtrVector resonanceDecayProducts;
        if(!analysis::FindDecayProducts(*final_state.resonance, resonanceDecay,resonanceDecayProducts)) {
            std::cout << "event id = " << event.eventId().eventId << std::endl;
//            genEvent.PrintChain(final_state.resonance);
//            throw std::runtime_error("Resonance does not decay into 2 taus");
            std::cerr << "Resonance does not decay into 2 taus\n";
            return false;
        }

        final_state.taus = resonanceDecayProducts;

        return true;
    }

    bool FindAnalysisFinalState_eTau(analysis::finalState::ETaujet& final_state)
    {
        final_state.electron = final_state.tau_jet = nullptr;
        for (const analysis::GenParticle* tau_MC : final_state.taus) {
            analysis::GenParticlePtrVector tauProducts;
            if (analysis::FindDecayProducts(*tau_MC, analysis::TauElectronDecay, tauProducts))
                final_state.electron = tauProducts.at(0);
            else if (!analysis::FindDecayProducts(*tau_MC, analysis::TauMuonicDecay, tauProducts))
                final_state.tau_jet = tau_MC;
        }

        if (!final_state.electron || !final_state.tau_jet) return false;
        return true;
    }

    bool FindAnalysisFinalState_muTau(analysis::finalState::MuTaujet& final_state)
    {
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

    bool FindAnalysisFinalState_tauTau(analysis::finalState::TaujetTaujet& final_state)
    {
        unsigned n_hadronic_taus = 0;

        if(final_state.taus.size() != 2)
            throw std::runtime_error("bad tautau MC final state");

        for (const analysis::GenParticle* tau_MC : final_state.taus) {
            analysis::GenParticlePtrVector TauProducts;
            if (!analysis::FindDecayProducts(*tau_MC, analysis::TauMuonicDecay, TauProducts)
                    && !analysis::FindDecayProducts(*tau_MC, analysis::TauElectronDecay, TauProducts))
                ++n_hadronic_taus;
        }

        if (n_hadronic_taus != 2) return false;
        if(final_state.taus.at(0)->momentum.Pt() > final_state.taus.at(1)->momentum.Pt()) {
            final_state.leading_tau_jet = final_state.taus.at(0);
            final_state.subleading_tau_jet = final_state.taus.at(1);
        } else {
            final_state.leading_tau_jet = final_state.taus.at(1);
            final_state.subleading_tau_jet = final_state.taus.at(0);
        }
        return true;
    }

    void FillSyncTree(ntuple::SyncTree& syncTree, const analysis::Candidate& higgs)
    {
        syncTree.run() = event.eventInfo().run;
        syncTree.lumi() = event.eventInfo().lumis;
        syncTree.evt() = event.eventInfo().EventId;
        for (unsigned n = 0; n < event.eventInfo().bunchCrossing.size(); ++n){
            if (event.eventInfo().bunchCrossing.at(n) == 0){
                syncTree.npu() = event.eventInfo().nPU.at(n); //only in-time PU
            }
        }
        syncTree.puweight() = weight;
        Double_t DMweight = 1;
        const analysis::Candidate* cand1 = higgs.finalStateDaughters.at(0);
        const analysis::Candidate* cand2 = higgs.finalStateDaughters.at(1);
        if (cand1->momentum.Pt() < cand2->momentum.Pt()){
            cand1 = higgs.finalStateDaughters.at(1);
            cand2 = higgs.finalStateDaughters.at(0);
        }
        const ntuple::Tau& leg1 = event.taus().at(cand1->index);
        const ntuple::Tau& leg2 = event.taus().at(cand2->index);
        if (leg1.decayMode == 0)
            DMweight *= 0.88;
        if (leg2.decayMode == 0)
            DMweight *= 0.88;
        syncTree.decaymodeweight() = DMweight;
        syncTree.mvis() = higgs.momentum.M();

        syncTree.pt_1() = cand1->momentum.Pt();
        syncTree.eta_1() = cand1->momentum.Eta();
        syncTree.phi_1() = cand1->momentum.Phi();
        syncTree.m_1() = cand1->momentum.M();
        syncTree.q_1() = leg1.charge;
        syncTree.iso_1() = leg1.byIsolationMVAraw;
        syncTree.passid_1() = 1;
        syncTree.passiso_1() = 1;
        syncTree.mt_1() = cand1->momentum.Mt();
        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = leg1.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        syncTree.againstElectronMVA3raw_1() = leg1.againstElectronMVA3raw;
        syncTree.againstElectronMVA3category_1() = leg1.againstElectronMVA3category;
        syncTree.byIsolationMVA2raw_1() = leg1.byIsolationMVA2raw;
        syncTree.againstMuonLoose_1() = leg1.againstMuonLoose;
        syncTree.againstMuonLoose2_1() = leg1.againstMuonLoose2;
        syncTree.againstMuonMedium2_1() = leg1.againstMuonMedium2;
        syncTree.againstMuonTight2_1() = leg1.againstMuonTight2;
        syncTree.againstElectronLooseMVA3_1() = leg1.againstElectronLooseMVA3;
        syncTree.againstElectronLoose_1() = leg1.againstElectronLoose;

        syncTree.pt_2() = cand2->momentum.Pt();
        syncTree.eta_2() = cand2->momentum.Eta();
        syncTree.phi_2() = cand2->momentum.Phi();
        syncTree.m_2() = cand2->momentum.M();
        syncTree.q_2() = leg2.charge;
        syncTree.iso_2() = leg2.byIsolationMVAraw;
        syncTree.passid_2() = 1;
        syncTree.passiso_2() = 1;
        syncTree.mt_2() = cand2->momentum.Mt();
        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = leg2.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        syncTree.againstElectronMVA3raw_2() = leg2.againstElectronMVA3raw;
        syncTree.againstElectronMVA3category_2() = leg2.againstElectronMVA3category;
        syncTree.byIsolationMVA2raw_2() = leg2.byIsolationMVA2raw;
        syncTree.againstMuonLoose_2() = leg2.againstMuonLoose;
        syncTree.againstMuonLoose2_2() = leg2.againstMuonLoose2;
        syncTree.againstMuonMedium2_2() = leg2.againstMuonMedium2;
        syncTree.againstMuonTight2_2() = leg2.againstMuonTight2;
        syncTree.againstElectronLooseMVA3_2() = leg2.againstElectronLooseMVA3;
        syncTree.againstElectronLoose_2() = leg2.againstElectronLoose;

        syncTree.met() = event.metPF().pt_uncorrected; //raw
        syncTree.metphi() = event.metPF().phi_uncorrected; //raw
        syncTree.mvamet() = event.metMVA().pt;
        syncTree.mvametphi() = event.metMVA().phi;
        syncTree.metcov00() = event.metPF().significanceMatrix.at(0);
        syncTree.metcov01() = event.metPF().significanceMatrix.at(1);
        syncTree.metcov10() = event.metPF().significanceMatrix.at(2);
        syncTree.metcov11() = event.metPF().significanceMatrix.at(3);
        syncTree.mvacov00() = event.metMVA().significanceMatrix.at(0);
        syncTree.mvacov01() = event.metMVA().significanceMatrix.at(1);
        syncTree.mvacov10() = event.metMVA().significanceMatrix.at(2);
        syncTree.mvacov11() = event.metMVA().significanceMatrix.at(3);

        syncTree.Fill();
    }

private:
    BaselineAnalyzerData anaData;
    ntuple::SyncTree etauTree, mutauTree, tautauTree;
};

