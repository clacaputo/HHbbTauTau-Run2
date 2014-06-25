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
    ENTRY_1D(double, Tau_Zele_deltaR)
};

class HetauBaseline_sync : public analysis::H_BaseAnalyzer {
public:
    HetauBaseline_sync(const std::string& inputFileName, const std::string& outputFileName,
                       const std::string& configFileName, const std::string& _prefix = "none",
                       size_t _maxNumberOfEvents = 0)
        : H_BaseAnalyzer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          anaData(*outputFile)
    {
        anaData.getOutputFile().cd();
    }

    virtual ~HetauBaseline_sync() override
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
        using namespace cuts::Htautau_Summer13::ETau;
        finalState::ETaujet eTau;
        if (useMCtruth && !FindAnalysisFinalState(eTau,false,true)) return;

        cuts::Cutter cut(&anaData.Selection("event"));
        cut(true, "total");

        cut(HaveTriggerFired(trigger::hltPaths), "trigger");

        const VertexVector vertices = CollectVertices();
        cut(vertices.size(), "vertex");
        primaryVertex = vertices.front();

        const auto z_electrons = CollectZelectrons();
//        cut(analysis::AllCandidatesHaveSameCharge(z_electrons), "no_electron_OSpair");
        const auto z_electron_candidates = FindCompatibleObjects(z_electrons, ZeeVeto::deltaR,Candidate::Z, "Z_e_e", 0);
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

        ApplyTauCorrections(eTau, false);

        const auto taus = CollectTaus();
        cut(taus.size(), "tau_cand");

        const auto higgses = FindCompatibleObjects(electrons, taus, DeltaR_betweenSignalObjects,
                                                   Candidate::Higgs, "H_e_tau", 0);
        cut(higgses.size(), "e_tau");

        const auto higgsTriggered = ApplyTriggerMatch(higgses, trigger::hltPaths,false);
        cut(higgsTriggered.size(), "trigger obj match");

//        const auto higgsZveto = ApplyZVeto(higgsTriggered, z_electrons);
//        cut(higgsZveto.size(), "Z veto with dR");

        const Candidate higgs = SelectSemiLeptonicHiggs(higgsTriggered);
//        if(!analysis::AllCandidatesHaveSameCharge(z_electrons)) {
//            const Candidate electron = higgs.GetDaughter(Candidate::Electron);
//            const Candidate tau = higgs.GetDaughter(Candidate::Tau);
//            double maxDeltaR = 0;
//            for(const Candidate& z_electron : z_electrons) {
//                if(electron.charge != z_electron.charge && tau.charge == z_electron.charge){
//                    const double DeltaR = tau.momentum.DeltaR(z_electron.momentum);
//                    anaData.Tau_Zele_deltaR().Fill(DeltaR);
//                    maxDeltaR = std::max(maxDeltaR,DeltaR);
//                }
//            }
//            cut(maxDeltaR < 0.03, "DR_eleTau");
//        }

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

        ApplyPostRecoilCorrections(higgs, eTau.resonance, jets.size());
        const double m_sv = CorrectMassBySVfit(higgs, postRecoilMET,1);

        CalculateFullEventWeight(higgs);

        FillSyncTree(higgs, m_sv, jets, filteredLooseJets, bjets, retagged_bjets, vertices);


    }

    virtual analysis::Candidate SelectElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                               root_ext::AnalyzerData& _anaData,
                                               const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::electronID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event.electrons().at(id);
        const analysis::Candidate electron(analysis::Candidate::Electron, id, object);

        cut(true, ">0 ele cand");
        cut(X(pt) > pt, "pt");
        const double eta = std::abs( X(eta) );
        cut(eta < eta_high, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        cut(X(missingHits) < missingHits, "missingHits");
        cut(X(hasMatchedConversion) == hasMatchedConversion, "conversion");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(ele_vertex,primaryVertex.position,electron.momentum); // same as dB
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[eta_index], "mva");

//        const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
//        cut(Y(haveTriggerMatch, 2, -0.5, 1.5), "triggerMatch");
        return analysis::Candidate(analysis::Candidate::Electron, id, object);
    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData,
                                          const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::tauID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Tau& object = correctedTaus.at(id);

        cut(true, ">0 tau cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(againstMuonLoose) > againstMuonLoose, "vs_mu_loose");
        const bool againstElectron =  againstElectronMediumMVA3_Custom(object);
        cut(Y(againstElectron), "vs_e_mediumMVA");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits) < byCombinedIsolationDeltaBetaCorrRaw3Hits, "looseIso3Hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
//        const TVector3 tau_vertex(object.vx, object.vy, object.vz);
//        const double dB_PV = (tau_vertex - primaryVertex.position).Perp();
//        cut(std::abs( Y(dB_PV) ) < dB, "dB");


//        const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
//        cut(Y(haveTriggerMatch, 2, -0.5, 1.5), "triggerMatch");

        return analysis::Candidate(analysis::Candidate::Tau, id, object);
    }

    bool againstElectronMediumMVA3_Custom(const ntuple::Tau& tau)
    {
        using namespace cuts::Htautau_Summer13::ETau::tauID;
        const int icut = std::round(tau.againstElectronMVA3category);
        if(icut < 0) return false;
		const size_t ucut = (size_t)icut;
        if(ucut >= againstElectronMediumMVA3_customValues.size()) return true;
        return tau.againstElectronMVA3raw > againstElectronMediumMVA3_customValues.at(ucut);
    }

    analysis::CandidateVector CollectZelectrons()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> analysis::Candidate
            { return SelectZelectron(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<analysis::Candidate>("z_electrons", base_selector, event.electrons().size());
    }

    virtual analysis::Candidate SelectZelectron(size_t id, cuts::ObjectSelector* objectSelector,
                                                root_ext::AnalyzerData& _anaData,
                                                const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::ETau::ZeeVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event.electrons().at(id);
        const analysis::Candidate electron(analysis::Candidate::Electron, id, object);

        cut(true, ">0 mu cand");
//        std::cout << "Z ee" <<  std::endl;
//        std::cout << "q = " << object.charge << std::endl;
//        std::cout << "pt = " << object.pt << std::endl;
        cut(X(pt) > pt, "pt");
//        std::cout << "eta = " << object.eta  << ", phi="<< object.phi << std::endl;
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
//        std::cout << "dz = " << DeltaZ << std::endl;
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(ele_vertex,primaryVertex.position,electron.momentum); // same as dB
//        std::cout << "d0 = " << d0_PV << std::endl;
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
//        std::cout << "pfRelIso = " << object.pfRelIso << std::endl;
        cut(X(pfRelIso) < pfRelIso, "pFRelIso");
        const size_t eta_index = eta <= barrel_eta_high ? barrel_index : endcap_index;
//        std::cout << "sigmaIEtaIEta = " << object.sigmaIEtaIEta << std::endl;
        cut(X(sigmaIEtaIEta) < sigma_ieta_ieta[eta_index], "sigmaIetaIeta");
//        std::cout << "deltaEtaTrkSC = " << object.deltaEtaTrkSC << std::endl;
        cut(X(deltaEtaTrkSC) < delta_eta[eta_index], "deltaEtaSC");
//        std::cout << "deltaPhiTrkSC = " << object.deltaPhiTrkSC << std::endl;
        cut(X(deltaPhiTrkSC) < delta_phi[eta_index], "deltaPhiSC");
//        std::cout << "eSuperClusterOverP = " << object.eSuperClusterOverP << std::endl;
        cut(X(eSuperClusterOverP) < HoverE[eta_index], "HoverE");

        return electron;
    }

    bool FindAnalysisFinalState(analysis::finalState::ETaujet& final_state, bool oneResonanceSelection, bool inclusive)
    {
        const bool base_result = H_BaseAnalyzer::FindAnalysisFinalState(final_state, oneResonanceSelection);
        if(inclusive || !base_result)
            return base_result;

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

    virtual void CalculateTriggerWeights(const analysis::Candidate& higgs) override
    {
        using namespace analysis::Htautau_Summer13::trigger::Run2012ABCD::ETau;
        const analysis::Candidate& ele = higgs.GetDaughter(analysis::Candidate::Electron);
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        triggerWeights = CalculateWeights(ele.momentum, tau.momentum);
    }

    virtual void CalculateIsoWeights(const analysis::Candidate& higgs) override
    {
        using namespace cuts::Htautau_Summer13::ETau::electronISOscaleFactor;
        IsoWeights.clear();
        const analysis::Candidate& ele = higgs.GetDaughter(analysis::Candidate::Electron);
        const double ele_pt = ele.momentum.Pt(), ele_eta = std::abs(ele.momentum.Eta());
        if(ele_pt < pt.at(0))
            throw std::runtime_error("No information about ISO. Electron pt is too small");
        const size_t pt_bin = ele_pt < pt.at(1) ? 0 : 1;
        if(ele_eta >= eta.at(1))
            throw std::runtime_error("No information about ISO. Electron eta is too big");
        const size_t eta_bin = ele_eta < eta.at(0) ? 0 : 1;
        const double scale = scaleFactors.at(pt_bin).at(eta_bin);
        // first ele, second tau
        IsoWeights.push_back(scale);
        IsoWeights.push_back(1);
    }

    virtual void CalculateIdWeights(const analysis::Candidate& higgs) override
    {
        using namespace cuts::Htautau_Summer13::ETau::electronIDscaleFactor;
        IDweights.clear();
        const analysis::Candidate& ele = higgs.GetDaughter(analysis::Candidate::Electron);
        const double ele_pt = ele.momentum.Pt(), ele_eta = std::abs(ele.momentum.Eta());
        if(ele_pt < pt.at(0))
            throw std::runtime_error("No information about ID. Electron pt is too small");
        const size_t pt_bin = ele_pt < pt.at(1) ? 0 : 1;
        if(ele_eta >= eta.at(1))
            throw std::runtime_error("No information about ID. Electron eta is too big");
        const size_t eta_bin = ele_eta < eta.at(0) ? 0 : 1;
        const double scale = scaleFactors.at(pt_bin).at(eta_bin);
        // first ele, second tau
        IDweights.push_back(scale);
        IDweights.push_back(1);
    }

    virtual void CalculateDMWeights(const analysis::Candidate& higgs) override
    {
        DMweights.clear();
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);
        const ntuple::Tau& tau_leg = correctedTaus.at(tau.index);
        const double weight = tau_leg.decayMode == ntuple::tau_id::kOneProng0PiZero ?
                                  cuts::Htautau_Summer13::tauCorrections::DecayModeWeight : 1;
        // first electron, second tau
        DMweights.push_back(1);
        DMweights.push_back(weight);
    }


    void FillSyncTree(const analysis::Candidate& higgs, double m_sv,
                      const analysis::CandidateVector& jets, const analysis::CandidateVector& jetsPt20,
                      const analysis::CandidateVector& bjets, const analysis::CandidateVector& retagged_bjets,
                      const analysis::VertexVector& vertices)
    {
        const analysis::Candidate& electron = higgs.GetDaughter(analysis::Candidate::Electron);
        const ntuple::Electron& ntuple_electron = event.electrons().at(electron.index);
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);

        H_BaseAnalyzer::
                FillSyncTree(higgs, m_sv, jets, jetsPt20, bjets, retagged_bjets, vertices, electron, tau);

        syncTree.iso_1() = ntuple_electron.pfRelIso;
        syncTree.mva_1() = ntuple_electron.mvaPOGNonTrig;
        //syncTree.passid_2();
        //syncTree.passiso_2();

        //syncTree.mva_2() = againstElectronMediumMVA3_Custom(ntuple_tau);

        syncTree.Fill();
    }

private:
    BaselineAnalyzerData anaData;
};

#include "METPUSubtraction/interface/GBRProjectDict.cxx"
