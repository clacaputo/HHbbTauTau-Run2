/*!
 * \file HetauFlatTreeProducer.C
 * \brief Generate flat-tree for H->tautau->e_taujet analysis using looser selection.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-07-11 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/BaseFlatTreeProducer.h"

class HHbbetau_FlatTreeProducer : public analysis::BaseFlatTreeProducer {
public:
    HHbbetau_FlatTreeProducer(const std::string& inputFileName, const std::string& outputFileName,
                              const std::string& configFileName, const std::string& _prefix = "none",
                              size_t _maxNumberOfEvents = 0,
                              std::shared_ptr<ntuple::FlatTree> _flatTree = std::shared_ptr<ntuple::FlatTree>())
        : BaseFlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents, _flatTree),
          anaData(*outputFile)
    {
        anaData.getOutputFile().cd();
    }

    virtual analysis::BaseAnalyzerData& GetAnaData() override { return anaData; }

    virtual void ProcessEvent(std::shared_ptr<const analysis::EventDescriptor> _event) override
    {
        BaseFlatTreeProducer::ProcessEvent(_event);
        using namespace analysis;
        using namespace cuts::Htautau_Summer13;
        using namespace cuts::Htautau_Summer13::ETau;
        finalState::ETaujet eTau;

        if (!FindAnalysisFinalState(eTau) && config.RequireSpecificFinalState()) return;

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

        const auto electrons = CollectSignalElectrons();
        cut(electrons.size(), "electron_cand");
        cut(electrons.size() == 1, "one_electron_cand");

        const auto electrons_bkg = CollectBackgroundElectrons();
        const bool have_bkg_electron = electrons_bkg.size() > 1 ||
                ( electrons_bkg.size() == 1 && electrons_bkg.front() != electrons.front() );
        cut(!have_bkg_electron, "no_bkg_electron");

        const auto allelectrons = CollectMuons();

        correctedTaus = config.ApplyTauESCorrection() ? ApplyTauCorrections(eTau.hadronic_taus,false) : event->taus();

        const auto alltaus = CollectTaus();
        cut(alltaus.size(), "tau_cand");


        const auto signaltaus = CollectSignalTaus() ;

        analysis::CandidateVector higgses ;

        // First check OS, isolated higgs candidates
        const auto higgses_sig = FindCompatibleObjects(electrons, signaltaus, DeltaR_betweenSignalObjects,
                                                       Candidate::Higgs, "H_e_tau",0);

        // If no OS candidate, keep any higgs-ish candidates for bkg estimation (don't cut on sign nor isolation)
        if (!higgses_sig.size())
        {
          const auto higgses_bkg = FindCompatibleObjects(allelectrons, alltaus, DeltaR_betweenSignalObjects,
                                                         Candidate::Higgs, "H_e_tau");
          higgses = higgses_bkg ;
        }
        else
        {
          higgses = higgses_sig ;
        }

        cut(higgses.size(), "e_tau");

        const auto higgsTriggered = ApplyTriggerMatch(higgses,trigger::hltPaths,false);
        cut(higgsTriggered.size(), "trigger obj match");


        const Candidate higgs = SelectSemiLeptonicHiggs(higgsTriggered);


        const ntuple::MET mvaMet = mvaMetProducer.ComputeMvaMet(higgs,event->pfCandidates(),
                                                                event->jets(),primaryVertex,
                                                                vertices,event->taus());

        correctedMET = config.ApplyTauESCorrection() ? ApplyTauCorrectionsToMVAMET(mvaMet, correctedTaus) : mvaMet;

        const auto looseJets = CollectLooseJets();

        const auto filteredLooseJets = FilterCompatibleObjects(looseJets, higgs,
                                                               cuts::Htautau_Summer13::jetID::deltaR_signalObjects);


        const auto jets = CollectJets(filteredLooseJets);
        const auto bjets = CollectBJets(filteredLooseJets, false);
        const auto retagged_bjets = CollectBJets(filteredLooseJets, config.isMC());


        postRecoilMET = ApplyRecoilCorrections(higgs, eTau.resonance, jets.size(), correctedMET);

        const double m_sv      = 1.;//CorrectMassBySVfit(higgs, postRecoilMET,1   );
        const double m_sv_Up   = 1.;//CorrectMassBySVfit(higgs, postRecoilMET,1.03);
        const double m_sv_Down = 1.;//CorrectMassBySVfit(higgs, postRecoilMET,0.97);

        CalculateFullEventWeight(higgs);

        const ntuple::MET pfMET = config.isMC() ? event->metPF() : mvaMetProducer.ComputePFMet(event->pfCandidates(), primaryVertex);


        FillFlatTree(higgs, m_sv, m_sv_Up, m_sv_Down, jets, filteredLooseJets, bjets, retagged_bjets, vertices, pfMET);
    }

protected:

    virtual analysis::Candidate SelectElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                               root_ext::AnalyzerData& _anaData,
                                               const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::electronID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event->electrons().at(id);
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
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[eta_index], "mva");
        //cut(X(pfRelIso) < cuts::skim::MuTau::pFRelIso , "pFRelIso");

//        const bool haveTriggerMatch = analysis::HaveTriggerMatched(object.matchedTriggerPaths, trigger::hltPaths);
//        cut(Y(haveTriggerMatch, 2, -0.5, 1.5), "triggerMatch");
        return analysis::Candidate(analysis::Candidate::Electron, id, object);
    }

    virtual analysis::Candidate SelectSignalElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                               root_ext::AnalyzerData& _anaData,
                                               const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::electronID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event->electrons().at(id);
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
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[eta_index], "mva");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");
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
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits) <
            cuts::skim::ETau::tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits, "looseIso3Hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");

        const analysis::Candidate tau(analysis::Candidate::Tau, id, object);

        return tau;
    }

    virtual analysis::Candidate SelectSignalTau(size_t id, cuts::ObjectSelector* objectSelector,
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

        return analysis::Candidate(analysis::Candidate::Tau, id, object);
    }

    virtual analysis::Candidate SelectBackgroundMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                                     root_ext::AnalyzerData& _anaData,
                                                     const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::muonVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Muon& object = event->muons().at(id);
        const analysis::Candidate muon(analysis::Candidate::Mu, id, object);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(mu_vertex,primaryVertex.position,muon.momentum);
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(isGlobalMuonPromptTight) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon) == isPFMuon, "PF");
        cut(X(nMatchedStations) > nMatched_Stations, "stations");
        cut(X(pixHits) > pixHits, "pix_hits");
        cut(X(trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");

        return muon;
    }

    virtual analysis::Candidate SelectBackgroundElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                                         root_ext::AnalyzerData& _anaData,
                                                         const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::electronVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event->electrons().at(id);
        const analysis::Candidate electron(analysis::Candidate::Electron, id, object);

        cut(true, ">0 ele cand");
        cut(X(pt) > pt, "pt");
        const double eta = std::abs( X(eta) );
        cut(eta < eta_high, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(ele_vertex,primaryVertex.position,electron.momentum); // same as dB
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        const size_t pt_index = object.pt < ref_pt ? 0 : 1;
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[pt_index][eta_index], "mva");
        cut(X(missingHits) < missingHits, "missingHits");
        cut(X(hasMatchedConversion) == hasMatchedConversion, "conversion");

        return electron;
    }

    analysis::CandidateVector CollectZelectrons()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> analysis::Candidate
            { return SelectZelectron(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<analysis::Candidate>("z_electrons", base_selector, event->electrons().size());
    }

    virtual analysis::Candidate SelectZelectron(size_t id, cuts::ObjectSelector* objectSelector,
                                                root_ext::AnalyzerData& _anaData,
                                                const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::ETau::ZeeVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event->electrons().at(id);
        const analysis::Candidate electron(analysis::Candidate::Electron, id, object);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(ele_vertex,primaryVertex.position,electron.momentum); // same as dB
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(pfRelIso) < pfRelIso, "pFRelIso");
        const size_t eta_index = std::abs(object.eta) <= barrel_eta_high ? barrel_index : endcap_index;
        cut(X(sigmaIEtaIEta) < sigma_ieta_ieta[eta_index], "sigmaIetaIeta");
        cut(X(deltaEtaTrkSC) < delta_eta[eta_index], "deltaEtaSC");
        cut(X(deltaPhiTrkSC) < delta_phi[eta_index], "deltaPhiSC");
        cut(X(eSuperClusterOverP) < HoverE[eta_index], "HoverE");

        return electron;
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

    bool FindAnalysisFinalState(analysis::finalState::ETaujet& final_state)
    {
        const bool base_result = BaseFlatTreeProducer::FindAnalysisFinalState(final_state);
        if(!base_result)
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

    void FillFlatTree(const analysis::Candidate& higgs, double m_sv,
                      double m_sv_Up, double m_sv_Down,
                      const analysis::CandidateVector& jets, const analysis::CandidateVector& jetsPt20,
                      const analysis::CandidateVector& bjets, const analysis::CandidateVector& retagged_bjets,
                      const analysis::VertexVector& vertices, const ntuple::MET& pfMET)
    {
        const analysis::Candidate& electron = higgs.GetDaughter(analysis::Candidate::Electron);
        const ntuple::Electron& ntuple_electron = event->electrons().at(electron.index);
        const analysis::Candidate& tau = higgs.GetDaughter(analysis::Candidate::Tau);

        BaseFlatTreeProducer::FillFlatTree(higgs, m_sv, m_sv_Up, m_sv_Down, jets, jetsPt20, bjets, retagged_bjets, vertices, electron, tau, pfMET);

        flatTree->channel() = static_cast<int>(ntuple::Channel::ETau);
        flatTree->pfRelIso_1() = ntuple_electron.pfRelIso;
        //flatTree->mva_1() = 0;
        //flatTree->passid_2();
        //flatTree->passiso_2();
        //flatTree->mva_2() = ntuple_tau.againstElectronLoose;
        flatTree->Fill();
    }

private:
    analysis::BaseAnalyzerData anaData;
};

#include "METPUSubtraction/interface/GBRProjectDict.cxx"
