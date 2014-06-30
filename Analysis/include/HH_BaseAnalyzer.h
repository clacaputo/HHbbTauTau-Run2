/*!
 * \file HH_BaseAnalyzer.h
 * \brief Definition of HH_BaseAnalyzer class which is the base class for all X->HH->bbTauTau analyzers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-05-07 created
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

#pragma once

#include "BaseAnalyzer.h"

namespace analysis {

class HH_BaseAnalyzerData : public BaseAnalyzerData {
public:
    HH_BaseAnalyzerData(TFile& outputFile) : BaseAnalyzerData(outputFile) {}

    ENTRY_1D(float, Resonance_Mass)
    ENTRY_1D(float, Resonance_Pt)
    ENTRY_1D(float, Resonance_Eta)
    ENTRY_1D(float, Resonance_Phi)
    ENTRY_1D(float, Bjets_Pt_MC)
    ENTRY_1D(float, Higgs_leptonic_MC_Pt)
    ENTRY_1D(float, Higgs_BB_MC_Pt)
    ENTRY_2D(float, DR_bjets_vs_HiggsPt_MC)
    ENTRY_2D(float, DR_Higgs_vs_ResonancePt_MC)
};

class HH_BaseAnalyzer : public BaseAnalyzer {
public:
    HH_BaseAnalyzer(const std::string& inputFileName, const std::string& outputFileName,
                 const std::string& _prefix = "none", size_t _maxNumberOfEvents = 0, bool _useMCtruth = false,
                 const std::string& reweightFileName = "none")
        : BaseAnalyzer(inputFileName, outputFileName, _prefix, _maxNumberOfEvents, _useMCtruth, reweightFileName) {}

protected:
    CandidateVector CollectBJets(double csv, const std::string& selection_label, bool signal = true)
    {
        const BaseSelector base_selector_signal = [&](unsigned id, cuts::ObjectSelector& _objectSelector,
                bool enabled, root_ext::AnalyzerData& _anaData)
                -> Candidate { return SelectBJet(id, _objectSelector, enabled, _anaData, csv, selection_label); };
        const BaseSelector base_selector_bkg = [&](unsigned id, cuts::ObjectSelector& _objectSelector,
                bool enabled, root_ext::AnalyzerData& _anaData)
                -> Candidate { return SelectBackgroundBJet(id, _objectSelector, enabled, _anaData); };
        const auto base_selector = signal ? base_selector_signal : base_selector_bkg;
        auto& objectSelector = signal ? GetAnaData().BJetSelection(selection_label)
                                          : GetAnaData().BJetSelectionBkg(selection_label);
        return CollectObjects<Candidate>(objectSelector, event->jets().size(), base_selector, "bjets");
    }

    Candidate SelectBJet(size_t id, cuts::ObjectSelector& objectSelector, bool enabled, root_ext::AnalyzerData& _anaData,
                         double csv, const std::string& _selection_label)
    {
        using namespace cuts::Htautau_Summer13::btag::signal;
        cuts::Cutter cut(objectSelector, enabled);
        const std::string selection_label = "bjet_" + _selection_label;
        const ntuple::Jet& object = event->jets().at(id);
        cut(true, ">0 b-jet cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(combinedSecondaryVertexBJetTags, 130, -11.0, 2.0) > csv, "CSV");

        return analysis::Candidate(analysis::Candidate::Bjet, id, object);
    }

    Candidate SelectBackgroundBJet(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                   root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::btag::veto;
        const std::string selection_label = "bjet_bkg";
        cuts::Cutter cut(objectSelector, enabled);

        const ntuple::Jet& object = event->jets().at(id);
        cut(true, ">0 b-jet cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(combinedSecondaryVertexBJetTags, 130, -11.0, 2.0) > CSV, "CSV");
        cut(X(passLooseID, 2, -0.5, 1.5) == passLooseID, "passLooseID");

        return analysis::Candidate(analysis::Candidate::Bjet, id, object);
    }

    virtual analysis::Candidate SelectBackgroundElectron(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                                         root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::electronID::veto;
        const std::string selection_label = "electron_bkg";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Electron& object = event->electrons().at(id);

        cut(true, ">0 ele cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        const double eta = std::abs( X(eta, 120, -6.0, 6.0) );
        cut(eta < eta_high && (eta < cuts::Htautau_Summer13::electronID::eta_CrackVeto_low ||
                               eta > cuts::Htautau_Summer13::electronID::eta_CrackVeto_high), "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");
        cut(X(pfRelIso, 1000, 0.0, 100.0) < pFRelIso, "pFRelIso");
        const size_t pt_index = object.pt < ref_pt ? 0 : 1;
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig, 300, -1.5, 1.5) > MVApogNonTrig[pt_index][eta_index], "mva");

        return analysis::Candidate(analysis::Candidate::Electron, id, object,object.charge);
    }

    virtual analysis::Candidate SelectBackgroundMuon(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                                     root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::muonID::veto;
        const std::string selection_label = "muon_bkg";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Muon& object = event->muons().at(id);

        cut(true, ">0 mu cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(isTightMuon, 2, -0.5, 1.5) == isTightMuon, "tight");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");
        cut(X(pfRelIso, 1000, 0.0, 100.0) < pFRelIso, "pFRelIso");

        return analysis::Candidate(analysis::Candidate::Mu, id, object,object.charge);
    }

    virtual analysis::Candidate SelectBackgroundTau(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                                    root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::tauID::veto;
        const std::string selection_label = "tau_bkg";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Tau& object = event->taus().at(id);

        cut(true, ">0 tau cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(decayModeFinding, 2, -0.5, 1.5) > decayModeFinding, "decay_mode");
        cut(X(byLooseCombinedIsolationDeltaBetaCorr3Hits, 2, -0.5, 1.5) >
            LooseCombinedIsolationDeltaBetaCorr3Hits, "looseIso3Hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");

        return analysis::Candidate(analysis::Candidate::Tau, id, object,object.charge);
    }

    CandidateVector ApplyVetos(const CandidateVector& Resonances, cuts::Cutter& cut)
    {
        const auto electrons_bkg = CollectElectrons(false);
        const auto resonances_noEle = FilterBackground(Resonances,electrons_bkg,
                              cuts::Htautau_Summer13::electronID::veto::deltaR_signalObjects,"resonances_noEle");
        cut(resonances_noEle.size(), "no_electrons");


        const auto muons_bkg = CollectMuons(false);
        const auto resonances_noMu = FilterBackground(resonances_noEle,muons_bkg,
                                 cuts::Htautau_Summer13::muonID::veto::deltaR_signalObjects,"resonances_noMu");
        cut(resonances_noMu.size(), "no_muons");


        const auto bjets_bkg = CollectBJets(cuts::Htautau_Summer13::btag::CSVL, "loose",false);
        const auto resonances_noBjets = FilterBackground(resonances_noMu,bjets_bkg,
                              cuts::Htautau_Summer13::btag::veto::deltaR_signalObjects,"resonances_noBjets");
        cut(resonances_noBjets.size(), "no_bjets");


        const auto taus_bkg = CollectTaus(false);
        const auto resonances_noTau = FilterBackground(resonances_noBjets,taus_bkg,
                                 cuts::Htautau_Summer13::tauID::veto::deltaR_signalObjects,"resonances_noTau");
        cut(resonances_noTau.size(), "no_taus");

        return resonances_noTau;
    }

    CandidateVector FilterBackground(const CandidateVector& candidates, const CandidateVector& backgroundCandidates,
                                            double minDeltaR, const std::string& hist_name)
    {
        CandidateVector result;
        for (const Candidate& candidate : candidates){
            if (FilterBackground(candidate, backgroundCandidates,minDeltaR)){
                result.push_back(candidate);
                GetAnaData().Mass(hist_name).Fill(candidate.momentum.M(),weight);
            }
        }
        GetAnaData().N_objects(hist_name).Fill(result.size(),weight);
        return result;
    }

    static bool FilterBackground(const Candidate& candidate, const CandidateVector& backgroundCandidates,
                                            double minDeltaR)
    {
        if(candidate.finalStateDaughters.size()) {
            for (const Candidate& bkg_candidate : backgroundCandidates) {
                bool hasMatchedDaughter = false;
                for(const Candidate& daughter : candidate.finalStateDaughters) {
                    if (bkg_candidate.momentum.DeltaR(daughter.momentum) <= minDeltaR) {
                        hasMatchedDaughter = true;
                        break;
                    }
                }
                if(!hasMatchedDaughter)
                    return false;
            }
            return true;
        }
        for (const Candidate& bkg_candidate : backgroundCandidates){
            if (bkg_candidate.momentum.DeltaR(candidate.momentum) > minDeltaR)
                return false;
        }
        return true;
    }

    void HistogramAfterFinalSelection(const CandidateVector& finalSelectionResonances)
    {
        for (const Candidate& resonance : finalSelectionResonances){
            for (const Candidate& finalStateDaughter : resonance.finalStateDaughters){
                if (finalStateDaughter.type == Candidate::Mu)
                    SelectMuon(finalStateDaughter.index,GetAnaData().MuonSelection(), false, anaDataFinalSelection);
                if (finalStateDaughter.type == Candidate::Tau)
                    SelectTau(finalStateDaughter.index,GetAnaData().TauSelection(), false, anaDataFinalSelection);
                if (finalStateDaughter.type == Candidate::Bjet)
                    SelectBJet(finalStateDaughter.index,GetAnaData().BJetSelection("loose"), false, anaDataFinalSelection,
                               cuts::Htautau_Summer13::btag::CSVL, "loose");
                if (finalStateDaughter.type == Candidate::Electron)
                    SelectElectron(finalStateDaughter.index,GetAnaData().ElectronSelection(), false, anaDataFinalSelection);
            }
//            GetAnaData().Mass("Final_H_tautau").Fill(resonance.daughters.at(0)->momentum.M(),weight);
//            GetAnaData().Mass("Final_H_bb").Fill(resonance.daughters.at(1)->momentum.M(),weight);
//            const Candidate corrected_h_tautua = CorrectMassBySVfit(*resonance.daughters.at(0), event->metMVA());
//            GetAnaData().Mass("Final_H_tautau_corr").Fill(corrected_h_tautua.momentum.M(), weight);
            GetAnaData().Mass("Final_H_tautau").Fill(resonance.momentum.M(), weight);
            genEvent.Initialize(event->genParticles());
            const auto higgsesMC = genEvent.GetParticles( { particles::Higgs } );
            //std::cerr << "N higgses" << higgsesMC.size() << std::endl;
            if(higgsesMC.size() != 1) {

                throw std::runtime_error("More then 1 MC higgs per event!");
            }
            const ntuple::MET correctedMET = ApplyPostRecoilCorrection(event->metMVA(), resonance.momentum,
                                                                       (*higgsesMC.begin())->momentum);
            //const Candidate corrected_h_tautua = CorrectMassBySVfit(resonance, event->metMVA());
            const Candidate corrected_h_tautua = CorrectMassBySVfit(resonance, correctedMET);
            GetAnaData().Mass("Final_H_tautau_corr").Fill(corrected_h_tautua.momentum.M(), weight);
        }
    }

    void FindAnalysisFinalState(finalState::bbTauTau& final_state)
    {
        static const particles::ParticleCodes resonanceCodes = { particles::Radion };
        static const particles::ParticleCodes resonanceDecay = { particles::Higgs, particles::Higgs };
        static const particles::ParticleCodes2D HiggsDecays = { { particles::b, particles::b },
                                                                { particles::tau, particles::tau } };

        genEvent.Initialize(event->genParticles());

        const GenParticleSet resonances = genEvent.GetParticles(resonanceCodes);
        if (resonances.size() != 1)
            throw std::runtime_error("not one resonance per event");

        final_state.resonance = *resonances.begin();

        GenParticlePtrVector HiggsBosons;
        if(!FindDecayProducts(*final_state.resonance, resonanceDecay,HiggsBosons))
            throw std::runtime_error("Resonance does not decay into 2 Higgs");

        GenParticleVector2D HiggsDecayProducts;
        GenParticleIndexVector HiggsIndexes;
        if(!FindDecayProducts2D(HiggsBosons,HiggsDecays,HiggsDecayProducts,HiggsIndexes))
            throw std::runtime_error("NOT HH -> bb tautau");

        final_state.b_jets = HiggsDecayProducts.at(0);
        final_state.taus = HiggsDecayProducts.at(1);

        final_state.Higgs_TauTau = HiggsBosons.at(HiggsIndexes.at(1));
        final_state.Higgs_BB = HiggsBosons.at(HiggsIndexes.at(0));

        GetAnaData().Resonance_Mass().Fill(final_state.resonance->momentum.M());
        GetAnaData().Resonance_Pt().Fill(final_state.resonance->momentum.Pt());
        GetAnaData().Resonance_Eta().Fill(final_state.resonance->momentum.Eta());
        GetAnaData().Resonance_Phi().Fill(final_state.resonance->momentum.Phi());

        GetAnaData().Bjets_Pt_MC().Fill(final_state.b_jets.at(0)->momentum.Pt());
        GetAnaData().Bjets_Pt_MC().Fill(final_state.b_jets.at(1)->momentum.Pt());
        GetAnaData().Higgs_leptonic_MC_Pt().Fill(final_state.Higgs_TauTau->momentum.Pt());
        GetAnaData().Higgs_BB_MC_Pt().Fill(final_state.Higgs_BB->momentum.Pt());
        GetAnaData().DR_bjets_vs_HiggsPt_MC().Fill(final_state.Higgs_BB->momentum.Pt(),
                                         final_state.b_jets.at(0)->momentum.DeltaR(final_state.b_jets.at(1)->momentum));
        GetAnaData().DR_Higgs_vs_ResonancePt_MC().Fill(final_state.resonance->momentum.Pt(),
                                         final_state.Higgs_TauTau->momentum.DeltaR(final_state.Higgs_BB->momentum));
    }
};
} // analysis
