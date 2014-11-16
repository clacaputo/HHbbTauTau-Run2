/*!
 * \file BaseFlatTreeProducer.h
 * \brief Definition of BaseFlatTreeProducer class, the base class for flat tree producers.
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

#pragma once

#include "AnalysisBase/include/FlatTree.h"

#include "BaseAnalyzer.h"

namespace analysis {

class BaseFlatTreeProducer : public BaseAnalyzer {
public:
    BaseFlatTreeProducer(const std::string& inputFileName, const std::string& outputFileName,
                         const std::string& configFileName, const std::string& _prefix = "none",
                         size_t _maxNumberOfEvents = 0,
                         std::shared_ptr<ntuple::FlatTree> _flatTree = std::shared_ptr<ntuple::FlatTree>())
        : BaseAnalyzer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          flatTree(_flatTree), writeFlatTree(!flatTree)
    {
        if(!flatTree)
            flatTree = std::shared_ptr<ntuple::FlatTree>(new ntuple::FlatTree("flatTree"));
    }

    virtual ~BaseFlatTreeProducer() override
    {
        if(writeFlatTree)
            flatTree->Write();
    }

    virtual void ProcessEvent(std::shared_ptr<const EventDescriptor> _event) override
    {
        using namespace analysis;
        using namespace cuts::Htautau_Summer13;

        BaseAnalyzer::ProcessEvent(_event);

        SelectionResults& selection = ApplyBaselineSelection();

        selection.svfitResults = sv_fit::FitWithUncertainties(*selection.higgs, selection.MET_with_recoil_corrections,
                                                              tauCorrections::energyUncertainty, true, true);

        selection.kinfitResults = RunKinematicFit(selection.bjets_all, *selection.higgs,
                                                  selection.MET_with_recoil_corrections, true, true);
        FillFlatTree(selection);
    }

protected:
    virtual void FillFlatTree(const SelectionResults& selection)
    {
        static const float default_value = ntuple::DefaultFloatFillValueForFlatTree();

        // Event
        flatTree->run() = event->eventInfo().run;
        flatTree->lumi() = event->eventInfo().lumis;
        flatTree->evt() = event->eventInfo().EventId;
        flatTree->eventType() = static_cast<int>(selection.eventType);

        flatTree->npv() = selection.vertices.size();
        if (config.ApplyPUreweight()){
            const size_t bxIndex = tools::find_index(event->eventInfo().bunchCrossing, 0);
            if(bxIndex >= event->eventInfo().bunchCrossing.size())
                throw std::runtime_error("in-time BX not found");
            flatTree->npu() = event->eventInfo().trueNInt.at(bxIndex);
        }

        // Weights
        flatTree->puweight()        = PUweight;
        flatTree->trigweight_1()    = triggerWeights.at(0);
        flatTree->trigweight_2()    = triggerWeights.at(1);
        flatTree->idweight_1()      = IDweights.at(0);
        flatTree->idweight_2()      = IDweights.at(1);
        flatTree->isoweight_1()     = IsoWeights.at(0);
        flatTree->isoweight_2()     = IsoWeights.at(1);
        flatTree->fakeweight_1()     = fakeWeights.at(0); // e -> tau fake rate
        flatTree->fakeweight_2()     = fakeWeights.at(1); // jet -> tau fake rate - default
        flatTree->weight()          = eventWeight;
        flatTree->embeddedWeight()  = config.isDYEmbeddedSample() ? event->genEvent().embeddedWeight : 1.;

        // HTT candidate
        flatTree->mvis() = selection.higgs->GetMomentum().M();
        flatTree->m_sv_vegas() = selection.svfitResults.fit_vegas.has_valid_mass
                ? selection.svfitResults.fit_vegas.mass : default_value;
        flatTree->m_sv_up_vegas() = selection.svfitResults.fit_vegas_up.has_valid_mass
                ? selection.svfitResults.fit_vegas_up.mass : default_value;
        flatTree->m_sv_down_vegas() = selection.svfitResults.fit_vegas_down.has_valid_mass
                ? selection.svfitResults.fit_vegas_down.mass : default_value;
        flatTree->m_sv_MC() = selection.svfitResults.fit_mc.has_valid_mass
                ? selection.svfitResults.fit_mc.mass : default_value;
        flatTree->pt_sv_MC() = selection.svfitResults.fit_mc.has_valid_momentum
                ? selection.svfitResults.fit_mc.momentum.Pt() : default_value;
        flatTree->m_sv_up_MC() = selection.svfitResults.fit_mc_up.has_valid_mass
                ? selection.svfitResults.fit_mc_up.mass : default_value;
        flatTree->pt_sv_up_MC() = selection.svfitResults.fit_mc_up.has_valid_momentum
                ? selection.svfitResults.fit_mc_up.momentum.Pt() : default_value;
        flatTree->m_sv_down_MC() = selection.svfitResults.fit_mc_down.has_valid_mass
                ? selection.svfitResults.fit_mc_down.mass : default_value;
        flatTree->pt_sv_down_MC() = selection.svfitResults.fit_mc_down.has_valid_momentum
                ? selection.svfitResults.fit_mc_down.momentum.Pt() : default_value;
        flatTree->eta_sv_MC() = selection.svfitResults.fit_mc.has_valid_momentum
                ? selection.svfitResults.fit_mc.momentum.Eta() : default_value;
        flatTree->phi_sv_MC() = selection.svfitResults.fit_mc.has_valid_momentum
                ? selection.svfitResults.fit_mc.momentum.Phi() : default_value;
        flatTree->eta_sv_up_MC() = selection.svfitResults.fit_mc_up.has_valid_momentum
                ? selection.svfitResults.fit_mc_up.momentum.Eta() : default_value;
        flatTree->phi_sv_up_MC() = selection.svfitResults.fit_mc_up.has_valid_momentum
                ? selection.svfitResults.fit_mc_up.momentum.Phi() : default_value;
        flatTree->eta_sv_down_MC() = selection.svfitResults.fit_mc_down.has_valid_momentum
                ? selection.svfitResults.fit_mc_down.momentum.Eta() : default_value;
        flatTree->phi_sv_down_MC() = selection.svfitResults.fit_mc_down.has_valid_momentum
                ? selection.svfitResults.fit_mc_down.momentum.Phi() : default_value;
        flatTree->DeltaR_leptons() = selection.GetLeg1()->GetMomentum().DeltaR(selection.GetLeg2()->GetMomentum()) ;
        flatTree->pt_tt()          = (selection.GetLeg1()->GetMomentum() + selection.GetLeg2()->GetMomentum()).Pt();

        // Kinematic fit
        for(const auto& fit_result_entry : selection.kinfitResults) {
            const kinematic_fit::four_body::FitResults& result_4body = fit_result_entry.second.fit_bb_tt;
            flatTree->kinfit_bb_tt_mass().push_back(result_4body.mass);
            flatTree->kinfit_bb_tt_convergence().push_back(result_4body.convergence);
            flatTree->kinfit_bb_tt_chi2().push_back(result_4body.chi2);
            flatTree->kinfit_bb_tt_pull_balance().push_back(result_4body.pull_balance);
            const kinematic_fit::four_body::FitResults& result_4body_tt_up = fit_result_entry.second.fit_bb_tt_up;
            flatTree->kinfit_bb_tt_up_mass().push_back(result_4body_tt_up.mass);
            flatTree->kinfit_bb_tt_up_convergence().push_back(result_4body_tt_up.convergence);
            flatTree->kinfit_bb_tt_up_chi2().push_back(result_4body_tt_up.chi2);
            flatTree->kinfit_bb_tt_up_pull_balance().push_back(result_4body_tt_up.pull_balance);
            const kinematic_fit::four_body::FitResults& result_4body_tt_down = fit_result_entry.second.fit_bb_tt_down;
            flatTree->kinfit_bb_tt_down_mass().push_back(result_4body_tt_down.mass);
            flatTree->kinfit_bb_tt_down_convergence().push_back(result_4body_tt_down.convergence);
            flatTree->kinfit_bb_tt_down_chi2().push_back(result_4body_tt_down.chi2);
            flatTree->kinfit_bb_tt_down_pull_balance().push_back(result_4body_tt_down.pull_balance);

//            const kinematic_fit::two_body::FitResults& result_2body = fit_result_entry.second.fit_bb;
//            TLorentzVector ditau_momentum;
//            ditau_momentum.SetPtEtaPhiM(svfitResults.fit_mc.pt, higgs.momentum.Eta(), higgs.momentum.Phi(),
//                                        svfitResults.fit_mc.mass);
//            const TLorentzVector combined_momentum = ditau_momentum + result_2body.bjet_momentums.at(0)
//                                                                    + result_2body.bjet_momentums.at(1);
//            flatTree->kinfit_bb_sv_mass().push_back(combined_momentum.M());
//            flatTree->kinfit_bb_convergence().push_back(result_2body.convergence);
//            flatTree->kinfit_bb_chi2().push_back(result_2body.chi2);
        }

        // Hhh generator info candidate
        if(selection.GetFinalStateMC().resonance) {
            const TLorentzVector& momentum = selection.GetFinalStateMC().resonance->momentum;
            flatTree->pt_resonance_MC()    = momentum.Pt()  ;
            flatTree->eta_resonance_MC()   = momentum.Eta() ;
            flatTree->phi_resonance_MC()   = momentum.Phi() ;
            flatTree->mass_resonance_MC()  = momentum.M()   ;
            flatTree->pdgId_resonance_MC() = selection.GetFinalStateMC().resonance->pdg.ToInteger();
        } else {
            flatTree->pt_resonance_MC()    = default_value  ;
            flatTree->eta_resonance_MC()   = default_value  ;
            flatTree->phi_resonance_MC()   = default_value  ;
            flatTree->mass_resonance_MC()  = default_value  ;
            flatTree->pdgId_resonance_MC() = particles::NONEXISTENT.RawCode() ;
        }

        if(selection.GetFinalStateMC().Higgs_TauTau) {
            const TLorentzVector& momentum = selection.GetFinalStateMC().Higgs_TauTau->momentum;
            flatTree->pt_Htt_MC()    = momentum.Pt()  ;
            flatTree->eta_Htt_MC()   = momentum.Eta() ;
            flatTree->phi_Htt_MC()   = momentum.Phi() ;
            flatTree->mass_Htt_MC()  = momentum.M()   ;
            flatTree->pdgId_Htt_MC() = selection.GetFinalStateMC().Higgs_TauTau->pdg.ToInteger();
        } else {
            flatTree->pt_Htt_MC()    = default_value  ;
            flatTree->eta_Htt_MC()   = default_value  ;
            flatTree->phi_Htt_MC()   = default_value  ;
            flatTree->mass_Htt_MC()  = default_value  ;
            flatTree->pdgId_Htt_MC() = particles::NONEXISTENT.RawCode() ;
        }

        if(selection.GetFinalStateMC().Higgs_BB) {
            const TLorentzVector& momentum = selection.GetFinalStateMC().Higgs_BB->momentum;
            flatTree->pt_Hbb_MC()    = momentum.Pt()  ;
            flatTree->eta_Hbb_MC()   = momentum.Eta() ;
            flatTree->phi_Hbb_MC()   = momentum.Phi() ;
            flatTree->mass_Hbb_MC()  = momentum.M()   ;
            flatTree->pdgId_Hbb_MC() = selection.GetFinalStateMC().Higgs_BB->pdg.ToInteger();
        } else {
            flatTree->pt_Hbb_MC()    = default_value  ;
            flatTree->eta_Hbb_MC()   = default_value  ;
            flatTree->phi_Hbb_MC()   = default_value  ;
            flatTree->mass_Hbb_MC()  = default_value  ;
            flatTree->pdgId_Hbb_MC() = particles::NONEXISTENT.RawCode()  ;
        }

        // needs to be filles with NUP!
        // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauTauAnalyzer.py#L51
        if (config.MaxTreeVersion() == 2)
            flatTree->n_extraJets_MC() = event->genEvent().nup;
        else
            flatTree->n_extraJets_MC() = default_value;

        // MET
        const TLorentzVector MET_momentum = MakeLorentzVectorPtEtaPhiM(selection.MET_with_recoil_corrections.pt, 0,
                                                                       selection.MET_with_recoil_corrections.phi, 0);
        flatTree->pt_tt_MET() = (selection.GetLeg1()->GetMomentum() + selection.GetLeg2()->GetMomentum()
                                 + MET_momentum).Pt();

        flatTree->met() = selection.pfMET.pt;
        flatTree->metphi() = selection.pfMET.phi;
        flatTree->mvamet() = MET_momentum.Pt();
        flatTree->mvametphi() = MET_momentum.Phi();
        //flatTree->pzetavis();
        //flatTree->pzetamiss();
        if(selection.pfMET.significanceMatrix.size()) {
            const TMatrixD metPFcov = ntuple::VectorToSignificanceMatrix(selection.pfMET.significanceMatrix);
            flatTree->metcov00() = metPFcov[0][0];
            flatTree->metcov01() = metPFcov[0][1];
            flatTree->metcov10() = metPFcov[1][0];
            flatTree->metcov11() = metPFcov[1][1];
        }
        const TMatrixD metMVAcov =
                ntuple::VectorToSignificanceMatrix(selection.MET_with_recoil_corrections.significanceMatrix);
        flatTree->mvacov00() = metMVAcov[0][0];
        flatTree->mvacov01() = metMVAcov[0][1];
        flatTree->mvacov10() = metMVAcov[1][0];
        flatTree->mvacov11() = metMVAcov[1][1];

        // Leg 1, lepton
        flatTree->pt_1()     = selection.GetLeg1()->GetMomentum().Pt();
        flatTree->phi_1()    = selection.GetLeg1()->GetMomentum().Phi();
        flatTree->eta_1()    = selection.GetLeg1()->GetMomentum().Eta();
        flatTree->m_1()      = selection.GetLeg1()->GetMomentum().M();
        flatTree->energy_1() = selection.GetLeg1()->GetMomentum().E();
        flatTree->q_1()      = selection.GetLeg1()->GetCharge();
        flatTree->mt_1()     = Calculate_MT(selection.GetLeg1()->GetMomentum(), MET_momentum.Pt(), MET_momentum.Phi());
        flatTree->d0_1()     = Calculate_dxy(selection.GetLeg1()->GetVertexPosition(), primaryVertex->GetPosition(),
                                             selection.GetLeg1()->GetMomentum());
        flatTree->dZ_1()     = selection.GetLeg1()->GetVertexPosition().Z() - primaryVertex->GetPosition().Z();

        // Leg 2, tau
        flatTree->pt_2()     = selection.GetLeg2()->GetMomentum().Pt();
        flatTree->phi_2()    = selection.GetLeg2()->GetMomentum().Phi();
        flatTree->eta_2()    = selection.GetLeg2()->GetMomentum().Eta();
        flatTree->m_2()      = selection.GetLeg2()->GetMomentum().M();
        flatTree->energy_2() = selection.GetLeg2()->GetMomentum().E();
        flatTree->q_2()      = selection.GetLeg2()->GetCharge();
        flatTree->mt_2()     = Calculate_MT(selection.GetLeg2()->GetMomentum(), MET_momentum.Pt(), MET_momentum.Phi());
        flatTree->d0_2()     = Calculate_dxy(selection.GetLeg2()->GetVertexPosition(), primaryVertex->GetPosition(),
                                             selection.GetLeg2()->GetMomentum());
        flatTree->dZ_2()     = selection.GetLeg2()->GetVertexPosition().Z() - primaryVertex->GetPosition().Z();

        // RM: for the three channels, mt, et, tt this leg is always a tau
        const ntuple::Tau& ntuple_tau_leg2 = selection.GetLeg2()->GetNtupleObject<ntuple::Tau>();
        flatTree->decayMode_2()                                = ntuple_tau_leg2.decayMode;
        flatTree->againstElectronLooseMVA_2() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau_leg2, 0);
        flatTree->againstElectronMediumMVA_2() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau_leg2, 1);
        flatTree->againstElectronTightMVA_2() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau_leg2, 2);
        flatTree->againstElectronVTightMVA_2() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau_leg2, 3);
        flatTree->againstElectronLoose_2()                     = ntuple_tau_leg2.againstElectronLoose  ;
        flatTree->againstElectronMedium_2()                    = ntuple_tau_leg2.againstElectronMedium ;
        flatTree->againstElectronTight_2()                     = ntuple_tau_leg2.againstElectronTight  ;
        flatTree->againstMuonLoose_2()                         = ntuple_tau_leg2.againstMuonLoose      ;
        flatTree->againstMuonMedium_2()                        = ntuple_tau_leg2.againstMuonMedium     ;
        flatTree->againstMuonTight_2()                         = ntuple_tau_leg2.againstMuonTight      ;
        flatTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = ntuple_tau_leg2.byCombinedIsolationDeltaBetaCorrRaw3Hits ;

        // Jets
        flatTree->njets()     = selection.jets.size();
        flatTree->njetspt20() = selection.jetsPt20.size();
        flatTree->nBjets()    = selection.bjets_all.size();
        flatTree->nBjets_retagged()    = selection.retagged_bjets.size();

        for (const CandidatePtr& jet : selection.bjets_all) {
            const ntuple::Jet& ntuple_jet = jet->GetNtupleObject<ntuple::Jet>();

            flatTree->pt_Bjets()      .push_back( jet->GetMomentum().Pt() );
            flatTree->eta_Bjets()     .push_back( jet->GetMomentum().Eta() );
            flatTree->phi_Bjets()     .push_back( jet->GetMomentum().Phi() );
            flatTree->energy_Bjets()  .push_back( jet->GetMomentum().E() );
            flatTree->chargedHadronEF_Bjets().push_back( ntuple_jet.chargedHadronEnergyFraction );
            flatTree->neutralHadronEF_Bjets()  .push_back( ntuple_jet.neutralHadronEnergyFraction );
            flatTree->photonEF_Bjets()         .push_back( ntuple_jet.photonEnergyFraction );
            flatTree->muonEF_Bjets()  .push_back( ntuple_jet.muonEnergyFraction );
            flatTree->electronEF_Bjets()  .push_back( ntuple_jet.electronEnergyFraction );
            flatTree->csv_Bjets()     .push_back( ntuple_jet.combinedSecondaryVertexBJetTags );
            // inspect the flavour of the gen jet
            const VisibleGenObjectVector matched_bjets_MC = FindMatchedObjects(jet->GetMomentum(),
                                                                               selection.GetFinalStateMC().b_jets,
                                                                               cuts::DeltaR_MC_Match);
            const bool isJet_MC_Bjet = matched_bjets_MC.size() != 0;
            const bool isJet_MC_Bjet_withLeptonicDecay = isJet_MC_Bjet
                    && matched_bjets_MC.at(0).finalStateChargedLeptons.size() != 0;
            flatTree->isBjet_MC_Bjet()                  .push_back( isJet_MC_Bjet );
            flatTree->isBjet_MC_Bjet_withLeptonicDecay().push_back( isJet_MC_Bjet_withLeptonicDecay );
        }

        flatTree->x_PV() = primaryVertex->GetPosition().x();
        flatTree->y_PV() = primaryVertex->GetPosition().y();
        flatTree->z_PV() = primaryVertex->GetPosition().z();
    }

protected:
    std::shared_ptr<ntuple::FlatTree> flatTree;
    bool writeFlatTree;
};
} // analysis
