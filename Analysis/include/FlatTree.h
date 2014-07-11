/*!
 * \file FlatTree.h
 * \brief Definiton of ntuple::FlatTree and ntuple::Flat classes.
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

#include "HHbbTauTau/TreeProduction/interface/SmartTree.h"

#define FLAT_DATA() \
    /* Event Variables */ \
    SIMPLE_VAR(Int_t, run) /* Run */ \
    SIMPLE_VAR(Int_t, lumi) /* Lumi */ \
    SIMPLE_VAR(Int_t, evt) /* Evt */ \
    /* First signal lepton :  muon for mu Tau, electron for e Tau, Leading (in pT) Tau for Tau Tau */ \
    SIMPLE_VAR(Int_t, type_1) /* pT */ \
    SIMPLE_VAR(Float_t, pt_1) /* pT */ \
    SIMPLE_VAR(Float_t, phi_1) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_1) /* Eta */ \
    SIMPLE_VAR(Float_t, m_1) /* Mass */ \
    SIMPLE_VAR(Float_t, energy_1) /* Mass */ \
    SIMPLE_VAR(Int_t, q_1) /* Charge */ \
    SIMPLE_VAR(Float_t, mt_1) /* mT of  first lepton wrt to MVA met */ \
    SIMPLE_VAR(Float_t, d0_1) /* d0 with respect to primary vertex */ \
    SIMPLE_VAR(Float_t, dZ_1) /* dZ with respect to primary vertex */ \
    SIMPLE_VAR(Bool_t, hasMatchedGenparticle_1) /* is matched with Gen Particle */ \
    /* Gen particle quantities of First signal lepton matched with the truth */ \
    SIMPLE_VAR(Int_t, pdgId_1_matched) /* pT */ \
    SIMPLE_VAR(Float_t, pt_1_matched) /* pT */ \
    SIMPLE_VAR(Float_t, phi_1_matched) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_1_matched) /* Eta */ \
    SIMPLE_VAR(Float_t, energy_1_matched) /* Energy */ \
    /* First lepton discriminators */ \
    SIMPLE_VAR(Int_t, decayMode_1) /* decayMode - tau tau channel*/ \
    SIMPLE_VAR(Float_t, iso_1) /* MVA iso for hadronic Tau, Delta Beta for muon and electron */ \
    SIMPLE_VAR(Float_t, mva_1) /* MVA id (when using electron) 0 otherwise */ \
    SIMPLE_VAR(Bool_t, passid_1) /* Whether it passes id  (not necessarily iso) */ \
    SIMPLE_VAR(Bool_t, passiso_1) /* Whether it passes iso (not necessarily id) */ \
    SIMPLE_VAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits_1) /*  */ \
    SIMPLE_VAR(Float_t, againstElectronMVA_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstElectronLoose_1) /*  */ \
    SIMPLE_VAR(Float_t, againstElectronMedium_1) /*  */ \
    SIMPLE_VAR(Float_t, againstElectronTight_1) /*  */ \
    SIMPLE_VAR(Float_t, againstMuonLoose_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonMedium_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonTight_1) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    /* Second lepton :  hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT)  Tau for Tau Tau */ \
    SIMPLE_VAR(Int_t, type_2) /* pT */ \
    SIMPLE_VAR(Float_t, pt_2) /* pT */ \
    SIMPLE_VAR(Float_t, phi_2) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_2) /* Eta */ \
    SIMPLE_VAR(Float_t, m_2) /* Mass */ \
    SIMPLE_VAR(Float_t, energy_2) /* Mass */ \
    SIMPLE_VAR(Int_t, q_2) /* Charge */ \
    SIMPLE_VAR(Float_t, mt_2) /* mT of  first lepton wrt to MVA met */ \
    SIMPLE_VAR(Float_t, d0_2) /* d0 with respect to primary vertex */ \
    SIMPLE_VAR(Float_t, dZ_2) /* dZ with respect to primary vertex */ \
    SIMPLE_VAR(Bool_t, hasMatchedGenparticle_2) /* is matched with Gen Particle */ \
    /* Gen particle quantities of second signal lepton matched with the truth */ \
    SIMPLE_VAR(Int_t, pdgId_2_matched) /* pT */ \
    SIMPLE_VAR(Float_t, pt_2_matched) /* pT */ \
    SIMPLE_VAR(Float_t, phi_2_matched) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_2_matched) /* Eta */ \
    SIMPLE_VAR(Float_t, energy_2_matched) /* Energy */ \
    /* Second lepton discriminators */ \
    SIMPLE_VAR(Int_t, decayMode_2) /* decayMode - tau tau channel*/ \
    SIMPLE_VAR(Float_t, iso_2) /* MVA iso for hadronic Tau, Delta Beta for muon and electron */ \
    SIMPLE_VAR(Float_t, mva_2) /* MVA id (when using electron) 0 otherwise */ \
    SIMPLE_VAR(Bool_t, passid_2) /* Whether it passes id  (not necessarily iso) */ \
    SIMPLE_VAR(Bool_t, passiso_2) /* Whether it passes iso (not necessarily id) */ \
    SIMPLE_VAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits_2) /*  */ \
    SIMPLE_VAR(Float_t, againstElectronMVA_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstElectronLoose_2) /*  */ \
    SIMPLE_VAR(Float_t, againstElectronMedium_2) /*  */ \
    SIMPLE_VAR(Float_t, againstElectronTight_2) /*  */ \
    SIMPLE_VAR(Float_t, againstMuonLoose_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonMedium_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    SIMPLE_VAR(Float_t, againstMuonTight_2) /* MVA iso for hadronic Tau, Delta Beta for muon */ \
    /* H_tautau variables */ \
    SIMPLE_VAR(Float_t, DeltaR_leptons) /* DeltaR between 2 legs of Higgs candidate */ \
    SIMPLE_VAR(Float_t, mvis) /* SV Fit using integration method */ \
    SIMPLE_VAR(Float_t, m_sv) /* SV Fit using integration method */ \
    SIMPLE_VAR(Float_t, pt_tt) /* pT of 2 legs w/o MVAMET*/ \
    SIMPLE_VAR(Float_t, pt_tt_MET) /* pT of 2 legs with MVAMET*/ \
    SIMPLE_VAR(Float_t, m_sv_Up) /* High Energy scale shape */ \
    SIMPLE_VAR(Float_t, m_sv_Down) /* Low Energy Scale Shape */ \
    /* Met related variables */ \
    SIMPLE_VAR(Float_t, met) /* pfmet */ \
    SIMPLE_VAR(Float_t, metphi) /* pfmet Phi */ \
    SIMPLE_VAR(Float_t, mvamet) /* mvamet */ \
    SIMPLE_VAR(Float_t, mvametphi) /* mvamet Phi */ \
    /* MET covariance matrices */ \
    SIMPLE_VAR(Float_t, metcov00) /* pf met covariance matrix 00 */ \
    SIMPLE_VAR(Float_t, metcov01) /* pf met covariance matrix 01 */ \
    SIMPLE_VAR(Float_t, metcov10) /* pf met covariance matrix 10 */ \
    SIMPLE_VAR(Float_t, metcov11) /* pf met covariance matrix 11 */ \
    /* MVAMet covariance matrices */ \
    SIMPLE_VAR(Float_t, mvacov00) /* mva met covariance matrix 00 */ \
    SIMPLE_VAR(Float_t, mvacov01) /* mva met covariance matrix 01 */ \
    SIMPLE_VAR(Float_t, mvacov10) /* mva met covariance matrix 10 */ \
    SIMPLE_VAR(Float_t, mvacov11) /* mva met covariance matrix 11 */ \
    /* Extra Leptons*/ \
    /* first 4 Muons sorted by pt*/ \
    VECTOR_VAR(Float_t, pt_muons) /* pt of first 4 muons */ \
    VECTOR_VAR(Float_t, eta_muons) /* eta of first 4 muons */ \
    VECTOR_VAR(Float_t, phi_muons) /* phi of first 4 muons */ \
    VECTOR_VAR(Float_t, mass_muons) /* mass of first 4 muons */ \
    VECTOR_VAR(Float_t, iso_muons) /* iso of first 4 muons */ \
    VECTOR_VAR(Int_t, charge_muons) /* q of first 4 muons */ \
    VECTOR_VAR(Bool_t, passId_muons) /* passId of first 4 muons */ \
    /* first 4 Electrons sorted by pt*/ \
    VECTOR_VAR(Float_t, pt_electrons) /* pt of first 4 electrons */ \
    VECTOR_VAR(Float_t, eta_electrons) /* eta of first 4 electrons */ \
    VECTOR_VAR(Float_t, phi_electrons) /* phi of first 4 electrons */ \
    VECTOR_VAR(Float_t, mass_electrons) /* mass of first 4 electrons */ \
    VECTOR_VAR(Float_t, iso_electrons) /* iso of first 4 electrons */ \
    VECTOR_VAR(Int_t, charge_electrons) /* q of first 4 electrons */ \
    VECTOR_VAR(Bool_t, passId_electrons) /* passId of first 4 electrons */ \
    /*useful info at gen level*/ \
    SIMPLE_VAR(Float_t, pt_resonance) /* pt resonance */ \
    SIMPLE_VAR(Float_t, eta_resonance) /* eta resonance */ \
    SIMPLE_VAR(Float_t, phi_resonance) /* phi resonance */ \
    SIMPLE_VAR(Float_t, mass_resonance) /* mass resonance */ \
    SIMPLE_VAR(Float_t, pt_Htt) /* pt Htt */ \
    SIMPLE_VAR(Float_t, eta_Htt) /* eta Htt */ \
    SIMPLE_VAR(Float_t, phi_Htt) /* phi Htt */ \
    SIMPLE_VAR(Float_t, mass_Htt) /* mass Htt */ \
    SIMPLE_VAR(Float_t, pt_Hbb) /* pt Hbb */ \
    SIMPLE_VAR(Float_t, eta_Hbb) /* eta Hbb */ \
    SIMPLE_VAR(Float_t, phi_Hbb) /* phi Hbb */ \
    SIMPLE_VAR(Float_t, mass_Hbb) /* mass Hbb */ \
    SIMPLE_VAR(Int_t, n_extraJets) /* number extra jets*/ \
    /* jets passing jet id ( pt > 20 ) sorted by pt*/ \
    SIMPLE_VAR(Int_t, njets) /*  */ \
    VECTOR_VAR(Float_t, pt_jets) /* pt of first 4 jets */ \
    VECTOR_VAR(Float_t, eta_jets) /* eta of first 4 jets */ \
    VECTOR_VAR(Float_t, phi_jets) /* phi of first 4 jets */ \
    VECTOR_VAR(Float_t, mass_jets) /* mass of first 4 jets */ \
    VECTOR_VAR(Float_t, csv_jets) /* csv of first 4 jets */ \
    VECTOR_VAR(Bool_t, isBjet) /* is a bjet */ \
    VECTOR_VAR(Bool_t, isB_leptonicDecay) /* is a b with leptonic decay */ \
    /* bjets passing jet id ( pt > 20 ) sorted by csv*/ \
    SIMPLE_VAR(Int_t, nBjets) /*  */ \
    VECTOR_VAR(Float_t, pt_Bjets) /* pt of first 4 Bjets */ \
    VECTOR_VAR(Float_t, eta_Bjets) /* eta of first 4 Bjets */ \
    VECTOR_VAR(Float_t, phi_Bjets) /* phi of first 4 Bjets */ \
    VECTOR_VAR(Float_t, mass_Bjets) /* mass of first 4 Bjets */ \
    VECTOR_VAR(Float_t, csv_Bjets) /* csv of first 4 Bjets */ \
    VECTOR_VAR(Bool_t, isBjet) /* is a bjet */ \
    VECTOR_VAR(Bool_t, isB_leptonicDecay) /* is a b with leptonic decay */ \
    /* vertices*/ \
    SIMPLE_VAR(Float_t, x_PV) /* x of PV*/ \
    SIMPLE_VAR(Float_t, y_PV) /* y of PV*/ \
    SIMPLE_VAR(Float_t, z_PV) /* z of PV*/ \
    SIMPLE_VAR(Int_t, npv) /* NPV */ \
    /* Event Weights */ \
    SIMPLE_VAR(Float_t, puweight) /* Pielup Weight */ \
    SIMPLE_VAR(Float_t, trigweight_1) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, trigweight_2) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, idweight_1) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, idweight_2) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, isoweight_1) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, isoweight_2) /* Effieiency Scale factor (all components multiplied in) */ \
    SIMPLE_VAR(Float_t, decayModeWeight) /* decay Mode Weight */ \
    SIMPLE_VAR(Float_t, embeddedWeight) /*  */ \
    SIMPLE_VAR(Float_t, mcweight) /* MC Weight (xs/nevents * additional wieght (ie pt weight for gghiggs)) */ \
    SIMPLE_VAR(Float_t, weight) /* mcweight*puweight*effweight */ \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Flat, FLAT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS(ntuple, FlatTree, FLAT_DATA, Flat, "flat", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_INITIALIZE(ntuple, FlatTree, FLAT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef FLAT_DATA

namespace ntuple {
inline double DefaultFillValueForFlatTree() { return -10000; }
}
