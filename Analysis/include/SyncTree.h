/*!
 * \file SyncTree.h
 * \brief Definiton of ntuple::SyncTree and ntuple::Sync classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-05 created
 */

#pragma once

#include "HHbbTauTau/TreeProduction/interface/SmartTree.h"

#define SYNC_DATA() \
    SIMPLE_VAR(Int_t, NUP) \
    SIMPLE_VAR(Int_t, run) \
    SIMPLE_VAR(Int_t, lumi) \
    SIMPLE_VAR(Int_t, evt) \
    SIMPLE_VAR(Int_t, npv) \
    VECTOR_VAR(Int_t, npu) \
    SIMPLE_VAR(Double_t, rho) \
    SIMPLE_VAR(Double_t, mcweight) \
    SIMPLE_VAR(Double_t, puweight) \
    SIMPLE_VAR(Double_t, effweight) \
    SIMPLE_VAR(Double_t, decaymodeweight) \
    SIMPLE_VAR(Double_t, weight) \
    SIMPLE_VAR(Double_t, m_sv) \
    SIMPLE_VAR(Double_t, mvis) \
    SIMPLE_VAR(Double_t, m_sv_Up) \
    SIMPLE_VAR(Double_t, m_sv_Down) \
    /*1st tau*/ \
    SIMPLE_VAR(Double_t, gen_pt_1) \
    SIMPLE_VAR(Double_t, pt_1) \
    SIMPLE_VAR(Double_t, phi_1) \
    SIMPLE_VAR(Double_t, eta_1) \
    SIMPLE_VAR(Double_t, m_1) \
    SIMPLE_VAR(Double_t, q_1) \
    SIMPLE_VAR(Double_t, iso_1) \
    SIMPLE_VAR(Double_t, antiele_1) \
    SIMPLE_VAR(Double_t, d0_1) \
    SIMPLE_VAR(Double_t, dZ_1) \
    SIMPLE_VAR(Int_t, passid_1) \
    SIMPLE_VAR(Int_t, passiso_1) \
    SIMPLE_VAR(Double_t, mt_1) \
    SIMPLE_VAR(Double_t, trigweight_1) \
    SIMPLE_VAR(Double_t, idweight_1) \
    SIMPLE_VAR(Double_t, isoweight_1) \
    SIMPLE_VAR(Double_t, byCombinedIsolationDeltaBetaCorrRaw3Hits_1) \
    SIMPLE_VAR(Double_t, againstElectronMVA3raw_1) \
    SIMPLE_VAR(Double_t, againstElectronMVA3category_1) \
    SIMPLE_VAR(Double_t, byIsolationMVA2raw_1) \
    SIMPLE_VAR(Double_t, againstMuonLoose_1) \
    SIMPLE_VAR(Double_t, againstMuonLoose2_1) \
    SIMPLE_VAR(Double_t, againstMuonMedium2_1) \
    SIMPLE_VAR(Double_t, againstMuonTight2_1) \
    SIMPLE_VAR(Double_t, againstElectronLooseMVA3_1) \
    SIMPLE_VAR(Double_t, againstElectronLoose_1) \
    SIMPLE_VAR(Double_t, againstElectronNewLooseMVA3_1) \
    /*2nd tau*/ \
    SIMPLE_VAR(Double_t, gen_pt_2) \
    SIMPLE_VAR(Double_t, pt_2) \
    SIMPLE_VAR(Double_t, phi_2) \
    SIMPLE_VAR(Double_t, eta_2) \
    SIMPLE_VAR(Double_t, m_2) \
    SIMPLE_VAR(Double_t, q_2) \
    SIMPLE_VAR(Double_t, iso_2) \
    SIMPLE_VAR(Double_t, antiele_2) \
    SIMPLE_VAR(Int_t, passid_2) \
    SIMPLE_VAR(Int_t, passiso_2) \
    SIMPLE_VAR(Double_t, mt_2) \
    SIMPLE_VAR(Double_t, trigweight_2) \
    SIMPLE_VAR(Double_t, idweight_2) \
    SIMPLE_VAR(Double_t, isoweight_2) \
    SIMPLE_VAR(Double_t, byCombinedIsolationDeltaBetaCorrRaw3Hits_2) \
    SIMPLE_VAR(Double_t, againstElectronMVA3raw_2) \
    SIMPLE_VAR(Double_t, againstElectronMVA3category_2) \
    SIMPLE_VAR(Double_t, byIsolationMVA2raw_2) \
    SIMPLE_VAR(Double_t, againstMuonLoose_2) \
    SIMPLE_VAR(Double_t, againstMuonLoose2_2) \
    SIMPLE_VAR(Double_t, againstMuonMedium2_2) \
    SIMPLE_VAR(Double_t, againstMuonTight2_2) \
    SIMPLE_VAR(Double_t, againstElectronLooseMVA3_2) \
    SIMPLE_VAR(Double_t, againstElectronLoose_2) \
    SIMPLE_VAR(Double_t, againstElectronNewLooseMVA3_2) \
    /*met*/ \
    SIMPLE_VAR(Double_t, met) \
    SIMPLE_VAR(Double_t, metphi) \
    SIMPLE_VAR(Double_t, mvamet) \
    SIMPLE_VAR(Double_t, mvametphi) \
    SIMPLE_VAR(Double_t, pzetavis) \
    SIMPLE_VAR(Double_t, pzetamiss) \
    SIMPLE_VAR(Double_t, metcov00) \
    SIMPLE_VAR(Double_t, metcov01) \
    SIMPLE_VAR(Double_t, metcov10) \
    SIMPLE_VAR(Double_t, metcov11) \
    SIMPLE_VAR(Double_t, mvacov00) \
    SIMPLE_VAR(Double_t, mvacov01) \
    SIMPLE_VAR(Double_t, mvacov10) \
    SIMPLE_VAR(Double_t, mvacov11) \
    /*jet*/ \
    SIMPLE_VAR(Double_t, jpt_1) \
    SIMPLE_VAR(Double_t, jeta_1) \
    SIMPLE_VAR(Double_t, jphi_1) \
    SIMPLE_VAR(Double_t, jptraw_1) \
    SIMPLE_VAR(Double_t, jptunc_1) \
    SIMPLE_VAR(Double_t, jmva_1) \
    SIMPLE_VAR(Int_t, jpass_1) \
    SIMPLE_VAR(Double_t, jpt_2) \
    SIMPLE_VAR(Double_t, jeta_2) \
    SIMPLE_VAR(Double_t, jphi_2) \
    SIMPLE_VAR(Double_t, jptraw_2) \
    SIMPLE_VAR(Double_t, jptunc_2) \
    SIMPLE_VAR(Double_t, jmva_2) \
    SIMPLE_VAR(Int_t, jpass_2) \
    SIMPLE_VAR(Double_t, bpt) \
    SIMPLE_VAR(Double_t, beta) \
    SIMPLE_VAR(Double_t, bphi) \
    SIMPLE_VAR(Double_t, mjj) \
    SIMPLE_VAR(Double_t, jdeta) \
    SIMPLE_VAR(Int_t, njetingap) \
    SIMPLE_VAR(Double_t, mva) \
    SIMPLE_VAR(Double_t, jdphi) \
    SIMPLE_VAR(Double_t, dijetpt) \
    SIMPLE_VAR(Double_t, dijetphi) \
    SIMPLE_VAR(Double_t, hdijetphi) \
    SIMPLE_VAR(Double_t, visjeteta) \
    SIMPLE_VAR(Double_t, ptvis) \
    SIMPLE_VAR(Int_t, nbtag) \
    SIMPLE_VAR(Int_t, njets) \
    SIMPLE_VAR(Int_t, njetspt20) \
    /*trigger*/ \
    SIMPLE_VAR(Int_t, l1TrigMatched_diTau) \
    SIMPLE_VAR(Int_t, l2TrigMatched_diTau) \
    SIMPLE_VAR(Int_t, l1TrigMatched_diTauJet) \
    SIMPLE_VAR(Int_t, l2TrigMatched_diTauJet) \
    SIMPLE_VAR(Int_t, jetTrigMatched_diTauJet) \
    SIMPLE_VAR(Double_t, triggerWeight_diTauJet) \
    SIMPLE_VAR(Double_t, triggerEffMC_diTauJet) \
    SIMPLE_VAR(Double_t, triggerEffData_diTauJet) \
    SIMPLE_VAR(Double_t, triggerEffMC_diTau) \
    SIMPLE_VAR(Double_t, triggerEffData_diTau) \
    SIMPLE_VAR(Double_t, triggerWeight_diTau) \
    /**/ \
    SIMPLE_VAR(Double_t, muon1Pt) \
    SIMPLE_VAR(Double_t, electron1Pt) \
    SIMPLE_VAR(Double_t, genfilter) \
    SIMPLE_VAR(Double_t, tauspin) \
    SIMPLE_VAR(Double_t, zmumusel) \
    SIMPLE_VAR(Double_t, muradcorr) \
    SIMPLE_VAR(Double_t, genTau2PtVsGenTau1Pt) \
    SIMPLE_VAR(Double_t, genTau2EtaVsGenTau1Eta) \
    SIMPLE_VAR(Double_t, genDiTauMassVsGenDiTauPt) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Sync, SYNC_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS(ntuple, SyncTree, SYNC_DATA, Sync, "sync", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_INITIALIZE(ntuple, SyncTree, SYNC_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef SYNC_DATA
