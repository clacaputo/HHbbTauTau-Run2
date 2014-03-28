/*!
 * \file Tau.h
 * \brief Definiton of ntuple::TauTree and ntuple::Tau classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-25 created
 */

#pragma once

#include "SmartTree.h"

#define TAU_DATA() \
    /* 4-momentum */ \
    SIMPLE_VAR(Float_t, eta, 0.0) \
    SIMPLE_VAR(Float_t, phi, 0.0) \
    SIMPLE_VAR(Float_t, pt, 0.0) \
    SIMPLE_VAR(Float_t, energy, 0.0) \
    /* Charge */ \
    SIMPLE_VAR(Int_t, charge, 0.0) \
    /* Leading particle pT */ \
    SIMPLE_VAR(Float_t, leadChargedParticlePt, 0.0) \
    SIMPLE_VAR(Float_t, leadNeutralParticleEt, 0.0) \
    SIMPLE_VAR(Float_t, leadParticleEt, 0.0) \
    /* Number of charged/neutral candidates and photons in different cones */ \
    SIMPLE_VAR(UInt_t, numChargedHadronsSignalCone, 0) \
    SIMPLE_VAR(UInt_t, numNeutralHadronsSignalCone, 0) \
    SIMPLE_VAR(UInt_t, numPhotonsSignalCone, 0) \
    SIMPLE_VAR(UInt_t, numParticlesSignalCone, 0) \
    SIMPLE_VAR(UInt_t, numChargedHadronsIsoCone, 0) \
    SIMPLE_VAR(UInt_t, numNeutralHadronsIsoCone, 0) \
    SIMPLE_VAR(UInt_t, numPhotonsIsoCone, 0) \
    SIMPLE_VAR(UInt_t, numParticlesIsoCone, 0) \
    SIMPLE_VAR(Float_t, ptSumPFChargedHadronsIsoCone, 0.0) \
    SIMPLE_VAR(Float_t, etSumPhotonsIsoCone, 0.0) \
    /* Charged hadron candidates */ \
    VECTOR_VAR(Float_t, ChHadCand_Pt) \
    VECTOR_VAR(Float_t, ChHadCand_Eta) \
    VECTOR_VAR(Float_t, ChHadCand_Phi) \
   /* MVA */ \
    SIMPLE_VAR(Float_t, leadPFCand_mva_e_pi, 0.0) \
    /* PFTau specific variables DataFormats/PatCandidates/interface/TauPFSpecific.h */ \
    /* see DataFormats/TauReco/interface/PFTau.h */ \
    SIMPLE_VAR(Float_t, emFraction, 0.0) \
    SIMPLE_VAR(Float_t, maximumHCALPFClusterEt, 0.0) \
    SIMPLE_VAR(Float_t, ecalStripSumEOverPLead, 0.0) \
    SIMPLE_VAR(Float_t, bremsRecoveryEOverPLead, 0.0) \
    SIMPLE_VAR(Float_t, hcalTotOverPLead, 0.0) \
    SIMPLE_VAR(Float_t, hcalMaxOverPLead, 0.0) \
    SIMPLE_VAR(Float_t, hcal3x3OverPLead, 0.0) \
    /* ET weighted eta-phi statistics starting from PFJet, see DataFormats/JetReco/interface/Jet.h */ \
    SIMPLE_VAR(Float_t, etaetaMoment, 0.0) \
    SIMPLE_VAR(Float_t, phiphiMoment, 0.0) \
    SIMPLE_VAR(Float_t, etaphiMoment, 0.0)  \
    /* Vertex */ \
    SIMPLE_VAR(Float_t, vx, 0.0) \
    SIMPLE_VAR(Float_t, vy, 0.0) \
    SIMPLE_VAR(Float_t, vz, 0.0) \
    /* tau discriminators */ \
    TAU_DISCRIMINATOR_DATA() \
    /**/

#define TAU_DISCRIMINATOR_DATA() \
    /* Discriminator by decay mode */ \
    SIMPLE_VAR(Float_t, decayModeFinding, 0.0) \
    /* Discriminators agains electron */ \
    SIMPLE_VAR(Float_t, againstElectronLoose, 0.0) \
    SIMPLE_VAR(Float_t, againstElectronMedium, 0.0) \
    SIMPLE_VAR(Float_t, againstElectronTight, 0.0) \
    SIMPLE_VAR(Float_t, againstElectronLooseMVA5, 0.0) \
    SIMPLE_VAR(Float_t, againstElectronMediumMVA5, 0.0) \
    SIMPLE_VAR(Float_t, againstElectronTightMVA5, 0.0) \
    SIMPLE_VAR(Float_t, againstElectronVTightMVA5, 0.0) \
    /* Discriminators agains muon */ \
    SIMPLE_VAR(Float_t, againstMuonLoose, 0.0) \
    SIMPLE_VAR(Float_t, againstMuonMedium, 0.0) \
    SIMPLE_VAR(Float_t, againstMuonTight, 0.0) \
    SIMPLE_VAR(Float_t, againstMuonLoose3, 0.0) \
    SIMPLE_VAR(Float_t, againstMuonTight3, 0.0) \
    SIMPLE_VAR(Float_t, againstMuonLooseMVA, 0.0) \
    SIMPLE_VAR(Float_t, againstMuonMediumMVA, 0.0) \
    SIMPLE_VAR(Float_t, againstMuonTightMVA, 0.0) \
    SIMPLE_VAR(Float_t, againstMuonMVAraw, 0.0) \
    /* Discriminators for isolation */ \
    SIMPLE_VAR(Float_t, byVLooseCombinedIsolationDeltaBetaCorr, 0.0) \
    SIMPLE_VAR(Float_t, byLooseCombinedIsolationDeltaBetaCorr, 0.0) \
    SIMPLE_VAR(Float_t, byMediumCombinedIsolationDeltaBetaCorr, 0.0) \
    SIMPLE_VAR(Float_t, byTightCombinedIsolationDeltaBetaCorr, 0.0) \
    SIMPLE_VAR(Float_t, byLooseCombinedIsolationDeltaBetaCorr3Hits, 0.0) \
    SIMPLE_VAR(Float_t, byMediumCombinedIsolationDeltaBetaCorr3Hits, 0.0) \
    SIMPLE_VAR(Float_t, byTightCombinedIsolationDeltaBetaCorr3Hits, 0.0) \
    SIMPLE_VAR(Float_t, byIsolationMVA3oldDMwoLTraw, 0.0) \
    SIMPLE_VAR(Float_t, byVLooseIsolationMVA3oldDMwoLT, 0.0) \
    SIMPLE_VAR(Float_t, byLooseIsolationMVA3oldDMwoLT, 0.0) \
    SIMPLE_VAR(Float_t, byMediumIsolationMVA3oldDMwoLT, 0.0) \
    SIMPLE_VAR(Float_t, byTightIsolationMVA3oldDMwoLT, 0.0) \
    SIMPLE_VAR(Float_t, byVTightIsolationMVA3oldDMwoLT, 0.0) \
    SIMPLE_VAR(Float_t, byVVTightIsolationMVA3oldDMwoLT, 0.0) \
    SIMPLE_VAR(Float_t, byIsolationMVA3oldDMwLTraw, 0.0) \
    SIMPLE_VAR(Float_t, byVLooseIsolationMVA3oldDMwLT, 0.0) \
    SIMPLE_VAR(Float_t, byLooseIsolationMVA3oldDMwLT, 0.0) \
    SIMPLE_VAR(Float_t, byMediumIsolationMVA3oldDMwLT, 0.0) \
    SIMPLE_VAR(Float_t, byTightIsolationMVA3oldDMwLT, 0.0) \
    SIMPLE_VAR(Float_t, byVTightIsolationMVA3oldDMwLT, 0.0) \
    SIMPLE_VAR(Float_t, byVVTightIsolationMVA3oldDMwLT, 0.0) \
    SIMPLE_VAR(Float_t, byIsolationMVA3newDMwoLTraw, 0.0) \
    SIMPLE_VAR(Float_t, byVLooseIsolationMVA3newDMwoLT, 0.0) \
    SIMPLE_VAR(Float_t, byLooseIsolationMVA3newDMwoLT, 0.0) \
    SIMPLE_VAR(Float_t, byMediumIsolationMVA3newDMwoLT, 0.0) \
    SIMPLE_VAR(Float_t, byTightIsolationMVA3newDMwoLT, 0.0) \
    SIMPLE_VAR(Float_t, byVTightIsolationMVA3newDMwoLT, 0.0) \
    SIMPLE_VAR(Float_t, byVVTightIsolationMVA3newDMwoLT, 0.0) \
    SIMPLE_VAR(Float_t, byIsolationMVA3newDMwLTraw, 0.0) \
    SIMPLE_VAR(Float_t, byVLooseIsolationMVA3newDMwLT, 0.0) \
    SIMPLE_VAR(Float_t, byLooseIsolationMVA3newDMwLT, 0.0) \
    SIMPLE_VAR(Float_t, byMediumIsolationMVA3newDMwLT, 0.0) \
    SIMPLE_VAR(Float_t, byTightIsolationMVA3newDMwLT, 0.0) \
    SIMPLE_VAR(Float_t, byVTightIsolationMVA3newDMwLT, 0.0) \
    SIMPLE_VAR(Float_t, byVVTightIsolationMVA3newDMwLT, 0.0) \
    /**/

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;
DATA_CLASS(ntuple, Tau, TAU_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) SIMPLE_DATA_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, TauTree, TAU_DATA, Tau, "taus")
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) AddSimpleBranch(#name, data.name);
#define VECTOR_VAR(type, name) AddVectorBranch(#name, data.name);
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, TauTree, TAU_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef TAU_DATA
