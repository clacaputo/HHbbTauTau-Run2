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
    SIMPLE_VAR(Float_t, eta, 0.0) \
    SIMPLE_VAR(Float_t, phi, 0.0) \
    SIMPLE_VAR(Float_t, pt, 0.0) \
    SIMPLE_VAR(Float_t, energy, 0.0) \
    SIMPLE_VAR(Int_t, charge, 0.0) \
    SIMPLE_VAR(Float_t, mass, 0.0) \
    /* Tracks */ \
    SIMPLE_VAR(Float_t, leadTrkPt, 0.0) \
    SIMPLE_VAR(Float_t, leadTrkPtError, 0.0) \
    SIMPLE_VAR(Float_t, leadTrkEta, 0.0) \
    SIMPLE_VAR(Float_t, leadTrkPhi, 0.0) \
    SIMPLE_VAR(Float_t, leadTrkCharge, 0.0) \
    SIMPLE_VAR(Float_t, leadTrkD0, 0.0) \
    SIMPLE_VAR(Float_t, leadTrkD0Error, 0.0) \
    SIMPLE_VAR(Float_t, leadTrkDz, 0.0) \
    SIMPLE_VAR(Float_t, leadTrkDzError, 0.0) \
    /* Vertex */ \
    SIMPLE_VAR(Int_t, vtxIndex, 0) \
    SIMPLE_VAR(Float_t, vtxDxy, 0.0) \
    SIMPLE_VAR(Float_t, vtxDz, 0.0) \
    /* Leading particle pT */ \
    SIMPLE_VAR(Float_t, leadChargedParticlePt, 0.0) \
    SIMPLE_VAR(Float_t, leadNeutralParticlePt, 0.0) \
    SIMPLE_VAR(Float_t, leadParticlePt, 0.0) \
    /* Number of charged/neutral candidates and photons in different cones */ \
    SIMPLE_VAR(Int_t, numChargedHadronsSignalCone, 0) \
    SIMPLE_VAR(Int_t, numNeutralHadronsSignalCone, 0) \
    SIMPLE_VAR(Int_t, numPhotonsSignalCone, 0) \
    SIMPLE_VAR(Int_t, numParticlesSignalCone, 0) \
    SIMPLE_VAR(Int_t, numChargedHadronsIsoCone, 0) \
    SIMPLE_VAR(Int_t, numNeutralHadronsIsoCone, 0) \
    SIMPLE_VAR(Int_t, numPhotonsIsoCone, 0) \
    SIMPLE_VAR(Int_t, numParticlesIsoCone, 0) \
    SIMPLE_VAR(Float_t, ptSumPFChargedHadronsIsoCone, 0.0) \
    SIMPLE_VAR(Float_t, ptSumPFNeutralHadronsIsoCone, 0.0) \
    SIMPLE_VAR(Float_t, ptSumPhotonsIsoCone, 0.0) \
    /* Signal cone */ \
    VECTOR_VAR(Float_t, sigChHadCandPt) \
    VECTOR_VAR(Float_t, sigChHadCandEta) \
    VECTOR_VAR(Float_t, sigChHadCandPhi) \
    VECTOR_VAR(Float_t, sigNeHadCandPt) \
    VECTOR_VAR(Float_t, sigNeHadCandEta) \
    VECTOR_VAR(Float_t, sigNeHadCandPhi) \
    VECTOR_VAR(Float_t, sigGammaCandPt) \
    VECTOR_VAR(Float_t, sigGammaCandEta) \
    VECTOR_VAR(Float_t, sigGammaCandPhi) \
    /* Isolation cone */ \
    VECTOR_VAR(Float_t, isoChHadCandPt) \
    VECTOR_VAR(Float_t, isoChHadCandEta) \
    VECTOR_VAR(Float_t, isoChHadCandPhi) \
    VECTOR_VAR(Float_t, isoNeHadCandPt) \
    VECTOR_VAR(Float_t, isoNeHadCandEta) \
    VECTOR_VAR(Float_t, isoNeHadCandPhi) \
    VECTOR_VAR(Float_t, isoGammaCandPt) \
    VECTOR_VAR(Float_t, isoGammaCandEta) \
    VECTOR_VAR(Float_t, isoGammaCandPhi) \
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
    /* Charged hadron candidates */ \
    SIMPLE_VAR(Float_t, ChHadCandPt3Prong_1track, 0.0) \
    SIMPLE_VAR(Float_t, ChHadCandPt3Prong_2track, 0.0) \
    SIMPLE_VAR(Float_t, ChHadCandPt3Prong_3track, 0.0) \
    SIMPLE_VAR(Float_t, ChHadCandEta3Prong_1track, 0.0) \
    SIMPLE_VAR(Float_t, ChHadCandEta3Prong_2track, 0.0) \
    SIMPLE_VAR(Float_t, ChHadCandEta3Prong_3track, 0.0) \
    SIMPLE_VAR(Float_t, ChHadCandPhi3Prong_1track, 0.0) \
    SIMPLE_VAR(Float_t, ChHadCandPhi3Prong_2track, 0.0) \
    SIMPLE_VAR(Float_t, ChHadCandPhi3Prong_3track, 0.0) \
    SIMPLE_VAR(Float_t, ChHadCandPt1Prong, 0.0) \
    SIMPLE_VAR(Float_t, ChHadCandEta1Prong, 0.0) \
    SIMPLE_VAR(Float_t, ChHadCandPhi1Prong, 0.0) \
    /* MVA */ \
    SIMPLE_VAR(Float_t, pfElectronMVA, 0.0) \
    /* kinematic variables for PFJet associated to PFTau */ \
    SIMPLE_VAR(Float_t, jetPt, 0.0) \
    SIMPLE_VAR(Float_t, jetEta, 0.0) \
    SIMPLE_VAR(Float_t, jetPhi, 0.0) \
    /* */ \
    SIMPLE_VAR(Float_t, emFraction, 0.0) \
    SIMPLE_VAR(Float_t, maximumHCALPFClusterEt, 0.0) \
    SIMPLE_VAR(Float_t, ecalStripSumEOverPLead, 0.0) \
    SIMPLE_VAR(Float_t, bremsRecoveryEOverPLead, 0.0) \
    SIMPLE_VAR(Float_t, hcalTotOverPLead, 0.0) \
    SIMPLE_VAR(Float_t, hcalMaxOverPLead, 0.0) \
    SIMPLE_VAR(Float_t, hcal3x3OverPLead, 0.0) \
    /* */ \
    SIMPLE_VAR(Float_t, etaetaMoment, 0.0) \
    SIMPLE_VAR(Float_t, phiphiMoment, 0.0) \
    SIMPLE_VAR(Float_t, etaphiMoment, 0.0) \
    /* */ \
    SIMPLE_VAR(Float_t, vx, 0.0) \
    SIMPLE_VAR(Float_t, vy, 0.0) \
    SIMPLE_VAR(Float_t, vz, 0.0) \
    /* */ \
    SIMPLE_VAR(Float_t, zvertex, 0.0) \
    SIMPLE_VAR(Float_t, ltsipt, 0.0) \
    SIMPLE_VAR(Int_t, NumChHad, 0) \
    /**/


#define SIMPLE_VAR(type, name, default_value) SIMPLE_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_TREE_BRANCH(type, name)

namespace ntuple {
class TauTree : public root_ext::SmartTree {
public:
    static const std::string& Name() { static const std::string name = "taus"; return name; }
    TauTree() : SmartTree(Name(), "/", false) {}
    TauTree(const std::string& prefix, TFile& file)
        : SmartTree(prefix + "/" + Name(), file) {}

    SIMPLE_TREE_BRANCH(UInt_t, EventId, 0) \
    TAU_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct Tau {
    TAU_DATA()
    Tau() {}
    inline Tau(TauTree& tree);
};
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) name = tree.name();
#define VECTOR_VAR(type, name) name = tree.name();

namespace ntuple {
inline Tau::Tau(TauTree& tree)
{
    TAU_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef TAU_DATA
