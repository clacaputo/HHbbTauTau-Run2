/*!
 * \file Muon.h
 * \brief Definiton of ntuple::MuonTree and ntuple::Muon classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-25 created
 */

#pragma once

#include "SmartTree.h"

#define MUON_DATA() \
    /* Muon ID */ \
    SIMPLE_VAR(Bool_t, isTrackerMuon, false) \
    SIMPLE_VAR(Bool_t, isPFMuon, false) \
    SIMPLE_VAR(Bool_t, isGlobalMuonPromptTight, false) \
    SIMPLE_VAR(Bool_t, isTightMuon, false) \
    SIMPLE_VAR(Bool_t, isHighPtMuon, false) \
    /* 4-momentum */ \
    SIMPLE_VAR(Float_t, eta, 0.0) \
    SIMPLE_VAR(Float_t, phi, 0.0) \
    SIMPLE_VAR(Float_t, pt, 0.0) \
    SIMPLE_VAR(Float_t, ptError, 0.0) \
    SIMPLE_VAR(Float_t, energy, 0.0) \
    SIMPLE_VAR(Int_t, charge, 0) \
    /* Track info */ \
    SIMPLE_VAR(Float_t, trkD0, 0.0) \
    SIMPLE_VAR(Float_t, trkD0Error, 0.0) \
    SIMPLE_VAR(Float_t, trkDz, 0.0) \
    SIMPLE_VAR(Float_t, trkDzError, 0.0) \
    SIMPLE_VAR(Float_t, globalChi2, 0.0) \
    SIMPLE_VAR(Bool_t, passID, 0) \
    SIMPLE_VAR(Int_t, pixHits, 0.0) \
    SIMPLE_VAR(Int_t, trkHits, 0.0) \
    SIMPLE_VAR(Int_t, muonHits, 0.0) \
    SIMPLE_VAR(Int_t, matches, 0.0) \
    SIMPLE_VAR(Int_t, trackerLayersWithMeasurement, 0.0) \
    /* Isolation */ \
    SIMPLE_VAR(Float_t, trkIso, 0.0) \
    SIMPLE_VAR(Float_t, ecalIso, 0.0) \
    SIMPLE_VAR(Float_t, hcalIso, 0.0) \
    SIMPLE_VAR(Float_t, hoIso, 0.0) \
    SIMPLE_VAR(Float_t, relIso, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIso, 0.0) \
    /* Vertex */ \
    SIMPLE_VAR(Float_t, vtxDist3D, 0.0) \
    SIMPLE_VAR(Int_t, vtxIndex, 0) \
    SIMPLE_VAR(Float_t, vtxDistZ, 0.0) \
    SIMPLE_VAR(Float_t, vx, 0.0) \
    SIMPLE_VAR(Float_t, vy, 0.0) \
    SIMPLE_VAR(Float_t, vz, 0.0)\
    /* PV2D*/ \
    SIMPLE_VAR(Float_t, dB, 0.0) \
    SIMPLE_VAR(Float_t, edB, 0.0) \
    /* PV3D */ \
    SIMPLE_VAR(Float_t, dB3d, 0.0) \
    SIMPLE_VAR(Float_t, edB3d, 0.0) \
    /* UW Recommendation*/ \
    SIMPLE_VAR(Bool_t, isAllArbitrated, false) \
    SIMPLE_VAR(Int_t, nChambers, 0) \
    SIMPLE_VAR(Int_t, nMatches, 0) \
    SIMPLE_VAR(Int_t, nMatchedStations, 0) \
    SIMPLE_VAR(UInt_t, stationMask, 0) \
    SIMPLE_VAR(UInt_t, stationGapMaskDistance, 0)\
    SIMPLE_VAR(UInt_t, stationGapMaskPull, 0) \
    SIMPLE_VAR(Int_t, muonID, 0) \
    /* MVA */ \
    SIMPLE_VAR(Float_t, idMVA, 0.0) \
    SIMPLE_VAR(Float_t, isoRingsMVA, 0.0) \
    SIMPLE_VAR(Float_t, isoRingsRadMVA, 0.0) \
    SIMPLE_VAR(Float_t, idIsoCombMVA, 0.0) \
    /* Iso variables*/ \
    SIMPLE_VAR(Float_t, pfRelIso03v1, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIso03v2, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIsoDB03v1, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIsoDB03v2, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIso04v1, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIso04v2, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIsoDB04v1, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIsoDB04v2, 0.0) \
    /**/

#define SIMPLE_VAR(type, name, default_value) SIMPLE_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_TREE_BRANCH(type, name)

namespace ntuple {
class MuonTree : public root_ext::SmartTree {
public:
    static const std::string& Name() { static const std::string name = "muons"; return name; }
    MuonTree() : SmartTree(Name(), "/", false) {}
    MuonTree(TFile& file)
        : SmartTree(Name(), file) {}

    SIMPLE_TREE_BRANCH(UInt_t, EventId, 0) \
    MUON_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct Muon {
    MUON_DATA()
    Muon() {}
    inline Muon(MuonTree& tree);
};

typedef std::vector<Muon> MuonVector;

} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) name = tree.name();
#define VECTOR_VAR(type, name) name = tree.name();

namespace ntuple {
inline Muon::Muon(MuonTree& tree)
{
    MUON_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef MUON_DATA
