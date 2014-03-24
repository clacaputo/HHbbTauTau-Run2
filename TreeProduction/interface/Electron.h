/*!
 * \file Electoron.h
 * \brief Definiton of ntuple::ElectronTree and ntuple::Electron classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-24 created
 */

#pragma once

#include "SmartTree.h"

#define ELECTRON_DATA() \
    SIMPLE_VAR(Float_t, eta, 0.0) \
    SIMPLE_VAR(Float_t, phi, 0.0) \
    SIMPLE_VAR(Float_t, pt, 0.0) \
    SIMPLE_VAR(Bool_t, ecalDriven, false) \
    SIMPLE_VAR(Bool_t, hasGsfTrack, false) \
    SIMPLE_VAR(Float_t, trackPt, 0.0) \
    SIMPLE_VAR(Float_t, trackPtError, 0.0) \
    SIMPLE_VAR(Float_t, energy, 0.0) \
    SIMPLE_VAR(Float_t, caloEnergy, 0.0) \
    SIMPLE_VAR(Float_t, caloEnergyError, 0.0) \
    SIMPLE_VAR(Int_t, charge, 0) \
    SIMPLE_VAR(Int_t, pixHits, 0) \
    SIMPLE_VAR(Int_t, trkHits, 0) \
    SIMPLE_VAR(Int_t, nValidHits, 0) \
    SIMPLE_VAR(Float_t, trkD0, 0.0) \
    SIMPLE_VAR(Float_t, trkD0Error, 0.0) \
    /* ID variables */ \
    SIMPLE_VAR(Float_t, hoe, 0.0) \
    SIMPLE_VAR(Float_t, hoeDepth1, 0.0) \
    SIMPLE_VAR(Float_t, eop, 0.0) \
    SIMPLE_VAR(Float_t, sigmaEtaEta, 0.0) \
    SIMPLE_VAR(Float_t, sigmaIEtaIEta, 0.0) \
    SIMPLE_VAR(Float_t, deltaPhiTrkSC, 0.0) \
    SIMPLE_VAR(Float_t, deltaEtaTrkSC, 0.0) \
    SIMPLE_VAR(Int_t, classif, 0) \
    SIMPLE_VAR(Float_t, e1x5overe5x5, 0.0) \
    SIMPLE_VAR(Float_t, e2x5overe5x5, 0.0) \
    /* Iso variables */ \
    SIMPLE_VAR(Float_t, isoEcal03, 0.0) \
    SIMPLE_VAR(Float_t, isoHcal03, 0.0) \
    SIMPLE_VAR(Float_t, isoTrk03, 0.0) \
    SIMPLE_VAR(Float_t, isoEcal04, 0.0) \
    SIMPLE_VAR(Float_t, isoHcal04, 0.0) \
    SIMPLE_VAR(Float_t, isoTrk04, 0.0) \
    SIMPLE_VAR(Float_t, isoRel03, 0.0) \
    SIMPLE_VAR(Float_t, isoRel04, 0.0) \
    /* Vertex */ \
    SIMPLE_VAR(Float_t, vx, 0.0) \
    SIMPLE_VAR(Float_t, vy, 0.0) \
    SIMPLE_VAR(Float_t, vz, 0.0) \
    /* SC associated with electron */ \
    SIMPLE_VAR(Float_t, scEn, 0.0) \
    SIMPLE_VAR(Float_t, scEta, 0.0) \
    SIMPLE_VAR(Float_t, scPhi, 0.0) \
    SIMPLE_VAR(Float_t, scET, 0.0) \
    SIMPLE_VAR(Float_t, scRawEnergy, 0.0) \
    /* Vertex association variables */ \
    SIMPLE_VAR(Float_t, vtxDist3D, 0.0) \
    SIMPLE_VAR(Int_t, vtxIndex, 0) \
    SIMPLE_VAR(Float_t, vtxDistZ, 0.0) \
    SIMPLE_VAR(Float_t, relIso, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIso, 0.0) \
    /* PFlow isolation variable */ \
    SIMPLE_VAR(Float_t, chargedHadronIso, 0.0) \
    SIMPLE_VAR(Float_t, neutralHadronIso, 0.0) \
    SIMPLE_VAR(Float_t, photonIso, 0.0) \
    /**/ \
    SIMPLE_VAR(Int_t, missingHits, 0) \
    SIMPLE_VAR(Float_t, dB, 0.0) \
    SIMPLE_VAR(Float_t, edB, 0.0) \
    SIMPLE_VAR(Float_t, dB3d, 0.0) \
    SIMPLE_VAR(Float_t, edB3d, 0.0) \
    SIMPLE_VAR(Float_t, scE1E9, 0.0) \
    SIMPLE_VAR(Float_t, scS4S1, 0.0) \
    SIMPLE_VAR(Float_t, sckOutOfTime, 0.0) \
    SIMPLE_VAR(Float_t, scEcalIso, 0.0) \
    SIMPLE_VAR(Float_t, scHEEPEcalIso, 0.0) \
    SIMPLE_VAR(Float_t, scHEEPTrkIso, 0.0) \
    SIMPLE_VAR(Int_t, nBrems, 0) \
    SIMPLE_VAR(Float_t, fbrem, 0.0) \
    SIMPLE_VAR(Float_t, dist_vec, 0.0) \
    SIMPLE_VAR(Float_t, dCotTheta, 0.0) \
    SIMPLE_VAR(Float_t, hasMatchedConv, 0.0) \
    SIMPLE_VAR(Float_t, mva, 0.0) \
    SIMPLE_VAR(Float_t, mvaPOGTrig, 0.0) \
    SIMPLE_VAR(Float_t, mvaPOGNonTrig, 0.0) \
    SIMPLE_VAR(Bool_t, mvaPreselection, false) \
    SIMPLE_VAR(Bool_t, isTriggerElectron, false) \
    SIMPLE_VAR(Float_t, isoMVA, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIso03v1, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIso03v2, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIsoDB03v1, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIsoDB03v2, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIsoDB03v3, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIso04v1, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIso04v2, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIsoDB04v1, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIsoDB04v2, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIsoDB04v3, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIso03, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIso04, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIsoDB03, 0.0) \
    SIMPLE_VAR(Float_t, pfRelIsoDB04, 0.0) \
    SIMPLE_VAR(Int_t, selbit, 0) \
    SIMPLE_VAR(Int_t, fidFlag, 0) \
    /**/

#define SIMPLE_VAR(type, name, default_value) SIMPLE_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_TREE_BRANCH(type, name)

namespace ntuple {
class ElectronTree : public root_ext::SmartTree {
public:
    static const std::string& Name() { static const std::string name = "electrons"; return name; }
    ElectronTree() : SmartTree(Name(), "/", false) {}
    ElectronTree(const std::string& prefix, TFile& file)
        : SmartTree(prefix + "/" + Name(), file) {}

    SIMPLE_TREE_BRANCH(UInt_t, EventId, 0) \
    ELECTRON_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct Electron {
    ELECTRON_DATA()
    Electron() {}
    inline Electron(ElectronTree& tree);
};
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) name = tree.name();
#define VECTOR_VAR(type, name) name = tree.name();

namespace ntuple {
inline Electron::Electron(ElectronTree& tree)
{
    ELECTRON_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef ELECTRON_DATA
