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
    /* 4-momentum */ \
    SIMPLE_VAR(Double_t, eta) \
    SIMPLE_VAR(Double_t, phi) \
    SIMPLE_VAR(Double_t, pt) \
    SIMPLE_VAR(Double_t, energy) \
    SIMPLE_VAR(Double_t, mass) \
    SIMPLE_VAR(Double_t, caloEnergy) \
    SIMPLE_VAR(Double_t, caloEnergyError) \
    /* Origin */ \
    SIMPLE_VAR(Bool_t, ecalDriven) \
    SIMPLE_VAR(Bool_t, hasGsfTrack) \
    /* Gsf based electrons */ \
    SIMPLE_VAR(Double_t, trackPt) \
    SIMPLE_VAR(Double_t, trackPtError) \
    SIMPLE_VAR(Int_t, charge) \
    SIMPLE_VAR(Int_t, pixHits) \
    SIMPLE_VAR(Int_t, trkHits) \
    SIMPLE_VAR(UInt_t, nValidHits) \
    SIMPLE_VAR(Double_t, trkD0) \
    SIMPLE_VAR(Double_t, trkD0Error) \
    SIMPLE_VAR(Int_t, missingHits) \
    /* ID variables */ \
    SIMPLE_VAR(Double_t, hcalOverEcal) \
    SIMPLE_VAR(Double_t, hcalDepth1OverEcal) \
    SIMPLE_VAR(Double_t, eSuperClusterOverP) \
    SIMPLE_VAR(Double_t, sigmaEtaEta) \
    SIMPLE_VAR(Double_t, sigmaIEtaIEta) \
    SIMPLE_VAR(Double_t, deltaPhiTrkSC) \
    SIMPLE_VAR(Double_t, deltaEtaTrkSC) \
    SIMPLE_VAR(Int_t, classification) \
    SIMPLE_VAR(Double_t, e1x5overe5x5) \
    SIMPLE_VAR(Double_t, e2x5overe5x5) \
    /* Iso variables */ \
    SIMPLE_VAR(Double_t, isoEcal03) \
    SIMPLE_VAR(Double_t, isoHcal03) \
    SIMPLE_VAR(Double_t, isoTrk03) \
    SIMPLE_VAR(Double_t, isoEcal04) \
    SIMPLE_VAR(Double_t, isoHcal04) \
    SIMPLE_VAR(Double_t, isoTrk04) \
    SIMPLE_VAR(Double_t, isoRel03) \
    SIMPLE_VAR(Double_t, isoRel04) \
    /* Vertex */ \
    SIMPLE_VAR(Double_t, vx) \
    SIMPLE_VAR(Double_t, vy) \
    SIMPLE_VAR(Double_t, vz) \
    /* SC associated with electron */ \
    SIMPLE_VAR(Double_t, scEn) \
    SIMPLE_VAR(Double_t, scEta) \
    SIMPLE_VAR(Double_t, scPhi) \
    SIMPLE_VAR(Double_t, scET) \
    SIMPLE_VAR(Double_t, scRawEnergy) \
    /* Vertex association variables */ \
    SIMPLE_VAR(Double_t, vtxDist3D) \
    SIMPLE_VAR(Int_t, vtxIndex) \
    SIMPLE_VAR(Double_t, vtxDistZ) \
    SIMPLE_VAR(Double_t, relIso) \
    SIMPLE_VAR(Double_t, pfRelIso) \
    /* PFlow isolation variable */ \
    SIMPLE_VAR(Double_t, chargedHadronIso) \
    SIMPLE_VAR(Double_t, neutralHadronIso) \
    SIMPLE_VAR(Double_t, photonIso) \
    /* IP information */ \
    SIMPLE_VAR(Double_t, dB) \
    SIMPLE_VAR(Double_t, edB) \
    SIMPLE_VAR(Double_t, dB3d) \
    SIMPLE_VAR(Double_t, edB3d) \
    SIMPLE_VAR(Int_t, nBrems) \
    SIMPLE_VAR(Double_t, fbrem) \
    SIMPLE_VAR(Double_t, dist_vec) \
    SIMPLE_VAR(Double_t, dCotTheta) \
    SIMPLE_VAR(Bool_t, hasMatchedConversion) \
    /* MVA */ \
    SIMPLE_VAR(Double_t, mva) \
    SIMPLE_VAR(Double_t, mvaPOGTrig) \
    SIMPLE_VAR(Double_t, mvaPOGNonTrig) \
    SIMPLE_VAR(Bool_t, mvaPreselection) \
    SIMPLE_VAR(Bool_t, isTriggerElectron) \
    SIMPLE_VAR(Double_t, isoMVA) \
    SIMPLE_VAR(Double_t, pfRelIso03v1) \
    SIMPLE_VAR(Double_t, pfRelIso03v2) \
    SIMPLE_VAR(Double_t, pfRelIsoDB03v1) \
    SIMPLE_VAR(Double_t, pfRelIsoDB03v2) \
    SIMPLE_VAR(Double_t, pfRelIsoDB03v3) \
    SIMPLE_VAR(Double_t, pfRelIso04v1) \
    SIMPLE_VAR(Double_t, pfRelIso04v2) \
    SIMPLE_VAR(Double_t, pfRelIsoDB04v1) \
    SIMPLE_VAR(Double_t, pfRelIsoDB04v2) \
    SIMPLE_VAR(Double_t, pfRelIsoDB04v3) \
    SIMPLE_VAR(Double_t, pfRelIso03) \
    SIMPLE_VAR(Double_t, pfRelIso04) \
    SIMPLE_VAR(Double_t, pfRelIsoDB03) \
    SIMPLE_VAR(Double_t, pfRelIsoDB04) \
    SIMPLE_VAR(UInt_t, fidFlag) \
    /* Trigger match information */ \
    VECTOR_VAR(std::string, matchedTriggerPaths)
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Electron, ELECTRON_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, ElectronTree, ELECTRON_DATA, Electron, "electrons", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, ElectronTree, ELECTRON_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef ELECTRON_DATA
