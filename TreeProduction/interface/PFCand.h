/*!
 * \file PFCand.h
 * \brief Definiton of ntuple::PFCandTree and ntuple::PFCand classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-30 created
 */

#pragma once

#include "SmartTree.h"

#define PFCANDIDATE_DATA() \
    /* 4-momentum */ \
    SIMPLE_VAR(Double_t, eta) \
    SIMPLE_VAR(Double_t, phi) \
    SIMPLE_VAR(Double_t, pt) \
    SIMPLE_VAR(Double_t, energy) \
    SIMPLE_VAR(Double_t, mass) \
    /* Charge */ \
    SIMPLE_VAR(Int_t, charge) \
    /* Vertex */ \
    SIMPLE_VAR(Double_t, vx) \
    SIMPLE_VAR(Double_t, vy) \
    SIMPLE_VAR(Double_t, vz) \
    /* track info */ \
    SIMPLE_VAR(Bool_t, haveTrackInfo) \
    SIMPLE_VAR(Double_t, trk_eta) \
    SIMPLE_VAR(Double_t, trk_phi) \
    SIMPLE_VAR(Double_t, trk_pt) \
    SIMPLE_VAR(Double_t, trk_vx) \
    SIMPLE_VAR(Double_t, trk_vy) \
    SIMPLE_VAR(Double_t, trk_vz) \
    /**/


#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, PFCand, PFCANDIDATE_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, PFCandTree, PFCANDIDATE_DATA, PFCand, "PFCand", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, PFCandTree, PFCANDIDATE_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef PFCANDIDATE_DATA

