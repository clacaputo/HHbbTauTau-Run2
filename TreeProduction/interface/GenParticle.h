/*!
 * \file GenParticle.h
 * \brief Definiton of ntuple::GenParticleTree and ntuple::GenParticle classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-24 created
 */

#pragma once

#include "SmartTree.h"

#define GEN_PARTICLE_DATA() \
    SIMPLE_VAR(Int_t, PdgId, 0) \
    SIMPLE_VAR(Int_t, Status, 0) \
    SIMPLE_VAR(Int_t, Charge, 0) \
    SIMPLE_VAR(Float_t, E, 0.0) \
    SIMPLE_VAR(Float_t, Px, 0.0) \
    SIMPLE_VAR(Float_t, Py, 0.0) \
    SIMPLE_VAR(Float_t, Pz, 0.0) \
    SIMPLE_VAR(Float_t, X, 0.0) \
    SIMPLE_VAR(Float_t, Y, 0.0) \
    SIMPLE_VAR(Float_t, Z, 0.0) \
    SIMPLE_VAR(UInt_t, Index, 0) \
    VECTOR_VAR(UInt_t, Mother_Indexes) \
    VECTOR_VAR(UInt_t, Daughter_Indexes) \
    /**/

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;
DATA_CLASS(ntuple, GenParticle, GEN_PARTICLE_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) SIMPLE_DATA_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, GenParticleTree, GEN_PARTICLE_DATA, GenParticle, "genParticles")
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) AddSimpleBranch(#name, data.name);
#define VECTOR_VAR(type, name) AddVectorBranch(#name, data.name);
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, GenParticleTree, GEN_PARTICLE_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef GEN_PARTICLE_DATA
