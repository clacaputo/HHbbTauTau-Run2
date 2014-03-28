/*!
 * \file Vertex.h
 * \brief Definiton of ntuple::VertexTree and ntuple::Vertex classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-25 created
 */

#pragma once

#include "SmartTree.h"

#define VERTEX_DATA() \
    SIMPLE_VAR(Float_t, x, 0.0) \
    SIMPLE_VAR(Float_t, y, 0.0) \
    SIMPLE_VAR(Float_t, z, 0.0) \
    SIMPLE_VAR(Float_t, xErr, 0.0) \
    SIMPLE_VAR(Float_t, yErr, 0.0) \
    SIMPLE_VAR(Float_t, zErr, 0.0) \
    SIMPLE_VAR(Float_t, rho, 0.0) \
    SIMPLE_VAR(Float_t, chi2, 0.0) \
    SIMPLE_VAR(Float_t, ndf, 0.0) \
    SIMPLE_VAR(UInt_t, ntracks, 0) \
    SIMPLE_VAR(UInt_t, ntracksw05, 0) \
    SIMPLE_VAR(Bool_t, isfake, false) \
    SIMPLE_VAR(Bool_t, isvalid, false) \
    SIMPLE_VAR(Float_t, sumPt, 0.0) \
    /**/

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;
DATA_CLASS(ntuple, Vertex, VERTEX_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) SIMPLE_DATA_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, VertexTree, VERTEX_DATA, Vertex, "vertices")
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) AddSimpleBranch(#name, data.name);
#define VECTOR_VAR(type, name) AddVectorBranch(#name, data.name);
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, VertexTree, VERTEX_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef VERTEX_DATA
