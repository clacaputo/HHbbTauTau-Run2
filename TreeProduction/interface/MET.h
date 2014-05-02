/*!
 * \file MET.h
 * \brief Definiton of ntuple::METTree and ntuple::MET classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-25 created
 */

#pragma once

#include "SmartTree.h"

#define MET_DATA() \
    SIMPLE_VAR(Float_t, pt) \
    SIMPLE_VAR(Float_t, phi) \
    SIMPLE_VAR(Float_t, eta) \
    SIMPLE_VAR(Float_t, sumEt) \
    SIMPLE_VAR(Float_t, pt_uncorrected) \
    SIMPLE_VAR(Float_t, phi_uncorrected) \
    SIMPLE_VAR(Float_t, sumEt_uncorrected) \
    VECTOR_VAR(Float_t, significanceMatrix) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, MET, MET_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, METTree, MET_DATA, MET, "METs", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, METTree, MET_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef MET_DATA
