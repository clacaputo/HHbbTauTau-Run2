/*!
 * \file GenMET.h
 * \brief Definiton of ntuple::GenMETTree and ntuple::GenMET classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-24 created
 */

#pragma once

#include "SmartTree.h"

#define GENMET_DATA() \
    SIMPLE_VAR(Double_t, met) \
    SIMPLE_VAR(Double_t, metphi) \
    SIMPLE_VAR(Double_t, sumet) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, GenMET, GENMET_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, GenMETTree, GENMET_DATA, GenMET, "genMETs", true)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, GenMETTree, GENMET_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef GENMET_DATA
