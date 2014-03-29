/*!
 * \file GenEvent.h
 * \brief Definiton of ntuple::GenEventTree and ntuple::GenEvent classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-24 created
 */

#pragma once

#include "SmartTree.h"

#define GENEVENT_DATA() \
    SIMPLE_VAR(UInt_t, processID) \
    SIMPLE_VAR(Float_t, ptHat) \
    VECTOR_VAR(Double_t, pdfWeights) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, GenEvent, GENEVENT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, GenEventTree, GENEVENT_DATA, GenEvent, "genEvents")
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, GenEventTree, GENEVENT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef GENEVENT_DATA
