/*!
 * \file TriggerObject.h
 * \brief Definiton of ntuple::TriggerObjectTree and ntuple::TriggerObject classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-25 created
 */

#pragma once

#include "SmartTree.h"

#define TRIGGER_OBJECT_DATA() \
    SIMPLE_VAR(Double_t, energy) \
    SIMPLE_VAR(Double_t, pt) \
    SIMPLE_VAR(Double_t, eta) \
    SIMPLE_VAR(Double_t, phi) \
    SIMPLE_VAR(Int_t, pdgId) \
    VECTOR_VAR(std::string, pathNames) \
    VECTOR_VAR(UInt_t, pathValues) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, TriggerObject, TRIGGER_OBJECT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, TriggerObjectTree, TRIGGER_OBJECT_DATA, TriggerObject, "triggerObjects", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, TriggerObjectTree, TRIGGER_OBJECT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef TRIGGER_OBJECT_DATA
