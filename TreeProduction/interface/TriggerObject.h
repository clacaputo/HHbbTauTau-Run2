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
    SIMPLE_VAR(Float_t, energy, 0.0) \
    SIMPLE_VAR(Float_t, pt, 0.0) \
    SIMPLE_VAR(Float_t, eta, 0.0) \
    SIMPLE_VAR(Float_t, phi, 0.0) \
    VECTOR_VAR(std::string, pathNames) \
    VECTOR_VAR(UInt_t, pathValues) \
    /**/

#define SIMPLE_VAR(type, name, default_value) SIMPLE_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_TREE_BRANCH(type, name)

namespace ntuple {
class TriggerObjectTree : public root_ext::SmartTree {
public:
    static const std::string& Name() { static const std::string name = "triggerObjects"; return name; }
    TriggerObjectTree() : SmartTree(Name(), "/", false) {}
    TriggerObjectTree(const std::string& prefix, TFile& file)
        : SmartTree(prefix + "/" + Name(), file) {}

    SIMPLE_TREE_BRANCH(UInt_t, EventId, 0) \
    TRIGGER_OBJECT_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct TriggerObject {
    TRIGGER_OBJECT_DATA()
    TriggerObject() {}
    inline TriggerObject(TriggerObjectTree& tree);
};
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) name = tree.name();
#define VECTOR_VAR(type, name) name = tree.name();

namespace ntuple {
inline TriggerObject::TriggerObject(TriggerObjectTree& tree)
{
    TRIGGER_OBJECT_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef TRIGGER_OBJECT_DATA
