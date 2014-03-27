/*!
 * \file Trigger.h
 * \brief Definiton of ntuple::TriggerTree and ntuple::Trigger classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-25 created
 */

#pragma once

#include "SmartTree.h"

#define TRIGGER_DATA() \
    VECTOR_VAR(Bool_t, l1physbits) \
    VECTOR_VAR(Bool_t, l1techbits) \
    VECTOR_VAR(std::string, hltpaths) \
    VECTOR_VAR(Bool_t, hltresults) \
    VECTOR_VAR(UInt_t, hltprescales) \
    /**/

#define SIMPLE_VAR(type, name, default_value) SIMPLE_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_TREE_BRANCH(type, name)

namespace ntuple {
class TriggerTree : public root_ext::SmartTree {
public:
    static const std::string& Name() { static const std::string name = "Triggers"; return name; }
    TriggerTree() : SmartTree(Name(), "/", false) {}
    TriggerTree(const std::string& prefix, TFile& file)
        : SmartTree(prefix + "/" + Name(), file) {}

    SIMPLE_TREE_BRANCH(UInt_t, EventId, 0) \
    TRIGGER_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct Trigger {
    TRIGGER_DATA()
    Trigger() {}
    inline Trigger(TriggerTree& tree);
};
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) name = tree.name();
#define VECTOR_VAR(type, name) name = tree.name();

namespace ntuple {
inline Trigger::Trigger(TriggerTree& tree)
{
    TRIGGER_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef TRIGGER_DATA
