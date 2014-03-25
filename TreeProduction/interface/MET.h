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
    SIMPLE_VAR(Float_t, met, 0.0) \
    SIMPLE_VAR(Float_t, metphi, 0.0) \
    SIMPLE_VAR(Float_t, sumet, 0.0) \
    SIMPLE_VAR(Float_t, metuncorr, 0.0) \
    SIMPLE_VAR(Float_t, metphiuncorr, 0.0) \
    SIMPLE_VAR(Float_t, sumetuncorr, 0.0) \
    /**/

#define SIMPLE_VAR(type, name, default_value) SIMPLE_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_TREE_BRANCH(type, name)

namespace ntuple {
class METTree : public root_ext::SmartTree {
public:
    static const std::string& Name() { static const std::string name = "METs"; return name; }
    METTree() : SmartTree(Name(), "/", false) {}
    METTree(const std::string& prefix, TFile& file)
        : SmartTree(prefix + "/" + Name(), file) {}

    SIMPLE_TREE_BRANCH(UInt_t, EventId, 0) \
    MET_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct MET {
    MET_DATA()
    MET() {}
    inline MET(METTree& tree);
};
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) name = tree.name();
#define VECTOR_VAR(type, name) name = tree.name();

namespace ntuple {
inline MET::MET(METTree& tree)
{
    MET_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef MET_DATA
