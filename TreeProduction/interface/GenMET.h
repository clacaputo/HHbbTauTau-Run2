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
    SIMPLE_VAR(Float_t, met, 0.0) \
    SIMPLE_VAR(Float_t, metphi, 0.0) \
    SIMPLE_VAR(Float_t, sumet, 0.0) \
    /**/

#define SIMPLE_VAR(type, name, default_value) SIMPLE_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_TREE_BRANCH(type, name)

namespace ntuple {
class GenMETTree : public root_ext::SmartTree {
public:
    static const std::string& Name() { static const std::string name = "genMETs"; return name; }
    GenMETTree() : SmartTree(Name(), "/", false) {}
    GenMETTree(TFile& file)
        : SmartTree(Name(), file) {}

    SIMPLE_TREE_BRANCH(UInt_t, EventId, 0) \
    GENMET_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct GenMET {
    GENMET_DATA()
    GenMET() {}
    inline GenMET(GenMETTree& tree);
};
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) name = tree.name();
#define VECTOR_VAR(type, name) name = tree.name();

namespace ntuple {
inline GenMET::GenMET(GenMETTree& tree)
{
    GENMET_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef GENMET_DATA
