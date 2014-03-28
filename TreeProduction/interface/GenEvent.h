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
    SIMPLE_VAR(UInt_t, processID, 0) \
    SIMPLE_VAR(Float_t, ptHat, 0.0) \
    VECTOR_VAR(Double_t, pdfWeights) \
    /**/

#define SIMPLE_VAR(type, name, default_value) SIMPLE_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_TREE_BRANCH(type, name)

namespace ntuple {
class GenEventTree : public root_ext::SmartTree {
public:
    static const std::string& Name() { static const std::string name = "genEvents"; return name; }
    GenEventTree() : SmartTree(Name(), "/", false) {}
    GenEventTree(TFile& file)
        : SmartTree(Name(), file) {}

    SIMPLE_TREE_BRANCH(UInt_t, EventId, 0) \
    GENEVENT_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct GenEvent {
    GENEVENT_DATA()
    GenEvent() {}
    inline GenEvent(GenEventTree& tree);
};
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) name = tree.name();
#define VECTOR_VAR(type, name) name = tree.name();

namespace ntuple {
inline GenEvent::GenEvent(GenEventTree& tree)
{
    GENEVENT_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef GENEVENT_DATA
