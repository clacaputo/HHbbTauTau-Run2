/*!
 * \file Event.h
 * \brief Definiton of ntuple::EventTree and ntuple::Event classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-24 created
 */

#pragma once

#include "SmartTree.h"

#define EVENT_DATA() \
    SIMPLE_VAR(UInt_t, run, 0) \
    SIMPLE_VAR(UInt_t, event, 0) \
    SIMPLE_VAR(UInt_t, lumis, 0) \
    SIMPLE_VAR(UInt_t, bunch, 0) \
    SIMPLE_VAR(UInt_t, orbit, 0) \
    SIMPLE_VAR(Double_t, time, 0.0) \
    SIMPLE_VAR(Bool_t, isdata, false) \
    SIMPLE_VAR(Bool_t, isPhysDeclared, false) \
    SIMPLE_VAR(Bool_t, isBPTX0, false) \
    SIMPLE_VAR(Bool_t, isBSCMinBias, false) \
    SIMPLE_VAR(Bool_t, isBSCBeamHalo, false) \
    SIMPLE_VAR(Bool_t, isPrimaryVertex, false) \
    SIMPLE_VAR(Bool_t, isBeamScraping, false) \
    SIMPLE_VAR(Bool_t, passHBHENoiseFilter, false) \
    VECTOR_VAR(Int_t, nPU) \
    VECTOR_VAR(Int_t, bunchCrossing) \
    VECTOR_VAR(Int_t, trueNInt) \
    /**/

#define SIMPLE_VAR(type, name, default_value) SIMPLE_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_TREE_BRANCH(type, name)

namespace ntuple {
class EventTree : public root_ext::SmartTree {
public:
    static const std::string& Name() { static const std::string name = "events"; return name; }
    EventTree() : SmartTree(Name()) {}
    EventTree(const std::string& prefix, TFile& file)
        : SmartTree(prefix + "/" + Name(), file) {}

    EVENT_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct Event {
    EVENT_DATA()
    Event() {}
    inline Electron(EventTree& tree);
};
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) name = tree.name();
#define VECTOR_VAR(type, name) name = tree.name();

namespace ntuple {
inline Event::Event(EventTree& tree)
{
    EVENT_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef EVENT_DATA
