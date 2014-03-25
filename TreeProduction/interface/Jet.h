/*!
 * \file Jet.h
 * \brief Definiton of ntuple::JetTree and ntuple::Jet classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-25 created
 */

#pragma once

#include "SmartTree.h"

#define JET_DATA() \
    SIMPLE_VAR(Float_t, eta, 0.0) \
    SIMPLE_VAR(Float_t, phi, 0.0) \
    SIMPLE_VAR(Float_t, pt, 0.0) \
    SIMPLE_VAR(Float_t, pt_raw, 0.0) \
    SIMPLE_VAR(Float_t, energy, 0.0) \
    SIMPLE_VAR(Float_t, energy_raw, 0.0) \
    SIMPLE_VAR(Float_t, jecUnc, 0.0) \
    SIMPLE_VAR(Float_t, resJEC, 0.0) \
    SIMPLE_VAR(Int_t, partonFlavour, 0) \
    SIMPLE_VAR(Float_t, puIdMVA, 0.0) \
    SIMPLE_VAR(Int_t, puIdFlag, 0) \
    SIMPLE_VAR(Int_t, puIdBits, 0) \
    SIMPLE_VAR(Float_t, chargedEmEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, chargedHadronEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, chargedMuEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, electronEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, muonEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, neutralEmEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, neutralHadronEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, photonEnergyFraction, 0.0) \
    SIMPLE_VAR(Int_t, chargedHadronMultiplicity, 0) \
    SIMPLE_VAR(Int_t, chargedMultiplicity, 0) \
    SIMPLE_VAR(Int_t, electronMultiplicity, 0) \
    SIMPLE_VAR(Int_t, muonMultiplicity, 0) \
    SIMPLE_VAR(Int_t, neutralHadronMultiplicity, 0) \
    SIMPLE_VAR(Int_t, neutralMultiplicity, 0) \
    SIMPLE_VAR(Int_t, photonMultiplicity, 0) \
    SIMPLE_VAR(Int_t, nConstituents, 0) \
    SIMPLE_VAR(Float_t, trackCountingHighEffBTag, 0.0) \
    SIMPLE_VAR(Float_t, trackCountingHighPurBTag, 0.0) \
    SIMPLE_VAR(Float_t, simpleSecondaryVertexHighEffBTag, 0.0) \
    SIMPLE_VAR(Float_t, simpleSecondaryVertexHighPurBTag, 0.0) \
    SIMPLE_VAR(Float_t, jetProbabilityBTag, 0.0) \
    SIMPLE_VAR(Float_t, jetBProbabilityBTag, 0.0) \
    SIMPLE_VAR(Float_t, combinedSecondaryVertexBTag, 0.0) \
    SIMPLE_VAR(Float_t, combinedSecondaryVertexMVABTag, 0.0) \
    SIMPLE_VAR(Float_t, combinedInclusiveSecondaryVertexBTag, 0.0) \
    SIMPLE_VAR(Float_t, combinedMVABTag, 0.0) \
    SIMPLE_VAR(Int_t, passLooseID, 0) \
    SIMPLE_VAR(Int_t, passTightID, 0) \
    SIMPLE_VAR(Int_t, selbit, 0) \
    /**/

#define SIMPLE_VAR(type, name, default_value) SIMPLE_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_TREE_BRANCH(type, name)

namespace ntuple {
class JetTree : public root_ext::SmartTree {
public:
    static const std::string& Name() { static const std::string name = "jets"; return name; }
    JetTree() : SmartTree(Name(), "/", false) {}
    JetTree(const std::string& prefix, TFile& file)
        : SmartTree(prefix + "/" + Name(), file) {}

    SIMPLE_TREE_BRANCH(UInt_t, EventId, 0) \
    JET_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct Jet {
    JET_DATA()
    Jet() {}
    inline Jet(JetTree& tree);
};
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) name = tree.name();
#define VECTOR_VAR(type, name) name = tree.name();

namespace ntuple {
inline Jet::Jet(JetTree& tree)
{
    JET_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef JET_DATA
