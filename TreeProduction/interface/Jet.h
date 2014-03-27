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
    /* 4-momentum */ \
    SIMPLE_VAR(Float_t, pt, 0.0) \
    SIMPLE_VAR(Float_t, pt_raw, 0.0) \
    SIMPLE_VAR(Float_t, eta, 0.0) \
    SIMPLE_VAR(Float_t, phi, 0.0) \
    SIMPLE_VAR(Float_t, energy, 0.0) \
    SIMPLE_VAR(Float_t, energy_raw, 0.0) \
    /* Jet energy correction */ \
    SIMPLE_VAR(Float_t, jecUnc, 0.0) \
    SIMPLE_VAR(Float_t, resJEC, 0.0) \
    SIMPLE_VAR(Int_t, partonFlavour, 0) \
    /* Jet identification in high PU environment */ \
    SIMPLE_VAR(Float_t, puIdMVA, 0.0) \
    SIMPLE_VAR(Int_t, puIdBits, 0) \
    /* Energy fraction information */ \
    SIMPLE_VAR(Float_t, chargedEmEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, chargedHadronEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, chargedMuEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, electronEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, muonEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, neutralEmEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, neutralHadronEnergyFraction, 0.0) \
    SIMPLE_VAR(Float_t, photonEnergyFraction, 0.0) \
    /* Multiplicity information */ \
    SIMPLE_VAR(Int_t, chargedHadronMultiplicity, 0) \
    SIMPLE_VAR(Int_t, chargedMultiplicity, 0) \
    SIMPLE_VAR(Int_t, electronMultiplicity, 0) \
    SIMPLE_VAR(Int_t, muonMultiplicity, 0) \
    SIMPLE_VAR(Int_t, neutralHadronMultiplicity, 0) \
    SIMPLE_VAR(Int_t, neutralMultiplicity, 0) \
    SIMPLE_VAR(Int_t, photonMultiplicity, 0) \
    SIMPLE_VAR(UInt_t, nConstituents, 0) \
    SIMPLE_VAR(Bool_t, passLooseID, 0) \
    SIMPLE_VAR(Bool_t, passTightID, 0) \
    /**/

#define B_TAG_DATA() \
    SIMPLE_VAR(Float_t, trackCountingHighEffBJetTags, 0.0) \
    SIMPLE_VAR(Float_t, trackCountingHighPurBJetTags, 0.0) \
    SIMPLE_VAR(Float_t, simpleSecondaryVertexHighEffBJetTags, 0.0) \
    SIMPLE_VAR(Float_t, simpleSecondaryVertexHighPurBJetTags, 0.0) \
    SIMPLE_VAR(Float_t, jetProbabilityBJetTags, 0.0) \
    SIMPLE_VAR(Float_t, jetBProbabilityBJetTags, 0.0) \
    SIMPLE_VAR(Float_t, combinedSecondaryVertexBJetTags, 0.0) \
    SIMPLE_VAR(Float_t, combinedSecondaryVertexMVABJetTags, 0.0) \
    SIMPLE_VAR(Float_t, combinedInclusiveSecondaryVertexBJetTags, 0.0) \
    SIMPLE_VAR(Float_t, combinedMVABJetTags, 0.0) \
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
    B_TAG_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct Jet {
    JET_DATA()
    B_TAG_DATA()
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
    B_TAG_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef JET_DATA
