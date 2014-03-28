/*!
 * \file Vertex.h
 * \brief Definiton of ntuple::VertexTree and ntuple::Vertex classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-25 created
 */

#pragma once

#include "SmartTree.h"

#define VERTEX_DATA() \
    SIMPLE_VAR(Float_t, x, 0.0) \
    SIMPLE_VAR(Float_t, y, 0.0) \
    SIMPLE_VAR(Float_t, z, 0.0) \
    SIMPLE_VAR(Float_t, xErr, 0.0) \
    SIMPLE_VAR(Float_t, yErr, 0.0) \
    SIMPLE_VAR(Float_t, zErr, 0.0) \
    SIMPLE_VAR(Float_t, rho, 0.0) \
    SIMPLE_VAR(Float_t, chi2, 0.0) \
    SIMPLE_VAR(Float_t, ndf, 0.0) \
    SIMPLE_VAR(UInt_t, ntracks, 0) \
    SIMPLE_VAR(UInt_t, ntracksw05, 0) \
    SIMPLE_VAR(Bool_t, isfake, false) \
    SIMPLE_VAR(Bool_t, isvalid, false) \
    SIMPLE_VAR(Float_t, sumPt, 0.0) \
    /**/

#define SIMPLE_VAR(type, name, default_value) SIMPLE_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_TREE_BRANCH(type, name)

namespace ntuple {
class VertexTree : public root_ext::SmartTree {
public:
    static const std::string& Name() { static const std::string name = "vertices"; return name; }
    VertexTree() : SmartTree(Name(), "/", false) {}
    VertexTree(TFile& file)
        : SmartTree(Name(), file) {}

    SIMPLE_TREE_BRANCH(UInt_t, EventId, 0) \
    VERTEX_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct Vertex {
    VERTEX_DATA()
    Vertex() {}
    inline Vertex(VertexTree& tree);
};

typedef std::vector<Vertex> VertexVector;

} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) name = tree.name();
#define VECTOR_VAR(type, name) name = tree.name();

namespace ntuple {
inline Vertex::Vertex(VertexTree& tree)
{
    VERTEX_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef VERTEX_DATA
