/*!
 * \file Track.h
 * \brief Definiton of ntuple::TrackTree and ntuple::Track classes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-25 created
 */

#pragma once

#include "SmartTree.h"

#define TRACK_DATA() \
    SIMPLE_VAR(Float_t, eta, 0.0) \
    SIMPLE_VAR(Float_t, etaError, 0.0) \
    SIMPLE_VAR(Float_t, theta, 0.0) \
    SIMPLE_VAR(Float_t, thetaError, 0.0) \
    SIMPLE_VAR(Float_t, phi, 0.0) \
    SIMPLE_VAR(Float_t, phiError, 0.0) \
    SIMPLE_VAR(Float_t, pt, 0.0) \
    SIMPLE_VAR(Float_t, ptError, 0.0) \
    SIMPLE_VAR(Float_t, qoverp, 0.0) \
    SIMPLE_VAR(Float_t, qoverpError, 0.0) \
    SIMPLE_VAR(Int_t, charge, 0) \
    /* */ \
    SIMPLE_VAR(UInt_t, nValidHits, 0) \
    SIMPLE_VAR(UInt_t, nLostHits, 0) \
    SIMPLE_VAR(Float_t, validFraction, 0.0) \
    SIMPLE_VAR(Int_t, nValidTrackerHits, 0) \
    SIMPLE_VAR(Int_t, nValidPixelHits, 0) \
    SIMPLE_VAR(Int_t, nValidStripHits, 0) \
    SIMPLE_VAR(Int_t, trackerLayersWithMeasurement, 0) \
    SIMPLE_VAR(Int_t, pixelLayersWithMeasurement, 0) \
    SIMPLE_VAR(Int_t, stripLayersWithMeasurement, 0) \
    /* */ \
    SIMPLE_VAR(Float_t, dxy, 0.0) \
    SIMPLE_VAR(Float_t, dxyError, 0.0) \
    SIMPLE_VAR(Float_t, dz, 0.0) \
    SIMPLE_VAR(Float_t, dzError, 0.0) \
    SIMPLE_VAR(Float_t, chi2, 0.0) \
    SIMPLE_VAR(Float_t, ndof, 0.0) \
    SIMPLE_VAR(Float_t, vx, 0.0) \
    SIMPLE_VAR(Float_t, vy, 0.0) \
    SIMPLE_VAR(Float_t, vz, 0.0) \
    /**/

#define SIMPLE_VAR(type, name, default_value) SIMPLE_TREE_BRANCH(type, name, default_value)
#define VECTOR_VAR(type, name) VECTOR_TREE_BRANCH(type, name)

namespace ntuple {
class TrackTree : public root_ext::SmartTree {
public:
    static const std::string& Name() { static const std::string name = "tracks"; return name; }
    TrackTree() : SmartTree(Name(), "/", false) {}
    TrackTree(const std::string& prefix, TFile& file)
        : SmartTree(prefix + "/" + Name(), file) {}

    SIMPLE_TREE_BRANCH(UInt_t, EventId, 0) \
    TRACK_DATA()
};
} // ntuple


#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) type name;
#define VECTOR_VAR(type, name) std::vector< type > name;

namespace ntuple {
struct Track {
    TRACK_DATA()
    Track() {}
    inline Track(TrackTree& tree);
};
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name, default_value) name = tree.name();
#define VECTOR_VAR(type, name) name = tree.name();

namespace ntuple {
inline Track::Track(TrackTree& tree)
{
    TRACK_DATA()
}
} // ntuple

#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef TRACK_DATA
