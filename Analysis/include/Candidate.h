/*!
 * \file Candidate.h
 * \brief Definition of Candidate class to store reconstructed object candidate.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-19 created
 */

#pragma once
#include <TLorentzVector.h>
#include <vector>

namespace analysis {

class Candidate;
typedef std::vector<Candidate> CandidateVector;
typedef std::vector<const Candidate*> CandidatePtrVector;

class Vertex;
typedef std::vector<Vertex> VertexVector;

class Candidate {
public:
    enum Type { Unknown, Mu, Electron, Tau, Bjet, Higgs, Resonance };

    Type type;
    size_t index;
    TLorentzVector momentum;
    CandidatePtrVector daughters;


    Candidate() : type(Unknown), index(0){}

    template<typename NtupleObject>
    Candidate(Type _type, size_t _index, const NtupleObject& ntupleObject) :
        type(_type), index(_index) {
        momentum.SetPtEtaPhiE(ntupleObject.pt, ntupleObject.eta, ntupleObject.phi, ntupleObject.energy);     
    }

    Candidate(Type _type, const Candidate& daughter1, const Candidate& daughter2) : type(_type), index(-1) {
        daughters.push_back(&daughter1);
        daughters.push_back(&daughter2);
        momentum = daughter1.momentum + daughter2.momentum;
    }

    bool operator < (const Candidate& other) const
    {
        return momentum.Pt() < other.momentum.Pt();
    }

    bool operator == (const Candidate& other) const
    {
        return type == other.type && index == other.index;
    }

    bool operator != (const Candidate& other) const
    {
        return type != other.type || index != other.index;
    }
};

class Vertex{
public:

    TVector3 position;
    size_t index;
    double sumPt;
    unsigned ndf;

    Vertex() : index(0), sumPt(0), ndf(0) {}

    template <typename NtupleObject>
    Vertex(size_t _index, const NtupleObject& ntupleObject) : index(_index), sumPt(ntupleObject.sumPt),
        ndf(ntupleObject.ndf)
    {
       position = TVector3(ntupleObject.x,ntupleObject.y,ntupleObject.z);
    }

    bool operator < (const Vertex& other) const
    {
        return sumPt < other.sumPt;
    }

};

} // analysis
