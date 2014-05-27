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
#include "EventDescriptor.h"
#include "Particles.h"

namespace analysis {

class Candidate;

typedef std::vector<Candidate> CandidateVector;

class Vertex;
typedef std::vector<Vertex> VertexVector;


class Candidate {
public:
    enum Type { Unknown, Mu, Electron, Tau, Jet, Bjet, Z, Higgs, Resonance };
    static const std::map<Type,particles::ParticleCode>& typeToPdgMap(){
        static const std::map<Type,particles::ParticleCode> map = {{Type::Mu, particles::mu},
                                                                      {Type::Electron, particles::e},
                                                                      {Type::Tau, particles::tau},
                                                                      {Type::Jet, particles::Jet}};
        return map;
    }
    Type type;
    size_t index;
    TLorentzVector momentum;
    CandidateVector daughters;
    CandidateVector finalStateDaughters;
    int charge;
    TVector3 vertexPosition;

    static int UnknownCharge() { return std::numeric_limits<int>::max(); }

    Candidate() : type(Unknown), index(0), charge(UnknownCharge()){}

    template<typename NtupleObject>
    Candidate(Type _type, size_t _index, const NtupleObject& ntupleObject) :
        type(_type), index(_index), charge(ntupleObject.charge),
        vertexPosition(ntupleObject.vx, ntupleObject.vy, ntupleObject.vz) {
        momentum.SetPtEtaPhiE(ntupleObject.pt, ntupleObject.eta, ntupleObject.phi, ntupleObject.energy);
    }

    Candidate(Type _type, size_t _index, const ntuple::Jet& ntupleObject) :
        type(_type), index(_index), charge(UnknownCharge()) {
        momentum.SetPtEtaPhiE(ntupleObject.pt, ntupleObject.eta, ntupleObject.phi, ntupleObject.energy);
    }

    Candidate(Type _type, size_t _index, const ntuple::Muon& ntupleObject) :
        type(_type), index(_index), charge(ntupleObject.charge) {
        momentum.SetPtEtaPhiM(ntupleObject.pt, ntupleObject.eta, ntupleObject.phi, 0.1057);
    }


    Candidate(Type _type, const Candidate& daughter1, const Candidate& daughter2) : type(_type),
        index(std::numeric_limits<size_t>::max()) {
        daughters.push_back(daughter1);
        daughters.push_back(daughter2);
        momentum = daughter1.momentum + daughter2.momentum;
        if (daughter1.charge == UnknownCharge() || daughter2.charge == UnknownCharge())
            charge = UnknownCharge();
        else charge = daughter1.charge + daughter2.charge;
        for(const Candidate& daughter : daughters) {
            if(daughter.finalStateDaughters.size()) {
                for(const Candidate& finalStateDaughter : daughter.finalStateDaughters)
                    finalStateDaughters.push_back(finalStateDaughter);
            } else
                finalStateDaughters.push_back(daughter);
        }
    }

    bool operator < (const Candidate& other) const
    {
        if(index == std::numeric_limits<size_t>::max() || other.index == std::numeric_limits<size_t>::max())
            return momentum.Pt() > other.momentum.Pt();
        return index < other.index;
    }

    bool operator == (const Candidate& other) const
    {
        return !(*this != other);
    }

    bool operator != (const Candidate& other) const
    {
        if(type != other.type || index != other.index || charge != other.charge ||
                daughters.size() != other.daughters.size()) return true;
        for (unsigned n = 0; n < daughters.size(); ++n){
            if (daughters.at(n) != other.daughters.at(n)) return true;
        }
        return false;
    }

    const Candidate& GetDaughter(Type daughterType) const
    {
        for(const Candidate& daughter : daughters) {
            if(daughter.type == daughterType)
                return daughter;
        }
        throw std::runtime_error("daughter with specified type not found.");
    }

    const Candidate& GetLeadingDaughter(Type expectedDaughterType = Unknown) const
    {
        if(!daughters.size())
            throw std::runtime_error("candidate has no daughters");
        if(daughters.size() != 2)
            throw std::runtime_error("candidate has too many daughters");
        const Candidate& leadingDaughter = daughters.at(0).momentum.Pt() > daughters.at(1).momentum.Pt()
                                         ? daughters.at(0) : daughters.at(1);
        if(expectedDaughterType != Unknown && leadingDaughter.type != expectedDaughterType)
            throw std::runtime_error("unexpected leading daughter type");
        return leadingDaughter;
    }

    const Candidate& GetSubleadingDaughter(Type expectedDaughterType = Unknown) const
    {
        if(!daughters.size())
            throw std::runtime_error("candidate has no daughters");
        if(daughters.size() != 2)
            throw std::runtime_error("candidate has too many daughters");
        const Candidate& subleadingDaughter = daughters.at(0).momentum.Pt() > daughters.at(1).momentum.Pt()
                                         ? daughters.at(1) : daughters.at(0);
        if(expectedDaughterType != Unknown && subleadingDaughter.type != expectedDaughterType)
            throw std::runtime_error("unexpected subleading daughter type");
        return subleadingDaughter;
    }

    const particles::ParticleCode& GetPdgId() const
    {
        return typeToPdgMap().at(type);
    }
};

class Vertex{
public:

    TVector3 position;
    size_t index;
    double sumPtSquared;
    unsigned ndf;

    Vertex() : index(0), sumPtSquared(0), ndf(0) {}

    template <typename NtupleObject>
    Vertex(size_t _index, const NtupleObject& ntupleObject) : index(_index), sumPtSquared(ntupleObject.sumPtSquared),
        ndf(ntupleObject.ndf)
    {
       position = TVector3(ntupleObject.x,ntupleObject.y,ntupleObject.z);
    }

    bool operator < (const Vertex& other) const
    {
        return index < other.index;
//        return sumPtSquared < other.sumPtSquared;
    }

};

} // analysis
