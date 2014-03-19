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

class Candidate {
  public:

    enum Type{
        Unknown, Mu, Electron, Tau, Bjet, Higgs, Resonance
    };

    Type type;
    int index;
    TLorentzVector momentum;
    CandidatePtrVector daughters;


    Candidate(): type(Unknown), index(-1){}
    Candidate(Type _type, int _index, const TLorentzVector& _momentum):
        type(_type), index(_index), momentum(_momentum) {}

    Candidate(Type _type, const Candidate& daughter1, const Candidate& daughter2) : type(_type), index(-1){
        daughters.push_back(&daughter1);
        daughters.push_back(&daughter2);
        momentum = daughter1.momentum + daughter2.momentum;
    }

    bool operator < (const Candidate& other) const
    {
        return momentum.Pt() < other.momentum.Pt();
    }


};




}
