/*!
 * \file iostream_operators.h
 * \brief Definition of iostream operators suitable for tau trigger analysis.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2013-05-10 created
 */

#pragma once

#include "TLorentzVector.h"
#include "Particles.h"


std::ostream& operator<<(std::ostream& s, const particles::Status& particleStatus){
    s << particles::NameProvider<particles::Status>::Name(particleStatus);
    return s;
}

std::ostream& operator<<(std::ostream& s, const particles::ParticleCode& code){
    s << code.Name();
    return s;
}

std::ostream& operator<<(std::ostream& s, const particles::PdgParticle& particle){
    s << particle.Name();
    return s;
}

std::ostream& operator<<(std::ostream& s, const TVector3& v){
    s << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return s;
}

std::ostream& operator<<(std::ostream& s, const TLorentzVector& v){
    s << "(" << v.E() << ", " << v.Px() << ", " << v.Py() << ", " << v.Pz() << ")";
    return s;
}



