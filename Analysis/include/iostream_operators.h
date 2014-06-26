/*!
 * \file iostream_operators.h
 * \brief Definition of iostream operators suitable for tau trigger analysis.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2013-05-10 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
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
    s << /*"(" <<*/ v.Pt() << " " << v.Eta() << " " << v.Phi() << " " << v.E() /*<< ")"*/;
    return s;
}



