/*!
 * \file GenParticle.h
 * \brief Definition of GenParticle class for analysis.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-17 created
 */

#pragma once

#include <boost/shared_ptr.hpp>

#include <TLorentzVector.h>

#include "../include/Particles.h"

namespace analysis {

class GenParticle;

typedef std::vector< boost::shared_ptr<GenParticle> > GenParticleVector;
typedef std::map<size_t, boost::shared_ptr<GenParticle> > GenParticleMap;

class GenParticle {
public:
    size_t index;
    particles::PdgParticle pdg;
    particles::Status status;
    TLorentzVector momentum;
    TVector3 direction;
    GenParticleVector mothers;
    GenParticleVector daughters;

public:
    template<typename EventTree>
    GenParticle(size_t _index, const EventTree& event)
        : index(_index)
    {
        if(event.nGenParticle < 0 || index >= static_cast<size_t>(event.nGenParticle))
            throw std::runtime_error("GenParticles index is out of range.");
        pdg = particles::PdgParticle(event.GenParticle_pdgId[index]);
        status = particles::NameProvider<particles::Status>::Convert(event.GenParticle_status[index]);
        momentum.SetXYZT(event.GenParticle_px[index], event.GenParticle_py[index], event.GenParticle_pz[index],
                         event.GenParticle_energy[index]);
        direction.SetPtEtaPhi(event.GenParticle_pt[index], event.GenParticle_eta[index], event.GenParticle_phi[index]);
    }

    template<typename EventTree>
    void Initialize(const EventTree& event, const GenParticleVector& particles)
    {
        mothers.clear();
        daughters.clear();
        const size_t n_mothers = event.GenParticle_motherIndices[index].size();
        const size_t n_daughters = event.GenParticle_daughtIndices[index].size();
        mothers.reserve(n_mothers);
        daughters.reserve(n_daughters);

        if( event.GenParticle_numMother[index] < 0
                || static_cast<size_t>(event.GenParticle_numMother[index]) != n_mothers ) {
            std::cerr << "index=" << index << "num=" << event.GenParticle_numMother[index]
                         << ",n=" << n_mothers << std::endl;

            throw std::runtime_error("Inconsistent GenParticles mother relations in the input event.");
        }
//        if( event.GenParticle_numDaught[index] < 0
//                || static_cast<size_t>(event.GenParticle_numDaught[index]) != n_daughters ) {
//            std::cerr << "index=" << index << "num=" << event.GenParticle_numDaught[index]
//                         << ",n=" << n_daughters << std::endl;
//            throw std::runtime_error("Inconsistent GenParticles daughter relations in the input event.");
//        }

        for(size_t n = 0; n < n_mothers; ++n)
            mothers.push_back(particles.at(event.GenParticle_motherIndices[index].at(n)));
//        for(size_t n = 0; n < n_daughters; ++n)
//            daughters.push_back(particles.at(event.GenParticle_daughtIndices[index].at(n)));
    }

};

class GenEvent {
public:
    GenParticleVector particles;

public:
    template<typename EventTree>
    GenEvent(const EventTree& event)
    {
        if(event.nGenParticle < 0)
            throw std::runtime_error("Number of GenParticles in the event is less than 0.");
        const size_t n_genParticles = static_cast<size_t>(event.nGenParticle);
        std::cout << "N_gen=" << n_genParticles << std::endl;
        particles.reserve(n_genParticles);
        for(size_t n = 0; n < n_genParticles; ++n)
            particles.push_back(boost::shared_ptr<GenParticle>(new GenParticle(n, event)));
        for(auto& particle : particles)
            particle->Initialize(event, particles);
    }
};

} // namespace analysis
