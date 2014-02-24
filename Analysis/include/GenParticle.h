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

typedef std::vector< GenParticle* > GenParticleVector;
typedef std::set< GenParticle* > GenParticleSet;
typedef std::map<size_t, boost::shared_ptr<GenParticle>> GenParticleMap;
typedef std::map<particles::ParticleCode, GenParticleSet> ParticleCodeMap;
typedef std::vector<particles::ParticleCode> ParticleCodes;

class GenParticle {
public:
    size_t GenParticle_index,index;
    particles::PdgParticle pdg;
    particles::Status status;
    TLorentzVector momentum;
    GenParticleVector mothers;
    GenParticleVector daughters;
    bool MissingMother;
    bool MissingDaughter;

public:
    template<typename EventTree>
    GenParticle(size_t n, const EventTree& event)
        : GenParticle_index(n), index(event.GenParticle_index[n]), MissingMother(false), MissingDaughter(false)
    {

        if(event.nGenParticle < 0 || GenParticle_index >= static_cast<size_t>(event.nGenParticle))
            throw std::runtime_error("GenParticles index is out of range.");
        pdg = particles::PdgParticle(event.GenParticle_pdgId[GenParticle_index]);
        status = particles::NameProvider<particles::Status>::Convert(event.GenParticle_status[GenParticle_index]);
        momentum.SetXYZT(event.GenParticle_px[GenParticle_index], event.GenParticle_py[GenParticle_index],
                         event.GenParticle_pz[GenParticle_index],event.GenParticle_energy[GenParticle_index]);
    }

    template<typename EventTree>
    void Initialize(const EventTree& event, const GenParticleMap& particles)
    {
        mothers.clear();
        daughters.clear();
        const size_t n_mothers = event.GenParticle_numMother[GenParticle_index];
        const size_t n_daughters = event.GenParticle_numDaught[GenParticle_index];
        mothers.reserve(n_mothers);
        daughters.reserve(n_daughters);


        bool noMother1 = false;
        bool noMother2 = false;

        if (event.GenParticle_motherIndex_1[GenParticle_index] != std::numeric_limits<unsigned>::max()){
            GenParticleMap::const_iterator particle_iter_1 = particles.find(event.GenParticle_motherIndex_1[GenParticle_index]);
            if (particle_iter_1 != particles.end()){
                GenParticle* mother = particle_iter_1->second.get();
                mothers.push_back(mother);
                //GenParticle* daughter = this;
                mother->daughters.push_back(this);
            }
            else noMother1 = true;
        }

        if (event.GenParticle_motherIndex_2[GenParticle_index] != std::numeric_limits<unsigned>::max()){
            GenParticleMap::const_iterator particle_iter_2 = particles.find(event.GenParticle_motherIndex_2[GenParticle_index]);
            if (particle_iter_2 != particles.end()){
                GenParticle* mother = particle_iter_2->second.get();
                mothers.push_back(mother);
                //GenParticle* daughter = this;
                mother->daughters.push_back(this);
            }
            else noMother2 = true;
        }

        if (noMother1 || noMother2) MissingMother = true;
    }

};

class GenEvent {
public:
    GenParticleMap genParticles;
    ParticleCodeMap particleCodeMap;
    GenParticleSet primaryParticles;

public:
    template<typename EventTree>
    GenEvent(const EventTree& event)
    {
        if(event.nGenParticle < 0)
            throw std::runtime_error("Number of GenParticles in the event is less than 0.");
        const size_t n_genParticles = static_cast<size_t>(event.nGenParticle);
        //std::cout << "N_gen=" << n_genParticles << std::endl;
        for(size_t n = 0; n < n_genParticles; ++n){
            auto particle = boost::shared_ptr<GenParticle>(new GenParticle(n, event));
            genParticles[particle->index]=particle;
        }
        for(GenParticleMap::value_type& particle : genParticles){
            particle.second->Initialize(event, genParticles);
            if (!particle.second->MissingMother && particle.second->mothers.size() !=
                    event.GenParticle_numMother[particle.second->GenParticle_index])
                throw std::runtime_error("inconsistent number of mothers");
            if (particle.second->mothers.size() == 0 && !particle.second->MissingMother){
                GenParticle* primaryParticle = particle.second.get();
                primaryParticles.insert(primaryParticle);
            }
        }
        //std::cout << "primary particles = " << primaryParticles.size() << std::endl;

        for(GenParticleMap::value_type& particle : genParticles){
            if (particle.second->daughters.size() !=
                    event.GenParticle_numDaught[particle.second->GenParticle_index])
                particle.second->MissingDaughter = true;
            if (particle.second->status == particles::FinalStateParticle ||
                    particle.second->status == particles::Decayed_or_fragmented)
                particleCodeMap[particle.second->pdg.Code].insert(particle.second.get());
        }

    }

    GenParticleSet GetParticles(const ParticleCodes& particleCodes) const {
        GenParticleSet results;
        for (const particles::ParticleCode& code : particleCodes){
            const ParticleCodeMap::const_iterator code_iter = particleCodeMap.find(code);
            if (code_iter == particleCodeMap.end())
                continue;
            for (GenParticle* particle : code_iter->second){
                results.insert(particle);
            }
        }
        return results;
    }
};

} // namespace analysis
