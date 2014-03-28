/*!
 * \file GenParticle.h
 * \brief Definition of GenParticle class for analysis.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-17 created
 */

#pragma once
#include <set>
#include <boost/shared_ptr.hpp>

#include <TLorentzVector.h>

#include "../include/Particles.h"
#include "../include/iostream_operators.h"
#include "TreeProduction/interface/GenParticle.h"


namespace analysis {

class GenParticle;

typedef std::vector< GenParticle* > GenParticleVector;
typedef std::vector< GenParticleVector > GenParticleVector2D;
typedef std::set< GenParticle* > GenParticleSet;
typedef std::map<size_t, boost::shared_ptr<GenParticle>> GenParticleMap;
typedef std::map<particles::ParticleCode, GenParticleSet> ParticleCodeMap;
typedef std::vector<particles::ParticleCode> ParticleCodes;
typedef std::vector<ParticleCodes> ParticleCodes2D;

typedef std::vector<size_t> GenParticleIndexVector;

class GenParticle {
public:
    size_t index;
    particles::PdgParticle pdg;
    particles::Status status;
    TLorentzVector momentum;
    GenParticleVector mothers;
    GenParticleVector daughters;


public:

    GenParticle(size_t n, const ntuple::GenParticleVector& genParticles)
        : index(n)
    {

        if(n >= genParticles.size())
            throw std::runtime_error("GenParticles index is out of range.");
        const ntuple::GenParticle& ntupleGenParticle = genParticles.at(n);
        pdg = particles::PdgParticle(ntupleGenParticle.PdgId);
        status = particles::NameProvider<particles::Status>::Convert(ntupleGenParticle.Status);
        momentum.SetXYZT(ntupleGenParticle.Px, ntupleGenParticle.Py,ntupleGenParticle.Pz,ntupleGenParticle.E);
    }


    void Initialize(const ntuple::GenParticleVector& genParticles, const GenParticleMap& particles)
    {
        const ntuple::GenParticle& ntupleGenParticle = genParticles.at(index);
        mothers.clear();
        daughters.clear();
        mothers.reserve(ntupleGenParticle.Mother_Indexes.size());
        daughters.reserve(ntupleGenParticle.Daughter_Indexes.size());
        for (unsigned motherIndex : ntupleGenParticle.Mother_Indexes){
            GenParticle* mother = particles.at(motherIndex).get();
            mothers.push_back(mother);
            mother->daughters.push_back(this);
        }
    }

};

class GenEvent {
public:
    GenParticleMap genParticles;
    ParticleCodeMap particleCodeMap;
    GenParticleSet primaryParticles;

public:

    GenEvent(const ntuple::GenParticleVector& ntupleGenParticles)
    {
        for(size_t n = 0; n < ntupleGenParticles.size(); ++n){
            auto particle = boost::shared_ptr<GenParticle>(new GenParticle(n, ntupleGenParticles));
            genParticles[particle->index]=particle;
        }
        for(GenParticleMap::value_type& particle : genParticles){
            particle.second->Initialize(ntupleGenParticles, genParticles);
            if (particle.second->mothers.size() == 0){
                GenParticle* primaryParticle = particle.second.get();
                primaryParticles.insert(primaryParticle);
            }
            if (particle.second->status == particles::Decayed_or_fragmented ||
                    particle.second->status == particles::FinalStateParticle)
                particleCodeMap[particle.second->pdg.Code].insert(particle.second.get());
        }    
    }

    GenParticleSet GetParticles(const ParticleCodes& particleCodes) const
    {
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

    void Print() const
    {
        for (const GenParticle* particle : primaryParticles) {
            PrintChain(particle);
        }
    }

    void PrintChain(const GenParticle* particle, unsigned iteration = 0) const
    {
        const particles::PdgParticle pdgParticle(particle->pdg);
        const particles::Status particleStatus = particles::NameProvider<particles::Status>::Convert(particle->status);
        const TLorentzVector genParticle_momentum = particle->momentum;
        for (unsigned n = 0; n < iteration; ++n)
            std::cout << "  ";
        std::cout << "index=" << particle->index << " name=" << pdgParticle << " status=" << particleStatus
                  <<  " pt= " << genParticle_momentum.Pt() <<"\n";
        for(unsigned n = 0; n < particle->daughters.size(); ++n) {
            const GenParticle* daughter = particle->daughters.at(n);
                PrintChain(daughter,iteration+1);
        }
    }
};

inline bool FindDecayProducts(const GenParticle& genParticle, const ParticleCodes& particleCodes,
                              GenParticleVector& decayProducts)
{
    if (genParticle.status != particles::Decayed_or_fragmented)
        throw std::runtime_error("particle type not supported");

    decayProducts.clear();
    const GenParticleVector* daughters = &genParticle.daughters;
    size_t expected_nDaughters = particleCodes.size();
    if (daughters->size() == 0){
        if(genParticle.mothers.size() != 1)
            throw std::runtime_error("cannot find decay products");
        const GenParticle& mother = *genParticle.mothers.front();
        if(mother.status != particles::HardInteractionProduct || mother.pdg != genParticle.pdg)
            throw std::runtime_error("cannot find decay products");
        daughters = &mother.daughters;
        ++expected_nDaughters;
    }
    if (daughters->size() != expected_nDaughters)
        return false;
    std::set<size_t> taken_daughters;
    for (const particles::ParticleCode& code : particleCodes){
        bool daughter_found = false;
        for (size_t n = 0; n < daughters->size(); ++n){
            if (taken_daughters.count(n))
                continue;
            if (code != daughters->at(n)->pdg.Code)
                continue;
            GenParticle* daughter = daughters->at(n);
            if (daughter->status == particles::HardInteractionProduct){
                bool grandDaughter_found = false;
                for (GenParticle* grandDaughter : daughter->daughters){
                    if (grandDaughter->pdg == daughter->pdg &&
                            grandDaughter->status == particles::Decayed_or_fragmented){
                        grandDaughter_found = true;
                        daughter = grandDaughter;
                        break;
                    }
                }
                if (!grandDaughter_found)
                    throw std::runtime_error("cannot find decay products");
            }

            decayProducts.push_back(daughter);
            taken_daughters.insert(n);
            daughter_found = true;
            break;
        }
        if (!daughter_found) return false;
    }
    return true;
}

inline bool FindDecayProducts2D(const GenParticleVector& genParticles, const ParticleCodes2D& particleCodes2D,
                                GenParticleVector2D& decayProducts2D, GenParticleIndexVector& indexes)
{
    std::set<size_t> taken_genParticles;
    if (genParticles.size() != particleCodes2D.size())
        throw std::runtime_error("mismatched vector size of particles");
    decayProducts2D.clear();
    for(const ParticleCodes& codes : particleCodes2D){
        bool particleFound = false;
        for (size_t n = 0; n < genParticles.size(); ++n ){
            if (taken_genParticles.count(n)) continue;
            const GenParticle& genParticle = *genParticles.at(n);
            GenParticleVector decayProducts;
            if (!FindDecayProducts(genParticle,codes,decayProducts)) continue;
            particleFound = true;
            decayProducts2D.push_back(decayProducts);
            taken_genParticles.insert(n);
            indexes.push_back(n);
        }
        if (!particleFound) return false;
    }
    return true;
}

} // analysis
