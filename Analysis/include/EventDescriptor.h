/*!
 * \file EventDescriptor.h
 * \brief Definition of EventDescriptor class that contains all ntuple data for one event.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-28 created
 */

#pragma once

#include "TreeProduction/interface/Event.h"
#include "TreeProduction/interface/Electron.h"
#include "TreeProduction/interface/Muon.h"
#include "TreeProduction/interface/Tau.h"
#include "TreeProduction/interface/Jet.h"
#include "TreeProduction/interface/Vertex.h"
#include "TreeProduction/interface/GenParticle.h"

namespace analysis{

class EventDescriptor{
public:
    ntuple::Event eventInfo;
    ntuple::ElectronVector electrons;
    ntuple::MuonVector muons;
    ntuple::TauVector taus;
    ntuple::JetVector jets;
    ntuple::VertexVector vertices;
    ntuple::GenParticleVector genParticles;

    void Clear(){
        electrons.clear();
        muons.clear();
        taus.clear();
        jets.clear();
        vertices.clear();
        genParticles.clear();
    }

    unsigned GetEventId() const{
        return eventInfo.EventId;
    }


};

}
