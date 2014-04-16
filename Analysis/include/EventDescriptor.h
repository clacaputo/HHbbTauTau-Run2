/*!
 * \file EventDescriptor.h
 * \brief Definition of EventDescriptor class that contains all ntuple data for one event.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-28 created
 */

#pragma once

#include <tuple>
#include "TreeProduction/interface/Event.h"
#include "TreeProduction/interface/Electron.h"
#include "TreeProduction/interface/Muon.h"
#include "TreeProduction/interface/Tau.h"
#include "TreeProduction/interface/Jet.h"
#include "TreeProduction/interface/Vertex.h"
#include "TreeProduction/interface/GenParticle.h"
#include "TreeProduction/interface/Trigger.h"
#include "EventId.h"

namespace analysis{

namespace detail {
typedef std::tuple< ntuple::Event,
                    ntuple::ElectronVector,
                    ntuple::MuonVector,
                    ntuple::TauVector,
                    ntuple::JetVector,
                    ntuple::VertexVector,
                    ntuple::GenParticleVector,
                    ntuple::TriggerVector> EventTuple;

template<typename ObjectType>
void ClearContainer(ObjectType& container) {}

template<typename ObjectType>
void ClearContainer(std::vector<ObjectType>& container) { container.clear(); }

template<size_t N = 0>
inline typename std::enable_if< N == std::tuple_size<EventTuple>::value >::type
ClearEventTuple(EventTuple& data) {}

template<size_t N = 0>
inline typename std::enable_if< (N < std::tuple_size<EventTuple>::value) >::type
ClearEventTuple(EventTuple& data)
{
    ClearContainer(std::get<N>(data));
    ClearEventTuple<N + 1>(data);
}

} // detail

class EventDescriptor{
public:
    EventId eventId() const { return EventId(eventInfo().run, eventInfo().lumis, eventInfo().EventId); }
    const ntuple::Event& eventInfo() const { return std::get<0>(_data); }
    const ntuple::ElectronVector& electrons() const { return std::get<1>(_data); }
    const ntuple::MuonVector& muons() const { return std::get<2>(_data); }
    const ntuple::TauVector& taus() const { return std::get<3>(_data); }
    const ntuple::JetVector& jets() const { return std::get<4>(_data); }
    const ntuple::VertexVector& vertices() const { return std::get<5>(_data); }
    const ntuple::GenParticleVector& genParticles() const { return std::get<6>(_data); }
    const ntuple::TriggerVector& triggers() const { return std::get<7>(_data); }

    detail::EventTuple& data() { return _data; }
    void Clear() { detail::ClearEventTuple(_data); }

private:
    detail::EventTuple _data;
};

}
