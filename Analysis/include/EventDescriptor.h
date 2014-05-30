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
#include "TreeProduction/interface/TriggerObject.h"
#include "TreeProduction/interface/MET.h"
#include "TreeProduction/interface/GenMET.h"
#include "TreeProduction/interface/PFCand.h"
#include "EventId.h"

namespace analysis{

namespace detail {
typedef std::tuple< ntuple::Event,
                    ntuple::ElectronVector,
                    ntuple::MuonVector,
                    ntuple::TauVector,
                    ntuple::PFCandVector,
                    ntuple::JetVector,
                    ntuple::VertexVector,
                    ntuple::GenParticleVector,
                    ntuple::TriggerVector,
                    ntuple::TriggerObjectVector,
                    ntuple::METVector,
                    ntuple::METVector,
                    ntuple::METVector,
//                    ntuple::METVector,
//                    ntuple::METVector,
//                    ntuple::METVector,
                    ntuple::GenMETVector > EventTuple;

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
    const ntuple::PFCandVector& pfCandidates() const { return std::get<4>(_data); }
    const ntuple::JetVector& jets() const { return std::get<5>(_data); }
    const ntuple::VertexVector& vertices() const { return std::get<6>(_data); }
    const ntuple::GenParticleVector& genParticles() const { return std::get<7>(_data); }
    const ntuple::TriggerVector& triggers() const { return std::get<8>(_data); }
    const ntuple::TriggerObjectVector& triggerObjects() const { return std::get<9>(_data); }
    const ntuple::MET& met() const {
        const ntuple::METVector& metVector = std::get<10>(_data);
        if (metVector.size() == 1)
            return metVector.at(0);
        throw std::runtime_error("more than 1 MET in the event");
    }
    const ntuple::MET& metPF() const {
        const ntuple::METVector& metVector = std::get<11>(_data);
        if (metVector.size() == 1)
            return metVector.at(0);
        throw std::runtime_error("more than 1 MET PF in the event");
    }
    const ntuple::MET& metTC() const {
        const ntuple::METVector& metVector = std::get<12>(_data);
        if (metVector.size() == 1)
            return metVector.at(0);
        throw std::runtime_error("more than 1 MET TC in the event");
    }
    const ntuple::GenMET& genMet() const {
        const ntuple::GenMETVector& genMetVector = std::get<13>(_data);
        if (genMetVector.size() == 1)
            return genMetVector.at(0);
        throw std::runtime_error("more than 1 gen MET in the event");
    }
//    const ntuple::MET& metMVAeTau() const {
//        const ntuple::METVector& metVector = std::get<14>(_data);
//        if (metVector.size() == 1)
//            return metVector.at(0);
//        throw std::runtime_error("more than 1 MET MVA in the event");
//    }
//    const ntuple::MET& metMVAmuTau() const {
//        const ntuple::METVector& metVector = std::get<15>(_data);
//        if (metVector.size() == 1)
//            return metVector.at(0);
//        throw std::runtime_error("more than 1 MET MVA in the event");
//    }
//    const ntuple::MET& metMVAtauTau() const {
//        const ntuple::METVector& metVector = std::get<16>(_data);
//        if (metVector.size() == 1)
//            return metVector.at(0);
//        throw std::runtime_error("more than 1 MET MVA in the event");
//    }

    detail::EventTuple& data() { return _data; }
    void Clear() { detail::ClearEventTuple(_data); }

private:
    detail::EventTuple _data;
};

}
