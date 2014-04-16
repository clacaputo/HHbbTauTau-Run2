/*!
 * \file EventId.h
 * \brief Definition of EventId class which represent unique CMS event identifier.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-04-16 created
 */

#pragma once

#include <limits>

namespace analysis {

struct EventId{
    unsigned runId;
    unsigned lumiBlock;
    unsigned eventId;
    static const EventId& Undef_event() {
        static const EventId undef_event;
        return undef_event;
    }

    EventId() : runId(std::numeric_limits<UInt_t>::max()), lumiBlock(std::numeric_limits<UInt_t>::max()),
                eventId(std::numeric_limits<UInt_t>::max()){}

    EventId(unsigned _runId, unsigned _lumiBlock, unsigned _eventId) : runId(_runId), lumiBlock(_lumiBlock),
                eventId(_eventId){}

    bool operator == (const EventId& other) const
    {
        return runId == other.runId && lumiBlock == other.lumiBlock && eventId == other.eventId;
    }

    bool operator != (const EventId& other) const
    {
        return runId != other.runId || lumiBlock != other.lumiBlock || eventId != other.eventId;
    }
};

} // analysis
