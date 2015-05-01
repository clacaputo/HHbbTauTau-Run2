/*!
 * \file EventId.h
 * \brief Definition of EventId class which represent unique CMS event identifier.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
 * \date 2014-04-16 created
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

#include <limits>
#include <iostream>

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
        return !(*this != other);
    }

    bool operator != (const EventId& other) const
    {
        return runId != other.runId || lumiBlock != other.lumiBlock || eventId != other.eventId;
    }

    bool operator < (const EventId& other) const
    {
        if(runId < other.runId) return true;
        if(runId > other.runId) return false;
        if(lumiBlock < other.lumiBlock) return true;
        if(lumiBlock > other.lumiBlock) return false;
        return eventId < other.eventId;
    }
};

} // analysis

inline std::ostream& operator <<(std::ostream& s, const analysis::EventId& event)
{
    s << "run = " << event.runId << ", lumi = " << event.lumiBlock << ", evt = " << event.eventId;
    return s;
}
