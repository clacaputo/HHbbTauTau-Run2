/*!
 * \file AnalysisTypes.h
 * \brief Common simple types for analysis purposes.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
 * \date 2014-10-09 created
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

#include "PhysicalValue.h"

namespace analysis {

enum class Channel { ETau = 0, MuTau = 1, TauTau = 2 };

enum class EventEnergyScale { Central = 0, TauUp = 1, TauDown = 2, JetUp = 3, JetDown = 4, BtagEfficiencyUp = 5,
                              BtagEfficiencyDown = 6 , BtagFakeUp = 7, BtagFakeDown = 8};

namespace detail {
const std::map<Channel, std::string> ChannelNameMap = {
    { Channel::ETau, "eTau" }, { Channel::MuTau, "muTau" }, { Channel::TauTau, "tauTau" }
};

const std::map<Channel, std::string> ChannelNameMapLatex = {
    { Channel::ETau, "e#tau_{h}" }, { Channel::MuTau, "#mu#tau_{h}" }, { Channel::TauTau, "#tau_{h}#tau_{h}" }
};

const std::map<EventEnergyScale, std::string> EventEnergyScaleNameMap = {
    { EventEnergyScale::Central, "Central" },
    { EventEnergyScale::TauUp, "TauUp" }, { EventEnergyScale::TauDown, "TauDown" },
    { EventEnergyScale::JetUp, "JetUp" }, { EventEnergyScale::JetDown, "JetDown" },
    { EventEnergyScale::BtagEfficiencyUp, "BtagEfficiencyUp" },
    { EventEnergyScale::BtagEfficiencyDown, "BtagEfficiencyDown" },
    { EventEnergyScale::BtagFakeUp, "BtagFakeUp" },
    { EventEnergyScale::BtagFakeDown, "BtagFakeDown" }
};

} // namespace detail

typedef std::set<EventEnergyScale> EventEnergyScaleSet;

const std::set<EventEnergyScale> AllEventEnergyScales = tools::collect_map_keys(detail::EventEnergyScaleNameMap);

std::ostream& operator<< (std::ostream& s, const Channel& c)
{
    s << detail::ChannelNameMap.at(c);
    return s;
}

std::istream& operator>> (std::istream& s, Channel& c)
{
    std::string name;
    s >> name;
    for(const auto& map_entry : detail::ChannelNameMap) {
        if(map_entry.second == name) {
            c = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown channel name '") << name << "'.";
}

std::ostream& operator<< (std::ostream& s, const EventEnergyScale& es)
{
    s << detail::EventEnergyScaleNameMap.at(es);
    return s;
}

std::istream& operator>> (std::istream& s, EventEnergyScale& es)
{
    std::string name;
    s >> name;
    for(const auto& map_entry : detail::EventEnergyScaleNameMap) {
        if(map_entry.second == name) {
            es = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown event energy scale '") << name << "'.";
}


} // namespace analysis
