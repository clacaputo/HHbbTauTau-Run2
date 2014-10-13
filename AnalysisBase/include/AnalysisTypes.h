/*!
 * \file AnalysisMath.h
 * \brief Common simple types for analysis purposes.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
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

#include <iostream>
#include <map>

#include "exception.h"

namespace analysis {

enum class Channel { ETau = 0, MuTau = 1, TauTau = 2 };

namespace detail {
std::map<Channel, std::string> ChannelNameMap = {
    { Channel::ETau, "eTau" }, { Channel::MuTau, "muTau" }, { Channel::TauTau, "tauTau" }
};
} // namespace detail

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


} // namespace analysis
