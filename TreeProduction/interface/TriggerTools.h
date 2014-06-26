/*!
 * \file TriggerTools.h
 * \brief Definiton of tools to work with embedded trigger information.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-04-29 created
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

#include "DataFormats/PatCandidates/interface/PATObject.h"

template<typename PatObject>
std::vector<std::string> CollectMatchedTriggerPaths(const PatObject& patObject)
{
    const pat::TriggerObjectStandAloneCollection& matchedTriggers = patObject.triggerObjectMatches();
    std::set<std::string> pathNames;
    for(const pat::TriggerObjectStandAlone& triggerObject : matchedTriggers) {
        const auto objectPathNames = triggerObject.pathNames();
        for(const std::string& pathName : objectPathNames)
            pathNames.insert(pathName);
    }
    return std::vector<std::string>(pathNames.begin(), pathNames.end());
}
