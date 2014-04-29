/*!
 * \file TriggerTools.h
 * \brief Definiton of tools to work with embedded trigger information.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-04-29 created
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
