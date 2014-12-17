/*!
 * \file EventSelectionStudy.C
 * \brief Study of possible issues in event selection.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-11-05 created
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

#include "Analysis/include/LightBaseFlatTreeAnalyzer.h"

class EventSelectionStudy : public analysis::LightBaseFlatTreeAnalyzer {
public:
    EventSelectionStudy(const std::string& _inputFileName, const std::string& _outputFileName)
         : LightBaseFlatTreeAnalyzer(_inputFileName,_outputFileName)
    {
        recalc_kinfit = false;
        GetOutputFile().cd();
    }

protected:

    virtual void AnalyzeEvent(const analysis::FlatEventInfo& eventInfo, analysis::EventCategory category) override
    {
        using analysis::EventCategory;
        if(DetermineEventRegion(eventInfo,category) != analysis::EventRegion::OS_Isolated) return;
        //if (!analysis::TwoJetsEventCategories_MediumBjets.count(category)) return;

        if(eventInfo.eventEnergyScale == analysis::EventEnergyScale::Central &&
                category == EventCategory::Inclusive &&
                eventInfo.eventType == ntuple::EventType::ZJ )
            std::cout << eventInfo.eventType << " " << eventInfo.event->evt << std::endl;
    }

};
