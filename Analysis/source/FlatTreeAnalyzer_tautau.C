/*!
 * \file FlatTreeAnalyzer_tautau.C
 * \brief Analyze flat-tree for Htautau analysis.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-09-03 created
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

#include "../include/BaseFlatTreeAnalyzer.h"

class FlatTreeAnalyzer_tautau : public analysis::BaseFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_tautau(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                         const std::string& outputFileName, const std::string& _signalName,
                         const std::string& _dataName)
        : BaseFlatTreeAnalyzer(source_cfg, hist_cfg, _inputPath, outputFileName, _signalName, _dataName)
    {
    }

protected:
    virtual EventType_QCD DetermineEventTypeForQCD(const ntuple::Flat& event) override
    {
        return EventType_QCD::Unknown;
    }

    virtual EventType_Wjets DetermineEventTypeForWjets(const ntuple::Flat& event) override
    {
        return EventType_Wjets::Unknown;
    }

    virtual EventCategory DetermineEventCategory(const ntuple::Flat& event) override
    {
        return EventCategory::TwoJets_OneBtag;
    }

    virtual void EstimateQCD() override
    {

    }

    virtual void EstimateWjets() override
    {

    }

};
