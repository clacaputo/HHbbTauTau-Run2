/*!
 * \file FlatTreeAnalyzer_etau.C
 * \brief Analyze flat-tree for etau channel for Htautau analysis.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
 * \author Rosamaria Venditti (INFN Bari, Bari University)
 * \author Claudio Caputo (INFN Bari, Bari University)
 * \date 2014-09-03 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>,
 *                Rosamaria Venditti,
 *                Claudio Caputo
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

#include "Analysis/include/SemileptonicFlatTreeAnalyzer.h"

class FlatTreeAnalyzer_etau : public analysis::SemileptonicFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_etau(const std::string& source_cfg, const std::string& _inputPath,
                          const std::string& outputFileName, const std::string& signal_list,
                          bool applyPostFitCorrections = false, bool saveFullOutput = false)
        : SemileptonicFlatTreeAnalyzer(analysis::DataCategoryCollection(source_cfg, signal_list, ChannelId()),
                                       _inputPath, outputFileName, applyPostFitCorrections, saveFullOutput)
    {
    }

protected:
    virtual analysis::Channel ChannelId() const override { return analysis::Channel::ETau; }

    virtual analysis::EventRegion DetermineEventRegion(const ntuple::Flat& event,
                                                       analysis::EventCategory eventCategory) override
    {
        using analysis::EventRegion;
        using namespace cuts::Htautau_Summer13::ETau;

        if(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || (event.pfRelIso_1 >= electronID::pFRelIso && !IsAntiIsolatedRegion(event))
                || (event.mt_1 >= electronID::mt && !IsHighMtRegion(event,eventCategory)) /*|| event.pt_2 <= 30*/ )
            return EventRegion::Unknown;

        const bool os = event.q_1 * event.q_2 == -1;
        const bool iso = event.pfRelIso_1 < electronID::pFRelIso;
        const bool low_mt = event.mt_1 < electronID::mt;


        if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_Iso_HighMt;
        if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_Iso_HighMt;
        if(os) return low_mt ? EventRegion::OS_AntiIsolated : EventRegion::OS_AntiIso_HighMt;
        return low_mt ? EventRegion::SS_AntiIsolated : EventRegion::SS_AntiIso_HighMt;
    }
};
