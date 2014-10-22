/*!
 * \file FlatTreeAnalyzer_etau.C
 * \brief Analyze flat-tree for etau channel for Htautau analysis.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
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

#include "Analysis/include/BaseFlatTreeAnalyzer.h"

class FlatAnalyzerData_etau : public analysis::FlatAnalyzerData {
public:
    TH1D_ENTRY(mt_1, 50, 0, 50)
};

class FlatTreeAnalyzer_etau : public analysis::BaseFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_etau(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                          const std::string& outputFileName, const std::string& signal_list)
        : BaseFlatTreeAnalyzer(source_cfg, hist_cfg, _inputPath, outputFileName, ChannelId(), signal_list)
    {
    }

protected:
    virtual analysis::Channel ChannelId() const override { return analysis::Channel::ETau; }

    virtual std::shared_ptr<analysis::FlatAnalyzerData> MakeAnaData() override
    {
        return std::shared_ptr<FlatAnalyzerData_etau>(new FlatAnalyzerData_etau());
    }

    virtual analysis::EventRegion DetermineEventRegion(const ntuple::Flat& event) override
    {
        using analysis::EventRegion;
        using namespace cuts::Htautau_Summer13::ETau;

        if(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || (event.mt_1 >= electronID::mt && event.mt_1 <= BackgroundEstimation::HighMtRegion) )
            return EventRegion::Unknown;

        const bool os = event.q_1 * event.q_2 == -1;
        const bool iso = event.pfRelIso_1 < electronID::pFRelIso;
        const bool low_mt = event.mt_1 < electronID::mt;

        if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_HighMt;
        if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_HighMt;
        return os ? EventRegion::OS_NotIsolated : EventRegion::SS_NotIsolated;
    }

    virtual bool FillHistograms(analysis::FlatAnalyzerData& _anaData, const analysis::FlatEventInfo& eventInfo,
                                double weight, bool fillAllHistograms) override
    {
        if(!BaseFlatTreeAnalyzer::FillHistograms(_anaData, eventInfo, weight, fillAllHistograms))
            return false;
        FlatAnalyzerData_etau& anaData = *dynamic_cast<FlatAnalyzerData_etau*>(&_anaData);
        const ntuple::Flat& event = *eventInfo.event;
        anaData.mt_1().Fill(event.mt_1, weight);
        return true;
    }
};
