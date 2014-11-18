/*!
 * \file FlatTreeAnalyzer_mutau.C
 * \brief Analyze flat-tree for mu-tau channel for HHbbtautau analysis.
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

class FlatAnalyzerData_mutau : public analysis::FlatAnalyzerData {
public:
    TH1D_ENTRY(mt_1, 50, 0, 50)

    virtual void Fill(const analysis::FlatEventInfo& eventInfo, double weight, bool fill_all,
                      analysis::EventEnergyScale eventEnergyScale) override
    {
        FlatAnalyzerData::Fill(eventInfo, weight, fill_all, eventEnergyScale);
        if(!fill_all) return;
        if (eventEnergyScale != analysis::EventEnergyScale::Central) return;
        const ntuple::Flat& event = *eventInfo.event;
        mt_1().Fill(event.mt_1, weight);
    }
};

class FlatTreeAnalyzer_mutau : public analysis::BaseFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_mutau(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                           const std::string& outputFileName, const std::string& signal_list)
         : BaseFlatTreeAnalyzer(source_cfg, hist_cfg, _inputPath, outputFileName, ChannelId(), signal_list)
    {
    }

protected:
    virtual analysis::Channel ChannelId() const override { return analysis::Channel::MuTau; }

    virtual std::shared_ptr<analysis::FlatAnalyzerData> MakeAnaData() override
    {
        return std::shared_ptr<FlatAnalyzerData_mutau>(new FlatAnalyzerData_mutau());
    }

    virtual analysis::EventRegion DetermineEventRegion(const ntuple::Flat& event) override
    {
        using analysis::EventRegion;
        using namespace cuts::Htautau_Summer13::MuTau;

        if(!event.againstMuonTight_2
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || (event.mt_1 >= muonID::mt && event.mt_1 <= BackgroundEstimation::HighMtRegion) /*|| event.pt_2 <= 30*/ )
            return EventRegion::Unknown;

        const bool os = event.q_1 * event.q_2 == -1;
        const bool iso = event.pfRelIso_1 < muonID::pFRelIso;
        const bool low_mt = event.mt_1 < muonID::mt;

        if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_Iso_HighMt;
        else if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_Iso_HighMt;
        else if(!iso && os) return low_mt ? EventRegion::OS_NotIsolated : EventRegion::OS_NotIso_HighMt;
        else return low_mt ? EventRegion::SS_NotIsolated : EventRegion::SS_NotIso_HighMt;
    }

    virtual bool PassMvaCut(const analysis::FlatEventInfo& eventInfo, analysis::EventCategory eventCategory) override
    {
        static const std::map<analysis::EventCategory, double> mva_BDT_cuts = {
            { analysis::EventCategory::TwoJets_ZeroBtag, -0.05 },
            { analysis::EventCategory::TwoJets_OneBtag, -0.05 },
            { analysis::EventCategory::TwoJets_TwoBtag, -0.05 }
        };
        if(!eventInfo.has_bjet_pair || !mva_BDT_cuts.count(eventCategory))
            return true;

        return eventInfo.mva_BDT > mva_BDT_cuts.at(eventCategory);
    }

    virtual analysis::PhysicalValue CalculateQCDScaleFactor(analysis::EventCategory /*eventCategory*/,
                                                            const std::string& /*hist_name*/) override
    {
        static const analysis::PhysicalValue sf(1.06, 0.001);
        return sf;
    }
};
