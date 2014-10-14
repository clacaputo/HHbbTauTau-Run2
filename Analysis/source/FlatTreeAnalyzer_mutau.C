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

class FlatTreeAnalyzer_mutau : public analysis::BaseFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_mutau(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                           const std::string& outputFileName, const std::string& signal_list, bool _WjetsData = false)
         : BaseFlatTreeAnalyzer(source_cfg, hist_cfg, _inputPath, outputFileName, ChannelId(), signal_list, _WjetsData)
    {
    }

protected:
    virtual analysis::Channel ChannelId() const override { return analysis::Channel::MuTau; }

    virtual analysis::EventRegion DetermineEventRegion(const ntuple::Flat& event) override
    {
        using analysis::EventRegion;
        using namespace cuts::Htautau_Summer13::MuTau;

        if(!event.againstMuonTight_2
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || (event.mt_1 >= muonID::mt && event.mt_1 <= BackgroundEstimation::HighMtRegion) )
            return EventRegion::Unknown;

        const bool os = event.q_1 * event.q_2 == -1;
        const bool iso = event.pfRelIso_1 < muonID::pFRelIso;
        const bool low_mt = event.mt_1 < muonID::mt;

        if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_HighMt;
        if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_HighMt;
        return os ? EventRegion::OS_NotIsolated : EventRegion::SS_NotIsolated;
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

    virtual double CalculateQCDScaleFactor(analysis::EventCategory eventCategory, const std::string& hist_name) override
    {
        using analysis::EventRegion;
        using analysis::DataCategoryType;

        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        auto hist_OSnotIso_data = GetHistogram(eventCategory, data.name, EventRegion::OS_NotIsolated, hist_name);
        auto hist_SSnotIso_data = GetHistogram(eventCategory, data.name, EventRegion::SS_NotIsolated, hist_name);
        if(!hist_OSnotIso_data || !hist_SSnotIso_data)
            throw analysis::exception("Unable to find histograms for QCD scale factor estimation");

        TH1D& hist_OSnotIso = CloneHistogram(eventCategory, qcd.name, EventRegion::OS_NotIsolated, *hist_OSnotIso_data);
        TH1D& hist_SSnotIso = CloneHistogram(eventCategory, qcd.name, EventRegion::SS_NotIsolated, *hist_SSnotIso_data);

        SubtractBackgroundHistograms(hist_OSnotIso, eventCategory, EventRegion::OS_NotIsolated, qcd.name);
        SubtractBackgroundHistograms(hist_SSnotIso, eventCategory, EventRegion::SS_NotIsolated, qcd.name);

        const double n_OSnotIso = hist_OSnotIso.Integral(0, hist_OSnotIso.GetNbinsX() + 1);
        const double n_SSnotIso = hist_SSnotIso.Integral(0, hist_SSnotIso.GetNbinsX() + 1);
        return n_OSnotIso / n_SSnotIso;
    }

    virtual void EstimateQCD(analysis::EventCategory eventCategory, const std::string& hist_name,
                             double scale_factor) override
    {
        using analysis::EventRegion;
        using analysis::DataCategoryType;

        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        auto hist_SSIso_data = GetHistogram(eventCategory, data.name, EventRegion::SS_Isolated, hist_name);
        if(!hist_SSIso_data) return;

        TH1D& histogram = CloneHistogram(eventCategory, qcd.name, EventRegion::OS_Isolated, *hist_SSIso_data);
        SubtractBackgroundHistograms(histogram, eventCategory, EventRegion::SS_Isolated, qcd.name);
        histogram.Scale(scale_factor);
    }

    virtual std::pair<double, double> CalculateWjetsScaleFactors(analysis::EventCategory eventCategory,
                                                                 const std::string& hist_name) override
    {
        using analysis::EventRegion;
        using analysis::DataCategoryType;

        const analysis::DataCategory& wjets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets);
        const analysis::DataCategory& wjets_mc = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets_MC);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        auto hist_OS_HighMt_mc = GetHistogram(eventCategory, wjets_mc.name, EventRegion::OS_HighMt, hist_name);
        auto hist_SS_HighMt_mc = GetHistogram(eventCategory, wjets_mc.name, EventRegion::SS_HighMt, hist_name);
        auto hist_OS_HighMt_data = GetHistogram(eventCategory, data.name, EventRegion::OS_HighMt, hist_name);
        auto hist_SS_HighMt_data = GetHistogram(eventCategory, data.name, EventRegion::SS_HighMt, hist_name);
        if(!hist_OS_HighMt_mc || !hist_SS_HighMt_mc || !hist_OS_HighMt_data || !hist_SS_HighMt_data)
            throw analysis::exception("Unable to find histograms for Wjet scale factors estimation");

        TH1D& hist_OS_HighMt = CloneHistogram(eventCategory, wjets.name, EventRegion::OS_HighMt, *hist_OS_HighMt_data);
        TH1D& hist_SS_HighMt = CloneHistogram(eventCategory, wjets.name, EventRegion::SS_HighMt, *hist_SS_HighMt_data);

        SubtractBackgroundHistograms(hist_OS_HighMt, eventCategory, EventRegion::OS_HighMt, wjets.name);
        SubtractBackgroundHistograms(hist_SS_HighMt, eventCategory, EventRegion::SS_HighMt, wjets.name);

        const double n_OS_HighMt = hist_OS_HighMt.Integral(0, hist_OS_HighMt.GetNbinsX() + 1);
        const double n_OS_HighMt_mc = hist_OS_HighMt_mc->Integral(0, hist_OS_HighMt_mc->GetNbinsX() + 1);
        const double n_SS_HighMt = hist_SS_HighMt.Integral(0, hist_SS_HighMt.GetNbinsX() + 1);
        const double n_SS_HighMt_mc = hist_SS_HighMt_mc->Integral(0, hist_SS_HighMt_mc->GetNbinsX() + 1);

        const double OS_ratio = n_OS_HighMt / n_OS_HighMt_mc;
        const double SS_ratio = n_SS_HighMt / n_SS_HighMt_mc;
        return std::pair<double, double>(OS_ratio, SS_ratio);
    }

    virtual void EstimateWjets(analysis::EventCategory eventCategory, const std::string& hist_name,
                               std::pair<double, double> scale_factors) override
    {
        using analysis::EventRegion;
        using analysis::DataCategoryType;

        const analysis::DataCategory& wjets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets);
        const analysis::DataCategory& wjets_mc = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets_MC);

        if(auto hist_OS_mc = GetHistogram(eventCategory, wjets_mc.name, EventRegion::OS_Isolated, hist_name)) {
            TH1D& hist_OS = CloneHistogram(eventCategory, wjets.name, EventRegion::OS_Isolated, *hist_OS_mc);
            hist_OS.Scale(scale_factors.first);
        }

        if(auto hist_SS_mc = GetHistogram(eventCategory, wjets_mc.name, EventRegion::SS_Isolated, hist_name)) {
            TH1D& hist_SS = CloneHistogram(eventCategory, wjets.name, EventRegion::SS_Isolated, *hist_SS_mc);
            hist_SS.Scale(scale_factors.second);
        }
    }
};
