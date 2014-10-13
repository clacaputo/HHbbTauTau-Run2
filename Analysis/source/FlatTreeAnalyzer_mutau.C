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

    virtual analysis::EventType_QCD DetermineEventTypeForQCD(const ntuple::Flat& event) override
    {
        using analysis::EventType_QCD;
        using namespace cuts::Htautau_Summer13::MuTau;

        if(!event.againstMuonTight_2
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || event.mt_1 >= muonID::mt)
            return EventType_QCD::Unknown;

        if(event.pfRelIso_1 < muonID::pFRelIso && (event.q_1 * event.q_2) == -1 )
            return EventType_QCD::OS_Isolated;

        if(event.pfRelIso_1 < muonID::pFRelIso && (event.q_1 * event.q_2) == +1 )
            return EventType_QCD::SS_Isolated;

        if(event.pfRelIso_1 >= muonID::pFRelIso && (event.q_1 * event.q_2) == -1 )
            return EventType_QCD::OS_NotIsolated;

        if(event.pfRelIso_1 >= muonID::pFRelIso && (event.q_1 * event.q_2) == +1 )
            return EventType_QCD::SS_NotIsolated;

        return EventType_QCD::Unknown;
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

    virtual double CalculateQCDScaleFactor(analysis::EventCategory /*eventCategory*/, AnaDataForDataCategory& anaData,
                                             const std::string& hist_name) override
    {
        using analysis::EventType_QCD;
        using analysis::DataCategoryType;

        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        auto hist_OSnotIso_data = anaData[data.name].QCD[EventType_QCD::OS_NotIsolated].GetPtr<TH1D>(hist_name);
        auto hist_SSnotIso_data = anaData[data.name].QCD[EventType_QCD::SS_NotIsolated].GetPtr<TH1D>(hist_name);
        if(!hist_OSnotIso_data || !hist_SSnotIso_data)
            throw analysis::exception("Unable to find histograms for QCD scale factor estimation");

        TH1D& hist_OSnotIso = anaData[qcd.name].QCD[EventType_QCD::OS_NotIsolated].Clone(*hist_OSnotIso_data);
        TH1D& hist_SSnotIso = anaData[qcd.name].QCD[EventType_QCD::SS_NotIsolated].Clone(*hist_SSnotIso_data);

        for (auto category : dataCategoryCollection.GetCategories(analysis::DataCategoryType::Background)) {
            if(category->types.count(DataCategoryType::Composit)) continue;

            if( TH1D* nonQCD_histIso = anaData[category->name].QCD[EventType_QCD::OS_NotIsolated].GetPtr<TH1D>(hist_name) )
                hist_OSnotIso.Add(nonQCD_histIso, -1);

            if( TH1D* nonQCD_histNotIso = anaData[category->name].QCD[EventType_QCD::SS_NotIsolated].GetPtr<TH1D>(hist_name) )
                hist_SSnotIso.Add(nonQCD_histNotIso, -1);
        }

        const double n_OSnotIso = hist_OSnotIso.Integral(0, hist_OSnotIso.GetNbinsX() + 1);
        const double n_SSnotIso = hist_SSnotIso.Integral(0, hist_SSnotIso.GetNbinsX() + 1);
        return n_OSnotIso / n_SSnotIso;
    }


    virtual void EstimateQCD(analysis::EventCategory /*eventCategory*/, AnaDataForDataCategory& anaData,
                             const std::string& hist_name, double scale_factor) override
    {
        using analysis::EventType_QCD;
        using analysis::DataCategoryType;

        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::Data);

        auto hist_SSIso_data = anaData[data.name].QCD[EventType_QCD::SS_Isolated].GetPtr<TH1D>(hist_name);
        if(!hist_SSIso_data) return;

        TH1D& histogram = anaData[qcd.name].QCD[EventType_QCD::OS_Isolated].Clone(*hist_SSIso_data);

        for (auto category : dataCategoryCollection.GetCategories(analysis::DataCategoryType::Background)) {
            if(category->types.count(DataCategoryType::Composit)) continue;

            if( TH1D* nonQCD_hist = anaData[category->name].QCD[EventType_QCD::SS_Isolated].GetPtr<TH1D>(hist_name) )
                histogram.Add(nonQCD_hist, -1);
        }

        histogram.Scale(scale_factor);
    }


    virtual analysis::EventType_Wjets DetermineEventTypeForWjets(const ntuple::Flat& event) override
    {
        using analysis::EventType_Wjets;
        using namespace cuts::Htautau_Summer13::MuTau;

        if(!event.againstMuonTight_2
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || event.pfRelIso_1 >= muonID::pFRelIso)
            return EventType_Wjets::Unknown;

        if(event.mt_1 < muonID::mt)
            return EventType_Wjets::Signal;
        if(event.mt_1 > BackgroundEstimation::HighMtRegion)
            return EventType_Wjets::HighMt;
        return EventType_Wjets::Unknown;
    }

    virtual double CalculateWjetsScaleFactor(analysis::EventCategory /*eventCategory*/, AnaDataForDataCategory& anaData,
                                             const std::string& hist_name) override
    {
        using analysis::EventType_Wjets;
        using analysis::DataCategoryType;

        const analysis::DataCategory& wjets_mc = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets_MC);

        TH1D* histWjetsHighMt = anaData[wjets_mc.name].Wjets[EventType_Wjets::HighMt].GetPtr<TH1D>(hist_name);
        if(!histWjetsHighMt)
            throw analysis::exception("histogram for Wjets High Mt category doesn't exist");

        TH1D* histWjetsSignal = anaData[wjets_mc.name].Wjets[EventType_Wjets::Signal].GetPtr<TH1D>(hist_name);
        if(!histWjetsSignal)
            throw analysis::exception("histogram for Wjets Signal category doesn't exist");

        const double n_wjets_signal = histWjetsSignal->Integral(0, histWjetsSignal->GetNbinsX() + 1);
        const double n_wjets_high_mt = histWjetsHighMt->Integral(0, histWjetsHighMt->GetNbinsX() + 1);
        return n_wjets_signal / n_wjets_high_mt;
    }

    virtual void EstimateWjets(analysis::EventCategory /*eventCategory*/, AnaDataForDataCategory& anaData,
                               const std::string& hist_name, double scale_factor) override
    {
        using analysis::EventType_Wjets;
        using analysis::DataCategoryType;

        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);
        const analysis::DataCategory& wjets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets);

        auto histData_HighMt = anaData[data.name].Wjets[EventType_Wjets::HighMt].GetPtr<TH1D>(hist_name);
        if(!histData_HighMt) return;
        TH1D& histogram = anaData[wjets.name].Signal().Clone(*histData_HighMt);

        for (auto category : dataCategoryCollection.GetCategories(DataCategoryType::Background)) {
            if(category->types.count(DataCategoryType::Composit)) continue;
            if(TH1D* nonWjets_hist = anaData[category->name].Wjets[EventType_Wjets::HighMt].GetPtr<TH1D>(hist_name))
                histogram.Add(nonWjets_hist,-1);
        }

        histogram.Scale(scale_factor);
    }
};
