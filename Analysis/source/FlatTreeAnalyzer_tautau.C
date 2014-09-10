/*!
 * \file FlatTreeAnalyzer_tautau.C
 * \brief Analyze flat-tree for tau-tau channel for HHbbtautau analysis.
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
                            const std::string& _dataName, bool _WjetsData = false, bool _isBlind=false)
        : BaseFlatTreeAnalyzer(source_cfg, hist_cfg, _inputPath, outputFileName, _signalName, _dataName, _WjetsData,
                               _isBlind)
    {
    }

protected:
    virtual const std::string& ChannelName() override
    {
        static const std::string channelName = "tauTau";
        return channelName;
    }

    virtual analysis::EventType_QCD DetermineEventTypeForQCD(const ntuple::Flat& event) override
    {
        using analysis::EventType_QCD;
        using namespace cuts::Htautau_Summer13::TauTau::tauID;
        using namespace cuts::Htautau_Summer13::TauTau::tauID::BackgroundEstimation;

        if (event.againstElectronLooseMVA_2 < againstElectronLooseMVA3)
            return EventType_QCD::Unknown;

        if (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                (event.q_1 * event.q_2 == -1))
            return EventType_QCD::OS_Isolated;

        if (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                (event.q_1 * event.q_2 == +1))
            return EventType_QCD::SS_Isolated;

        if (((event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits &&
              event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < Isolation_upperLimit) ||
                 (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                  event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < Isolation_upperLimit)) &&
                 (event.q_1 * event.q_2 == -1))
            return EventType_QCD::OS_NotIsolated;

        if (((event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits &&
              event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < Isolation_upperLimit) ||
                 (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                  event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < Isolation_upperLimit)) &&
                 (event.q_1 * event.q_2 == +1))
            return EventType_QCD::SS_NotIsolated;

        return EventType_QCD::Unknown;
    }

    virtual void EstimateQCD(analysis::EventCategory eventCategory, AnaDataForDataCategory& anaData,
                             const analysis::HistogramDescriptor& hist) override
    {
        using analysis::EventType_QCD;

        const analysis::DataCategory& qcd = FindCategory("QCD");
        const analysis::DataCategory& data = FindCategory("DATA");

        auto histogram_data = anaData[data.name].QCD[EventType_QCD::OS_NotIsolated].GetPtr<TH1D>(hist.name);
        auto hist_SSIso_data = anaData[data.name].QCD[EventType_QCD::SS_Isolated].GetPtr<TH1D>(hist.name);
        auto hist_SSnotIso_data = anaData[data.name].QCD[EventType_QCD::SS_NotIsolated].GetPtr<TH1D>(hist.name);
        if(!histogram_data || !hist_SSIso_data || !hist_SSnotIso_data) return;

        TH1D& histogram = anaData[qcd.name].QCD[EventType_QCD::OS_Isolated].Clone(*histogram_data);
        TH1D& hist_SSIso = anaData[qcd.name].QCD[EventType_QCD::SS_Isolated].Clone(*hist_SSIso_data);
        TH1D& hist_SSnotIso = anaData[qcd.name].QCD[EventType_QCD::SS_NotIsolated].Clone(*hist_SSnotIso_data);

        for (const analysis::DataCategory& category : categories) {
            if(category.IsData() || category.IsSignal() || category.name == qcd.name
                    || category.IsForLimitsOnly()) continue;

            if( TH1D* nonQCD_hist = anaData[category.name].QCD[EventType_QCD::OS_NotIsolated].GetPtr<TH1D>(hist.name) )
                histogram.Add(nonQCD_hist, -1);

            if( TH1D* nonQCD_histIso = anaData[category.name].QCD[EventType_QCD::SS_Isolated].GetPtr<TH1D>(hist.name) )
                hist_SSIso.Add(nonQCD_histIso, -1);

            if( TH1D* nonQCD_histNotIso = anaData[category.name].QCD[EventType_QCD::SS_NotIsolated].GetPtr<TH1D>(hist.name) )
                hist_SSnotIso.Add(nonQCD_histNotIso, -1);
        }

        const double ratio = hist_SSIso.Integral()/hist_SSnotIso.Integral();
        std::cout << eventCategory << " SS_Iso/SS_NotIso = " << ratio << std::endl;
        histogram.Scale(ratio);
    }

    virtual analysis::EventType_Wjets DetermineEventTypeForWjets(const ntuple::Flat& event) override
    {
        using analysis::EventType_Wjets;
        return EventType_Wjets::Unknown;
    }

    virtual void EstimateWjets(analysis::EventCategory eventCategory, AnaDataForDataCategory& anaData,
                               const analysis::HistogramDescriptor& hist) override
    {

    }

};
