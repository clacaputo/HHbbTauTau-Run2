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
#include "../include/Htautau_Summer13.h"
#include "../include/exception.h"

class FlatTreeAnalyzer_tautau : public analysis::BaseFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_tautau(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                         const std::string& outputFileName, const std::string& _signalName,
                         const std::string& _dataName)
        : BaseFlatTreeAnalyzer(source_cfg, hist_cfg, _inputPath, outputFileName, _signalName, _dataName)
    {
    }

protected:
    virtual analysis::EventType_QCD DetermineEventTypeForQCD(const ntuple::Flat& event) override
    {
        using analysis::EventType_QCD;

        if (event.againstElectronLooseMVA_2 > 0.5){
            // OS - isolated
            if (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < 1.5 &&
                    event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < 1.5 &&
                    (event.q_1 * event.q_2 == -1))
                return EventType_QCD::OS_Isolated;
            else if (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < 1.5 &&
                    event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < 1.5 &&
                    (event.q_1 * event.q_2 == +1)) // SS - isolated
                return EventType_QCD::SS_Isolated;
            else if ((event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= 1.5 ||
                     event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= 1.5) &&
                     (event.q_1 * event.q_2 == -1)) // OS - not isolated
                return EventType_QCD::OS_NotIsolated;
            else if ((event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= 1.5 ||
                     event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= 1.5) &&
                     (event.q_1 * event.q_2 == +1)) // SS - not isolated
                return EventType_QCD::SS_NotIsolated;
            else
                return EventType_QCD::Unknown;
        }
        else
            return EventType_QCD::Unknown;
    }

    virtual analysis::EventType_Wjets DetermineEventTypeForWjets(const ntuple::Flat& event) override
    {
        using analysis::EventType_Wjets;
        return EventType_Wjets::Unknown;
    }

    virtual analysis::EventCategory DetermineEventCategory(const ntuple::Flat& event) override
    {
        using analysis::EventCategory;

        std::vector<Float_t> goodCVSvalues;
        for (unsigned i = 0; i < event.eta_Bjets.size(); ++i){
            if ( std::abs(event.eta_Bjets.at(i)) >= cuts::Htautau_Summer13::btag::eta) continue;
            goodCVSvalues.push_back(event.csv_Bjets.at(i));
        }

        if (goodCVSvalues.size() < 2) return EventCategory::Unknown;

        if (goodCVSvalues.at(0) <= cuts::Htautau_Summer13::btag::CSVM )
            return EventCategory::TwoJets_ZeroBtag;
        else if ( goodCVSvalues.at(1) <= cuts::Htautau_Summer13::btag::CSVM )
            return EventCategory::TwoJets_OneBtag;
        else
            return EventCategory::TwoJets_TwoBtag;
    }

    virtual void EstimateQCD() override
    {
        analysis::DataCategory qcd;
        qcd.name = "QCD";
        qcd.title = "QCD";
        qcd.color = kPink;
        for (auto& fullAnaDataEntry : fullAnaData){
            const EventCategory eventCategory = fullAnaDataEntry.first;
            AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
            for (const analysis::HistogramDescriptor& hist : histograms){
                const analysis::DataCategory& data = FindCategory("DATA");
                if(!anaData[data.name].QCD[EventType_QCD::OS_NotIsolated].Contains(hist.name)) continue;
                anaData[qcd.name].QCD[EventType_QCD::OS_Isolated].Get<TH1D>(hist.name) =
                        anaData[data.name].QCD[EventType_QCD::OS_NotIsolated].Get<TH1D>(hist.name);
                TH1D& histogram = anaData[qcd.name].QCD[EventType_QCD::OS_Isolated].Get<TH1D>(hist.name);
                for (const analysis::DataCategory& category : categories){
                    const bool isData = category.name.find("DATA") != std::string::npos;
                    const bool isSignal = category.name.find("SIGNAL") != std::string::npos;
                    if (isData || isSignal) continue;
                    if(!anaData[category.name].QCD[EventType_QCD::OS_NotIsolated].Contains(hist.name)) continue;
                    TH1D& nonQCD_hist = anaData[category.name].QCD[EventType_QCD::OS_NotIsolated].Get<TH1D>(hist.name);
                    histogram.Add(&nonQCD_hist,-1);
                }
                histogram.Scale(1.06);
            }
        }
        categories.push_back(qcd);
    }

    const analysis::DataCategory& FindCategory(const std::string& prefix) const
    {
        for (const analysis::DataCategory& category : categories){
            const bool foundPrefix = category.name.find(prefix) != std::string::npos;
            if (foundPrefix)
                return category;
        }
        throw analysis::exception("category not found : ") << prefix ;
    }

    virtual void EstimateWjets() override
    {

    }

};
