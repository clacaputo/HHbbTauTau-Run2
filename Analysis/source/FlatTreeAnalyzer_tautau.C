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
                            const std::string& _dataName)
        : BaseFlatTreeAnalyzer(source_cfg, hist_cfg, _inputPath, outputFileName, _signalName, _dataName)
    {
    }

protected:
    virtual analysis::EventType_QCD DetermineEventTypeForQCD(const ntuple::Flat& event) override
    {
        using analysis::EventType_QCD;
        using namespace cuts::Htautau_Summer13::TauTau::tauID;

        if (event.againstElectronLooseMVA_2 < againstElectronLooseMVA3)
            return EventType_QCD::Unknown;

        // OS - isolated
        if (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                (event.q_1 * event.q_2 == -1))
            return EventType_QCD::OS_Isolated;

        if (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                (event.q_1 * event.q_2 == +1)) // SS - isolated
            return EventType_QCD::SS_Isolated;

        if ((event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits ||
                 event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits) &&
                 (event.q_1 * event.q_2 == -1)) // OS - not isolated
            return EventType_QCD::OS_NotIsolated;

        if ((event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits ||
                 event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits) &&
                 (event.q_1 * event.q_2 == +1)) // SS - not isolated
            return EventType_QCD::SS_NotIsolated;

        return EventType_QCD::Unknown;
    }

    virtual void EstimateQCD() override
    {
        using analysis::EventType_QCD;
        const analysis::DataCategory& qcd = FindCategory("QCD");
        for (auto& fullAnaDataEntry : fullAnaData){
            AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
            for (const analysis::HistogramDescriptor& hist : histograms){
                const analysis::DataCategory& data = FindCategory("DATA");
                if(!anaData[data.name].QCD[EventType_QCD::SS_Isolated].Contains(hist.name)) continue;
                if(!anaData[data.name].QCD[EventType_QCD::SS_NotIsolated].Contains(hist.name)) continue;
                if(!anaData[data.name].QCD[EventType_QCD::OS_NotIsolated].Contains(hist.name)) continue;
                TH1D& histogram = anaData[qcd.name].QCD[EventType_QCD::OS_Isolated].Clone(
                            anaData[data.name].QCD[EventType_QCD::OS_NotIsolated].Get<TH1D>(hist.name));
                TH1D& hist_SSIso = anaData[qcd.name].QCD[EventType_QCD::SS_Isolated].Clone(
                            anaData[data.name].QCD[EventType_QCD::SS_Isolated].Get<TH1D>(hist.name));
                TH1D& hist_SSnotIso = anaData[qcd.name].QCD[EventType_QCD::SS_NotIsolated].Clone(
                            anaData[data.name].QCD[EventType_QCD::SS_NotIsolated].Get<TH1D>(hist.name));
                for (const analysis::DataCategory& category : categories){
                    const bool isData = category.name.find("DATA") != std::string::npos;
                    const bool isSignal = category.name.find("SIGNAL") != std::string::npos;
                    if (isData || isSignal) continue;
                    if(anaData[category.name].QCD[EventType_QCD::OS_NotIsolated].Contains(hist.name)) {
                        TH1D& nonQCD_hist = anaData[category.name].QCD[EventType_QCD::OS_NotIsolated].Get<TH1D>(hist.name);
                        histogram.Add(&nonQCD_hist,-1);
                    }
                    if(anaData[category.name].QCD[EventType_QCD::SS_Isolated].Contains(hist.name)) {
                        TH1D& nonQCD_histIso = anaData[category.name].QCD[EventType_QCD::SS_Isolated].Get<TH1D>(hist.name);
                        hist_SSIso.Add(&nonQCD_histIso,-1);
                    }
                    if(anaData[category.name].QCD[EventType_QCD::SS_NotIsolated].Contains(hist.name)) {
                        TH1D& nonQCD_histNotIso =
                                anaData[category.name].QCD[EventType_QCD::SS_NotIsolated].Get<TH1D>(hist.name);
                        hist_SSnotIso.Add(&nonQCD_histNotIso,-1);
                    }
                }
                const double ratio = hist_SSIso.Integral()/hist_SSnotIso.Integral();
                std::cout << fullAnaDataEntry.first << " SS_Iso/SS_NotIso = " << ratio << std::endl;
                histogram.Scale(ratio);
            }
        }
    }

    virtual analysis::EventType_Wjets DetermineEventTypeForWjets(const ntuple::Flat& event) override
    {
        using analysis::EventType_Wjets;
        return EventType_Wjets::Unknown;
    }

    virtual void EstimateWjets() override
    {

    }

};
