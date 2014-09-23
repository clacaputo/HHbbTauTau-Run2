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
                           const std::string& outputFileName, const std::string& _signalName,
                           const std::string& _dataName, const std::string& _mvaXMLpath, bool _WjetsData = false,
                           bool _isBlind=false)
         : BaseFlatTreeAnalyzer(source_cfg, hist_cfg, _inputPath, outputFileName, _signalName, _dataName, _mvaXMLpath,
                                _WjetsData, _isBlind)
    {
    }

protected:
    virtual const std::string& ChannelName() override
    {
        static const std::string channelName = "muTau";
        return channelName;
    }

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

    virtual void EstimateQCD(analysis::EventCategory /*eventCategory*/, AnaDataForDataCategory& anaData,
                             const analysis::HistogramDescriptor& hist) override
    {
        static const double scale_factor = 1.06;

        using analysis::EventType_QCD;
        const analysis::DataCategory& qcd = FindCategory("QCD");

        const analysis::DataCategory& data = FindCategory("DATA");
        root_ext::SmartHistogram<TH1D>* hist_data;
        if(!(hist_data = anaData[data.name].QCD[EventType_QCD::SS_Isolated].GetPtr<TH1D>(hist.name))) return;
        TH1D& histogram = anaData[qcd.name].QCD[EventType_QCD::OS_Isolated].Clone(*hist_data);
        for (const analysis::DataCategory& category : categories){
            if (category.IsData() || category.IsSignal() || category.name == qcd.name
                    || category.IsForLimitsOnly()) continue;
            TH1D* nonQCD_hist;
            if(!(nonQCD_hist = anaData[category.name].QCD[EventType_QCD::SS_Isolated].GetPtr<TH1D>(hist.name))) continue;
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

    virtual void EstimateWjets(analysis::EventCategory eventCategory, AnaDataForDataCategory& anaData,
                               const analysis::HistogramDescriptor& hist) override
    {
        using analysis::EventType_Wjets;
        using analysis::EventType_QCD;

        //data
        const analysis::DataCategory& data = FindCategory("DATA");
        TH1D* histData_HighMt;
        if(!(histData_HighMt = anaData[data.name].Wjets[EventType_Wjets::HighMt].GetPtr<TH1D>(hist.name))) return;

        //MC wjets
        const analysis::DataCategory& wjets = FindCategory("REFERENCE Wjets");
        TH1D* histWjetsHighMt;
        if(!(histWjetsHighMt = anaData[wjets.name].Wjets[EventType_Wjets::HighMt].GetPtr<TH1D>(hist.name)))
            throw analysis::exception("histogram for Wjets High Mt category doesn't exist");

        TH1D* histWjetsSignal;
        if(!(histWjetsSignal = anaData[wjets.name].Wjets[EventType_Wjets::Signal].GetPtr<TH1D>(hist.name)))
            throw analysis::exception("histogram for Wjets Signal category doesn't exist");

        for (const analysis::DataCategory& category : categories){
            if (category.IsData() || category.IsSignal() || category.IsReference() || category.IsVirtual()
                    || category.IsForLimitsOnly()) continue;
            TH1D* nonWjets_hist;
            if(!(nonWjets_hist = anaData[category.name].Wjets[EventType_Wjets::HighMt].GetPtr<TH1D>(hist.name))) continue;
            histData_HighMt->Add(nonWjets_hist,-1);
        }

        const double ratio = histWjetsSignal->Integral()/histWjetsHighMt->Integral();
        std::cout << eventCategory << " Wjets_signal/Wjets_HighMt = " << ratio << std::endl;
        histData_HighMt->Scale(ratio);

        const analysis::DataCategory& ewk = FindCategory("EWK");
        TH1D* histEWK;
        if(!(histEWK = anaData[ewk.name].QCD[EventType_QCD::OS_Isolated].GetPtr<TH1D>(hist.name)))
            throw analysis::exception("histogram not found: ") << hist.name;
        histEWK->Add(histData_HighMt);
    }

};
