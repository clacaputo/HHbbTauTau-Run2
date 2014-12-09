/*!
 * \file LightBaseFlatTreeAnalyzer.h
 * \brief Definition of LightBaseFlatTreeAnalyzer class, the base class for separate studies on flat trees.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-10-08 created
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

#pragma once

#ifndef __APPLE__
#define override
#endif

#include <iostream>
#include <cmath>
#include <set>
#include <list>
#include <locale>
#include <numeric>
#include <algorithm>


#include <TColor.h>
#include <TLorentzVector.h>


#include "AnalysisBase/include/AnalyzerData.h"
#include "AnalysisBase/include/FlatEventInfo.h"
#include "AnalysisBase/include/AnalysisMath.h"
#include "AnalysisBase/include/exception.h"
#include "AnalysisBase/include/Particles.h"
#include "PrintTools/include/RootPrintToPdf.h"

#include "MVASelections/include/MvaReader.h"

#include "Htautau_Summer13.h"
#include "AnalysisCategories.h"

namespace analysis {

class LightBaseFlatTreeAnalyzer {
public:
    LightBaseFlatTreeAnalyzer(const std::string& inputFileName, const std::string& outputFileName)
        : inputFile(new TFile(inputFileName.c_str(),"READ")),
          outputFile(new TFile(outputFileName.c_str(), "RECREATE")),
          flatTree(new ntuple::FlatTree(*inputFile, "flatTree")), recalc_kinfit(true)
    {
        TH1::SetDefaultSumw2();
    }

    virtual ~LightBaseFlatTreeAnalyzer() {}

    void Run()
    {
        using namespace cuts::Htautau_Summer13::btag;
        for(Long64_t current_entry = 0; current_entry < flatTree->GetEntries(); ++current_entry) {
            flatTree->GetEntry(current_entry);
            const ntuple::Flat& event = flatTree->data;
            const EventCategoryVector eventCategories =
                    DetermineEventCategories(event.csv_Bjets, event.nBjets_retagged, CSVL, CSVM,true);
            FlatEventInfo eventInfo(event, FlatEventInfo::BjetPair(0, 1), recalc_kinfit);
            for (EventCategory eventCategory : eventCategories)
                AnalyzeEvent(eventInfo, eventCategory);
        }
        EndOfRun();
    }

protected:
    virtual void AnalyzeEvent(const FlatEventInfo& eventInfo, EventCategory eventCategory) = 0;
    virtual void EndOfRun(){}


    bool IsHighMtRegion(const ntuple::Flat& event, analysis::EventCategory eventCategory) const
    {
        using namespace cuts;
        if (eventCategory == analysis::EventCategory::TwoJets_TwoBtag)
            return event.mt_1 > WjetsBackgroundEstimation::HighMtRegion_low &&
                    event.mt_1 < WjetsBackgroundEstimation::HighMtRegion_high;
        else
            return event.mt_1 > WjetsBackgroundEstimation::HighMtRegion;
    }

    bool IsAntiIsolatedRegion(const ntuple::Flat& event) const
    {
        using namespace cuts;
        return event.pfRelIso_1 > IsolationRegionForLeptonicChannel::isolation_low &&
                event.pfRelIso_1 < IsolationRegionForLeptonicChannel::isolation_high;
    }

    TFile& GetOutputFile() { return *outputFile; }

    analysis::EventRegion DetermineEventRegion(const FlatEventInfo& eventInfo, analysis::EventCategory category) const
    {
        using analysis::EventRegion;
        const ntuple::Flat& event = *eventInfo.event;
        if (eventInfo.channel == Channel::MuTau){
            using namespace cuts::Htautau_Summer13::MuTau;

            if(!event.againstMuonTight_2
                    || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                    || (event.pfRelIso_1 >= muonID::pFRelIso && !IsAntiIsolatedRegion(event))
                    || (event.mt_1 >= muonID::mt && !IsHighMtRegion(event,category)))
                return EventRegion::Unknown;

            const bool os = event.q_1 * event.q_2 == -1;
            const bool iso = event.pfRelIso_1 < muonID::pFRelIso;
            const bool low_mt = event.mt_1 < muonID::mt;


            if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_Iso_HighMt;
            if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_Iso_HighMt;
            if(os) return low_mt ? EventRegion::OS_AntiIsolated : EventRegion::OS_AntiIso_HighMt;
            return low_mt ? EventRegion::SS_AntiIsolated : EventRegion::SS_AntiIso_HighMt;
        }
        if (eventInfo.channel == Channel::ETau){

            using namespace cuts::Htautau_Summer13::ETau;

            if(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                    || (event.pfRelIso_1 >= electronID::pFRelIso && !IsAntiIsolatedRegion(event))
                    || (event.mt_1 >= electronID::mt && !IsHighMtRegion(event,category)) /*|| event.pt_2 <= 30*/ )
                return EventRegion::Unknown;

            const bool os = event.q_1 * event.q_2 == -1;
            const bool iso = event.pfRelIso_1 < electronID::pFRelIso;
            const bool low_mt = event.mt_1 < electronID::mt;


            if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_Iso_HighMt;
            if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_Iso_HighMt;
            if(os) return low_mt ? EventRegion::OS_AntiIsolated : EventRegion::OS_AntiIso_HighMt;
            return low_mt ? EventRegion::SS_AntiIsolated : EventRegion::SS_AntiIso_HighMt;
        }
        if (eventInfo.channel == Channel::TauTau){

            using namespace cuts::Htautau_Summer13::TauTau::tauID;

            if(!event.againstElectronLooseMVA_2
                    || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= BackgroundEstimation::Isolation_upperLimit
                    || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= BackgroundEstimation::Isolation_upperLimit
                    || (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits
                        && event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits))
                return EventRegion::Unknown;

            const bool os = event.q_1 * event.q_2 == -1;
            const bool iso = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                             event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < byCombinedIsolationDeltaBetaCorrRaw3Hits;

            if(iso) return os ? EventRegion::OS_Isolated : EventRegion::SS_Isolated;
            return os ? EventRegion::OS_AntiIsolated : EventRegion::SS_AntiIsolated;
        }
        throw analysis::exception("unsupported channel ") << eventInfo.channel;
    }

private:
    std::shared_ptr<TFile> inputFile, outputFile;
    std::shared_ptr<ntuple::FlatTree> flatTree;


protected:
    bool recalc_kinfit;
};

} // namespace analysis
