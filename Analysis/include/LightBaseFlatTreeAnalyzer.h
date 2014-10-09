/*!
 * \file BaseFlatTreeAnalyzer.h
 * \brief Definition of BaseFlatTreeAnalyzer class, the base class for flat tree analyzers.
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
//#include "KinFit.h"

#include "MVASelections/include/MvaReader.h"

#include "Htautau_Summer13.h"
#include "AnalysisCategories.h"

namespace analysis {

class LightFlatAnalyzerData : public root_ext::AnalyzerData {
public:
    LightFlatAnalyzerData(TFile& outputFile) : AnalyzerData(outputFile) {}

    TH1D_ENTRY(MET, 35, 0, 350)

};

class LightBaseFlatTreeAnalyzer {
public:


    LightBaseFlatTreeAnalyzer(const std::string& _inputFileName, const std::string& _outputFileName)
        : outputFile(new TFile(_outputFileName.c_str(), "RECREATE")),anaData(*outputFile)
    {
        TH1::SetDefaultSumw2();
        TFile* inputFile = new TFile(_inputFileName.c_str(),"READ");
        flatTree = std::shared_ptr<ntuple::FlatTree>(new ntuple::FlatTree(*inputFile, "flatTree"));

    }

    virtual ~LightBaseFlatTreeAnalyzer() {}

    void Run()
    {
        for(size_t current_entry = 0; current_entry < flatTree->GetEntries(); ++current_entry) {
            flatTree->GetEntry(current_entry);
            const ntuple::Flat& event = flatTree->data;
//            std::cout << "event = " << event.evt << std::endl;
            const EventCategoryVector eventCategories = DetermineEventCategories(event);
            FlatEventInfo eventInfo(event, FlatEventInfo::BjetPair(0, 1));
            //std::cout << "MET: " << eventInfo.MET.Pt() << std::endl;
            //anaData.MET().Fill(eventInfo.MET.Pt());
            for (auto eventCategory : eventCategories){
                FullSelection(eventInfo, eventCategory);
            }
        }

    }



protected:
    virtual LightFlatAnalyzerData& GetAnaData() = 0;
    virtual const std::string& ChannelName() = 0;
    virtual void FullSelection(const FlatEventInfo& eventInfo, const EventCategory& eventCategory) = 0;

    virtual EventCategoryVector DetermineEventCategories(const ntuple::Flat& event)
    {
        EventCategoryVector categories;
        categories.push_back(EventCategory::Inclusive);

        if (event.csv_Bjets.size() == 1) {
            if (event.csv_Bjets.at(0) <= cuts::Htautau_Summer13::btag::CSVT )
                categories.push_back(EventCategory::OneJet_ZeroBtag);
            else
                categories.push_back(EventCategory::OneJet_OneBtag);
        }

        if (event.csv_Bjets.size() >= 2) {
            if (event.csv_Bjets.at(0) <= cuts::Htautau_Summer13::btag::CSVM )
                categories.push_back(EventCategory::TwoJets_ZeroBtag);
            else if ( event.csv_Bjets.at(1) <= cuts::Htautau_Summer13::btag::CSVM )
                categories.push_back(EventCategory::TwoJets_OneBtag);
            else
                categories.push_back(EventCategory::TwoJets_TwoBtag);
        }
        return categories;
    }

    std::vector<size_t> SortBjetsByPt(const FlatEventInfo& eventInfo)
    {
        std::vector<size_t> bjet_indexes(eventInfo.bjet_momentums.size());
        std::iota(bjet_indexes.begin(), bjet_indexes.end(), 0);

        const auto bjetsSelector = [&] (const size_t& first, const size_t& second) -> bool
        {
            const TLorentzVector first_bjet = eventInfo.bjet_momentums.at(first);
            const TLorentzVector second_bjet = eventInfo.bjet_momentums.at(second);
            return first_bjet.Pt() > second_bjet.Pt();
        };

        std::sort(bjet_indexes.begin(), bjet_indexes.end(), bjetsSelector);
        return bjet_indexes;
    }

    std::vector<size_t> SortBjetsByChiSquare(const FlatEventInfo& eventInfo)
    {
        std::vector<size_t> bjetsPair_indexes(eventInfo.event->kinfit_bb_tt_chi2.size());
        std::iota(bjetsPair_indexes.begin(), bjetsPair_indexes.end(), 0);

        const auto bjetsSelector = [&] (const size_t& first, const size_t& second) -> bool
        {
            const Float_t first_chi2 = eventInfo.event->kinfit_bb_tt_chi2.at(first);
            const Float_t second_chi2 = eventInfo.event->kinfit_bb_tt_chi2.at(second);
            return first_chi2 > second_chi2;
        };

        std::sort(bjetsPair_indexes.begin(), bjetsPair_indexes.end(), bjetsSelector);
        return bjetsPair_indexes;
    }

    unsigned MatchedBjetsByCSV(const FlatEventInfo& eventInfo)
    {
        unsigned matchedBjets = 0;
        if (eventInfo.event->isBjet_MC_Bjet.at(0))
            ++matchedBjets;
        if (eventInfo.event->isBjet_MC_Bjet.at(1))
            ++matchedBjets;
        return matchedBjets;
    }

    unsigned MatchedBjetsByPt(const FlatEventInfo& eventInfo, const std::vector<size_t>& indexes)
    {
        unsigned matchedBjets = 0;
        const size_t first_index = indexes.at(0);
        const size_t second_index = indexes.at(1);
        if (eventInfo.event->isBjet_MC_Bjet.at(first_index))
            ++matchedBjets;
        if (eventInfo.event->isBjet_MC_Bjet.at(second_index))
            ++matchedBjets;
        return matchedBjets;
    }

    unsigned MatchedBjetsByChi2(const FlatEventInfo& eventInfo, const std::vector<size_t>& indexes)
    {
        unsigned matchedBjets = 0;
        for (unsigned i = 0; i < indexes.size(); ++i){
//            FlatEventInfo::BjetPair bjetPair =
//                    FlatEventInfo::CombinationIndexToPair(indexes.at(i),eventInfo.event->energy_Bjets.size());
//            std::cout << "bjet pair(" << bjetPair.first << "," << bjetPair.second << ") - match bjet size = " <<
//                      eventInfo.event->isBjet_MC_Bjet.size() << std::endl;
        }
//        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.first))
//            ++matchedBjets;
//        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.second))
//            ++matchedBjets;

        return matchedBjets;
    }


private:


protected:
    std::shared_ptr<TFile> outputFile;
    std::shared_ptr<ntuple::FlatTree> flatTree;
    LightFlatAnalyzerData anaData;

};

} // namespace analysis
