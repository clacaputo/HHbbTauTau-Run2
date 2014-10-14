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

        if (eventInfo.event->kinfit_bb_tt_chi2.size() != eventInfo.event->kinfit_bb_tt_convergence.size())
            std::cout << "aaa - size chi2 = " << eventInfo.event->kinfit_bb_tt_chi2.size() <<
                         "- size convergence = " << eventInfo.event->kinfit_bb_tt_convergence.size() << std::endl;


        const auto bjetsSelector = [&] (const size_t& first, const size_t& second) -> bool
        {
            if (eventInfo.event->kinfit_bb_tt_convergence.at(first) > 0
                    && eventInfo.event->kinfit_bb_tt_convergence.at(second) <= 0)
                return true;
            if (eventInfo.event->kinfit_bb_tt_convergence.at(first) <= 0
                    && eventInfo.event->kinfit_bb_tt_convergence.at(second) > 0)
                return false;
            Float_t first_chi2 = 0;
            Float_t second_chi2 = 0;
            if ((eventInfo.event->kinfit_bb_tt_convergence.at(first) > 0
                    && eventInfo.event->kinfit_bb_tt_convergence.at(second) > 0) ||
                    (eventInfo.event->kinfit_bb_tt_convergence.at(first) <= 0
                     && eventInfo.event->kinfit_bb_tt_convergence.at(second) <= 0)){
                first_chi2 = eventInfo.event->kinfit_bb_tt_chi2.at(first);
                second_chi2 = eventInfo.event->kinfit_bb_tt_chi2.at(second);
            }
            return first_chi2 < second_chi2;
        };

        std::sort(bjetsPair_indexes.begin(), bjetsPair_indexes.end(), bjetsSelector);
        return bjetsPair_indexes;
    }

    std::vector<size_t> SortBjetsByMassPair(const FlatEventInfo& eventInfo)
    {
        std::vector<size_t> bjetsPair_indexes(eventInfo.event->kinfit_bb_tt_chi2.size());
        std::iota(bjetsPair_indexes.begin(), bjetsPair_indexes.end(), 0);

        const auto bjetsSelector = [&] (const size_t& first, const size_t& second) -> bool
        {
            const FlatEventInfo::BjetPair first_bjetPair =
                    FlatEventInfo::CombinationIndexToPair(first, eventInfo.event->energy_Bjets.size());
            const FlatEventInfo::BjetPair second_bjetPair =
                    FlatEventInfo::CombinationIndexToPair(second, eventInfo.event->energy_Bjets.size());
            const TLorentzVector b1_firstPair = eventInfo.bjet_momentums.at(first_bjetPair.first);
            const TLorentzVector b2_firstPair = eventInfo.bjet_momentums.at(first_bjetPair.second);
            const TLorentzVector firstMbb = b1_firstPair + b2_firstPair;
            const TLorentzVector b1_secondPair = eventInfo.bjet_momentums.at(second_bjetPair.first);
            const TLorentzVector b2_secondPair = eventInfo.bjet_momentums.at(second_bjetPair.second);
            const TLorentzVector secondMbb = b1_secondPair + b2_secondPair;
            return std::abs(firstMbb.M() - 125) < std::abs(secondMbb.M() - 125);
        };

        std::sort(bjetsPair_indexes.begin(), bjetsPair_indexes.end(), bjetsSelector);
        return bjetsPair_indexes;
    }

    analysis::FlatEventInfo::BjetPair GetBjetPairByMass(const TLorentzVector& Hbb_CSV, const FlatEventInfo& eventInfo)
    {
        analysis::FlatEventInfo::BjetPair pairCombined;
        if ((Hbb_CSV.M() < 60 || Hbb_CSV.M() > 160) && eventInfo.event->energy_Bjets.size() > 2){
            pairCombined.first = 0;
            pairCombined.second = 2;
            TLorentzVector mBB = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(2);
            if ((mBB.M() < 60 || mBB.M() > 160) && eventInfo.event->energy_Bjets.size() > 3){
                pairCombined.first = 0;
                pairCombined.second = 3;
            }
            else {
                pairCombined.first = 0;
                pairCombined.second = 2;
            }
        }
        else {
            pairCombined.first = 0;
            pairCombined.second = 1;
        }
        return pairCombined;
    }

    analysis::FlatEventInfo::BjetPair GetBjetPairByMass(const TLorentzVector& Hbb_CSV, const FlatEventInfo& eventInfo)
    {
        analysis::FlatEventInfo::BjetPair pairCombined;
        if ((Hbb_CSV.M() < 60 || Hbb_CSV.M() > 160) && eventInfo.event->energy_Bjets.size() > 2){
            pairCombined.first = 0;
            pairCombined.second = 2;
        }
        else {
            pairCombined.first = 0;
            pairCombined.second = 1;
        }
        return pairCombined;
    }

    analysis::FlatEventInfo::BjetPair GetBjetPairByBestMass(const FlatEventInfo& eventInfo)
    {
        analysis::FlatEventInfo::BjetPair candidatePair(0,1);
        if (eventInfo.event->energy_Bjets.size() > 3){
            const TLorentzVector bbPair01 = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(1);
            const TLorentzVector bbPair02 = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(2);
            const TLorentzVector bbPair03 = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(3);
            if (std::abs(bbPair01.M() - 110) < std::abs(bbPair02.M() - 110) &&
                    std::abs(bbPair01.M() - 110) < std::abs(bbPair03.M() - 110)){
                candidatePair.first = 0;
                candidatePair.second = 1;
            }
            else if (std::abs(bbPair02.M() - 110) < std::abs(bbPair01.M() - 110) &&
                     std::abs(bbPair02.M() - 110) < std::abs(bbPair03.M() - 110)){
                candidatePair.first = 0;
                candidatePair.second = 2;
            }
            else if (std::abs(bbPair03.M() - 110) < std::abs(bbPair01.M() - 110) &&
                     std::abs(bbPair03.M() - 110) < std::abs(bbPair02.M() - 110)){
                candidatePair.first = 0;
                candidatePair.second = 3;
            }
        }
        return candidatePair;
    }

    analysis::FlatEventInfo::BjetPair GetBjetPairByBestMass(const FlatEventInfo& eventInfo)
    {
        analysis::FlatEventInfo::BjetPair candidatePair(0,1);
        if (eventInfo.event->energy_Bjets.size() > 2){
            const TLorentzVector bbPair01 = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(1);
            const TLorentzVector bbPair02 = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(2);

            if (std::abs(bbPair01.M() - 110) < std::abs(bbPair02.M() - 110)){
                candidatePair.first = 0;
                candidatePair.second = 1;
            }
            else {
                candidatePair.first = 0;
                candidatePair.second = 2;
            }

        }
        return candidatePair;
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
        //std::cout << "pt pair (" << first_index << "," << second_index << ")" << std::endl;
        if (eventInfo.event->isBjet_MC_Bjet.at(first_index))
            ++matchedBjets;
        if (eventInfo.event->isBjet_MC_Bjet.at(second_index))
            ++matchedBjets;
        return matchedBjets;
    }

    unsigned MatchedBjetsByChi2(const FlatEventInfo& eventInfo, const std::vector<size_t>& indexes)
    {
        unsigned matchedBjets = 0;

        FlatEventInfo::BjetPair bjetPair =
                FlatEventInfo::CombinationIndexToPair(indexes.at(0),eventInfo.event->energy_Bjets.size());
        //std::cout << "bjets chi2 pair (" << bjetPair.first << "," << bjetPair.second << ")" << std::endl;
        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.first))
            ++matchedBjets;
        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.second))
            ++matchedBjets;

        return matchedBjets;
    }

    unsigned MatchedBjetsByMassPair(const FlatEventInfo& eventInfo, const std::vector<size_t>& indexes)
    {
        unsigned matchedBjets = 0;

        FlatEventInfo::BjetPair bjetPair =
                FlatEventInfo::CombinationIndexToPair(indexes.at(0),eventInfo.event->energy_Bjets.size());
        //std::cout << "bjets mass pair (" << bjetPair.first << "," << bjetPair.second << ")" << std::endl;
        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.first))
            ++matchedBjets;
        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.second))
            ++matchedBjets;

        return matchedBjets;
    }

    unsigned MatchedCombinedCSVandMASS(const FlatEventInfo& eventInfo, const FlatEventInfo::BjetPair& bjetPair)
    {
        unsigned matchedBjets = 0;
        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.first))
            ++matchedBjets;
        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.second))
            ++matchedBjets;
        return matchedBjets;
    }


private:


protected:
    std::shared_ptr<TFile> outputFile;
    std::shared_ptr<ntuple::FlatTree> flatTree;
    LightFlatAnalyzerData anaData;

};

} // namespace analysis
