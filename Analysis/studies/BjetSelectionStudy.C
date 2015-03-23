/*!
 * \file BjetSelectionStudy.C
 * \brief Study of different posibilities to select signal b-jets.
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

#include "Analysis/include/LightBaseFlatTreeAnalyzer.h"
#include "Analysis/include/FlatAnalyzerData.h"

class BjetSelectionStudyData : public root_ext::AnalyzerData {
public:
    BjetSelectionStudyData(std::shared_ptr<TFile> outputFile) : root_ext::AnalyzerData(outputFile) {}

    TH1D_ENTRY(MET, 35, 0, 350)
    TH1D_ENTRY(totalMatchedBjets, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjets, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjetsByPt, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjetsByChi2, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjets_byMassPair, 3, -0.5, 2.5)
    TH1D_ENTRY(Hbb_trueM, 30, 0, 300)
    TH1D_ENTRY(Hbb_CSV_fail, 35, 0, 350)
    TH1D_ENTRY(Hbb_CSV_good, 35, 0, 350)
    TH1D_ENTRY(Chi2_CSVfail, 10, 0, 50)
    TH1D_ENTRY(Chi2_CSVgood, 10, 0, 50)
    TH1D_ENTRY(matchedBjets_combinedCSVandMASS, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjets_CSVandBestMassPair, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjets_combinedCSVandMASS_2, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjets_CSVandBestMassPair_2, 3, -0.5, 2.5)
    TH2D_ENTRY(csv_b1_vs_ptb1, 20, 0, 200 ,25, 0.2, 1.2)
};

class BjetSelectionStudy : public analysis::LightBaseFlatTreeAnalyzer {
public:
    BjetSelectionStudy(const std::string& _inputFileName, const std::string& _outputFileName)
         : LightBaseFlatTreeAnalyzer(_inputFileName,_outputFileName), anaData(GetOutputFile())
    {
        recalc_kinfit = true;
        do_retag = false;
    }

protected:
    virtual PairSelectionMap SelectBjetPairs(const ntuple::Flat& event) override
    {
        PairSelectionMap pairMap;
        pairMap["CSV"] = DefaultPair();
        pairMap["Pt"] = SelectBestPtPair(event);
        pairMap["Chi2"] = SelectBestChi2Pair(event);
        return pairMap;
    }

    virtual void AnalyzeEvent(const analysis::FlatEventInfo& eventInfo, analysis::EventCategory category,
                              const std::string& selectionLabel) override
    {
        using analysis::EventCategory;
        static const std::set<EventCategory> categoriesToProcess = {
            EventCategory::TwoJets_Inclusive, EventCategory::TwoJets_ZeroBtag,
            EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag
        };

        if(DetermineEventRegion(eventInfo, category) != analysis::EventRegion::OS_Isolated) return;
        if (!categoriesToProcess.count(category)) return;

        std::ostringstream ss_label;
        ss_label << selectionLabel << "_" << category;
        const std::string label = ss_label.str();

        anaData.MET(label).Fill(eventInfo.MET.Pt());

        anaData.csv_b1_vs_ptb1(label).Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Pt(),
                                           eventInfo.event->csv_Bjets.at(eventInfo.selected_bjets.first));

        size_t totalMatchedBjets = 0;
        std::vector<size_t> indexes;
        TLorentzVector Hbb_true;

        for (size_t k = 0; k < eventInfo.event->energy_Bjets.size(); ++k) {
            if (eventInfo.event->isBjet_MC_Bjet.at(k)){
                Hbb_true += eventInfo.bjet_momentums.at(k);
                ++totalMatchedBjets;
                indexes.push_back(k);
            }
        }
        if(totalMatchedBjets > 2)
            throw analysis::exception("Too many matched b-jets.");
        anaData.totalMatchedBjets(label).Fill(totalMatchedBjets);

        if (totalMatchedBjets < 2) return;

        anaData.Hbb_trueM(label).Fill(Hbb_true.M());
        const size_t matchedBjets = N_SelectedMatchedBjets(*eventInfo.event, eventInfo.selected_bjets);
        anaData.matchedBjets(label).Fill(matchedBjets);
//        TLorentzVector Hbb_CSV = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(1);
//        if (matchedBjets_byCSV < 2){
//            //std::cout << "csv fail (" << indexes.at(0) << "," << indexes.at(1) << ")" << std::endl;
//            anaData.Hbb_CSV_fail(category).Fill(Hbb_CSV.M());
//            analysis::FlatEventInfo::BjetPair pair(0,1);
//            const size_t indexFromPair = eventInfo.CombinationPairToIndex(pair,eventInfo.event->energy_Bjets.size());
//            if (eventInfo.event->kinfit_bb_tt_convergence.at(indexFromPair) > 0)
//                anaData.Chi2_CSVfail(category).Fill(eventInfo.event->kinfit_bb_tt_chi2.at(indexFromPair));
//            else
//                anaData.Chi2_CSVfail(category).Fill(1000);
//        }
//        else{
//            anaData.Hbb_CSV_good(category).Fill(Hbb_CSV.M());
//            analysis::FlatEventInfo::BjetPair pair(0,1);
//            const size_t indexFromPair = eventInfo.CombinationPairToIndex(pair,eventInfo.event->energy_Bjets.size());
//            if (eventInfo.event->kinfit_bb_tt_convergence.at(indexFromPair) > 0)
//                anaData.Chi2_CSVgood(category).Fill(eventInfo.event->kinfit_bb_tt_chi2.at(indexFromPair));
//            else
//                anaData.Chi2_CSVgood(category).Fill(1000);
//        }

//        std::vector<size_t> bjetsIndexes_byPt = SortBjetsByPt(eventInfo);
//        unsigned matchedBjets_byPt = MatchedBjetsByPt(eventInfo,bjetsIndexes_byPt);
//        anaData.matchedBjetsByPt(category).Fill(matchedBjets_byPt);


//        std::vector<size_t> bjetsCoupleIndexes_byChi2 = SortBjetsByChiSquare(eventInfo);
//        unsigned matchedBjets_byChi2 = MatchedBjetsByChi2(eventInfo,bjetsCoupleIndexes_byChi2);
//        anaData.matchedBjetsByChi2(category).Fill(matchedBjets_byChi2);

//        std::vector<size_t> bjetsCoupleIndexes_byMassPair = SortBjetsByMassPair(eventInfo);
//        unsigned matchedBjets_byMassPair = MatchedBjetsByMassPair(eventInfo,bjetsCoupleIndexes_byMassPair);
//        anaData.matchedBjets_byMassPair(category).Fill(matchedBjets_byMassPair);


//        analysis::FlatEventInfo::BjetPair pairCombined = GetBjetPairByMass(Hbb_CSV,eventInfo);
//        unsigned matchedBjets_combinedCSVandMASS = MatchedCombinedCSVandMASS(eventInfo,pairCombined);
//        anaData.matchedBjets_combinedCSVandMASS(category).Fill(matchedBjets_combinedCSVandMASS);

//        if (matchedBjets_combinedCSVandMASS < 2){
//            std::cout << "combined fail (" << indexes.at(0) << "," << indexes.at(1) << ")" << std::endl;
//        }

//        std::cout << "bjets: " <<eventInfo.event->energy_Bjets.size() << std::endl;
//        analysis::FlatEventInfo::BjetPair candidatePair = GetBjetPairByBestMass(eventInfo);
//        std::cout << "event: " << eventInfo.event->evt << ", bestPair(" << candidatePair.first << "," <<
//                     candidatePair.second << ")" << std::endl;
//        unsigned matchedBjets_CSVandBestMassPair = MatchedCombinedCSVandMASS(eventInfo,candidatePair);
//        anaData.matchedBjets_CSVandBestMassPair(category).Fill(matchedBjets_CSVandBestMassPair);

//        //only first 2 couples
//        analysis::FlatEventInfo::BjetPair pairCombined_2 = GetBjetPairByMass_2(Hbb_CSV,eventInfo);
//        unsigned matchedBjets_combinedCSVandMASS_2 = MatchedCombinedCSVandMASS(eventInfo,pairCombined_2);
//        anaData.matchedBjets_combinedCSVandMASS_2(category).Fill(matchedBjets_combinedCSVandMASS_2);

//        analysis::FlatEventInfo::BjetPair candidatePair_2 = GetBjetPairByBestMass_2(eventInfo);
//        unsigned matchedBjets_CSVandBestMassPair_2 = MatchedCombinedCSVandMASS(eventInfo,candidatePair_2);
//        anaData.matchedBjets_CSVandBestMassPair_2(category).Fill(matchedBjets_CSVandBestMassPair_2);
    }

private:
    static size_t N_SelectedMatchedBjets(const ntuple::Flat& event, const analysis::FlatEventInfo::BjetPair& bjet_pair)
    {
        size_t matchedBjets = 0;
        const std::vector<size_t> bjet_indexes = { bjet_pair.first, bjet_pair.second };
        for(size_t index : bjet_indexes) {
            if(index < event.isBjet_MC_Bjet.size() && event.isBjet_MC_Bjet.at(index))
                ++matchedBjets;
        }
        return matchedBjets;
    }

    template<typename Comparator>
    static std::vector<size_t> SortBjets(size_t n_bjets, const Comparator& comparator)
    {
        std::vector<size_t> bjet_indexes(n_bjets);
        std::iota(bjet_indexes.begin(), bjet_indexes.end(), 0);
        std::sort(bjet_indexes.begin(), bjet_indexes.end(), comparator);
        return bjet_indexes;
    }

    static analysis::FlatEventInfo::BjetPair DefaultPair() { return analysis::FlatEventInfo::BjetPair(0, 1); }

    template<typename Comparator>
    static analysis::FlatEventInfo::BjetPair SelectBestPair(const ntuple::Flat& event, const Comparator& comparator)
    {
        if(event.pt_Bjets.size() < 2) return DefaultPair();
        const auto bjet_indexes = SortBjets(event.pt_Bjets.size(), comparator);
        return analysis::FlatEventInfo::BjetPair(bjet_indexes.at(0), bjet_indexes.at(1));
    }


    analysis::FlatEventInfo::BjetPair SelectBestPtPair(const ntuple::Flat& event)
    {
        const auto comparator = [&] (size_t first, size_t second) -> bool
        {
            return event.pt_Bjets.at(first) > event.pt_Bjets.at(second);
        };

        return SelectBestPair(event, comparator);
    }

    analysis::FlatEventInfo::BjetPair SelectBestChi2Pair(const ntuple::Flat& event)
    {
        if(event.pt_Bjets.size() < 2) return DefaultPair();

        const auto comparator = [&] (size_t first, size_t second) -> bool
        {
            return event.pt_Bjets.at(first) > event.pt_Bjets.at(second);
        };

        return SelectBestPair(event, comparator);
    }


//    std::vector<size_t> SortBjetsByChiSquare(const analysis::FlatEventInfo& eventInfo)
//    {
//        std::vector<size_t> bjetsPair_indexes(eventInfo.event->kinfit_bb_tt_chi2.size());
//        std::iota(bjetsPair_indexes.begin(), bjetsPair_indexes.end(), 0);

//        if (eventInfo.event->kinfit_bb_tt_chi2.size() != eventInfo.event->kinfit_bb_tt_convergence.size())
//            std::cout << "aaa - size chi2 = " << eventInfo.event->kinfit_bb_tt_chi2.size() <<
//                         "- size convergence = " << eventInfo.event->kinfit_bb_tt_convergence.size() << std::endl;


//        const auto bjetsSelector = [&] (const size_t& first, const size_t& second) -> bool
//        {
//            if (eventInfo.event->kinfit_bb_tt_convergence.at(first) > 0
//                    && eventInfo.event->kinfit_bb_tt_convergence.at(second) <= 0)
//                return true;
//            if (eventInfo.event->kinfit_bb_tt_convergence.at(first) <= 0
//                    && eventInfo.event->kinfit_bb_tt_convergence.at(second) > 0)
//                return false;
//            Float_t first_chi2 = 0;
//            Float_t second_chi2 = 0;
//            if ((eventInfo.event->kinfit_bb_tt_convergence.at(first) > 0
//                    && eventInfo.event->kinfit_bb_tt_convergence.at(second) > 0) ||
//                    (eventInfo.event->kinfit_bb_tt_convergence.at(first) <= 0
//                     && eventInfo.event->kinfit_bb_tt_convergence.at(second) <= 0)){
//                first_chi2 = eventInfo.event->kinfit_bb_tt_chi2.at(first);
//                second_chi2 = eventInfo.event->kinfit_bb_tt_chi2.at(second);
//            }
//            return first_chi2 < second_chi2;
//        };

//        std::sort(bjetsPair_indexes.begin(), bjetsPair_indexes.end(), bjetsSelector);
//        return bjetsPair_indexes;
//    }

//    std::vector<size_t> SortBjetsByMassPair(const analysis::FlatEventInfo& eventInfo)
//    {
//        std::vector<size_t> bjetsPair_indexes(eventInfo.event->kinfit_bb_tt_chi2.size());
//        std::iota(bjetsPair_indexes.begin(), bjetsPair_indexes.end(), 0);

//        const auto bjetsSelector = [&] (const size_t& first, const size_t& second) -> bool
//        {
//            const analysis::FlatEventInfo::BjetPair first_bjetPair =
//                    analysis::FlatEventInfo::CombinationIndexToPair(first, eventInfo.event->energy_Bjets.size());
//            const analysis::FlatEventInfo::BjetPair second_bjetPair =
//                    analysis::FlatEventInfo::CombinationIndexToPair(second, eventInfo.event->energy_Bjets.size());
//            const TLorentzVector b1_firstPair = eventInfo.bjet_momentums.at(first_bjetPair.first);
//            const TLorentzVector b2_firstPair = eventInfo.bjet_momentums.at(first_bjetPair.second);
//            const TLorentzVector firstMbb = b1_firstPair + b2_firstPair;
//            const TLorentzVector b1_secondPair = eventInfo.bjet_momentums.at(second_bjetPair.first);
//            const TLorentzVector b2_secondPair = eventInfo.bjet_momentums.at(second_bjetPair.second);
//            const TLorentzVector secondMbb = b1_secondPair + b2_secondPair;
//            return std::abs(firstMbb.M() - 110) < std::abs(secondMbb.M() - 110);
//        };

//        std::sort(bjetsPair_indexes.begin(), bjetsPair_indexes.end(), bjetsSelector);
//        return bjetsPair_indexes;
//    }

//    analysis::FlatEventInfo::BjetPair GetBjetPairByMass(const TLorentzVector& Hbb_CSV,
//                                                        const analysis::FlatEventInfo& eventInfo)
//    {
//        analysis::FlatEventInfo::BjetPair pairCombined;
//        if ((Hbb_CSV.M() < 60 || Hbb_CSV.M() > 160) && eventInfo.event->energy_Bjets.size() > 2){
//            pairCombined.first = 0;
//            pairCombined.second = 2;
//            TLorentzVector mBB = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(2);
//            if ((mBB.M() < 60 || mBB.M() > 160) && eventInfo.event->energy_Bjets.size() > 3){
//                pairCombined.first = 0;
//                pairCombined.second = 3;
//            }
//            else {
//                pairCombined.first = 0;
//                pairCombined.second = 2;
//            }
//        }
//        else {
//            pairCombined.first = 0;
//            pairCombined.second = 1;
//        }
//        return pairCombined;
//    }

//    analysis::FlatEventInfo::BjetPair GetBjetPairByMass_2(const TLorentzVector& Hbb_CSV,
//                                                        const analysis::FlatEventInfo& eventInfo)
//    {
//        analysis::FlatEventInfo::BjetPair pairCombined;
//        if ((Hbb_CSV.M() < 60 || Hbb_CSV.M() > 160) && eventInfo.event->energy_Bjets.size() > 2){
//            pairCombined.first = 0;
//            pairCombined.second = 2;
//        }
//        else {
//            pairCombined.first = 0;
//            pairCombined.second = 1;
//        }
//        return pairCombined;
//    }

//    analysis::FlatEventInfo::BjetPair GetBjetPairByBestMass(const analysis::FlatEventInfo& eventInfo)
//    {
//        analysis::FlatEventInfo::BjetPair candidatePair(0,1);
//        if (eventInfo.event->energy_Bjets.size() > 3){
//            const TLorentzVector bbPair01 = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(1);
//            const TLorentzVector bbPair02 = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(2);
//            const TLorentzVector bbPair03 = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(3);
//            if (std::abs(bbPair01.M() - 110) < std::abs(bbPair02.M() - 110) &&
//                    std::abs(bbPair01.M() - 110) < std::abs(bbPair03.M() - 110)){
//                candidatePair.first = 0;
//                candidatePair.second = 1;
//            }
//            else if (std::abs(bbPair02.M() - 110) < std::abs(bbPair01.M() - 110) &&
//                     std::abs(bbPair02.M() - 110) < std::abs(bbPair03.M() - 110)){
//                candidatePair.first = 0;
//                candidatePair.second = 2;
//            }
//            else if (std::abs(bbPair03.M() - 110) < std::abs(bbPair01.M() - 110) &&
//                     std::abs(bbPair03.M() - 110) < std::abs(bbPair02.M() - 110)){
//                candidatePair.first = 0;
//                candidatePair.second = 3;
//            }
//        }
//        return candidatePair;
//    }

//    analysis::FlatEventInfo::BjetPair GetBjetPairByBestMass_2(const analysis::FlatEventInfo& eventInfo)
//    {
//        analysis::FlatEventInfo::BjetPair candidatePair(0,1);
//        if (eventInfo.event->energy_Bjets.size() > 2){
//            const TLorentzVector bbPair01 = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(1);
//            const TLorentzVector bbPair02 = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(2);

//            if (std::abs(bbPair01.M() - 110) < std::abs(bbPair02.M() - 110)){
//                candidatePair.first = 0;
//                candidatePair.second = 1;
//            }
//            else {
//                candidatePair.first = 0;
//                candidatePair.second = 2;
//            }
//        }
//        return candidatePair;
//    }

//    unsigned MatchedBjetsByPt(const analysis::FlatEventInfo& eventInfo, const std::vector<size_t>& indexes)
//    {
//        unsigned matchedBjets = 0;
//        const size_t first_index = indexes.at(0);
//        const size_t second_index = indexes.at(1);
//        //std::cout << "pt pair (" << first_index << "," << second_index << ")" << std::endl;
//        if (eventInfo.event->isBjet_MC_Bjet.at(first_index))
//            ++matchedBjets;
//        if (eventInfo.event->isBjet_MC_Bjet.at(second_index))
//            ++matchedBjets;
//        return matchedBjets;
//    }

//    unsigned MatchedBjetsByChi2(const analysis::FlatEventInfo& eventInfo, const std::vector<size_t>& indexes)
//    {
//        unsigned matchedBjets = 0;

//        analysis::FlatEventInfo::BjetPair bjetPair =
//                analysis::FlatEventInfo::CombinationIndexToPair(indexes.at(0),eventInfo.event->energy_Bjets.size());
//        //std::cout << "bjets chi2 pair (" << bjetPair.first << "," << bjetPair.second << ")" << std::endl;
//        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.first))
//            ++matchedBjets;
//        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.second))
//            ++matchedBjets;

//        return matchedBjets;
//    }

//    unsigned MatchedBjetsByMassPair(const analysis::FlatEventInfo& eventInfo, const std::vector<size_t>& indexes)
//    {
//        unsigned matchedBjets = 0;

//        analysis::FlatEventInfo::BjetPair bjetPair =
//                analysis::FlatEventInfo::CombinationIndexToPair(indexes.at(0),eventInfo.event->energy_Bjets.size());
//        //std::cout << "bjets mass pair (" << bjetPair.first << "," << bjetPair.second << ")" << std::endl;
//        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.first))
//            ++matchedBjets;
//        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.second))
//            ++matchedBjets;

//        return matchedBjets;
//    }

//    unsigned MatchedCombinedCSVandMASS(const analysis::FlatEventInfo& eventInfo,
//                                       const analysis::FlatEventInfo::BjetPair& bjetPair)
//    {
//        unsigned matchedBjets = 0;
//        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.first))
//            ++matchedBjets;
//        if (eventInfo.event->isBjet_MC_Bjet.at(bjetPair.second))
//            ++matchedBjets;
//        return matchedBjets;
//    }

private:
    BjetSelectionStudyData anaData;
};
