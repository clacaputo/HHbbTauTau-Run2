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
    BjetSelectionStudyData(TFile& outputFile) : root_ext::AnalyzerData(outputFile) {}

    TH1D_ENTRY(MET, 35, 0, 350)
    TH1D_ENTRY(matchedBjets, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjetsByCSV, 3, -0.5, 2.5)
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
        anaData.getOutputFile().cd();
    }

protected:
    //virtual analysis::FlatAnalyzerData& GetAnaData() override { return anaData; }

    virtual void AnalyzeEvent(const analysis::FlatEventInfo& eventInfo, analysis::EventCategory category) override
    {
        using analysis::EventCategory;
        if(!PassSelection(eventInfo)) return;
        anaData.MET(category).Fill(eventInfo.MET.Pt());
        if (category != EventCategory::TwoJets_ZeroBtag && category != EventCategory::TwoJets_OneBtag
                && category != EventCategory::TwoJets_TwoBtag) return;
        anaData.csv_b1_vs_ptb1(category).Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Pt(),
                                              eventInfo.event->csv_Bjets.at(eventInfo.selected_bjets.first));

        unsigned matchedBjets = 0;
        std::vector<unsigned> indexes;
        TLorentzVector Hbb_true;

        for (unsigned k = 0; k < eventInfo.event->energy_Bjets.size(); ++k){

            if (eventInfo.event->isBjet_MC_Bjet.at(k)){
                Hbb_true += eventInfo.bjet_momentums.at(k);
                ++matchedBjets;
                indexes.push_back(k);
            }
        }
        anaData.matchedBjets(category).Fill(matchedBjets);

//        if (matchedBjets < 2) return;
//        anaData.Hbb_trueM(category).Fill(Hbb_true.M());
//        unsigned matchedBjets_byCSV = MatchedBjetsByCSV(eventInfo);
//        anaData.matchedBjetsByCSV(category).Fill(matchedBjets_byCSV);
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

    bool PassSelection(const analysis::FlatEventInfo& eventInfo)
    {
        using analysis::EventRegion;
        const ntuple::Flat& event = *eventInfo.event;
        if (eventInfo.channel == analysis::Channel::MuTau){
            using namespace cuts::Htautau_Summer13::MuTau;

            return !(!event.againstMuonTight_2
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || event.pfRelIso_1 >= muonID::pFRelIso || event.q_1 * event.q_2 == +1);
        }
        if (eventInfo.channel == analysis::Channel::ETau){
            using namespace cuts::Htautau_Summer13::ETau;

            return !(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || event.pfRelIso_1 >= electronID::pFRelIso || event.q_1 * event.q_2 == +1);
        }
        if (eventInfo.channel == analysis::Channel::TauTau){
            using namespace cuts::Htautau_Summer13::TauTau::tauID;

            if(!event.againstElectronLooseMVA_2
                    || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= BackgroundEstimation::Isolation_upperLimit
                    || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= BackgroundEstimation::Isolation_upperLimit
                    || (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits
                        && event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits))
                return false;

            const bool os = event.q_1 * event.q_2 == -1;
            const bool iso = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                             event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < byCombinedIsolationDeltaBetaCorrRaw3Hits;

            if(eventInfo.bjet_momentums.size() < 2 || !eventInfo.fitResults.has_valid_mass)
                return false;

            using namespace cuts::massWindow;
            const double mass_tautau = eventInfo.event->m_sv_MC;
            const double mass_bb = eventInfo.Hbb.M();
            const bool inside_mass_window = mass_tautau > m_tautau_low && mass_tautau < m_tautau_high
                    && mass_bb > m_bb_low && mass_bb < m_bb_high;

            return os && iso && inside_mass_window;
        }
        throw analysis::exception("unsupported channel ") << eventInfo.channel;
    }

//    std::vector<size_t> SortBjetsByPt(const analysis::FlatEventInfo& eventInfo)
//    {
//        std::vector<size_t> bjet_indexes(eventInfo.bjet_momentums.size());
//        std::iota(bjet_indexes.begin(), bjet_indexes.end(), 0);

//        const auto bjetsSelector = [&] (const size_t& first, const size_t& second) -> bool
//        {
//            const TLorentzVector first_bjet = eventInfo.bjet_momentums.at(first);
//            const TLorentzVector second_bjet = eventInfo.bjet_momentums.at(second);
//            return first_bjet.Pt() > second_bjet.Pt();
//        };

//        std::sort(bjet_indexes.begin(), bjet_indexes.end(), bjetsSelector);
//        return bjet_indexes;
//    }

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

//    unsigned MatchedBjetsByCSV(const analysis::FlatEventInfo& eventInfo)
//    {
//        unsigned matchedBjets = 0;
//        if (eventInfo.event->isBjet_MC_Bjet.at(0))
//            ++matchedBjets;
//        if (eventInfo.event->isBjet_MC_Bjet.at(1))
//            ++matchedBjets;
//        return matchedBjets;
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
