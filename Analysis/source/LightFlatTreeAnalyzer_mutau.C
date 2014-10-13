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

#include "Analysis/include/LightBaseFlatTreeAnalyzer.h"

class LightFlatAnalyzerDataMuTau : public analysis::LightFlatAnalyzerData {
public:
    LightFlatAnalyzerDataMuTau(TFile& outputFile) : LightFlatAnalyzerData(outputFile) {}


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


};

class LightFlatTreeAnalyzer_mutau : public analysis::LightBaseFlatTreeAnalyzer {
public:
    LightFlatTreeAnalyzer_mutau(const std::string& _inputFileName, const std::string& _outputFileName)
         : LightBaseFlatTreeAnalyzer(_inputFileName,_outputFileName), anaData(*outputFile)
    {
        anaData.getOutputFile().cd();
    }

protected:

    virtual const std::string& ChannelName() override
    {
        static const std::string channelName = "muTau";
        return channelName;
    }

    virtual analysis::LightFlatAnalyzerData& GetAnaData() override { return anaData; }

    virtual void FullSelection(const analysis::FlatEventInfo& eventInfo, const analysis::EventCategory& eventCategory) override
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        if(!eventInfo.event->againstMuonTight_2
                || eventInfo.event->byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || eventInfo.event->mt_1 >= muonID::mt || eventInfo.event->pfRelIso_1 >= muonID::pFRelIso ||
                (eventInfo.event->q_1 * eventInfo.event->q_2) == +1)
            return;
        const std::string category = analysis::eventCategoryMapName.at(eventCategory);
        anaData.MET(category).Fill(eventInfo.MET.Pt());
        if (eventCategory == analysis::EventCategory::Inclusive || eventCategory == analysis::EventCategory::OneJet_OneBtag ||
                eventCategory == analysis::EventCategory::OneJet_ZeroBtag) return;

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

        if (matchedBjets < 2) return;
        anaData.Hbb_trueM(category).Fill(Hbb_true.M());
        unsigned matchedBjets_byCSV = MatchedBjetsByCSV(eventInfo);
        anaData.matchedBjetsByCSV(category).Fill(matchedBjets_byCSV);
        TLorentzVector Hbb_CSV = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(1);
        if (matchedBjets_byCSV < 2){
            //std::cout << "csv fail (" << indexes.at(0) << "," << indexes.at(1) << ")" << std::endl;
            anaData.Hbb_CSV_fail(category).Fill(Hbb_CSV.M());
            analysis::FlatEventInfo::BjetPair pair(0,1);
            const size_t indexFromPair = eventInfo.CombinationPairToIndex(pair,eventInfo.event->energy_Bjets.size());
            if (eventInfo.event->kinfit_bb_tt_convergence.at(indexFromPair) > 0)
                anaData.Chi2_CSVfail(category).Fill(eventInfo.event->kinfit_bb_tt_chi2.at(indexFromPair));
            else
                anaData.Chi2_CSVfail(category).Fill(1000);
        }
        else{
            anaData.Hbb_CSV_good(category).Fill(Hbb_CSV.M());
            analysis::FlatEventInfo::BjetPair pair(0,1);
            const size_t indexFromPair = eventInfo.CombinationPairToIndex(pair,eventInfo.event->energy_Bjets.size());
            if (eventInfo.event->kinfit_bb_tt_convergence.at(indexFromPair) > 0)
                anaData.Chi2_CSVgood(category).Fill(eventInfo.event->kinfit_bb_tt_chi2.at(indexFromPair));
            else
                anaData.Chi2_CSVgood(category).Fill(1000);
        }

        std::vector<size_t> bjetsIndexes_byPt = SortBjetsByPt(eventInfo);
        unsigned matchedBjets_byPt = MatchedBjetsByPt(eventInfo,bjetsIndexes_byPt);
        anaData.matchedBjetsByPt(category).Fill(matchedBjets_byPt);


        std::vector<size_t> bjetsCoupleIndexes_byChi2 = SortBjetsByChiSquare(eventInfo);
        unsigned matchedBjets_byChi2 = MatchedBjetsByChi2(eventInfo,bjetsCoupleIndexes_byChi2);
        anaData.matchedBjetsByChi2(category).Fill(matchedBjets_byChi2);

        std::vector<size_t> bjetsCoupleIndexes_byMassPair = SortBjetsByMassPair(eventInfo);
        unsigned matchedBjets_byMassPair = MatchedBjetsByMassPair(eventInfo,bjetsCoupleIndexes_byMassPair);
        anaData.matchedBjets_byMassPair(category).Fill(matchedBjets_byMassPair);


        analysis::FlatEventInfo::BjetPair pairCombined = GetBjetPairByMass(Hbb_CSV,eventInfo);
        unsigned matchedBjets_combinedCSVandMASS = MatchedCombinedCSVandMASS(eventInfo,pairCombined);
        anaData.matchedBjets_combinedCSVandMASS(category).Fill(matchedBjets_combinedCSVandMASS);

//        if (matchedBjets_combinedCSVandMASS < 2){
//            std::cout << "combined fail (" << indexes.at(0) << "," << indexes.at(1) << ")" << std::endl;
//        }

        std::cout << "bjets: " <<eventInfo.event->energy_Bjets.size() << std::endl;
        analysis::FlatEventInfo::BjetPair candidatePair = GetBjetPairByBestMass(eventInfo);
        std::cout << "event: " << eventInfo.event->evt << ", bestPair(" << candidatePair.first << "," <<
                     candidatePair.second << ")" << std::endl;
        unsigned matchedBjets_CSVandBestMassPair = MatchedCombinedCSVandMASS(eventInfo,candidatePair);
        anaData.matchedBjets_CSVandBestMassPair(category).Fill(matchedBjets_CSVandBestMassPair);


    }

private:
    LightFlatAnalyzerDataMuTau anaData;




};
