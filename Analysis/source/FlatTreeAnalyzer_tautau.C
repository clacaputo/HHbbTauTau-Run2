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

#include "Analysis/include/BaseFlatTreeAnalyzer.h"

class FlatAnalyzerData_tautau : public analysis::FlatAnalyzerData {
public:
    TH1D_ENTRY_EX(mt_1, 20, 0, 200, "M_{T}[GeV]", "Events", false, 1.1)
    TH1D_ENTRY_EX(iso_tau1, 100, 0, 10, "Iso#tau_{1}", "Events", false, 1)
    TH1D_ENTRY_EX(iso_tau2, 100, 0, 10, "Iso#tau_{2}", "Events", false, 1)

    virtual void Fill(const analysis::FlatEventInfo& eventInfo, double weight,
                      analysis::EventEnergyScale eventEnergyScale, analysis::EventSubCategory subCategory) override
    {
        FlatAnalyzerData::Fill(eventInfo, weight, eventEnergyScale, subCategory);
        if (eventEnergyScale != analysis::EventEnergyScale::Central) return;
        const ntuple::Flat& event = *eventInfo.event;
        const std::string key = HistogramSuffix(subCategory, eventEnergyScale);
        mt_1(key).Fill(event.mt_1, weight);
        iso_tau1(key).Fill(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1,weight);
        iso_tau2(key).Fill(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2,weight);
    }

    virtual const std::vector<double>& M_ttbb_Bins() const override
    {
        static const std::vector<double> bins = { 200, 250, 280, 310, 340, 370, 400, 500, 600, 700 };
        return bins;
    }
};

class FlatAnalyzerData_tautau_2tag : public FlatAnalyzerData_tautau {
public:
    virtual const std::vector<double>& M_tt_Bins() const override
    {
        static const std::vector<double> bins = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 250, 300, 350 };
        return bins;
    }
};

class FlatTreeAnalyzer_tautau : public analysis::BaseFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_tautau(const std::string& source_cfg, const std::string& _inputPath,
                            const std::string& outputFileName, const std::string& signal_list)
          : BaseFlatTreeAnalyzer(analysis::DataCategoryCollection(source_cfg, signal_list, ChannelId()), _inputPath,
                                 outputFileName)
    {
    }

protected:
    virtual analysis::Channel ChannelId() const override { return analysis::Channel::TauTau; }

    virtual std::shared_ptr<analysis::FlatAnalyzerData> MakeAnaData(analysis::EventCategory eventCategory) override
    {
        if(eventCategory == analysis::EventCategory::TwoJets_TwoBtag
                || eventCategory == analysis::EventCategory::TwoJets_TwoLooseBtag)
            return std::shared_ptr<analysis::FlatAnalyzerData>(new FlatAnalyzerData_tautau_2tag());
        return std::shared_ptr<analysis::FlatAnalyzerData>(new FlatAnalyzerData_tautau());
    }

    virtual analysis::EventRegion DetermineEventRegion(const ntuple::Flat& event,
                                                       analysis::EventCategory /*eventCategory*/) override
    {
        using analysis::EventRegion;
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

    virtual analysis::PhysicalValue CalculateQCDYield(const std::string& hist_name,
                                                      analysis::EventCategory eventCategory,
                                                      std::ostream& s_out) override
    {
        using analysis::EventCategory;
        static const analysis::EventCategorySet categories =
            { EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag };
        analysis::EventCategory refEventCategory = eventCategory;
        if(categories.count(eventCategory))
            refEventCategory = analysis::MediumToLoose_EventCategoryMap.at(eventCategory);

        //loose
//        const analysis::PhysicalValue yield_OSAntiIso =
//                CalculateYieldsForQCD(hist_name, refEventCategory, analysis::EventRegion::OS_AntiIsolated, s_out);
        //medium
        const analysis::PhysicalValue yield_OSAntiIso =
                CalculateYieldsForQCD(hist_name, eventCategory, analysis::EventRegion::OS_AntiIsolated, s_out);

        const analysis::PhysicalValue yield_SSIso =
                CalculateYieldsForQCD(hist_name, refEventCategory, analysis::EventRegion::SS_Isolated, s_out);

        const analysis::PhysicalValue yield_SSAntiIso =
                CalculateYieldsForQCD(hist_name, refEventCategory, analysis::EventRegion::SS_AntiIsolated, s_out);


        const auto iso_antiIso_sf = yield_SSIso / yield_SSAntiIso;
        s_out << eventCategory << ": QCD SF Iso / AntiIso = " << iso_antiIso_sf << ". \n";
//        if(refEventCategory == eventCategory)
        return yield_OSAntiIso * iso_antiIso_sf;

//        const analysis::PhysicalValue yield_SSAntiIso_Medium =
//                CalculateYieldsForQCD(hist_name, eventCategory, analysis::EventRegion::SS_AntiIsolated, s_out);

//        const auto medium_loose_sf = yield_SSAntiIso_Medium / yield_SSAntiIso;
//        s_out << eventCategory << ": QCD SF Medium b-tag / Loose b-tag = " << medium_loose_sf << ".\n";

//        return yield_OSAntiIso * iso_antiIso_sf * medium_loose_sf;
    }

    virtual void EstimateQCD(const std::string& hist_name, analysis::EventCategory eventCategory,
                             const analysis::PhysicalValue& yield) override
    {
        using analysis::EventCategory;
        static const analysis::EventCategorySet categories =
            { EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag };
        analysis::EventCategory refEventCategory = eventCategory;
        if(categories.count(eventCategory))
            refEventCategory = analysis::MediumToLoose_EventCategoryMap.at(eventCategory);

        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::Data);

        auto hist_shape_data = GetHistogram(refEventCategory, data.name, analysis::EventRegion::OS_AntiIsolated, hist_name);
        if(!hist_shape_data) return;

        TH1D& histogram = CloneHistogram(eventCategory, qcd.name, analysis::EventRegion::OS_Isolated, *hist_shape_data);
        std::string debug_info, negative_bins_info;
        SubtractBackgroundHistograms(histogram, refEventCategory, analysis::EventRegion::OS_AntiIsolated, qcd.name,
                                     debug_info, negative_bins_info);
        if(negative_bins_info.size())
            std::cerr << negative_bins_info;
        analysis::RenormalizeHistogram(histogram, yield, true);
    }

    virtual PhysicalValueMap CalculateWjetsYields(analysis::EventCategory eventCategory, const std::string& hist_name,
                                                  bool fullEstimate) override
    {
        using analysis::EventCategory;

        const analysis::DataCategoryPtrSet& wjets_mc_categories =
                dataCategoryCollection.GetCategories(analysis::DataCategoryType::WJets_MC);


        static const analysis::EventCategorySet NotEstimated_categories = {
            EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_OneLooseBtag,
            EventCategory::TwoJets_TwoLooseBtag, EventCategory::TwoJets_AtLeastOneBtag,
            EventCategory::TwoJets_AtLeastOneLooseBtag
        };

        PhysicalValueMap valueMap;
        if (fullEstimate || !NotEstimated_categories.count(eventCategory)) {
            for (analysis::EventRegion eventRegion : analysis::QcdRegions)
                valueMap[eventRegion] = CalculateFullIntegral(eventCategory, eventRegion, hist_name, wjets_mc_categories);
        }
        return valueMap;
    }

    virtual void CreateHistogramForZTT(analysis::EventCategory eventCategory, const std::string& hist_name,
                               const PhysicalValueMap& ztt_yield_map, bool useEmbedded) override
    {
        using analysis::DataCategory;
        using analysis::DataCategoryType;
        using analysis::EventCategory;

        const DataCategory& embedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded);
        const DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);
        const DataCategory& ZTT = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT);
        const DataCategory& ZTT_L = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_L);
        const DataCategory& TTembedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::TT_Embedded);

        static const analysis::EventCategorySet categoriesToRelax =
            { EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag };

        EventCategory shapeEventCategory = eventCategory;
        if (categoriesToRelax.count(eventCategory))
            shapeEventCategory = analysis::MediumToLoose_EventCategoryMap.at(eventCategory);

        for(const auto& eventRegionKey : ztt_yield_map) {
            const analysis::EventRegion eventRegion = eventRegionKey.first;
            const analysis::PhysicalValue& ztt_yield = eventRegionKey.second;
            auto ztt_l_hist = GetHistogram(eventCategory, ZTT_L.name, eventRegion, hist_name);
            const std::string embeddedName = useEmbedded && eventRegion == analysis::EventRegion::OS_Isolated ?
                                             embedded.name : ZTT_MC.name;

            auto embedded_hist = GetHistogram(shapeEventCategory, embeddedName, eventRegion, hist_name);
            auto TTembedded_hist = GetHistogram(shapeEventCategory, TTembedded.name, eventRegion, hist_name);

            if (embedded_hist){
                TH1D& ztt_hist = CloneHistogram(eventCategory, ZTT.name, eventRegion, *embedded_hist);
                if (TTembedded_hist && useEmbedded)
                    ztt_hist.Add(TTembedded_hist, -1);
                RenormalizeHistogram(ztt_hist, ztt_yield, true);
                if (ztt_l_hist)
                    ztt_hist.Add(ztt_l_hist);                
            }
            if (!embedded_hist && ztt_l_hist)
                CloneHistogram(eventCategory, ZTT.name, eventRegion, *ztt_l_hist);
        }
    }

//    virtual void CreateHistogramForZcategory(analysis::EventCategory eventCategory, const std::string& hist_name) override
//    {
//        const std::map<analysis::DataCategoryType, analysis::DataCategoryType> z_type_category_map = {
//            { analysis::DataCategoryType::ZL_MC, analysis::DataCategoryType::ZL },
//            { analysis::DataCategoryType::ZJ_MC, analysis::DataCategoryType::ZJ }
//        };

//        for (const auto& z_category : z_type_category_map){
//            const analysis::DataCategory& originalZcategory = dataCategoryCollection.GetUniqueCategory(z_category.first);
//            const analysis::DataCategory& newZcategory = dataCategoryCollection.GetUniqueCategory(z_category.second);

//            for(analysis::EventRegion eventRegion : analysis::AllEventRegions) {
//                auto z_hist_shape = GetHistogram(eventCategory, originalZcategory.name, eventRegion, hist_name);
//                if (z_hist_shape)
//                    CloneHistogram(eventCategory, newZcategory.name, eventRegion, *z_hist_shape);
//            }
//        }
//    }


    virtual void CreateHistogramForVVcategory(analysis::EventCategory eventCategory,const std::string& hist_name) override
    {
        const std::map<analysis::DataCategoryType, analysis::DataCategoryType> diboson_category_map = {
            { analysis::DataCategoryType::DiBoson_MC, analysis::DataCategoryType::DiBoson }
        };

        for (const auto& diboson_category : diboson_category_map){
            const analysis::DataCategory& originalVVcategory = dataCategoryCollection.GetUniqueCategory(diboson_category.first);
            const analysis::DataCategory& newVVcategory = dataCategoryCollection.GetUniqueCategory(diboson_category.second);

            PhysicalValueMap valueMap;

            for(analysis::EventRegion eventRegion : analysis::AllEventRegions) {
                auto vv_hist_yield = GetHistogram(eventCategory, originalVVcategory.name, eventRegion, hist_name);
                if (vv_hist_yield)
                    valueMap[eventRegion] = analysis::Integral(*vv_hist_yield,true);
            }

            static const analysis::EventCategorySet categoriesToRelax =
                { analysis::EventCategory::TwoJets_OneBtag, analysis::EventCategory::TwoJets_TwoBtag,
                  analysis::EventCategory::TwoJets_AtLeastOneBtag };
            const analysis::EventCategory shapeEventCategory = categoriesToRelax.count(eventCategory)
                    ? analysis::MediumToLoose_EventCategoryMap.at(eventCategory) : eventCategory;

            for(const auto& yield_iter : valueMap) {
                const analysis::EventRegion eventRegion = yield_iter.first;
                const analysis::PhysicalValue& yield = yield_iter.second;
                auto vv_hist_shape = GetHistogram(shapeEventCategory, originalVVcategory.name, eventRegion, hist_name);
                if (vv_hist_shape){
                    TH1D& vv_hist = CloneHistogram(eventCategory, newVVcategory.name, eventRegion, *vv_hist_shape);
                    RenormalizeHistogram(vv_hist,yield,true);
                }
            }
        }
    }

};
